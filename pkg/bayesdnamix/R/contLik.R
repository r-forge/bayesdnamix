
#' @title contLik
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description contLik evaluates the marginal likelihood of the STR DNA mixture given some assumed model by integrate out parameters.
#' @details The procedure are doing numerical integration to approximate the marginal probability by integrate over noisance parameters. Mixture proportions have flat prior.
#' 
#' The peak heights are scaled between [0,1] such that the peak heights are not accounted for in the models.
#' 
#' Function calls procedure in c++ by using the package Armadillo
#'
#' @param nC Number of contributors in model.
#' @param mixData Evidence object with list elements adata[[i]] and hdata[[i]]. Each element has a loci-list with list-element 'i' storing qualitative data in 'adata' and quantitative data in 'hdata'.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects with list element [[s]]$adata[[i]]. The list element has reference-list with list-element 's' having a loci-list adata with list-element 'i storing qualitative data.
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model.
#' @param model A integer for specification of model. See details for more information.
#' @param pTau Prior function for tau-parameter. Flat prior is default.
#' @param taumax Maximum range of tau-parameter. Default is 100.
#' @param maxeval Maxumum number of evaluations in the interale function.
#' @param threshT The analytical threshold given. Used when considering allele drop-outs.
#' @param relaxq A parameter for model 6 which allows for relaxation of variance.
#' @return lik Marginalized likelihood of the hypothesis (model) given observed evidence.
#' @keywords continuous, Bayesian models

contLik = function(nC,mixData,popFreq,refData=NULL,condOrder=NULL,model=1,pTau=function(x) { return(1)},taumax=100, maxeval=10000,threshT=50,relaxq=0.2){
 require(cubature)

##########HELPFUNCTIONS############
  #Function that get table of allele-counts from Q
  getPtab = function(Q) {
   g1 = c(t(replicate(2,1:length(Q))))
   g2 = unlist(strsplit(Q,''))
   gdat = rep(1,length(Q)*2)
   Ptab = tapply(X=gdat, INDEX=list(g1, g2), FUN=sum)
   Ptab[is.na(Ptab)] = 0 #must insert 0 manually
   return(t(Ptab))
  }
  #get population genotypes:
  getGlist <- function(popFreq) {
   locs <- names(popFreq)
   nL <- length(locs)
   Glist <- list()
   for(i in 1:nL) {
    G = t(as.matrix(expand.grid(rep(list(as.numeric(names(popFreq[[i]])),as.numeric(names(popFreq[[i]])) ))))) #one genotype per column
    keep = G[2,]>=G[1,] #unique genotypes 
    G <- G[,keep]  #store genotypes
    G <- matrix(as.character(G),nrow=2) #make string names again
    tmpP = t(as.matrix(expand.grid(rep(list(as.numeric(popFreq[[i]]),as.numeric(popFreq[[i]]) )))))
    Gprob = exp(colSums(log(tmpP[,keep]))) #get allele probs
    ishet = G[1,]!=G[2,]
    Gprob[ishet] = 2*Gprob[ishet] #multiply with two to get heterozygote prob
    Glist[[locs[i]]] <- list(G=G,Gprob=Gprob)
   }
   return(Glist)
  }

  #function which returns combinations
  getQlist <- function(mixData,refData,condOrder) {
   Qcombs <- QcombsN <- list()
   nL <- length(mixData$adata)  
   for(i in 1:nL) {
    Alist <- mixData$adata[[i]] #==names(popFreq[[i]])?
    contr_combs <- Qcombs[[i]] <- getContrCombs(Alist,nC,symmetry=TRUE,refs=refData[[i]],condOrder,complete=TRUE) 
    ascii <- sort(unique(unlist(strsplit(contr_combs,"")))) 
    nCombs <- nrow(contr_combs)
    QcombsN[[i]] <- list()
    for(j in 1:nCombs) {  #for each genotype combination
     Qi = contr_combs[j,]
     cc <- lapply(strsplit(Qi,""),function(x) {
      tmp <- Alist[ascii%in%x]
      if(length(tmp)==1) tmp = rep(tmp,2)
      return(tmp)
     })
     QcombsN[[i]][[j]] <- t(matrix(unlist(cc),nrow=2))
     colnames(QcombsN[[i]][[j]]) <- paste("a",1:2,sep="")
     rownames(QcombsN[[i]][[j]]) <- paste("c",1:nC,sep="")
    }#end j
   } #end i
   return(list(Qcombs=Qcombs,QcombsN=QcombsN))
  } #end function

 #predefine all possible Q in Xlist,Olist. Also calculate genoProb
 getCombLists <- function(Qlist,Glist,nC,nAi,sY,condOrder) {
  maxQicomb = max(sapply(Qlist$Qcombs,nrow))
  pG <- matrix(nrow=nL,ncol=maxQicomb)
  Xlist <- Plist <- matrix(list(),nrow=nL,ncol=maxQicomb)
  popG <- list()
  for(i in 1:nL) {
   for(j in 1:nrow(Qlist$Qcombs[[i]])) {  #for each genotype combination
      Qi = Qlist$Qcombs[[i]][j,]
      Tmp <- Qlist$QcombsN[[i]][[j]]
      pG[i,j] <- 1
      for(k in 1:nC) {
        if(!any(condOrder==k)) { #if no restriction
         a <- Glist[[i]]$G[1,]==Tmp[k,1] & Glist[[i]]$G[2,]==Tmp[k,2]
         b <- Glist[[i]]$G[1,]==Tmp[k,2] & Glist[[i]]$G[2,]==Tmp[k,1]
         pG[i,j] <- pG[i,j]*Glist[[i]]$Gprob[a | b] #genotype probabilities
        }
      }
      #factorial(length(Qi))/prod( factorial(table(Qi)) ) #permutation factor:
      Pi = matrix(getPtab(Qi),ncol=nC,nrow=nAi[i]) #assign matrix
      Pitilde <- matrix(Pi[,-nC],ncol=nC-1) #(n_i X (nC-1)) 
      Pic <- as.matrix(Pi[,nC]) 
      Xlist[[i,j]] <- (Pitilde - Pic%*%t(matrix(1,nrow=nC-1)))
      Xlist[[i,j]] <- cbind(Xlist[[i,j]],Pic) #include offset vector 
      Plist[[i,j]] <- Pi
   } #end comb:j
  } #end loci:i
  return(list(Xlist=Xlist,Plist=Plist,pG=pG))
 }

##########END HELPFUNCTIONS############

##########Start function here:############
  #Insert missing allele into response data:
  #taken from deconvolve()
   locinames <- names(popFreq)
   nL <- length(locinames)

   if (!is.null(condOrder)) {
    for (i in 1:nL) {
      for (k in 1:length(condOrder)) {
        if (condOrder[k] == 0 | length(refData[[i]][[k]]) == 0)   next
          aref = refData[[i]][[k]]
          anew = aref[!aref %in% mixData$adata[[i]]]
          if (length(anew) > 0) {
           mixData$adata[[i]] = c(mixData$adata[[i]],anew)
           mixData$hdata[[i]] = c(mixData$hdata[[i]],rep(threshT, length(anew)))
           #print(paste("WARNING: At locus ",locinames[i],", the allele(s) ", paste(anew, collapse = "/", sep = ""), " was added  with threshold height ", threshT, sep = ""))
          }
        }
     }
  }
  #need to check if number of unique reference samples extends the model
 if(!is.null(refData) && !is.null(condOrder)) {
  nR = sum(condOrder>0)
  for(i in 1:nL) {
   refA = numeric()
   for(k in 1:length(condOrder)) {
    if(condOrder[k]==0) next
    refA = c(refA, refData[[i]][[k]] )
   }
   #get number of references on locus i:
   refA = unique(refA) #unique refererence alleles (uncselected are not counted)
   Ai = mixData$adata[[i]]
   leftA = Ai[!Ai%in%refA] #undescribed alleles for unknown
   if(length(leftA)>2*(nC-nR)) { 
    msg <- paste('For locus ',locinames[i],', number of unique allele left (after restriction) is ',length(leftA), ', while number of unknowns are ',nC-nR,'. Specify more unknowns or change reference conditioning.',sep='')
    stop(msg)
   }
  }
 }
 nA = unlist(lapply(mixData$adata,length)) #number of alleles of selected loci
 if(max(nA)>(2*nC)) {
    msg <- paste('Max alleles in a locus is ',max(nA),'. You must specify a greater number of contributors',sep='')
    stop(msg)
 }
 #END taken from deconvolve()
  for(i in 1:nL) mixData$hdata[[i]] <- mixData$hdata[[i]]/sum(mixData$hdata[[i]]) #scaled
  Qlist <-  getQlist(mixData,refData,condOrder)
  Glist <- getGlist(popFreq)
  sY = sapply(mixData$hdata,sum)
  lists <- getCombLists(Qlist,Glist,nC,nA,sY,condOrder) 
  nQi <- sapply(Qlist$Qcombs,nrow) #number of combs
  cdfX <- cumsum( c(0,rep(nA*nC,nQi)))
  allX <- unlist(t(lists$Xlist)) #should be (Cx1) matrix
  if(model==6)   allX <- unlist(t(lists$Plist)) #using P-matrix instead!!
  allY <- unlist(mixData$hdata)
  pG <- c(t(lists$pG))
  pG <- pG[!is.na(pG)]
  #get cumulative indices
  CnQ <- cumsum(c(0,nQi)) #start indices as 0
  CnA <- cumsum(c(0,nA))  #start indices as 0
  CnA2 <- cumsum(c(0,nA^2)) #start indices as 0
  PE = 1  

 likYtheta <- function(theta) {   #call c++- function: length(theta)=nC
  val <- .C("contlikC",as.numeric(PE),as.numeric(theta),as.integer(model),as.integer(nC),as.integer(nL),as.integer(nA), as.integer(nQi),as.numeric(pG),as.integer(CnQ),as.numeric(allY), as.integer(CnA),as.numeric(allX),as.integer(cdfX), as.numeric(sY),PACKAGE="bayesdnamix")[[1]]
  return(val*pTau(theta[nC]))
 }

# if(model==6) { #q=0.2 used as default
#  wrapperM6 <- function(th) likYtheta(c(th[1:nC-1],relaxq,th[nC]))
#  lik<- adaptIntegrate(wrapperM6, lowerLimit = rep(0,nC), upperLimit = c(rep(1,nC-1),taumax),maxEval = maxeval )[[1]]
# }

 #END taken from deconvolve()
 if(model==0) { #binary model contains 0 paramaters.
  pTau=function(x) { return(1)} #be sure of no scaling
  lik <- likYtheta(c(rep(1,nC))) #just send a default value ones
 } else { #C param
  lik <- adaptIntegrate(likYtheta, lowerLimit = rep(0,nC), upperLimit = c(rep(1,nC-1),taumax),maxEval = maxeval )[[1]]
 }


 return(lik)
}

