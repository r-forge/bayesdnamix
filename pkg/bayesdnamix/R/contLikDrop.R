
#' @title contLikDrop
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description contLikDrop evaluates the marginal likelihood of the STR DNA mixture given some assumed model by integrate out parameters for given allele drop-out probability for each contributors.
#' @details The procedure are doing numerical integration to approximate the marginal probability by integrate over noisance parameters. Mixture proportions have flat prior.
#' 
#' The user may specify probability of drop-out for each contributors. 
#' 
#' The peak heights are scaled between [0,1] such that the peak heights are not accounted for in the models.
#' 
#' Function calls procedure in c++ by using the package Armadillo
#'
#' Models: model={'0'-unit weights,'1'-mixsep}
#'
#' @param nC Number of contributors in model.
#' @param mixData Evidence object with list elements adata[[i]] and hdata[[i]]. Each element has a loci-list with list-element 'i' storing qualitative data in 'adata' and quantitative data in 'hdata'.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects with list element [[s]]$adata[[i]]. The list element has reference-list with list-element 's' having a loci-list adata with list-element 'i storing qualitative data.
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. condOrder=-1 means the reference is known-non contributor!
#' @param prD A vector of allele drop-out probabilities (p_hom=p_het^2) for each contributors. Default is 0.
#' @param prC A numeric for allele drop-in probability. Default is 0.
#' @param model A integer for specification of model. See details for more information.
#' @param pTau Prior function for tau-parameter. Flat prior is default.
#' @param taumax Maximum range of tau-parameter. Default is 100.
#' @param maxeval Maxumum number of evaluations in the interale function.
#' @param threshT The analytical threshold given. Used when considering allele drop-outs.
#' @param fst is the coancestry coeffecient. Default is 0.
#' @return lik Marginalized likelihood of the hypothesis (model) given observed evidence.
#' @keywords continuous, Bayesian models

contLikDrop = function(nC,mixData,popFreq,refData=NULL,condOrder=NULL,prD=NULL,prC=0,model=1,pTau=function(x) { return(1)},taumax=100, maxeval=10000,threshT=50,fst=0){
 require(cubature)
 if(is.null(prD)) prD <- rep(0,nC) #no dropout are default
 #call on the simpler function if all 
 if(all(c(prD,prC,fst)==0)) contLik(nC,mixData,popFreq,refData,condOrder,model,pTau,taumax, maxeval,threshT)

 nA = unlist(lapply(mixData$adata,length)) #number of alleles of selected loci
 if(max(nA)>(2*nC)) {
  stop(paste('Max alleles in a locus is ',max(nA),'. You must specify a greater number of contributors',sep=''))
 }
 nL <- length(nA) #number of loci in data
 locnames <- toupper(names(mixData$adata))
 missind <- !locnames%in%toupper(names(popFreq))  #indice of missing loci
 if(any(missind)) {
  stop(paste('Missing locus (',locnames[missind],') in popFreq.',sep=''))
 }

 #convertion of values in popFreq, mixData and Glist$G:
 #loci-order follows as in mixData: "locnames". Rearrange names:
 names(popFreq) <- toupper(names((popFreq))) #make invariant to capital
 popFreq <- popFreq[locnames] #order popFreq to mixData-order

 #get population genotypes:
 getGlist <- function(popFreq) {
  locs <- names(popFreq)
  Glist <- list()
  for(i in 1:length(locs)) {
   G = t(as.matrix(expand.grid(rep(list(as.numeric(names(popFreq[[i]])),as.numeric(names(popFreq[[i]])) ))))) #one genotype per column
   keep = G[2,]>=G[1,] #unique genotypes 
   G <- G[,keep]  #store genotypes
   G <- matrix(as.character(G),nrow=2) #make string names again
   tmpP = t(as.matrix(expand.grid(rep(list(as.numeric(popFreq[[i]]),as.numeric(popFreq[[i]]) )))))
   Gprob = exp(colSums(log(tmpP[,keep]))) #get allele probs
   ishet = G[1,]!=G[2,]
   Gprob[ishet] = 2*Gprob[ishet] #multiply with two to get heterozygote prob
   Glist[[locs[i]]] <- list(G=t(G),Gprob=Gprob)
  }
  return(Glist)
 }

 #Fix references: Assign condition to condM-matrix
 Glist <- getGlist(popFreq) #get population genotype information
 condM <- matrix(-1,nrow=nL,ncol=nC) #default is no references (=-1)
 #assign references to condM-matrix by values of Glist
 if(!is.null(refData) && !is.null(condOrder) && any(condOrder>0)) {
  for(i in 1:nL) {
   subRef <- refData[[locnames[i]]]
   Gset <- Glist[[locnames[i]]]$G #genotype combinations
   if(length(subRef)==0)  stop(paste('Missing locus (',locnames[i],') in refData.',sep=''))
   for(k in 1:length(subRef)) { #for each reference
    if(condOrder[k]>0) {
     Gind1 <- subRef[[k]][1]==Gset[,1] & subRef[[k]][2]==Gset[,2]
     Gind2 <- subRef[[k]][2]==Gset[,1] & subRef[[k]][1]==Gset[,2]
     condM[i,condOrder[k]] = which(Gind1 | Gind2) - 1 #subtract with one since we work from 0-indice
    }
   }
  }
 }

 #fix known sample-information:
 mkvec <- numeric()
 nkval <- rep(0,nL)
 for(i in 1:nL) {
  tmp <- rep(0, length(popFreq[[i]]))
  if(!is.null(condOrder) & !is.null(refData)) {
   for(k in 1:length(condOrder)) {
    if(condOrder[k]!=0) { 
     ind <- which( names(popFreq[[i]])%in%refData[[locnames[i]]][[k]] )
     tmp[ind] = (length(ind)==1) + 1
    }
   }
  }
  nkval[i] <- sum(tmp) #number of sampled (for each loci)
  mkvec <- c(mkvec,tmp) 
 }

 for(i in 1:nL) {
   anames <- names(popFreq[[i]]) #old names
   anames2 <- 0:(length(popFreq[[i]])-1) #new names
   names(popFreq[[i]]) <- anames2
   for(j in 1:nA[i]) { #for each allele
    mixData$adata[[i]][j] <- anames2[anames==mixData$adata[[i]][j]] #update name
   }
  }

 #fix genotypes:
 Glist <- getGlist(popFreq) #get population genotype information
 Gprob <- lapply(Glist,function(x) return(x$Gprob))
 Gset <- lapply(Glist,function(x) return(x$G))


 #Fix input-variables for C-function
 PE <- 1
 CnA <- c(0,cumsum(nA))
 allA <- as.integer(unlist(mixData$adata))
 allY <- as.numeric(unlist(mixData$hdata))
 sY <- rep(1,nL) #we use scaled models
 nG <- sapply(Gprob,length) #number of genotype combinations
 CnG <- c(0,cumsum(nG))
 CnG2 <- c(0,cumsum(nG*2)) #note: 2 columns for each genotype!!
 pG <- unlist(Gprob) #vectorize over all loci
 Gvec <- as.integer(rbind(unlist(Gset))) #vectorize a big matrix (loci are put chronologic)
 condRef <- c(condM) #vectorized over all loci
 nAall <- sapply(popFreq,length) #Number of population-alleles on each loci
 CnAall <- c(0,cumsum(nAall)) #cumulative number of alleles
 pA <- unlist(popFreq) #need each allele probability for drop-in probabilities
 t0 <- threshT #peak height imputed threshold


 likYtheta <- function(theta) {   #call c++- function: length(theta)=nC
  val <- .C("contlikdropC",as.numeric(PE),as.numeric(theta),as.integer(model),as.integer(nC),as.integer(nL),as.integer(nA), as.numeric(allY),as.integer(allA),as.integer(CnA),as.numeric(sY),as.integer(nAall),as.integer(CnAall),as.integer(Gvec),as.integer(nG),as.integer(CnG),as.integer(CnG2),as.numeric(pG),as.numeric(pA),as.numeric(prD), as.numeric(prC), as.integer(condRef),as.numeric(t0), as.numeric(fst),as.integer(mkvec),as.integer(nkval),PACKAGE="bayesdnamix")[[1]]
  return(val*pTau(theta[nC]))
 }

 #END taken from deconvolve()
 if(model==0) { #binary model contains 0 paramaters.
  pTau=function(x) { return(1)} #be sure of no scaling
  lik <- likYtheta(c(rep(1,nC))) #just send a default value ones
 } else { #C param
  lik <- adaptIntegrate(likYtheta, lowerLimit = rep(0,nC), upperLimit = c(rep(1,nC-1),taumax),maxEval = maxeval )[[1]]
 }
 return(lik)
}

