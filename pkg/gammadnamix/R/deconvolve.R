#' @title deconvolve
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description deconvolve ranks the set of the most conditional posterior probability of genotypes the STR DNA mixture given a fitted model under a hypothesis.
#' @details The procedure calculates the likelihood for each single locus. Then it combines the most probable genotypes from each loci to produce a ranked list of deconvolved profiles.
#' 
#' Function calls procedure in c++ by using the package Armadillo and Boost.
#'
#' @param mlefit Fitted object using contLikMLE function.
#' @param alpha Required sum of the listed posterior probabilities.
#' @param maxlist The ranked deconvolved profile list will not exceed this number (used to avoid endless search).
#' @param unknownonly A boolean whether table should only contain the unknown genotypes or both known and unknown genotypes.
#' @return ret A list(table1,rankG,pG) where rankG is the ranked genotypes with corresponding probabilities in pG. table1 is a formated version of these two.
#' @export
#' @references Cowell,R.G. et.al. (2014). Analysis of forensic DNA mixtures with artefacts. Applied Statistics, 64(1),1-32.
#' @keywords deconvolution
deconvolve = function(mlefit,alpha=0.95,maxlist=1000,unknownonly=TRUE){
 theta <- mlefit$fit$thetahat #condition on mle parameter
 model <- mlefit$model #take out assumed model with given data
 locs <- names(model$popFreq)
 nL <- length(locs)
 nC <- model$nC
 np <- length(theta)#number of unknown parameters
 loglikYtheta <- function() {   #call c++- function: length(theta)=nC+1
   Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(model$prC), ret$condRef,as.numeric(model$threshT),as.numeric(model$fst),ret$mkvec,ret$nkval,as.numeric(model$lambda),as.integer(0),PACKAGE="gammadnamix")[[1]]
   return(Cval + log(model$pXi(theta[ret$nC+2]))) #weight with prior of tau and 
 }
 if(!is.null(model$xi)) {
  theta <- c(theta,model$xi)
  loglikYtheta <- function() {   #call c++- function: length(theta)=nC
   Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(model$prC), ret$condRef,as.numeric(model$threshT),as.numeric(model$fst),ret$mkvec,ret$nkval,as.numeric(model$lambda),as.integer(0),PACKAGE="gammadnamix")[[1]]
   return(Cval)
  }
 }

 #Using information in ret to try out different genotypes:

 #Step 1) Calculate L(E|g,thetahat) for each marker
 dlist <- list()
 GClist <- list()
 for(loc in locs) {
  samples <- lapply(model$samples,function(x) x[loc])
  ret <- prepareC(nC=nC,samples,popFreq=model$popFreq[loc],refData=model$refData[loc],condOrder=model$condOrder,knownRef=model$knownRef)
  uind <- which(ret$condRef==-1) #unknown genotype indices
  nU <- length(uind) #number of unknowns
  if(nU==0) stop("There was no unknown genotype profiles to estimate. The evaluation will not be done!")
  ret$nK <- nC  #number of known is equal number of contributors
  Gset <- matrix(ret$Gvec,ncol=2) #genotype possibilities
  Glist <- list()
  for(k in 1:nU) {
   Glist[[k]] <- 1:nrow(Gset)
  }
  combGind <- expand.grid(Glist) #get all combinations
  combGind <- as.matrix(combGind,nrow=nrow(combGind))
  dvec <- rep(NA,nrow(combGind))
  #calculate for each genotypes:
  for(gind in 1:nrow(combGind)) { #for each possible genotypes:
   ret$condRef[uind] <- as.integer(combGind[gind,] - 1) #genotypes to consider
   dvec[gind] <- loglikYtheta() 
   dvec[gind] <- dvec[gind] + sum(log(ret$pG[combGind[gind,]])) #add genotype probability as well
  }
  isOK <- !is.infinite(dvec)
  combGind <- combGind[isOK,]
  dvec <- dvec[isOK]
  rank <- order(dvec,decreasing=TRUE)
  dlist[[loc]] <- dvec[rank] 

  if(is.null(dim(combGind))) { #threat the case of one unknown
   GClist[[loc]] <- as.matrix(combGind[rank]) 
  } else {
   GClist[[loc]] <- combGind[rank,]
  }
 }

 #Step 2) Combine markers to create full profiles (this is it's own function):
 rankGlist <- combineRank(dlist,loghdval=mlefit$fit$loglik,alpha=alpha,maxsearch=maxlist)
 pG <- rankGlist$pG
 Gset  <- rankGlist$rankG

 kvec <- 1:nC
 if(unknownonly) kvec <- uind  

 #Step 3) Convert rank-list to list with allele-names
 Glist <- getGlist(model$popFreq) #get genotype list with genotypes and corresponding frequencies
 deconvlist <- list()
 for(i in 1:nL) { #convert names for each locus
  rankgeno <- GClist[[locs[i]]][Gset[,i],]
  if(is.null(dim(rankgeno))) rankgeno <- as.matrix(rankgeno) #make matrix again
  rankgenos <- numeric()
  for(k in 1:nC) { #for each contributor
   if(k%in%uind) { #if unknown contributor
    geno <- Glist[[locs[i]]]$G[rankgeno[,which(k==uind)],] #get allele named genotype
   } else if(!unknownonly) { #if known contributors in addition(they are given in reference) 
    geno <- sort(model$refData[[locs[i]]][[which(model$condOrder==k)]])
    geno <- matrix( rep(geno,nrow(rankgeno)),ncol=2,byrow=TRUE)
   } else {
    next #skip to next contributor
   }
   geno <- paste0(geno[,1],"/",geno[,2])
   rankgenos <- cbind(rankgenos,geno)
  }
  colnames(rankgenos) <- paste0("g",kvec)
  deconvlist[[locs[i]]] <- rankgenos 
 }

 #Step4) Create table layouts:
 table1 <- numeric()
 for(loc in locs) { #convert names for each locus
  table1 <- cbind(table1,deconvlist[[loc]])
 }
 table1 <-  cbind(table1,pG)
 colnames(table1) <- c(paste0(c(t(replicate(length(kvec),locs))),"_g",kvec),"posterior")
 return(list(table1=table1,rankG=deconvlist,pG=pG))
} #end function
