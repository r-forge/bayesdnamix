
#' @title contLikMLE
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description contLikMLE optimizes the likelihood of the STR DNA mixture given some assumed a bayesian model.
#' @details The procedure are doing numerical optimization to approximate the marginal probabilit over noisance parameters. Mixture proportions have flat prior.
#' 
#' Function calls procedure in c++ by using the package Armadillo and Boost.
#'
#'
#' @param nC Number of contributors in model.
#' @param mixData Evidence object with list elements adata[[i]] and hdata[[i]]. Each element has a loci-list with list-element 'i' storing qualitative data in 'adata' and quantitative data in 'hdata'.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects with list element [[s]]$adata[[i]]. The list element has reference-list with list-element 's' having a loci-list adata with list-element 'i storing qualitative data.
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. condOrder=-1 means the reference is known-non contributor!
#' @param xi A numeric giving stutter-ratio if it is known. Default is NULL, meaning it is integrated out.
#' @param prC A numeric for allele drop-in probability. Default is 0.
#' @param pRhoTau Prior function for (rho,tau)-parameters. Flat prior is default.
#' @param rhotaumax Maximum range of rho and tau-parameter. Default is c(1000,1000).
#' @param nDone Maxumum number of random evaluations nlm-optimizing routing. Default is 1.
#' @param threshT The detection threshold given. Used when considering probability of allele drop-outs.
#' @param fst is the coancestry coeffecient. Default is 0.
#' @param lambda Parameter in modeled peak height shifted exponential model. Default is 0.
#' @param pXi Prior function for xi-parameter (stutter). Flat prior on [0,1] is default.
#' @param M Number of evaluations for importance sampling.
#' @return ret A list(MLE,Sigma,CI,loglik,margL) with Maximixed likelihood elements for hypothesis (model) given observed evidence.
#' @keywords continuous, Bayesian models

contLikMLE = function(nC,mixData,popFreq,refData=NULL,condOrder=NULL,xi=NULL,prC=0,pRhoTau=function(x) {1},rhotaumax=c(1000,1000),nDone=1,threshT=50,fst=0,lambda=0,pXi=function(x) {1},M=NULL ){
 nA = unlist(lapply(mixData$adata,length)) #number of alleles of selected loci
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

 #decode allele-names to index names + Stutter-preparation
 allAbpind <- as.numeric()
 for(i in 1:nL) {
   anames <- names(popFreq[[i]]) #old names
   #find corresponding index of 1-ahead pb-twin
   numAnames <- as.numeric(anames)
   BPind <- rep(0,length(numAnames)) #init 0 if allele doesn't have any bp-twin
   for(j in 1:length(numAnames)) {
    bpind <- which((numAnames[j]-1)==numAnames) #store index of what allele it stutters to
    if(length(bpind)>0) BPind[j] <- bpind #assign index-placement
   }
   allAbpind <- c(allAbpind,BPind)
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


 #Fix input-variables for Cpp-function
 logPE <- 1
 CnA <- c(0,cumsum(nA))
 allA <- as.integer(unlist(mixData$adata))
 allY <- as.numeric(unlist(mixData$hdata))
 sY <- sapply(mixData$hdata,sum)
 nG <- sapply(Gprob,length) #number of genotype combinations
 CnG <- c(0,cumsum(nG))
 CnG2 <- c(0,cumsum(nG*2)) #note: 2 columns for each genotype!!
 pG <- unlist(Gprob) #vectorize over all loci
 Gvec <- as.integer(rbind(unlist(Gset))) #vectorize a big matrix (loci are put chronologic)
 condRef <- c(condM) #vectorized over all loci
 nAall <- sapply(popFreq,length) #Number of population-alleles on each loci
 CnAall <- c(0,cumsum(nAall)) #cumulative number of alleles
 pA <- unlist(popFreq) #need each allele probability for drop-in probabilities



#Optimize MLE:
 #Two cases: Integrate over Stutter or Stutter known
 negloglikYtheta <- function(theta) {   #call c++- function: length(theta)=nC+1
  if(any(theta<lower | theta>upper)) return(Inf)
  if(sum(theta[1:(nC-1)])>1) return(Inf)
  Cval  <- .C("loglikgammaC",as.numeric(logPE),as.numeric(theta),as.integer(nC),as.integer(nL),as.integer(nA), as.numeric(allY),as.integer(allA),as.integer(CnA),as.integer(allAbpind),as.integer(nAall),as.integer(CnAall),as.integer(Gvec),as.integer(nG),as.integer(CnG),as.integer(CnG2),as.numeric(pG),as.numeric(pA), as.numeric(prC), as.integer(condRef),as.numeric(threshT),as.numeric(fst),as.integer(mkvec),as.integer(nkval),as.numeric(lambda),PACKAGE="gammadnamix")[[1]]
  loglik <- Cval + log(pRhoTau(theta[c(nC,nC+1)])) + log(pXi(theta[nC+2])) #weight with prior of tau and 
  return(-loglik) #weight with prior of tau and stutter.
 } 
 negloglikYthetaS <- function(theta2) {   #call c++- function: length(theta)=nC
  theta <- c(theta2,xi) #stutter-parameter added as known
  if(any(theta<lower | theta>upper)) return(Inf)
  if(sum(theta[1:(nC-1)])>1) return(Inf)
  Cval  <- .C("loglikgammaC",as.numeric(logPE),as.numeric(theta),as.integer(nC),as.integer(nL),as.integer(nA), as.numeric(allY),as.integer(allA),as.integer(CnA),as.integer(allAbpind),as.integer(nAall),as.integer(CnAall),as.integer(Gvec),as.integer(nG),as.integer(CnG),as.integer(CnG2),as.numeric(pG),as.numeric(pA), as.numeric(prC), as.integer(condRef),as.numeric(threshT),as.numeric(fst),as.integer(mkvec),as.integer(nkval),as.numeric(lambda),PACKAGE="gammadnamix")[[1]]
  loglik <- Cval + log(pRhoTau(theta[c(nC,nC+1)])) + log(pXi(xi)) #weight with prior of tau and stutter.
  return(-loglik)
 }

 if(is.null(xi)) { #if stutter model
  lower <- rep(0,nC+2)
  upper = c(rep(1,nC-1),rhotaumax,1)
  minL <- Inf
  nOK <- 0 #number of times for reaching optimum
  while(nOK<nDone) {
   tryCatch( {
    foo <- nlm(negloglikYtheta , p=runif(length(lower),lower,upper))
    if(foo$min<minL) {
     minL <- foo$min
     topfoo <- foo
     if(topfoo$code==1 && topfoo$iterations>2) nOK=nOK+1
    }
   },error=function(e) e) #end trycatch 
  } #end loop
  foo2 <- nlm(negloglikYtheta , p=topfoo$est ,hessian=TRUE)
 } else { #if stutter ratio assumed known known
  lower <- rep(0,nC+1)
  upper = c(rep(1,nC-1),rhotaumax)
  minL <- Inf
  nOK <- 0 #number of times for reaching optimum
  while(nOK<nDone) {
   tryCatch( {
    foo <- nlm(negloglikYthetaS , p=runif(length(lower),lower,upper))
    if(foo$min<minL) {
     minL <- foo$min
     topfoo <- foo
     if(topfoo$code==1 && topfoo$iterations>2) nOK=nOK+1
    }
   },error=function(e) e) #end trycatch 
  } #end loop
  foo2 <- nlm(negloglikYthetaS , p=topfoo$est ,hessian=TRUE)
 }
 Sigma <- solve(foo2$hessian)
 dev <- qnorm(0.975)*sqrt(diag(Sigma))
 tab <- cbind( foo2$est - dev ,foo2$est , foo2$est + dev)
 colnames(tab) <- c("2.5%","MLE","97.5%")
 if(!is.null(xi)) rownames(tab) <- c(paste0("omega",1:(nC-1)),"rho","tau")
 if(is.null(xi)) rownames(tab) <- c(paste0("omega",1:(nC-1)),"rho","tau","xi")

 #I) Trying to do Importance samplgin using Normal(MLEhat,Jhat)
 margL <- NULL
 if(!is.null(M)) {
  logdmvnorm <- function(X,mean,cholC) { #function taken from mvtnorm-package
   p <- nrow(cholC)
   tmp <- backsolve(cholC,t(X)-mean,p,transpose=TRUE)
   rss <- colSums(tmp^2)
   logretval <- -sum(log(diag(cholC))) - 0.5*p*log(2*pi) - 0.5*rss
   return(logretval)
  }
  p <- nrow(Sigma) #dimension
  C <- chol(Sigma)
  X <- t( t(C)%*%matrix(rnorm(p*M),ncol=M,nrow=p) + foo2$est )
  logpX <- logdmvnorm(X=X,mean=foo2$est,cholC=C) #get distr-values
  indrm <- numeric()
  for(i in 1:p) indrm <- c(indrm, which( X[,i]<lower[i] | X[,i]>upper[i] ) )
  if(length(indrm)>0) X <- X[-indrm,] 
  #go through each proposal
  M <- nrow(X) #update counter
  bigsum <- 0
  if(is.null(xi)) {
   for(m in 1:M) bigsum <- bigsum + exp( -negloglikYtheta(X[m,]) - logpX[m] ) 
  } else {
   for(m in 1:M) bigsum <- bigsum + exp( -negloglikYthetaS(X[m,]) - logpX[m] ) 
  }
  margL <- bigsum/M #estimated marginal likelihood
 }
 ret <- list(MLE=foo2$est,Sigma=Sigma,CI=tab,loglik=-minL,margL=margL)
 return(ret)
} #end function

