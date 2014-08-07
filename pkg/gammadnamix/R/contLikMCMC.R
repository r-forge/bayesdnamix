
#' @title contLikMCMC
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description contLikMCMC simulates from the posterior distribution for a bayesian STR DNA mixture model.
#' @details The procedure are doing MCMC to approximate the marginal probability over noisance parameters. Mixture proportions have flat prior.
#' 
#' If no initial values or covariance matrix has been provided to the function, a call to the MLE function is applied.
#' The Metropolis Hastings routine uses a Multivariate Normal distribution with mean 0 and covariance as delta multiplied with the inverse negative hessian with MLE inserted as transistion kernel.
#' Function calls procedure in c++ by using the package Armadillo and Boost.
#'
#' @param nC Number of contributors in model.
#' @param mixData Evidence object with list elements adata[[i]] and hdata[[i]]. Each element has a loci-list with list-element 'i' storing qualitative data in 'adata' and quantitative data in 'hdata'.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects with list element [[s]]$adata[[i]]. The list element has reference-list with list-element 's' having a loci-list adata with list-element 'i storing qualitative data.
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. condOrder=-1 means the reference is known-non contributor!
#' @param knownRef Specify known references from refData (index). For instance knownRef=(1,2) means that reference 1 and 2 is known allele samples in the hypothesis. This is affected by fst-correction.
#' @param xi A numeric giving stutter-ratio if it is known. Default is NULL, meaning it is integrated out.
#' @param prC A numeric for allele drop-in probability. Default is 0.
#' @param method Selected MCMC-routine for calculate marginal Likelihood (1=Importance sampling with normal approximation, 2=Metropolis Hastings with "Gelfand and Dey" method).
#' @param nDone Required number of attained optimization for random start points. Default is 1.
#' @param threshT The detection threshold given. Used when considering probability of allele drop-outs.
#' @param fst is the coancestry coeffecient. Default is 0.
#' @param lambda Parameter in modeled peak height shifted exponential model. Default is 0.
#' @param musigmamax Upper boundary of rho-tau-parameters. Default (10000,1)
#' @param pXi Prior function for xi-parameter (stutter). Flat prior on [0,1] is default.
#' @param M Number of samples in the MCMC-sampling.
#' @param mlefit Fitted object using contLikMLE
#' @param delta A numerical parameter to scale with the covariance function Sigma. Default is 2. Should be higher to obtain lower acception rate.
#' @param isPhi Boolean whether we are working with reparameterizated parameters
#' @return ret A list (margL,posttheta,postlogL,logpX,accrat) where margL is Marginalized likelihood for hypothesis (model) given observed evidence, posttheta is the posterior samples from a MC routine, postlogL is sampled log-likelihood values, accrat is ratio of accepted samples.
#' @export 
#' @references Craiu,R.V. and Rosenthal, J.S. (2014). Bayesian Computation Via Markov Chain Monte Carlo. Annu. Rev. Stat. Appl., 1,179-201.
#' @keywords continuous, Bayesian models, MCMC, Metropolis Hastings, Marginalized Likelihood estimation
#' @examples
#' \dontrun{
#' threshT <- 50 #threshold
#' load(system.file("mcData.Rdata", package = "gammadnamix"))
#' mixData <- data$samples$MC15
#' retlist <- Qassignate(mixData,popFreq=data$popFreq,refData=data$refData) #Q-assignation
#' popFreq <- retlist$popFreq #get updated popFreq
#' refData <- retlist$refData #get updated popFreq
#' #Hypothesis:
#' nC <- 2 #number of contributors
#' condhp <- condhd <- rep(0,length(refData)) 
#' condhp [1] <- 1 #condition on first one in refdata
#' M <- 1e2 #number of samples from posterior
#' delta = 10  #Step parameter for covariance of the proposal used in the Metropolis Hastings routine
#' #without stutter:
#' hpD <- contLikMCMC(nC,mixData,popFreq,refData=refData,condOrder=condhp,xi=0,threshT=threshT,M=M,delta=delta)
#' validMCMC(hpD$posttheta)
#' hdD <- contLikMCMC(nC,mixData,popFreq,refData=refData,condOrder=condhd,xi=0,threshT=threshT,M=M,delta=delta)
#' validMCMC(hdD$posttheta)
#' print(hpD$margL/hdD$margL) #estimated LR
#' #with stutter:
#' hpS <- contLikMCMC(nC,mixData,popFreq,refData=refData,condOrder=condhp,threshT=threshT,M=M,delta=delta,ximax=0.2)
#' validMCMC(hpS$posttheta)
#' hdS <- contLikMCMC(nC,mixData,popFreq,refData=refData,condOrder=condhd,threshT=threshT,M=M,delta=delta,ximax=0.2)
#' validMCMC(hdS$posttheta)
#' print(hpS$margL/hdS$margL) #estimated LR
#' }

contLikMCMC = function(nC,mixData,popFreq,refData=NULL,condOrder=NULL,knownRef=NULL,xi=NULL,prC=0,method=2,nDone=1,threshT=50,fst=0,lambda=0,musigmamax=c(10000,1),pXi=function(x)1,M=1e4,mlefit=NULL,delta=2,isPhi=FALSE){
 #Phi should be faster than theta since it doesn't check limits
 #Optimize with MLE if theta0,Sigma not given
 if(is.null(mlefit)) { #if any is missing:
  mlefit <- contLikMLE(nC,mixData,popFreq,refData=refData,condOrder=condOrder,knownRef=knownRef,xi=xi,prC=prC,threshT=threshT,fst=fst,lambda=lambda,pXi=pXi,nDone=3)
 }
 loglik0 <- mlefit$loglik0
 if(isPhi) {
  theta0 <- mlefit$theta0
  Sigma0 <- mlefit$Sigma0
 } else {
  theta0 <- mlefit$theta1
  Sigma0 <- mlefit$Sigma1
 }
 if(!all(length(theta0)%in%dim(Sigma0))) stop("Length of theta0 and dimension of Sigma was not the same!")
 ret <- prepareC(nC,mixData,popFreq,refData,condOrder,knownRef)

 if(isPhi) {
  if(is.null(xi)) {
    loglikYtheta <- function(theta) {   #call c++- function: length(theta)=nC+1
     Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nA,ret$allY,ret$allA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(prC), ret$condRef,as.numeric(threshT),as.numeric(fst),ret$mkvec,ret$nkval,as.numeric(lambda),as.integer(1),PACKAGE="gammadnamix")[[1]]
     loglik <- Cval + log(pXi(1/(1+exp(-theta[nC+2])))) #weight with prior of tau and 
     return(loglik) #weight with prior of tau and stutter.
    }
  } else {  
    loglikYtheta <- function(theta2) {   #call c++- function: length(theta)=nC
     theta <- c(theta2,xi) #stutter-parameter added as known
     Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nA,ret$allY,ret$allA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(prC), ret$condRef,as.numeric(threshT),as.numeric(fst),ret$mkvec,ret$nkval,as.numeric(lambda),as.integer(1),PACKAGE="gammadnamix")[[1]]
     return(Cval)
    }
  }
 } else {
  if(is.null(xi)) {
    loglikYtheta <- function(theta) {   #call c++- function: length(theta)=nC+1
     if(any(theta<lower | theta>upper)) return(-Inf)
     Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nA,ret$allY,ret$allA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(prC), ret$condRef,as.numeric(threshT),as.numeric(fst),ret$mkvec,ret$nkval,as.numeric(lambda),as.integer(0),PACKAGE="gammadnamix")[[1]]
     loglik <- Cval + log(pXi(theta[nC+2])) #weight with prior of tau and 
     return(loglik) #weight with prior of tau and stutter.
    }
  } else {  
    loglikYtheta <- function(theta2) {   #call c++- function: length(theta)=nC
     if(any(theta2<lower | theta2>upper)) return(-Inf)
     theta <- c(theta2,xi) #stutter-parameter added as known
     Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nA,ret$allY,ret$allA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(prC), ret$condRef,as.numeric(threshT),as.numeric(fst),ret$mkvec,ret$nkval,as.numeric(lambda),as.integer(0),PACKAGE="gammadnamix")[[1]]
     return(Cval)
    }
  }
 }
 np <- nC + 1 + sum(is.null(xi)) #number of unknown parameters
 C <- chol(delta*Sigma0) #scale variance with a factor 2: ensures broad posterior
 X <- t( t(C)%*%matrix(rnorm(np*M),ncol=M,nrow=np)) #proposal values

 if(!isPhi) {
  lower <- rep(0,nC+1)
  upper <- c(rep(1,nC-1),musigmamax)
  if(is.null(xi)) {
   lower <- c(lower,0)
   upper <- c(upper,1)
  }
 }

 logdmvnorm <- function(X,mean,cholC) { #function taken from mvtnorm-package
   p <- nrow(cholC)
   tmp <- backsolve(cholC,t(X)-mean,p,transpose=TRUE)
   rss <- colSums(tmp^2)
   logretval <- -sum(log(diag(cholC))) - 0.5*p*log(2*pi) - 0.5*rss
   return(logretval)
 }
 #v1: Importance sampling using Normal(theta0,delta*Sigma)
 if(method==1) {
  print("Removed")
  return(NULL)

 #v2: MCMC by Gelfand and Dey (1994), using h() = Normal(theta0,delta*Sigma)
 } else if(method==2) { #Simulate variable at-the time: Two blocks
   rlist <- list()
   rlist[[1]] <- 1:(nC-1)
   rlist[[2]] <- nC:(nC+1)
   if(is.null(xi))  rlist[[3]] <- np
   nB <- length(rlist) #number of blocks
   M2 <- nB*M+1
   posttheta <- matrix(NA,ncol=np,nrow=M2) #accepted theta
   postlogL <- rep(NA,M2) #accepted theta
   posttheta[1,] <- theta0
   postlogL[1] <- loglik0 #loglikYtheta(theta0) #get start-likelihood   
   U <- runif(M2) #random numbers
   m <- 2 #counter for samples
   m2 <- 1 #counter for proposal
   nacc <- 0  
   while(m<=M2) {
    for(r in 1:nB ) { #for each blocks
     range <- rlist[[r]]
     posttheta[m,] <-  posttheta[m-1,] #proposed theta
     posttheta[m,range] <- X[m2,range] + posttheta[m,range]
     postlogL[m] <- loglikYtheta(posttheta[m,])
     pr <- exp(postlogL[m]- postlogL[m-1]) #acceptance rate
     if(U[m]>pr) { #if not accepted, i.e. random pr too large
      posttheta[m,] <-  posttheta[m-1,]
      postlogL[m] <- postlogL[m-1 ]
     } else {
      nacc <- nacc + 1
     }
#     print(posttheta[m,])
#     print(1-sum(posttheta[m,1:(nC-1)]))
     m <- m + 1 #update counter
    } #end for each blocks
    m2 <- m2 +1 #update proposal counter
   } #end while not done
  accrat <- nacc/M2 #acceptance ratio
  logpX <- logdmvnorm(posttheta,mean=theta0,cholC=chol(Sigma0)) #insert with Normal-approx of post-theta
  #plot(postlogL,ty="l")
  #plot(logpX,ty="l")
  margL <- 1/mean(exp(logpX - postlogL)) #estimated marginal likelihood
 }
 nU <- nC-ret$nK #number of unknowns
 if(nU>1) { #if more than 1 unknown 
  margL <- factorial(nU)*margL #get correct ML adjusting for symmetry
 }

 return(list(margL=margL,posttheta=posttheta,postlogL=postlogL,logpX=logpX,accrat=accrat))
} #end function

