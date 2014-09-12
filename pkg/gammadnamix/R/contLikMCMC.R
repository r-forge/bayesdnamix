
#' @title contLikMCMC
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description contLikMCMC simulates from the posterior distribution for a bayesian STR DNA mixture model.
#' @details The procedure are doing MCMC to approximate the marginal probability over noisance parameters. Mixture proportions have flat prior.
#' 
#' If no initial values or covariance matrix has been provided to the function, a call to the MLE function is applied.
#' The Metropolis Hastings routine uses a Multivariate Normal distribution with mean 0 and covariance as delta multiplied with the inverse negative hessian with MLE inserted as transistion kernel.
#' Function calls procedure in c++ by using the package Armadillo and Boost.
#' Marginalized likelihood is estimated using Metropolis Hastings with the "Gelfand and Dey" method.
#'
#' @param mlefit Fitted object using contLikMLE
#' @param method Selected MCMC-routine for calculate marginal Likelihood (1=Importance sampling with normal approximation, 2=.
#' @param musigmamax Upper boundary of rho-tau-parameters. Default (10000,1)
#' @param M Number of samples in the MCMC-sampling.
#' @param delta A numerical parameter to scale with the covariance function Sigma. Default is 2. Should be higher to obtain lower acception rate.
#' @param isPhi Boolean whether we are working with reparameterizated parameters
#' @return ret A list (margL,posttheta,postlogL,logpX,accrat) where margL is Marginalized likelihood for hypothesis (model) given observed evidence, posttheta is the posterior samples from a MC routine, postlogL is sampled log-likelihood values, accrat is ratio of accepted samples.
#' @export 
#' @references Craiu,R.V. and Rosenthal, J.S. (2014). Bayesian Computation Via Markov Chain Monte Carlo. Annu. Rev. Stat. Appl., 1,179-201.
#' @keywords continuous, BayesianModels, MCMC, MetropolisHastings, MarginalizedLikelihoodEstimation

contLikMCMC = function(mlefit,musigmamax=c(10000,1),M=1e4,delta=2,isPhi=FALSE) {
 #A mlefit object returned from contLikMLE is required to do MCMC!
 #Phi should be faster than theta since it doesn't check limits

 loglik0 <-  mlefit$fit$loglik #get maximized likelihood
 model <- mlefit$model
 if(isPhi) {
  theta0 <- mlefit$fit$phihat
  Sigma0 <- mlefit$fit$phiSigma
 } else {
  theta0 <- mlefit$fit$thetahat
  Sigma0 <- mlefit$fit$thetaSigma
 }
 varnames <- names(theta0) #variable names
 if(!all(length(theta0)%in%dim(Sigma0))) stop("Length of theta0 and dimension of Sigma was not the same!")
 ret <- prepareC(model$nC,model$samples,model$popFreq,model$refData,model$condOrder,model$knownRef)

 if(isPhi) {
  if(is.null(model$xi)) {
    loglikYtheta <- function(theta) {   #call c++- function: length(theta)=nC+1
     Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(model$prC), ret$condRef,as.numeric(model$threshT),as.numeric(model$fst),ret$mkvec,ret$nkval,as.numeric(model$lambda),as.integer(1),PACKAGE="gammadnamix")[[1]]
     loglik <- Cval + log(model$pXi(theta[nC+2])) #weight with prior of tau and 
     return(loglik) #weight with prior of tau and stutter.
    }
  } else {  
     loglikYtheta <- function(theta2) {   #call c++- function: length(theta)=nC
     theta <- c(theta2,xi) #stutter-parameter added as known
     Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(model$prC), ret$condRef,as.numeric(model$threshT),as.numeric(model$fst),ret$mkvec,ret$nkval,as.numeric(model$lambda),as.integer(1),PACKAGE="gammadnamix")[[1]]
     return(Cval)
    }
  }
 } else {
  if(is.null(model$xi)) {
    loglikYtheta <- function(theta) {   #call c++- function: length(theta)=nC+1
     if(any(theta<lower) || any(theta>upper) ) return(-Inf)
     Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(model$prC), ret$condRef,as.numeric(model$threshT),as.numeric(model$fst),ret$mkvec,ret$nkval,as.numeric(model$lambda),as.integer(0),PACKAGE="gammadnamix")[[1]]
     loglik <- Cval + log(model$pXi(theta[nC+2])) #weight with prior of tau and 
     return(loglik) #weight with prior of tau and stutter.
    }
  } else {  
    loglikYtheta <- function(theta2) {   #call c++- function: length(theta)=nC
     if(any(theta2<lower) || any(theta2>upper)) return(-Inf)
     theta <- c(theta2,model$xi) #stutter-parameter added as known
     Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(model$prC), ret$condRef,as.numeric(model$threshT),as.numeric(model$fst),ret$mkvec,ret$nkval,as.numeric(model$lambda),as.integer(0),PACKAGE="gammadnamix")[[1]]
     return(Cval)
    }
  }
 }
 np <- nC + 1 + sum(is.null(model$xi)) #number of unknown parameters
 C <- chol(delta*Sigma0) #scale variance with a factor 2: ensures broad posterior
 X <- t( t(C)%*%matrix(rnorm(np*M),ncol=M,nrow=np)) #proposal values

 if(!isPhi) {
  lower <- rep(0,nC+1)
  upper <- c(rep(1,nC-1),musigmamax)
  if(is.null(model$xi)) {
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
 #removed:  Importance sampling using Normal(theta0,delta*Sigma)
 if(1) { #MCMC by Gelfand and Dey (1994), using h() = Normal(theta0,delta*Sigma)
   rlist <- list()
   if(nC>1) rlist[[length(rlist)+1]] <- 1:(nC-1)
   rlist[[length(rlist)+1]] <- nC:(nC+1)
   if(is.null(model$xi))  rlist[[length(rlist)+1]] <- np
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
 colnames(posttheta) <- varnames  #save variable names
 return(list(margL=margL,posttheta=posttheta,postlogL=postlogL,logpX=logpX,accrat=accrat,MLE=theta0,Sigma=Sigma0))
} #end function

