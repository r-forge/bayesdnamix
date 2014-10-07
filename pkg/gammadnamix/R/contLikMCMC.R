
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
#' @param niter Number of samples in the MCMC-sampling.
#' @param delta A numerical parameter to scale with the covariance function Sigma. Default is 2. Should be higher to obtain lower acception rate.
#' @return ret A list (margL,posttheta,postlogL,logpX,accrat) where margL is Marginalized likelihood for hypothesis (model) given observed evidence, posttheta is the posterior samples from a MC routine, postlogL is sampled log-likelihood values, accrat is ratio of accepted samples.
#' @export 
#' @references Craiu,R.V. and Rosenthal, J.S. (2014). Bayesian Computation Via Markov Chain Monte Carlo. Annu. Rev. Stat. Appl., 1,179-201.
#' @keywords continuous, BayesianModels, MCMC, MetropolisHastings, MarginalizedLikelihoodEstimation

contLikMCMC = function(mlefit,niter=1e4,delta=2) {
 #A mlefit object returned from contLikMLE is required to do MCMC!
 loglik0 <-  mlefit$fit$loglik #get maximized likelihood
 model <- mlefit$model
 phi0 <- mlefit$fit$phihat
 Sigma0 <- mlefit$fit$phiSigma
 varnames <- names(mlefit$fit$thetahat) #variable names
 if(!all(length(phi0)%in%dim(Sigma0))) stop("Length of phi0 and dimension of Sigma was not the same!")
 ret <- prepareC(model$nC,model$samples,model$popFreq,model$refData,model$condOrder,model$knownRef)
 nC <- ret$nC

 if(is.null(model$xi)) {
   loglikYphi <- function(phi) {   #call c++- function: length(phi)=nC+1
    Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(phi),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(model$prC), ret$condRef,as.numeric(model$threshT),as.numeric(model$fst),ret$mkvec,ret$nkval,as.numeric(model$lambda),as.integer(1),PACKAGE="gammadnamix")[[1]]
    loglik <- Cval + log(model$pXi(phi[ret$nC+2])) #weight with prior of tau and 
    return(loglik) #weight with prior of tau and stutter.
   }
 } else {  
   loglikYphi <- function(phi2) {   #call c++- function: length(phi)=nC
    phi <- c(phi2,xi) #stutter-parameter added as known
    Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(phi),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(model$prC), ret$condRef,as.numeric(model$threshT),as.numeric(model$fst),ret$mkvec,ret$nkval,as.numeric(model$lambda),as.integer(1),PACKAGE="gammadnamix")[[1]]
    return(Cval)
   }
 }
 
 np <- nC + 1 + sum(is.null(model$xi)) #number of unknown parameters
 C <- chol(delta*Sigma0) #scale variance with a factor 2: ensures broad posterior
 X <- t( t(C)%*%matrix(rnorm(np*niter),ncol=niter,nrow=np)) #proposal values

 logdmvnorm <- function(X,mean,cholC) { #function taken from mvtnorm-package
   p <- nrow(cholC)
   tmp <- backsolve(cholC,t(X)-mean,p,transpose=TRUE)
   rss <- colSums(tmp^2)
   logretval <- -sum(log(diag(cholC))) - 0.5*p*log(2*pi) - 0.5*rss
   return(logretval)
 }
 #removed:  Importance sampling using Normal(phi0,delta*Sigma)
 if(1) { #MCMC by Gelfand and Dey (1994), using h() = Normal(phi0,delta*Sigma)
   rlist <- list()
   if(nC>1) rlist[[length(rlist)+1]] <- 1:(nC-1)
   rlist[[length(rlist)+1]] <- nC:(nC+1)
   if(is.null(model$xi))  rlist[[length(rlist)+1]] <- np
   nB <- length(rlist) #number of blocks
   M2 <- nB*niter+1
   postphi <- matrix(NA,ncol=np,nrow=M2) #accepted phi
   postlogL <- rep(NA,M2) #accepted phi
   postphi[1,] <- phi0
   postlogL[1] <- loglik0 #loglikYphi(phi0) #get start-likelihood   
   U <- runif(M2) #random numbers
   m <- 2 #counter for samples
   m2 <- 1 #counter for proposal
   nacc <- 0  
   while(m<=M2) {
    for(r in 1:nB ) { #for each blocks
     range <- rlist[[r]]
     postphi[m,] <-  postphi[m-1,] #proposed phi
     postphi[m,range] <- X[m2,range] + postphi[m,range]
     postlogL[m] <- loglikYphi(postphi[m,])
     pr <- exp(postlogL[m]- postlogL[m-1]) #acceptance rate
     if(U[m]>pr) { #if not accepted, i.e. random pr too large
      postphi[m,] <-  postphi[m-1,]
      postlogL[m] <- postlogL[m-1 ]
     } else {
      nacc <- nacc + 1
     }
#     print(postphi[m,])
#     print(1-sum(postphi[m,1:(nC-1)]))
     m <- m + 1 #update counter
    } #end for each blocks
    m2 <- m2 +1 #update proposal counter
   } #end while not done
  accrat <- nacc/M2 #acceptance ratio
  logpX <- logdmvnorm(postphi,mean=phi0,cholC=chol(Sigma0)) #insert with Normal-approx of post-phi
  #plot(postlogL,ty="l")
  #plot(logpX,ty="l")
  margL <- 1/mean(exp(logpX - postlogL)) #estimated marginal likelihood
 }
 nU <- nC-ret$nK #number of unknowns
 if(nU>1) { #if more than 1 unknown 
  margL <- factorial(nU)*margL #get correct ML adjusting for symmetry
 } #end method
 colnames(postphi) <- varnames  #save variable names

 #transfer back to theta: 
 if(nC>1) {
  postphi[,1:(nC-1)] <- 1/(1+exp(-postphi[,1:(nC-1)]))
  if(nC>2) { #need to transfer back
   for(i in 2:(nC-1)) {
    if(i==2) postphi[,i] <- postphi[,i]*(1-postphi[,1:(i-1)])
    if(i>2) postphi[,i] <- postphi[,i]*( 1 - rowSums(postphi[,1:(i-1)]) )
   }
  }
 }
 postphi[,nC:(nC+1)] <- exp(postphi[,nC:(nC+1)])  #inverse-log
 if(is.null(model$xi))  postphi[,nC+2] <- 1/(1+exp(-postphi[,nC+2])) 
 return(list(margL=margL,posttheta=postphi,postlogL=postlogL,logpX=logpX,accrat=accrat,MLE=mlefit$fit$thetahat,Sigma=mlefit$fit$thetaSigma))
} #end function

