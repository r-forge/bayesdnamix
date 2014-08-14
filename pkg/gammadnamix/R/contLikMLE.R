
#' @title contLikMLE
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description contLikMLE optimizes the likelihood of the STR DNA mixture given some assumed a bayesian model.
#' @details The procedure are doing numerical optimization to approximate the marginal probabilit over noisance parameters. Mixture proportions have flat prior.
#' 
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
#' @param nDone Maxumum number of random evaluations nlm-optimizing routing. Default is 1.
#' @param threshT The detection threshold given. Used when considering probability of allele drop-outs.
#' @param fst is the coancestry coeffecient. Default is 0.
#' @param lambda Parameter in modeled peak height shifted exponential model. Default is 0.
#' @param pXi Prior function for xi-parameter (stutter). Flat prior on [0,1] is default.
#' @param delta Standard deviation of normal distribution when drawing random startpoints. Default is 10.
#' @return ret A list(MLE,Sigma,Sigma2,CI,loglik) with Maximixed likelihood elements for hypothesis (model) given observed evidence.
#' @export
#' @references Cowell,R.G. et.al. (2014). Analysis of forensic DNA mixtures with artefacts. Applied Statistics, 64(1),1-32.
#' @keywords continuous model, Maximum Likelihood Estimation
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
#' nDone <- 3 #require 3 optimizations
#' hpD <- contLikMLE(nC,mixData,popFreq,refData=refData,condOrder=condhp,xi=0,threshT=threshT,nDone=nDone )
#' hdD <- contLikMLE(nC,mixData,popFreq,refData=refData,condOrder=condhd,xi=0,threshT=threshT,nDone=nDone )
#' print(exp(hpD$loglik - hdD$loglik)) #estimated LR=21554.24
#' #with stutter, but stuttarratio at most 0.2:
#' hpS <- contLikMLE(nC,mixData,popFreq,refData=refData,condOrder=condhp,threshT=threshT,nDone=nDone,ximax=0.2)
#' hdS <- contLikMLE(nC,mixData,popFreq,refData=refData,condOrder=condhd,threshT=threshT,nDone=nDone,ximax=0.2)
#' print(exp(hpS$loglik - hdS$loglik)) #estimated LR=21554.4
#' }

contLikMLE = function(nC,mixData,popFreq,refData=NULL,condOrder=NULL,knownRef=NULL,xi=NULL,prC=0,nDone=1,threshT=50,fst=0,lambda=0,pXi=function(x)1,delta=10){
 ret <- prepareC(nC,mixData,popFreq,refData,condOrder,knownRef)

 if(is.null(xi)) {
  negloglikYtheta <- function(theta) {   #call c++- function: length(theta)=nC+1
   Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nA,ret$allY,ret$allA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(prC), ret$condRef,as.numeric(threshT),as.numeric(fst),ret$mkvec,ret$nkval,as.numeric(lambda),as.integer(1),PACKAGE="gammadnamix")[[1]]
   loglik <- Cval + log(pXi(1/(1+exp(-theta[nC+2])))) #weight with prior of tau and 
   return(-loglik) #weight with prior of tau and stutter.
  }
 } else {  
  negloglikYtheta <- function(theta2) {   #call c++- function: length(theta)=nC
   theta <- c(theta2,xi) #stutter-parameter added as known
   Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nA,ret$allY,ret$allA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(prC), ret$condRef,as.numeric(threshT),as.numeric(fst),ret$mkvec,ret$nkval,as.numeric(lambda),as.integer(1),PACKAGE="gammadnamix")[[1]]
   return(-Cval)
  }
 }
 np <- nC + 1 + sum(is.null(xi)) #number of unknown parameters
 maxL <- -Inf #
 nOK <- 0 #number of times for reaching optimum
 while(nOK<nDone) {
   p0 <- rnorm(np,sd=delta)
   if(is.infinite(negloglikYtheta(p0))) next #skip to next if still INF
   tryCatch( {
      foo <- nlm(f=negloglikYtheta, p=p0,hessian=TRUE)
      Sigma <- solve(foo$hessian)
      if(any(diag(Sigma)<=0)) next; #not a local optimum
      if(foo$code==1 && foo$iterations>2) {
       likval <- -foo$min
       nOK=nOK+1 #it was accepted as an optimum
       if(likval>maxL) {
        maxL <- likval #maximized likelihood
        maxTheta <- foo$est #set as topfoo     
        maxSigma <- Sigma 
       }
      }
    },error=function(e) e) #end trycatch 
  } #end while loop
 #transfer back: 
 mx <- numeric()
 if(nC>1) {
  mx <- 1/(1+exp(-maxTheta[1:(nC-1)]))
  if(nC>2) { #need to transfer back
   for(i in 2:(nC-1)) {
    mx[i] <- mx[i]*(1-sum(mx[1:(i-1)]))
   }
  }
 }
 musigma <- exp(maxTheta[nC:(nC+1)]) #inverse-log
 mle <- c(mx,musigma) #last index is removed. This could again be a known contributor
 mle2 <- c(mx,1-sum(mx),musigma) #last index is removed. This could again be a known contributor
 if(is.null(xi)) {
  tmp <- 1/(1+exp(-maxTheta[nC+2]))
  mle <- c(mle,tmp) #add xi to parameters
  mle2 <- c(mle2,tmp) #add xi to parameters
 }

 alpha <- 0.05
 #CI:
 SD0 <- sqrt(diag(maxSigma))
 CL <- maxTheta + qnorm(alpha/2)*SD0
 CU <- maxTheta + qnorm(1-alpha/2)*SD0 
 CI0 <- maxTheta + cbind(CL,maxTheta ,CU)

 #Delta-method to Sigma matrix
 Jacob <- function(phi,mle) { #Jabobian matrix
  J <- matrix(0,length(phi),length(phi))
  if(nC>1) {
   DmDm <- matrix(0,nC-1,nC-1)
   for(i in 1:(nC-1)) {
    mitmp <- 1/(1+exp(-phi[i])) #temporary variable
    for(j in 1:i) {
     if(j==i) {
       DmDm[i,i] <- exp(-phi[i])*mitmp^2
       if(i>1) DmDm[i,i] <- DmDm[i,i]*(1-sum(mle[1:(j-1)])) #note using mle(theta) here!
     } else {
       DmDm[i,j] <- -mitmp*(sum( DmDm[1:i,j] ))
     } #end case
    } #end for each col j (xj)
   } #end for each row i (fi)
  J[1:(nC-1),1:(nC-1)] <- DmDm
  }
  for(i in nC:(nC+1)) J[i,i] <- exp(phi[i])
  if(is.null(xi)) {
   tmp <- exp(-phi[nC+2])
   J[nC+2,nC+2] <- tmp*(1+tmp)^(-2)
  }
  return(J)
 } #end jacobian
 J <- Jacob(maxTheta,mle)
 Sigma <- (t(J)%*%maxSigma%*%J) #this is correct covariance of mle. Observed hessian is used
 #Sigma <- (J%*%maxSigma%*%t(J))/sqrt(np) #this formula is wrong!

 #get extended Sigma (all parameters)
 Sigma2 <- matrix(NA,nrow=np+1,ncol=np+1) #extended covariance matrix also including mx[nC]
 Sigma2[nC:np+1,nC:np+1] <- Sigma[nC:np,nC:np] 
 Sigma2[nC:np+1,nC:np+1] <- Sigma[nC:np,nC:np] 
 if(nC>1) {
  Sigma2[nC:np+1,1:(nC-1)] <- Sigma[nC:np,1:(nC-1)] 
  Sigma2[1:(nC-1),1:(nC-1)] <- Sigma[1:(nC-1),1:(nC-1)] 
  Sigma2[1:(nC-1),nC:np+1] <- Sigma[1:(nC-1),nC:np] 
  Sigma2[nC,nC] <- sum(Sigma[1:(nC-1),1:(nC-1)])
  for(k in (1:(np+1))[-nC]) {
   Sigma2[nC,k] <- Sigma2[k,nC] <- -sum(Sigma[1:(nC-1),k-sum(k>nC)]) 
  }
 } else {
  Sigma2[1,1:(np+1)] <- Sigma2[1:(np+1),1] <- 0 #no uncertainty
 }
 SD <- sqrt(diag(Sigma2))
 CL <- mle2 + qnorm(alpha/2)*SD 
 CU <- mle2 + qnorm(1-alpha/2)*SD
 CI <- cbind(CL,mle2,CU)
 ret <- list(theta0=maxTheta,Sigma0=maxSigma,loglik0=maxL,CI0=CI0,theta1=mle,Sigma1=Sigma,MLE=mle2,SIGMA=Sigma2,CI=CI)
 return(ret)
} #end function

