
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
#' @param musigmamax Maximum range of mu and sigma-parameter. Default is c(10000,1000).
#' @param nDone Maxumum number of random evaluations nlm-optimizing routing. Default is 1.
#' @param threshT The detection threshold given. Used when considering probability of allele drop-outs.
#' @param fst is the coancestry coeffecient. Default is 0.
#' @param lambda Parameter in modeled peak height shifted exponential model. Default is 0.
#' @param pXi Prior function for xi-parameter (stutter). Flat prior on [0,1] is default.
#' @param ximax Maximum range of xi-parameter. Default is 1.
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

contLikMLE = function(nC,mixData,popFreq,refData=NULL,condOrder=NULL,knownRef=NULL,xi=NULL,prC=0,musigmamax =c(10000,1000),nDone=1,threshT=50,fst=0,lambda=0,pXi=function(x)1,ximax=1){
 require(Rsolnp)
 require(numDeriv)
 ret <- prepareC(nC,mixData,popFreq,refData,condOrder,knownRef)
 if(is.null(condOrder)) condOrder <- rep(0,nC) #insert condorder if missing
 unRange <- (1:nC)[!(1:nC)%in%condOrder] #get the range of the unknowns
 nU <- length(unRange)

#Optimize MLE:
 #Two cases: Integrate over Stutter or Stutter known
 if(is.null(xi)) {
  negloglikYtheta <- function(theta) {   #call c++- function: length(theta)=nC+1
   Cval  <- .C("loglikgammaC",ret$logPE,as.numeric(theta),ret$nC,ret$nL,ret$nA,ret$allY,ret$allA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(prC), ret$condRef,as.numeric(threshT),as.numeric(fst),ret$mkvec,ret$nkval,as.numeric(lambda),PACKAGE="gammadnamix")[[1]]
   loglik <- Cval + log(pXi(theta[nC+2])) #weight with prior of tau and 
   return(-loglik) #weight with prior of tau and stutter.
  }
 } else {  
  negloglikYtheta <- function(theta2) {   #call c++- function: length(theta)=nC
   theta <- c(theta2,xi) #stutter-parameter added as known
   Cval  <- .C("loglikgammaC",ret$logPE,as.numeric(theta),ret$nC,ret$nL,ret$nA,ret$allY,ret$allA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(prC), ret$condRef,as.numeric(threshT),as.numeric(fst),ret$mkvec,ret$nkval,as.numeric(lambda),PACKAGE="gammadnamix")[[1]]
   loglik <- Cval + log(pXi(xi)) #weight with prior of tau and stutter.
   return(-loglik)
  }
 }

 lower <- lower0 <- rep(0,nC+1)
 upper <- upper0 <- c(rep(1,nC-1),musigmamax)
 if(is.null(xi)) { #if stutter model
  lower <- c(lower0,0)
  upper = c(upper0,ximax)
 }
 np <- length(lower) #number of unknown parameters
 mixsum <- function(x) sum(x[1:(nC-1)]) #avoid that the sum exeeds 1
 maxL <- -Inf #
 nOK <- 0 #number of times for reaching optimum
 while(nOK<nDone) {
   p0 <- runif(length(lower0),lower0,upper0)
   if(is.null(xi)) {
     xi0 <- rexp(1,20)
     if(xi0>ximax) next
     p0  <- c(p0,xi0) #decay start-value of stutter
   }
   if(is.infinite(negloglikYtheta(p0))) next #skip to next if still INF
   tryCatch( {
      foo <- solnp(pars=p0, fun=negloglikYtheta, ineqfun = mixsum ,ineqLB=0, ineqUB=1,LB = lower, UB = upper,control=list(trace=0))
      if(foo$convergence==0 && foo$nfuneval>2) {
       likval <- -negloglikYtheta(foo$pars) #get maximized log-likelihood
       if(is.infinite(likval)) next
       hess <- hessian(negloglikYtheta,foo$pars) #get the hessian
       if(any(is.infinite(hess) | is.nan(hess))) next
       if(any(diag(solve(hess))<=0)) next; #not a local optimum
       nOK=nOK+1 #it was accepted as an optimum
       if(likval>maxL) {
        maxL <- likval #maximized likelihood
        topfoo <- foo #set as topfoo
       }
      }
   },error=function(e) e) #end trycatch 
 } #end while loop
 mle <- topfoo$pars #get MLE
 mx <- c(mle[1:(nC-1)],1-sum(mle[1:(nC-1)])) #get all mix-props
 if(nU>1) {
  mx[unRange] <- sort(mx[unRange],decreasing=TRUE) #sort unknowns in decreasing mix prop
  mle <- c(mx[-nC],mle[nC:np]) #last index is removed. This could again be a known contributor
 }
 Sigma <- solve(hessian(negloglikYtheta,mle))
 mle2 <- c(mx,mle[nC:np])
 Sigma2 <- matrix(NA,nrow=np+1,ncol=np+1) #extended covariance matrix also including mx[nC]
 Sigma2[1:(nC-1),1:(nC-1)] <- Sigma[1:(nC-1),1:(nC-1)] 
 Sigma2[nC:np+1,nC:np+1] <- Sigma[nC:np,nC:np] 
 Sigma2[nC:np+1,nC:np+1] <- Sigma[nC:np,nC:np] 
 Sigma2[nC:np+1,1:(nC-1)] <- Sigma[nC:np,1:(nC-1)] 
 Sigma2[1:(nC-1),nC:np+1] <- Sigma[1:(nC-1),nC:np] 
 Sigma2[nC,nC] <- sum(Sigma[1:(nC-1),1:(nC-1)])
 for(k in (1:(np+1))[-nC]) {
  Sigma2[nC,k] <- Sigma2[k,nC] <- -sum(Sigma[1:(nC-1),k-sum(k>nC)]) 
 }
 dev <- qnorm(0.975)*sqrt(diag(Sigma2))  #constrained part is first
 tab <- cbind(mle2 - dev ,mle2, mle2 + dev)
 colnames(tab) <- c("2.5%","MLE","97.5%")
 if(!is.null(xi)) rownames(tab) <- c(paste0("mx",1:nC),"mu","sigma")
 if(is.null(xi)) rownames(tab) <- c(paste0("mx",1:nC),"mu","sigma","xi")
 ret <- list(MLE=mle,Sigma=Sigma,Sigma2=Sigma2,CI=tab,loglik=maxL)
 return(ret)
} #end function

