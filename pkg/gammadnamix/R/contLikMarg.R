
#' @title contLikMarg
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description contLikMarg marginalized the likelihood of the STR DNA mixture given some assumed a bayesian model.
#' @details The procedure are doing numerical integration to approximate the marginal probability over noisance parameters. Mixture proportions have flat prior.
#' 
#' Function calls procedure in c++ by using the package Armadillo and Boost.
#' Using "cubature" method (cubature) or "divonne" method (R2Cuba). "cubature" gives more accuracy but require more evaluations than "divonne".
#' "cuhre","suave" and "vegas" not recommended by experience.
#'
#' @param nC Number of contributors in model.
#' @param mixData Evidence object with list elements adata[[i]] and hdata[[i]]. Each element has a loci-list with list-element 'i' storing qualitative data in 'adata' and quantitative data in 'hdata'.
#' @param popFreq A list of allele frequencies for a given population.
#' @param lower Lower bounds of parameters
#' @param upper Upper bounds of parameters
#' @param refData Reference objects with list element [[s]]$adata[[i]]. The list element has reference-list with list-element 's' having a loci-list adata with list-element 'i storing qualitative data.
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. condOrder=-1 means the reference is known-non contributor!
#' @param knownRef Specify known references from refData (index). For instance knownRef=(1,2) means that reference 1 and 2 is known allele samples in the hypothesis. This is affected by fst-correction.
#' @param xi A numeric giving stutter-ratio if it is known. Default is NULL, meaning it is integrated out.
#' @param prC A numeric for allele drop-in probability. Default is 0.
#' @param reltol Required relative tolerance error of evaluations in integration routine. Default is 0.001.
#' @param threshT The detection threshold given. Used when considering probability of allele drop-outs.
#' @param fst is the coancestry coeffecient. Default is 0.
#' @param lambda Parameter in modeled peak height shifted exponential model. Default is 0.
#' @param pXi Prior function for xi-parameter (stutter). Flat prior on [0,1] is default.
#' @param isPhi Boolean whether we are working with reparameterizated parameters
#' @return ret A list(margL,deviation,nEvals) where margL is Marginalized likelihood for hypothesis (model) given observed evidence, deviation is the confidence-interval of margL, nEvals is number of evaluations.
#' @export 
#' @references Hahn,T. (2005). CUBA - a library for multidimensional numerical integration. Computer Physics Communications, 168(2),78-95.
#' @keywords continuous, Bayesian models, Marginalized Likelihood estimation
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
#' #without stutter:
#' hpD <- contLikMarg(nC,mixData,popFreq,refData=refData,condOrder=condhp,xi=0,threshT=threshT)
#' hdD <- contLikMarg(nC,mixData,popFreq,refData=refData,condOrder=condhd,xi=0,threshT=threshT)
#' print(hpD$margL/hdD$margL) #estimated LR=8670.289
#' #with stutter:
#' hpS <- contLikMarg(nC,mixData,popFreq,refData=refData,condOrder=condhp,threshT=threshT,ximax=0.2)
#' hdS <- contLikMarg(nC,mixData,popFreq,refData=refData,condOrder=condhd,threshT=threshT,ximax=0.2)
#' print(hpS$margL/hdS$margL) #estimated LR=8645.057
#' }

contLikMarg = function(nC,mixData,popFreq,lower,upper,refData=NULL,condOrder=NULL,knownRef=NULL,xi=NULL,prC=0,reltol=0.001,threshT=50,fst=0,lambda=0,pXi=function(x)1,isPhi=FALSE){
 require(cubature) 
 if(length(lower)!=length(upper)) stop("Length of integral limits differs")
 
 ret <- prepareC(nC,mixData,popFreq,refData,condOrder,knownRef)

 if(isPhi) {
  if(is.null(xi)) {
    likYtheta <- function(theta) {   #call c++- function: length(theta)=nC+1
     Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nA,ret$allY,ret$allA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(prC), ret$condRef,as.numeric(threshT),as.numeric(fst),ret$mkvec,ret$nkval,as.numeric(lambda),as.integer(1),PACKAGE="gammadnamix")[[1]]
     loglik <- Cval + log(pXi(theta[nC+2])) #weight with prior of tau and 
     return(exp(loglik)) #weight with prior of tau and stutter.
    }
  } else {  
    likYtheta <- function(theta2) {   #call c++- function: length(theta)=nC
     theta <- c(theta2,xi) #stutter-parameter added as known
     Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nA,ret$allY,ret$allA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(prC), ret$condRef,as.numeric(threshT),as.numeric(fst),ret$mkvec,ret$nkval,as.numeric(lambda),as.integer(1),PACKAGE="gammadnamix")[[1]]
     return(exp(Cval))
    }
  }
 } else {
  if(is.null(xi)) {
    likYtheta <- function(theta) {   #call c++- function: length(theta)=nC+1
     Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nA,ret$allY,ret$allA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(prC), ret$condRef,as.numeric(threshT),as.numeric(fst),ret$mkvec,ret$nkval,as.numeric(lambda),as.integer(0),PACKAGE="gammadnamix")[[1]]
     loglik <- Cval + log(pXi(theta[nC+2])) #weight with prior of tau and 
     return(exp(loglik)) #weight with prior of tau and stutter.
    }
  } else {  
    likYtheta <- function(theta2) {   #call c++- function: length(theta)=nC
     if(theta2<lower || theta2>upper) return(-Inf)
     theta <- c(theta2,xi) #stutter-parameter added as known
     Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nA,ret$allY,ret$allA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(prC), ret$condRef,as.numeric(threshT),as.numeric(fst),ret$mkvec,ret$nkval,as.numeric(lambda),as.integer(0),PACKAGE="gammadnamix")[[1]]
     return(exp(Cval))
    }
  }
 }
 np <- nC + 1 + sum(is.null(xi)) #number of unknown parameters
 nU <- nC-ret$nK #number of unknowns
 if(!isPhi) {
  bisectMx <- (nC==2 && nU==2) #if exact 2 unknown contributors in the hypothesis
  if(bisectMx) lower[1] <- 0.5 #restrict outcome of mixture proportions
 }
 foo <- adaptIntegrate(likYtheta, lowerLimit = lower , upperLimit = upper , tol = reltol)
 val <- foo$integral
 dev <- val + c(-1,1)*foo$error
 nEvals <- foo[[3]]

 if(nU>1) { #if more than 1 unknown 
  comb <- factorial(nU)
  val <- comb*val
  dev <- comb*dev
 }

 return(list(margL=val,deviation=dev,nEvals=nEvals))
}

