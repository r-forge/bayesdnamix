#' @title logLiki
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description logLiki calculates the likelihood of each marker of the STR DNA mixture given a theta
#' @details
#' Function calls procedure in c++ by using the package Armadillo and Boost.
#'
#' @param mlefit Fitted object using contLikMLE
#' @return ret A vector with log-likelihood-values for each locus for given model
#' @export

#Get value of likelihood for each marker given theta
logLiki <- function(mlefit){
 theta2 <- theta <- mlefit$fit$thetahat #condition on mle parameter
 model <- mlefit$model #take out assumed model with given data
 np <- length(theta) #number of unknown parameters
 locs <- names(model$popFreq)
 nL <- length(locs)

 nodeg <- is.null(model$kit) #check for degradation
 if(nodeg) theta <- c(theta[1:(model$nC+1)],1) #insert beta variable equal 1
 if(!is.null(model$xi)) {
  theta <- c(theta,as.numeric(model$xi))
 } else {
  if(nodeg)  theta <- c(theta,theta2[np]) #insert fitted xi last again
 }

 if(any(is.na(theta))) {
  logLi <- rep(NA,nL)
 } else {
  logLi <- numeric()
  if(!is.null(model$xi)) {
   loglikYtheta <- function() {   #call c++- function: length(theta)=nC
    Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(model$prC), ret$condRef,as.numeric(model$threshT),as.numeric(model$fst),ret$mkvec,ret$nkval,as.numeric(model$lambda),ret$bp,as.integer(0),PACKAGE="gammadnamix")[[1]]
    return(Cval)
   }
  } else {
   loglikYtheta <- function() {   #call c++- function: length(theta)=nC+1
    Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(model$prC), ret$condRef,as.numeric(model$threshT),as.numeric(model$fst),ret$mkvec,ret$nkval,as.numeric(model$lambda),ret$bp,as.integer(0),PACKAGE="gammadnamix")[[1]]
    return(Cval + log(model$pXi(theta2[np]))) #weight with prior of xi 
   }
  }
  for(loc in locs) { #traverse for each locus
   samples <- lapply(model$samples,function(x) x[loc])
   ret <- prepareC(nC=model$nC,samples,popFreq=model$popFreq[loc],refData=model$refData[loc],condOrder=model$condOrder,knownRef=model$knownRef,kit=model$kit)
   logLi <- c(logLi,loglikYtheta())
  }
 } #end case
 names(logLi) <- locs
 return(logLi)
}



