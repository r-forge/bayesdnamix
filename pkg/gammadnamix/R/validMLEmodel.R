#' @title validMLEmodel
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description validMLEmodel makes model validation whether the observed peak heights fits the maximum likelihood fitted gamma distribution.
#' @details The cumulative probability of the observed allele peaks are calculated and compared with a uniform distribution.
#' Function calls procedure in c++ by using the package Armadillo and Boost.
#'
#' @param mlefit Fitted object using contLikMLE
#' @return ret A vector for each marker with cumulative probabilities
#' @export
validMLEmodel <- function(mlefit) {
 theta2 <- theta <- mlefit$fit$thetahat #condition on mle parameter
 np <- length(theta) #number of unknown parameters
 model <- mlefit$model #take out assumed model with given data
 nC <- model$nC #number of contributors
 nodeg <- is.null(model$kit) #check for degradation
 if(nodeg) theta <- c(theta[1:(nC+1)],1) #insert beta variable equal 1
 if(!is.null(model$xi)) {
  theta <- c(theta,as.numeric(model$xi))
 } else {
  if(nodeg)  theta <- c(theta,theta2[np]) #insert fitted xi last again
 }
 locs <- names(model$popFreq)
 nL <- length(locs)
 mvec <- mlefit$fit$thetahat2[1:nC]
 mu <- theta[nC]
 sigma <- theta[nC+1]

 
if(is.null(model$xi)) { #stutter is unknown
  likYtheta <- function(yval) {   #call c++- function: length(theta)=nC+2
    ret$obsY[j] <- yval
    Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(model$prC), ret$condRef,as.numeric(model$threshT),as.numeric(model$fst),ret$mkvec,ret$nkval,as.numeric(model$lambda),as.numeric(ret$bp),as.integer(0),PACKAGE="gammadnamix")[[1]]
    return(exp(Cval + log(model$pXi(theta[ret$nC+3])))) #weight with prior of tau and 
  }
} else { #stutter is known. call c++- function: length(theta)=nC+1
  likYtheta <- function(yval) {   
    ret$obsY[j] <- yval
    Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(model$prC), ret$condRef,as.numeric(model$threshT),as.numeric(model$fst),ret$mkvec,ret$nkval,as.numeric(model$lambda),as.numeric(ret$bp),as.integer(0),PACKAGE="gammadnamix")[[1]]
    return(exp(Cval))
  }
}

  alpha <- 0.001 
  alpha2 <- alpha/sum(sapply(model$samples,function(x) sapply(x,function(y) length(y$adata)))) #"bonferroni outlier"
  maxYobs <- max(sapply(model$samples,function(x) sapply(x,function(y) max(y$hdata)))) #max observation
  maxYexp <- qgamma(1-alpha2,1/sigma^2,scale=2*mu*sigma^2) #max observation in theory
  maxY <- ceiling(max(maxYobs,maxYexp)) #get max observed
  minY <- model$threshT

  cumprobi <- numeric()
  for(loc in locs) { #traverse for each locus
   samples <- lapply(model$samples,function(x) x[loc])
   ret <- prepareC(nC=model$nC,samples,popFreq=model$popFreq[loc],refData=model$refData[loc],condOrder=model$condOrder,knownRef=model$knownRef,kit=model$kit)
   Yupper <- ret$obsY  #observed peak heights is upper limit in integral
   n <- length(Yupper) #number of observed peak heights 
   for(j in 1:n) {
    ret$obsY <- Yupper #reset observations
    num <- integrate(Vectorize(likYtheta),lower=minY,upper=Yupper[j])[[1]]
    denom <- integrate(Vectorize(likYtheta),lower=minY,upper=maxY)[[1]]
    val <- num/denom
    cumprobi <- c(cumprobi,val) #get cumulative probability
   }
  }
  N <- length(cumprobi)
  cumunif <-  punif((1:N)-0.5,0,N)

  #Goodness of fit test
  #pval <- ks.test(cumprobi, "punif")$p.value
  #txt <- paste0("p-value from Goodness-of-fit test = ",format(pval,digits=3))
  #print(txt)
  qqplot(cumunif,cumprobi,xlim=0:1,ylim=0:1,main="PP-plot between fitted model and theoretical model",xlab="Expected: Unif(0,1)",ylab="Observed: (Pr(Yj<=yj|Y_{-j}<=y_{-j},Yj>=thresh,model))")
  abline(0,1)
  #mtext(txt)
  return(cumunif)
}



