#' @title validMLEmodel
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description validMLEmodel makes model validation whether the observed peak heights fits the maximum likelihood fitted gamma distribution.
#' @details The cumulative probability of the observed allele peaks are calculated and compared with a uniform distribution.
#' Function calls procedure in c++ by using the package Armadillo and Boost.
#'
#' @param mlefit Fitted object using contLikMLE
#' @param plotYdistr A boolean whether to plot the model fitted peak height density for each contributors
#' @param plotDdistr A boolean whether to plot the model fitted heterozygote dropout probabilities as a function of amount of DNA(mu) for each contributors
#' @return ret A vector for each marker with cumulative probabilities
#' @export

validMLEmodel <- function(mlefit,plotYdistr=FALSE,plotDdistr=FALSE) {
 theta <- mlefit$fit$thetahat #condition on mle parameter
 np <- length(theta) #number of unknown parameters
 model <- mlefit$model #take out assumed model with given data
 xi <- model$xi
 if(!is.null(xi)) theta <- c(theta,as.numeric(xi))
 locs <- names(model$popFreq)
 nL <- length(locs)
 nC <- model$nC #number of contributors
 mvec <- mlefit$fit$thetahat2[1:nC]
 mu <- theta[nC]
 sigma <- theta[nC+1]
 rho <- sigma^(-2)
 tau <- mu*sigma^2
 
#Show assumed peak height distribution for a non-stuttered single heterozygote allele:
if(plotYdistr) {
 rfu <- model$threshT: qgamma(0.99,rho,scale=tau) 
 maxY <- max(dgamma(rfu ,rho*min(mvec),scale=tau))
 plot(0,0,xlim=range(rfu),ylim=c(0,maxY),ty="n",xlab="Peak height Y",ylab="f(Y)",main="MLE fitted het. peak height density for each contributor")
 abline(v=model$threshT,col="gray")
 for(m in 1:length(mvec)) {
  lines(rfu,dgamma(rfu ,rho*mvec[m],scale=tau),ty="l",col=m)
  abline(h=0)
 }
 legend("topright",legend=paste("Contributor ",1:nC),col=1:nC,lty=1)
 x11()
}

#Show assumed dropout distribution as a function of amount of dna mu:
if(plotDdistr) {
  threshT <- mlefit$model$threshT 
  muvec2 <- seq(1,2*mu,l=100)
  plot(0,0,xlim=c(0,max(muvec2 )),ylim=c(0,1),ty="n",xlab=paste0("Amount of DNA (mu)"),main=paste0("Het. dropout distribution for each contributors"),ylab="Probability of dropout")
  for(k in 1:nC)  lines(muvec2,pgamma(threshT, rho*mvec[k],scale=muvec2*sigma^2),col=k)
  legend("topright",legend=paste0("Contributor ",1:nC,"(mx=",round(mvec,2),")") ,col=1:nC,lty=1)
  x11()
}

if(is.null(xi)) { #stutter is unknown
  likYtheta <- function(yval) {   #call c++- function: length(theta)=nC+2
    ret$obsY[j] <- yval
    Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(model$prC), ret$condRef,as.numeric(model$threshT),as.numeric(model$fst),ret$mkvec,ret$nkval,as.numeric(model$lambda),as.integer(0),PACKAGE="gammadnamix")[[1]]
    return(exp(Cval + log(model$pXi(theta[ret$nC+2])))) #weight with prior of tau and 
  }
} else { #stutter is known. call c++- function: length(theta)=nC+1
  likYtheta <- function(yval) {   
    ret$obsY[j] <- yval
    Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(model$prC), ret$condRef,as.numeric(model$threshT),as.numeric(model$fst),ret$mkvec,ret$nkval,as.numeric(model$lambda),as.integer(0),PACKAGE="gammadnamix")[[1]]
    return(exp(Cval))
  }
}
  alpha <- 0.01
  alpha2 <- alpha/sum(sapply(model$samples,function(x) sapply(x,function(y) length(y$adata)))) #"bonferroni outlier"
  maxYobs <- max(sapply(model$samples,function(x) sapply(x,function(y) max(y$hdata)))) #max observation
  maxYexp <- qgamma(1-alpha2,2*nC*rho,scale=tau) #max observation in theory
  maxY <- ceiling(max(maxYobs,maxYexp)) #get max observed
  minY <- model$threshT

  cumprobi <- numeric()
  for(loc in locs) { #traverse for each locus
   samples <- lapply(model$samples,function(x) x[loc])
   ret <- prepareC(nC=model$nC,samples,popFreq=model$popFreq[loc],refData=model$refData[loc],condOrder=model$condOrder,knownRef=model$knownRef)
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
  pval <- ks.test(cumprobi, "punif")$p.value
  txt <- paste0("p-value from Goodness-of-fit test = ",format(pval,digits=3))
  print(txt)
  qqplot(cumunif,cumprobi,xlim=0:1,ylim=0:1,main="PP-plot between fitted model and theoretical model",xlab="Expected: Unif(0,1)",ylab="Observed: (Pr(Yj<=yj|Y_{-j}<=y_{-j},Yj>=thresh,model))")
  abline(0,1)
  mtext(txt)
 
  return(pval)
}



