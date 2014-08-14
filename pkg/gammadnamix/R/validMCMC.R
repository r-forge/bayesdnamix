#' @title validMCMC
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description Validates aposteriori samples from MCMC method
#' @details This function takes samples from the MCMC as given in a matrix and shows the aposterior functions.
#' @param X A(Mxp) matrix with MCMC samples.
#' @param trace Boolean for whether showing trace of samples.
#' @param acf Boolean for whether showing autocorrelation function of samples.
#' @param mlefit Fitted object using contLikMLE
#' @param isPhi Boolean whether we are working with reparameterizited parameters
#' @export
validMCMC <- function(X,trace=TRUE,acf=TRUE,mlefit=NULL,isPhi=FALSE) {
 p <- ncol(X)
 par(mfrow=c(p,1+sum(c(trace,acf)) ),mar = c(1,1,1,1), mgp = c(0,0.2,0))
 for(i in 1:p) {
  if(!isPhi) {
   txt <- paste0("theta",i)
  } else {
   txt <- paste0("phi",i)
  }
  dens <- density(X[,i])
  plot(dens$x,dens$y,ty="l",main=txt,xlab="",ylab="")
  if(!is.null(mlefit)) {
   if(isPhi) lines(dens$x,dnorm(dens$x,mlefit$theta0[i],sqrt(mlefit$Sigma0[i,i])),col=2,lty=2,ylab="",xlab="")
   if(!isPhi) lines(dens$x,dnorm(dens$x,mlefit$theta1[i],sqrt(mlefit$Sigma1[i,i])),col=2,lty=2,ylab="",xlab="")
  }
  if(trace) plot(X[,i],ty="l",ylab="",xlab="")
  if(acf) acf(X[,i],lag.max=200,ylab="",xlab="")
 }
 par(mfrow=c(1,1)) 
}

