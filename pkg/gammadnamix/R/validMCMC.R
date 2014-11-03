#' @title validMCMC
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description Validates aposteriori samples from MCMC method
#' @details This function takes samples from the MCMC as given in a matrix and shows the aposterior functions.
#' @param mcmcfit A object returned by contLikMCMC
#' @param trace Boolean for whether showing trace of samples.
#' @param acf Boolean for whether showing autocorrelation function of samples.
#' @export
validMCMC <- function(mcmcfit,trace=TRUE,acf=TRUE) {
 txt <- colnames(mcmcfit$posttheta)
 p <- length(txt)
 par(mfrow=c(p,1+sum(c(trace,acf)) ),mar = c(1,1,1,1), mgp = c(0,0.2,0))
 for(i in 1:p) {
  if( grepl("mx",txt[i])) {
   dens <- density(mcmcfit$posttheta[,i],from=0,to=1)
  } else if("mu"==txt[i] || "sigma"==txt[i] ) {
   dens <- density(mcmcfit$posttheta[,i],from=0)
  } else {
   dens <- density(mcmcfit$posttheta[,i])
  }
  xrange <- range(mcmcfit$posttheta[,i])
  mled <-dnorm(dens$x,mcmcfit$MLE[i],sqrt(mcmcfit$Sigma[i,i])) #density of lazy bayes
  plot(dens$x,dens$y,ty="l",main=txt[i],xlab="",ylab="",ylim=c(0,max(mled,dens$y)),xlim=xrange )
  lines(dens$x,mled,col=2,lty=2,ylab="",xlab="")
  if(trace) plot(mcmcfit$posttheta[,i],ty="l",ylab="",xlab="")
  if(acf) acf(mcmcfit$posttheta[,i],lag.max=200,ylab="",xlab="")
 }
 par(mfrow=c(1,1))
}

