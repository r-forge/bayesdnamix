#' @title validMCMC
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description Validates aposteriori samples from MCMC method
#' @details This function takes samples from the MCMC as given in a matrix and shows the aposterior functions.
#' @param X A(Mxp) matrix with MCMC samples.
#' @param trace Boolean for whether showing trace of samples.
#' @param acf Boolean for whether showing autocorrelation function of samples.
#' @export

validMCMC <- function(X,trace=TRUE,acf=TRUE) {
 p <- ncol(X)
 par(mfrow=c(p,1+sum(c(trace,acf)) ))
 for(i in 1:p) {
  if(i<nC) txt <- paste0("Mx",i)
  if(i==nC) txt <- paste0("mu (DNA amount)")
  if(i==nC+1) txt <- paste0("sigma (imbalancy)")
  if(i==nC+2) txt <- paste0("xi (stutter ratio)")
  plot(density(X[,i],from=0),ty="l",main=txt)
  if(trace) plot(X[,i],ty="l")
  if(acf) acf(X[,i],lag.max=200)
 }
 par(mfrow=c(1,1))
}
