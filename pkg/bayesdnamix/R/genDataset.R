#' @title genDataset
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description genDataset samples a random mixture with Normal[+] peak heights.
#' @export
#' @details Simple function for generating a random peak height mixture distributed as normal positive truncated.
#' @param nC Number of contributors in model.
#' @param popFreq A list of allele frequencies for a given population.
#' @param mu Expected peak heights for an allele of a contributor.
#' @param sd Standard deviation of peak heights for an allele of a contributor.
#' @return List with elements trueMx,mixData,refData

genDataset = function(nC,popFreq,mu=10000,sd=1000) {
  require(gtools)
  nL<-length(popFreq)
  Mx=rdirichlet(1,rep(1,nC))  #true Mx for contributors
  refData <- list() 
  mixData <- list(adata=list(),hdata=list()) 
  for(i in 1:nL) {
   refData[[i]] <- list()
   mixA = numeric()
   mixH = numeric()
   for(s in 1:nC) {
    refData[[i]][[s]] = sample(names(popFreq[[i]]),size=2,prob=popFreq[[i]],replace=TRUE)
    mixA = c(mixA, refData[[i]][[s]] )
    mixH = c(mixH, qnorm(runif(2,0.5,1))*sd+(mu*Mx[s])) #generate with N+(mu,sd)
#    mixH = c(mixH, abs(rnorm(2,mu*Mx[s],sd))) #generate with N+(mu,sd)
   }
   agg=aggregate( mixH,by=list(mixA),sum)
   mixData$adata[[i]] = agg[,1]
   mixData$hdata[[i]] = agg[,2]
  }
  locs <- names(popFreq)
  if(is.null(locs)) locs = paste("Loci",1:length(popFreq),sep="")
  names(mixData$adata) <- locs
  names(mixData$hdata) <- locs
  names(refData) <- locs
  return(list(trueMx=Mx,mixData=mixData,refData=refData))
}


