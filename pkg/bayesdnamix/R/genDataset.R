#' @title genDataset
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description genDataset samples a random mixture with Normal[+] peak heights.
#' @export
#' @details Simple function for generating a random peak height mixture distributed as normal positive truncated.
#' @param nC Number of contributors in model.
#' @param popFreq A list of allele frequencies for a given population.
#' @param mu Expected peak heights for an allele of a contributor.
#' @param sd Standard deviation of peak heights for an allele of a contributor.
#' @param PrD Fixed Allele drop-out probability for each contributors. 
#' @param sorted Boolean for wheter sorting the contributors with respect to decreasingly mixture proportions.
#' @return List with elements Mx,mixData,refData

genDataset = function(nC,popFreq,mu=10000,sd=1000,PrD = 0,sorted=FALSE) {
  require(gtools)
  nL<-length(popFreq)
  Mx=rdirichlet(1,rep(1,nC))  #true Mx for contributors
  if(sorted)  Mx  <- sort(Mx,decreasing=TRUE)
  refData <- list() 
  mixData <- list(adata=list(),hdata=list()) 
  if(length(PrD)==1) PrD <- rep(PrD,nC)
  if(length(PrD)!=nC) stop("Wrong size of PrD given")
  for(i in 1:nL) {
   refData[[i]] <- list()
   mixA = numeric()
   mixH = numeric()
   for(s in 1:nC) {
    Asim <- refData[[i]][[s]] <-  sample(names(popFreq[[i]]),size=2,prob=popFreq[[i]],replace=TRUE)
    Hsim <- qnorm(runif(2,0.5,1))*sd+(mu*Mx[s])
#    Hsim <-  abs(rnorm(2,mu*Mx[s],sd))
    droped <- sample(0:1,2,replace=TRUE,prob=c(1-PrD[s],PrD[s]))
    mixA = c(mixA, Asim[droped==0]) #keep not droped
    mixH = c(mixH, Hsim[droped==0]) #keep not droped
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
  return(list(Mx=Mx,mixData=mixData,refData=refData))
}


