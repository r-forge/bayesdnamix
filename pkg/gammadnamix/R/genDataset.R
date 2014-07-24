#' @title genDataset
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description genDataset samples a random mixture peak heights given as gamma(rho*sum(h_k),tau), with h_k as peak height of k-te contributor.
#' @details Simple function for generating a random peak height mixture distributed as normal positive truncated.
#' @param nC Number of contributors in model.
#' @param popFreq A list of allele frequencies for a given population.
#' @param mu Expected peak heights for a het. single contributor allele
#' @param sd Standard deviation of peak heights for a het. single contributor allele
#' @param sorted Boolean for wheter sorting the contributors with respect to decreasingly mixture proportions.
#' @param threshT Required allele peak height in mixture
#' @param refData A list of given reference profiles given as refData[[i]][[s]]. Default is random from population. 
#' @param Mx A vector of known mixture proportions. Default is random uniform.
#' @return List with elements theta,mixData,refData,nDrop where theta is the true parameters of the model and nDrop number of dropped alleles for each loci.
#' @export
#' @examples
#' \dontrun{
#' load(system.file("mcData.Rdata", package = "gammadnamix"))
#' mixdata <- genDataset(nC=2,popFreq=data$popFreq,mu=1000,sd=100,threshT=50) #generate a random mixture
#' }
genDataset = function(nC,popFreq,mu=5000,sd=400,sorted=FALSE,threshT=50,refData=NULL,Mx=NULL) {
  #mu = total expected sum peak height
  #sd = standard deviations
  nL<-length(popFreq)
  if(is.null(Mx)) {
   Mx <- rgamma(nC,1)
   Mx=Mx/sum(Mx) #rdirichlet(1,rep(1,nC))  #simulate Mx for contributors
  } else {
   if(length(Mx)!=nC) stop("Length of Mx not equal nC!")
  }
  if(sorted)  Mx  <- sort(Mx,decreasing=TRUE)
  Mx <- c(Mx)
  if(is.null(refData )) refData <- list() 
  mixData <- list(adata=list(),hdata=list()) 
  nDrop <- rep(NA,nL)

  #convert to gamma-parameters
  rho <- (mu/sd)^2
  tau <- sd^2/mu 

  for(i in 1:nL) {
   if( length(refData) < i) { #if no refData given
    refData[[i]] <- list()
    mixA = numeric()
    mixH = numeric()
    for(s in 1:nC) {
     Asim <- refData[[i]][[s]] <-  sample(names(popFreq[[i]]),size=2,prob=popFreq[[i]],replace=TRUE)
     mixA = c(mixA,Asim) #keep not droped
    }
   } else {
    mixA <- unlist(refData[[i]]) #vectorize
   } 
   mixH = c(t(replicate(2,Mx)))
   agg=aggregate(mixH,by=list(mixA),sum) #aggregate contributions
   mixH <- rgamma(length(agg$x),shape=rho*agg$x,scale=tau) #shape/scale given
   dropped <- mixH<threshT 
   nDrop[i] <- sum(dropped) #number of dropped alleles on each loci
   mixData$adata[[i]] = agg$Group.1[!dropped]
   mixData$hdata[[i]] = mixH[!dropped]
  }
  locs <- names(popFreq)
  if(is.null(locs)) locs = paste("Loci",1:length(popFreq),sep="")
  names(mixData$adata) <- locs
  names(mixData$hdata) <- locs
  names(refData) <- locs
  return(list(theta=list(Mx=Mx,rho=rho,tau=tau),mixData=mixData,refData=refData,nDrop=nDrop))
}


