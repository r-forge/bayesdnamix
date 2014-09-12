#' @title prepareC 
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description prepareC is used in the functions contLikMLE,contLikMarg and contLikMCMC to prepare input to C-call
#' @details The function builds the data input to the C-code
#' @param nC Number of contributors in model.
#' @param samples A List with samples which for each samples has locus-list elements with list elements adata and hdata. 'adata' is a qualitative (allele) data vector and 'hdata' is a quantitative (peak heights) data vector.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects with list element [[s]][[i]]. The list element has reference-list with list-element 's' having a loci-list adata with list-element 'i storing qualitative data.
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. condOrder=-1 means the reference is known-non contributor!
#' @param knownRef Specify known references from refData (index). For instance knownRef=(1,2) means that reference 1 and 2 is known allele samples in the hypothesis. This is affected by fst-correction.
#' @return ret A list of data input to call the C-code with
#' @export 

prepareC = function(nC,samples,popFreq,refData=NULL,condOrder=NULL,knownRef=NULL){
 #Note: Supports Invariant order of markers and caseletters!!
 #Supports replicates in the samples. 
 #should support to have zero contribution in some markers (all dropped out)-> all y=0
 #Evaluates only intercept of mixture-markers and freq-markers
 nS <- length(samples) #number of replicates
 #test 1) Require same markers in all replicates. 
 locs <- lapply(samples,function(x) toupper(names(x)))  #marker names for each replicates
 nLs <- sapply(locs,length)
 nL <- unique(nLs)
 if(length(nL)>1) stop("Number of markers in the replicates was not the same")
 locs <- unique(unlist(locs)) #get unique markernames
 if(nL!=length(locs)) stop("Different markers was specified in the replicates")
 names(popFreq) <- toupper(names(popFreq)) #toupper case!
 locs <- locs[locs%in%names(popFreq)] #take intercept with 
 #print( paste0("Evaluated loci: ", paste0(locs,collapse=",") ) )

 if(is.null(condOrder)) {
  condOrder <- rep(0,nC) #insert condorder if missing
  nK <- 0 #number of known contributors
 } else { #check that they are unique and in right order
  tmp <- condOrder[condOrder>0]
  nK <- length(tmp) #number of known contributors
  if( nK!=length(unique(tmp)) ) stop("Specify unique positions!")
  if( any( sort(tmp,decreasing=FALSE)!=(1:nK)) ) stop("Please condition references starting from 1. position")
 }
 
 #convertion of values in popFreq, mixData and Glist$G:
 #loci-order follows as in mixData: "locs". Rearrange names:
 popFreq <- popFreq[locs] #order popFreq to mixData-order

 #Fix references: Assign condition to condM-matrix
 Glist <- getGlist(popFreq) #get population genotype information
 condM <- matrix(-1,nrow=nL,ncol=nC) #default is no references (=-1)
 #assign references to condM-matrix by values of Glist
 if(!is.null(refData) && !is.null(condOrder) && any(condOrder>0)) {
  names(refData) <- toupper(names(refData)) #toupper case!
  for(i in 1:nL) {
   subRef <- refData[[locs[i]]] #take out relevant reference
   if(length(subRef)==0)  stop(paste('Missing locus (',locs[i],') in refData.',sep=''))
   Gset <- Glist[[locs[i]]]$G #genotype combinations
   for(k in 1:length(subRef)) { #for each reference
    if(condOrder[k]>0) {
     Gind1 <- subRef[[k]][1]==Gset[,1] & subRef[[k]][2]==Gset[,2]
     Gind2 <- subRef[[k]][2]==Gset[,1] & subRef[[k]][1]==Gset[,2]
     condM[i,condOrder[k]] = which(Gind1 | Gind2) - 1 #subtract with one since we work from 0-indice
    }
   }
  }
 }
 #fix known sample-information:
 mkvec <- numeric() #number of times each alleles are sampled
 nkval <- numeric() #total number of sampled alleles in each marker
 for(loc in locs) {
  tmp <- rep(0, length(popFreq[[loc]]))
  if(!is.null(condOrder) & !is.null(refData)) {
   for(k in 1:length(condOrder)) {
    if(condOrder[k]!=0) { 
     ind <- which( names(popFreq[[loc]])%in%refData[[loc]][[k]] )
     tmp[ind] = (length(ind)==1) + 1 #add twice sampled if homozygote
    }
   }
  }
  if(!is.null(knownRef)) { #known non-contributors
   for(k in knownRef) {
    ind <- which( names(popFreq[[loc]])%in%refData[[loc]][[k]] )
    tmp[ind] = (length(ind)==1) + 1 #add twice sampled if homozygote
   }
  }
  nkval <- c(nkval,sum(tmp)) #number of sampled (for each loci)
  mkvec <- c(mkvec,tmp) 
 }

 #decode allele-names to index names + Stutter-preparation
 allAbpind <- as.numeric()
 for(loc in locs) { #for each marker
   anames <- names(popFreq[[loc]]) #old names
   #find corresponding index of 1-ahead pb-twin
   numAnames <- as.numeric(anames) #convert from string to numbers
   BPind <- rep(0,length(numAnames)) #init 0 if allele doesn't have any bp-twin
   for(j in 1:length(numAnames)) {
    bpind <- which((numAnames[j]-1)==numAnames) #store index of what allele it stutters to
    if(length(bpind)>0) BPind[j] <- bpind #assign index-placement
   }
   allAbpind <- c(allAbpind,BPind)
   anames2 <- 0:(length(popFreq[[loc]])-1) #new names
   names(popFreq[[loc]]) <- anames2 #update names in popFreq
   #go through each observed alleles in samples and update names
   for(s in 1:nS) { #for each sample
    obsA <- samples[[s]][[loc]]$adata  #observed alllees
    if(length(obsA)>0) { #if atleast 1 observed
     for(j in 1:length(obsA)) samples[[s]][[loc]]$adata[j] <- anames2[anames==obsA[j]] #update allele-name
    } #dont do anything if non-observed!
   } #end for each samples
 } #end for each marker
 #fix genotypes:
 Glist <- getGlist(popFreq) #get population genotype information
 Gprob <- lapply(Glist,function(x) return(x$Gprob))
 Gset <- lapply(Glist,function(x) return(x$G))

 #take into account for replicates here!
 #count and vectorize alleles in mixtures here:
 nA <- obsA <- obsY <- numeric()
 for(loc in locs) { #for each marker
  for(s in 1:nS) { #for each replicates
   nA <- c(nA, length(samples[[s]][[loc]]$adata)) #count number of alleles 
   obsA <- c(obsA, samples[[s]][[loc]]$adata) #count number of alleles 
   obsY <- c(obsY, samples[[s]][[loc]]$hdata) #count number of alleles 
  }
 }
 CnA <- c(0,cumsum(nA)) #cumulative number of alleles for each marker

 #dependent on popFreq and ref-conditions
 nG <- sapply(Gprob,length) #number of genotype combinations
 CnG <- c(0,cumsum(nG))
 CnG2 <- c(0,cumsum(nG*2)) #note: 2 columns for each genotype!!
 pG <- unlist(Gprob) #vectorize over all loci
 Gvec <- as.integer(rbind(unlist(Gset))) #vectorize a big matrix (loci are put chronologic)
 condRef <- c(condM) #vectorized over all loci
 nAall <- sapply(popFreq,length) #Number of population-alleles on each loci
 CnAall <- c(0,cumsum(nAall)) #cumulative number of alleles
 pA <- unlist(popFreq) #need each allele probability for drop-in probabilities

 retlist <- list(nC=as.integer(nC),nK=as.integer(nK),nL=as.integer(nL),nA=as.integer(nA), obsY=as.numeric(obsY),obsA=as.integer(obsA),CnA=as.integer(CnA),allAbpind=as.integer(allAbpind),nAall=as.integer(nAall),CnAall=as.integer(CnAall),Gvec=as.integer(Gvec),nG=as.integer(nG),CnG=as.integer(CnG),CnG2=as.integer(CnG2),pG=as.numeric(pG),pA=as.numeric(pA), condRef=as.integer(condRef),mkvec=as.integer(mkvec),nkval=as.integer(nkval),nS=as.integer(nS))
 return(retlist)
} #end function

