#' @title prepareC 
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description prepareC is used in the functions contLikMLE,contLikMarg and contLikMCMC to prepare input to C-call
#' @details The function builds the data input to the C-code
#' @param nC Number of contributors in model.
#' @param mixData Evidence object with list elements adata[[i]] and hdata[[i]]. Each element has a loci-list with list-element 'i' storing qualitative data in 'adata' and quantitative data in 'hdata'.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects with list element [[s]]$adata[[i]]. The list element has reference-list with list-element 's' having a loci-list adata with list-element 'i storing qualitative data.
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. condOrder=-1 means the reference is known-non contributor!
#' @param knownRef Specify known references from refData (index). For instance knownRef=(1,2) means that reference 1 and 2 is known allele samples in the hypothesis. This is affected by fst-correction.
#' @return ret A list of data input to call the C-code with
#' @export 

prepareC = function(nC,mixData,popFreq,refData=NULL,condOrder=NULL,knownRef=NULL){
 nA = unlist(lapply(mixData$adata,length)) #number of alleles of selected loci
 nL <- length(nA) #number of loci in data
 locnames <- toupper(names(mixData$adata))
 missind <- !locnames%in%toupper(names(popFreq))  #indice of missing loci
 if(any(missind)) {
  stop(paste('Missing locus (',locnames[missind],') in popFreq.',sep=''))
 }
 if(is.null(condOrder)) {
  condOrder <- rep(0,nC) #insert condorder if missing
  nK <- 0
 } else { #check that they are unique and in right order
  tmp <- condOrder[condOrder>0]
  nK <- length(tmp)
  if( nK!=length(unique(tmp)) ) stop("Specify unique positions!")
  if( any( sort(tmp,decreasing=FALSE)!=(1:nK)) ) stop("Please condition references starting from 1. position")
  if( nK==nC ) stop("Not implemented yet!")
 }
 
 #convertion of values in popFreq, mixData and Glist$G:
 #loci-order follows as in mixData: "locnames". Rearrange names:
 names(popFreq) <- toupper(names((popFreq))) #make invariant to capital
 popFreq <- popFreq[locnames] #order popFreq to mixData-order

 #get population genotypes:
 getGlist <- function(popFreq) {
  locs <- names(popFreq)
  Glist <- list()
  for(i in 1:length(locs)) {
   G = t(as.matrix(expand.grid(rep(list(as.numeric(names(popFreq[[i]])),as.numeric(names(popFreq[[i]])) ))))) #one genotype per column
   keep = G[2,]>=G[1,] #unique genotypes 
   G <- G[,keep]  #store genotypes
   G <- matrix(as.character(G),nrow=2) #make string names again
   tmpP = t(as.matrix(expand.grid(rep(list(as.numeric(popFreq[[i]]),as.numeric(popFreq[[i]]) )))))
   Gprob = exp(colSums(log(tmpP[,keep]))) #get allele probs
   ishet = G[1,]!=G[2,]
   Gprob[ishet] = 2*Gprob[ishet] #multiply with two to get heterozygote prob
   Glist[[locs[i]]] <- list(G=t(G),Gprob=Gprob)
  }
  return(Glist)
 }

 #Fix references: Assign condition to condM-matrix
 Glist <- getGlist(popFreq) #get population genotype information
 condM <- matrix(-1,nrow=nL,ncol=nC) #default is no references (=-1)
 #assign references to condM-matrix by values of Glist
 if(!is.null(refData) && !is.null(condOrder) && any(condOrder>0)) {
  for(i in 1:nL) {
   subRef <- refData[[locnames[i]]]
   Gset <- Glist[[locnames[i]]]$G #genotype combinations
   if(length(subRef)==0)  stop(paste('Missing locus (',locnames[i],') in refData.',sep=''))
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
 nkval <- rep(0,nL) #total number of sampled alleles in each marker
 for(i in 1:nL) {
  tmp <- rep(0, length(popFreq[[i]]))
  if(!is.null(condOrder) & !is.null(refData)) {
   for(k in 1:length(condOrder)) {
    if(condOrder[k]!=0) { 
     ind <- which( names(popFreq[[i]])%in%refData[[locnames[i]]][[k]] )
     tmp[ind] = (length(ind)==1) + 1
    }
   }
  }
  if(!is.null(knownRef)) { #known non-contributors
   for(k in knownRef) {
    ind <- which( names(popFreq[[i]])%in%refData[[locnames[i]]][[k]] )
    tmp[ind] = (length(ind)==1) + 1 #add twice sampled if homozygote
   }
  }
  nkval[i] <- sum(tmp) #number of sampled (for each loci)
  mkvec <- c(mkvec,tmp) 
 }

 #decode allele-names to index names + Stutter-preparation
 allAbpind <- as.numeric()
 for(i in 1:nL) {
   anames <- names(popFreq[[i]]) #old names
   #find corresponding index of 1-ahead pb-twin
   numAnames <- as.numeric(anames)
   BPind <- rep(0,length(numAnames)) #init 0 if allele doesn't have any bp-twin
   for(j in 1:length(numAnames)) {
    bpind <- which((numAnames[j]-1)==numAnames) #store index of what allele it stutters to
    if(length(bpind)>0) BPind[j] <- bpind #assign index-placement
   }
   allAbpind <- c(allAbpind,BPind)
   anames2 <- 0:(length(popFreq[[i]])-1) #new names
   names(popFreq[[i]]) <- anames2
   for(j in 1:nA[i]) { #for each allele
    mixData$adata[[i]][j] <- anames2[anames==mixData$adata[[i]][j]] #update name
   }
  }
 #fix genotypes:
 Glist <- getGlist(popFreq) #get population genotype information
 Gprob <- lapply(Glist,function(x) return(x$Gprob))
 Gset <- lapply(Glist,function(x) return(x$G))

 #Fix input-variables for Cpp-function
 CnA <- c(0,cumsum(nA))
 allA <- unlist(mixData$adata)
 allY <- unlist(mixData$hdata)
 nG <- sapply(Gprob,length) #number of genotype combinations
 CnG <- c(0,cumsum(nG))
 CnG2 <- c(0,cumsum(nG*2)) #note: 2 columns for each genotype!!
 pG <- unlist(Gprob) #vectorize over all loci
 Gvec <- as.integer(rbind(unlist(Gset))) #vectorize a big matrix (loci are put chronologic)
 condRef <- c(condM) #vectorized over all loci
 nAall <- sapply(popFreq,length) #Number of population-alleles on each loci
 CnAall <- c(0,cumsum(nAall)) #cumulative number of alleles
 pA <- unlist(popFreq) #need each allele probability for drop-in probabilities

 retlist <- list(nC=as.integer(nC),nK=as.integer(nK),nL=as.integer(nL),nA=as.integer(nA), allY=as.numeric(allY),allA=as.integer(allA),CnA=as.integer(CnA),allAbpind=as.integer(allAbpind),nAall=as.integer(nAall),CnAall=as.integer(CnAall),Gvec=as.integer(Gvec),nG=as.integer(nG),CnG=as.integer(CnG),CnG2=as.integer(CnG2),pG=as.numeric(pG),pA=as.numeric(pA), condRef=as.integer(condRef),mkvec=as.integer(mkvec),nkval=as.integer(nkval) )
 return(retlist)
} #end function

