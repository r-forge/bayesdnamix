#' @title getContrCombs
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description getContrCombs is a helpfunction to specify all possible combined genotypes for a give hypothesis.
#' @export
#' @details Note: If a reference is unselected, it as numeric() in its loci
#' 
#' The allele-names are coded with a ascii-symbol to have threat more than 10-alleles in a loci.
#' 
#' @param Alist Vector of possible alleles 
#' @param nC Number of contributors in model.
#' @param symmetry Boolean of whether all combinations should be considered (permutation/combination)
#' @param refs List of contributors [[s]] where list-element 's' is the reference index with genotype vector ("A1","A2"). #' @param condOrder Specify conditioning references in refs (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model.
#' @param complete Boolean for requiring the set of combinations to fill the Alist. 
#' @return A matrix where each row represents a genotype combination, while each column correspond to a contributor


getContrCombs <- function(Alist,nC,symmetry=FALSE,refs=NULL,condOrder=NULL,complete=TRUE) {
 #Alist is Possible alleles.
 #symmetry=true gives all combinations
 #ref[[k]]=c("A1","A2") is a list of alleles to condition on.
 #condOrder gives order of references in system. It must be same length as number of references in refs.
 #complete: Boolean for requiring the set of combinations to fill the Alist.
 require(gtools)

 #helpfunctions:
 toAscii = function(x) { #takes vector with integers 1 to 62
  asciRange = c(48:57,65:90,97:122) #max range of 62. This is type constant
  return( sapply( as.raw(asciRange[x]) , FUN = rawToChar) )
 }
 #function taking in a set, number of samples and boolean as ordered
 getPerm = function(set,n,ordered) {
  if(!ordered) perm = combinations(length(set),n,v=set,repeats.allowed=TRUE)
  if(ordered) perm = permutations(length(set),n,v=set,repeats.allowed=TRUE)
  return(perm)
 }
 #method to get only allele with all contributors:
 removeCombsNotInA = function(combs,A) {
  temp = rep(TRUE,dim(combs)[1])
  bool = temp
  for(i in 1:length(temp)) temp[i] = paste(combs[i,],collapse='|')
  for(i in 1:length(A))   bool =   bool & grepl(A[i],temp)
  return(combs[bool,])
 }
 #function that goes through perm-set of alleles and get genoypes using paste:
 permGenotype = function(permset) {
  permGen = numeric(0)
  for(i in 1:nrow(permset)) {
   permGen = c(permGen,paste(permset[i,],collapse=''))
  }
  return(permGen)
 }
 Ai = toAscii(1:length(Alist))
 singleCombs = getPerm(Ai,2,ordered=FALSE) #always two alleles pr. ind.
 if(any(Alist%in%c('X','Y'))) {
  singleCombs = singleCombs[!rowSums(singleCombs==1)==2,]  #if amel: REMOVE YY-combination!
  if(is.null(dim(singleCombs))) singleCombs <- matrix(singleCombs,nrow=1)
 }
 genSet = permGenotype(singleCombs) #collapse genotypes
 #the alleles in ref must be a subset of alleles in Alist
 ncond = rep(FALSE,length(condOrder)) #boolean of vector wheter the condition is fine
 if(!is.null(refs)) { #if references given
  nRefs = length(refs) #number of references selected from import-frame
  condCC = rep(NA,nRefs) #vector of conditional combinations 
  for(k in 1:nRefs) { #for each reference we want to keep only combinations where these are included
   #check if all ref-alleles are in data A:
   #add the reference if both alleles are included:
   if(condOrder[k]==0) next #no condidtioning on the k-te reference
   condCC[k] = '' #init as no string
   if(length(refs[[k]])==0) next #continue if no reference there
   if(!all(refs[[k]]%in%Alist)) next #require both alleles to be in Alist
   CC = toAscii (which(Alist%in%refs[[k]])) #get combinations 
   if(length(CC)==1) condCC[k] = paste(rep(CC,2),collapse='') #if hom
   if(length(CC)==2) condCC[k] = paste(sort(CC),collapse='') #if het, its sorted
   ncond[k] = TRUE #the conditional reference was accepted to use!
  }
 }
 nUnk <- nC-sum(ncond) #number of unknowns
 if(nUnk==0) multiCombs <- getPerm(condCC[ncond],nC,ordered=symmetry)
 if(nUnk>0) multiCombs = getPerm(genSet,nUnk,ordered=symmetry) #get unknown CC
 if(is.null(dim(multiCombs)))  return(t(as.matrix(multiCombs))) #return one combination
 #combine unknown CC with reference CC:
 #reference CC are ordered as in condOrder.
 if(!is.null(refs) & sum(ncond)>0) { #important to use the ncond here as check here!
  finalcombs = matrix(NA,ncol=nC,nrow=nrow(multiCombs))
  for(k in 1:length(ncond)) { #fill each condition
   if(!ncond[k]) next #skip ref if not included
   finalcombs[,condOrder[k]] = condCC[k] #or else insert conditional profile
  }
  finalcombs[,-condOrder[ncond]] = multiCombs #set in the randoms for those not condidtioned on
  multiCombs = unique(finalcombs) #NEW:added unique here
 }
 if(complete) multiCombs = removeCombsNotInA(multiCombs,Ai)
 if(is.null(dim(multiCombs)))  return(t(as.matrix(multiCombs)))
 return(multiCombs)
} 

