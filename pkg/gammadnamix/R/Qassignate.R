#' @title Qassignate
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description Q-assignation
#' @details Assignes non-shown alleles as one single allele.
#' @param mixData Evidence object with list elements adata[[i]] and hdata[[i]]. Each element has a loci-list with list-element 'i' storing qualitative data in 'adata' and quantitative data in 'hdata'.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects with list element [[s]]$adata[[i]]. The list element has reference-list with list-element 's' having a loci-list adata with list-element 'i storing qualitative data.
#' @return ret A list(popFreq,refData) with Q-assignated alleles. 
#' @export

Qassignate <- function(mixData,popFreq,refData=NULL) {
 popFreq2 <- popFreq
 refData2 <- refData
 for(i in 1:length(popFreq)) { #make Q-assignation for each loci
  tmp <- popFreq[[i]][names(popFreq[[i]])%in%mixData$adata[[i]]]
  tmp <- c(tmp,1-sum(tmp))
  names(tmp)[length(tmp)] <- "99"
  popFreq2[[i]] <- tmp
  if(!is.null(refData)) { #insert 99 as default allele of missing refs
   newP <- names(popFreq2[[i]]) 
   if(!all(unlist(refData[[i]])%in%newP)) { #there was some missing alleles
    for(k in 1:length(refData[[i]])) {
     refData2[[i]][[k]][!refData[[i]][[k]]%in%newP] <- "99" #insert missing
    }
   }
  }
 }
 return(list(popFreq=popFreq2,refData=refData2))
}
