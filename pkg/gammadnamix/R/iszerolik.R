#' @title iszerolik 
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description iszerolik Determine if the likelihood becomes zero
#' @details The function is determining whether likelihood vil be 0 (TRUE or FALSE) assuming zero drop-in 
#' @param evid A vector with allele names
#' @param ref A vector with conditioned alleles in given hypothesis
#' @param nU Number of unknown individuals in given hypothesis
#' @return TRUE/FALSE A boolean whether the likelihood vil be zero.
#' @export 

iszerolik <- function(evid,ref,nU,xi=0) {
   Ei2 <- evid[!evid%in%ref] #set of unknown alleles (not explained by ref0) 
   if(length(Ei2)<=(2*nU)) return(FALSE) #unknown contributors explains E
   if(!is.null(xi) && xi==0) return(TRUE) #too many unexplained alleles
   #case of assuming stutters:
   Ei3 <- Ei2[!as.character(Ei2)%in%as.character(as.numeric(ref)-1)] #set of unknown alleles (after explained by being stutter from ref0)
   Ei4 <- Ei3[as.character(as.numeric(Ei3)-1)%in%Ei3] #set of unknown alleles (after explaining possible stutter from unknown)
   Ei5 <- unique(c(Ei4,Ei3[!Ei3%in%as.character(as.numeric(Ei3)-1)])) #Alleles to explain removed stutter 
   if(length(Ei5)>(2*nU)) {
    return(TRUE)  #too many unexplained alleles
   } else {
    return(FALSE) #the alleles can be explained
   }
}
