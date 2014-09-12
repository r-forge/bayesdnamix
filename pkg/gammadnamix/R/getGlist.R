#' @title getGlist 
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description getGlist Returns a list of genotypes with corresponding probabilities for a given allele frequency-list.
#' @details The function returns the list of all possible observed genotypes, with corresponding Hardy-Weinberg assumed probabilities. The allele-names in popFreq needs to be numbers.
#' @param popFreq A list of allele frequencies for a given population. Each named element in the list must be a allele-named vector with allele frequencies. 
#' @return Glist A list with genotypes and genotype probabilities for each locus.
#' @export 

 getGlist <- function(popFreq) {
  locs <- names(popFreq)
  Glist <- list()
  for (i in 1:length(locs)) {
   G = t(as.matrix(expand.grid(rep(list(as.numeric(names(popFreq[[i]])),as.numeric(names(popFreq[[i]])))))))
   keep = G[2, ] >= G[1, ]
   G <- G[, keep]
   G <- matrix(as.character(G), nrow = 2)
   tmpP = t(as.matrix(expand.grid(rep(list(as.numeric(popFreq[[i]]),as.numeric(popFreq[[i]]))))))
   Gprob = exp(colSums(log(tmpP[, keep])))
   ishet = G[1, ] != G[2, ]
   Gprob[ishet] = 2 * Gprob[ishet]
   Glist[[locs[i]]] <- list(G = t(G), Gprob = Gprob)
  }
  return(Glist)
 }
