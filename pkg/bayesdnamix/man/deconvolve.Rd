\name{deconvolve}
\alias{deconvolve}
\title{deconvolve}
\usage{
  deconvolve(nC, mixData, refData = NULL, condOrder = NULL,
    locsel_Mix = NULL, locsel_Ref = NULL, eps = 100,
    zeroMx = FALSE, threshT = 50, verbose = FALSE)
}
\arguments{
  \item{nC}{Number of contributors in model.}

  \item{mixData}{Evidence object with list elements
  adata[[i]] and hdata[[i]]. Each element has a loci-list
  with list-element 'i' storing qualitative data in 'adata'
  and quantitative data in 'hdata'.}

  \item{refData}{Reference objects with list element
  [[i]][[s]] where list-element 's' is the reference index
  and the list-element 'i is the loci index where the
  qualitative data is stored as a length two vector.}

  \item{condOrder}{Specify conditioning references from
  refData (must be consistent order). For instance
  condOrder=(0,2,1,0) means that we restrict the model such
  that Ref2 and Ref3 are respectively conditioned as 2.
  contributor and 1. contributor in the model.}

  \item{locsel_Mix}{Boolean-vector with Selected loci in
  mixData to deconvolve. locsel_Mix=NULL; accepts all
  loci.}

  \item{locsel_Ref}{Boolean-matrix for specifying
  conditional loci (row) for each reference
  (column).locsel_Ref=NULL; accepts all loci.}

  \item{eps}{Number of best combinations to keep during the
  search.}

  \item{zeroMx}{boolean of allowing zero mixture proportion
  as an estimate for any contributors.}

  \item{threshT}{Imputet quantitative value when
  conditioned reference alleles are non-observed.}

  \item{verbose}{Boolean for whether in-process information
  should be printed}
}
\value{
  Optimized deconvolution model object.
  \item{simpleList}{Table of loci independent
  optimizations. Uses independent covariance structure.}
  \item{pList}{Resultlist of optimized combinations,
  mixture proportions and error-distances (MD).}
  \item{locinames}{Name of loci in mixData}
  \item{result1}{Tabled optimized results in format 1.}
  \item{result2}{Tabled optimized results in format 2.}
  \item{data}{All data used as input in analysis.}
  \item{options}{Input parameters used in analysis.}
}
\description{
  deconvolve is a linear deconvolution procedure for STR
  DNA mixtures.
}
\details{
  The procedure optimizes the mixture proportion
  simultaneous with combined genotypes by assuming the STR
  response variation as normal distributed. The criterion
  for optimization is the error distance Mahalanobis
  Distance (MD) between the fitting model and observed
  responses.

  Conditioning on referenced genotypes is possible.
  Selection of conditioned loci for each of the references
  may be specified. Unobserved alleles from references will
  be imputed as observed alleles with the input threshold
  as the quantitative information. Non-selected or empty
  observed loci will return NA as genotype combinations and
  not treated in model.

  The search strategy is called keepElite which optimizes
  over all loci simultaniously by storing the 'eps' best
  fitted combinations during the search. The function also
  returns the optimized marginal result (each loci
  optimized).

  The covariance structures taking all loci into account
  assumes a compound symmetry structure which takes the
  number of alleles and peak heights into account (this
  ensures 'proportion of variance').

  The user may choose whether combinations giving zero
  mixture propotion (gives overfitting model) for any
  contributors are accepted.
}
\author{
  Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
}
\references{
  Tvedebrink,T, et.al.(2012). Identifying contributors of
  DNA mixtures by means of quantitative information of STR
  typing. Journal of Computational Biology, 19(7),887-902.
}
\keyword{deconvolution,}
\keyword{optimization}

