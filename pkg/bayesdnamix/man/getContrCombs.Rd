\name{getContrCombs}
\alias{getContrCombs}
\title{getContrCombs}
\usage{
  getContrCombs(Alist, nC, symmetry = FALSE, refs = NULL,
    condOrder = NULL, complete = TRUE)
}
\arguments{
  \item{Alist}{Vector of possible alleles}

  \item{nC}{Number of contributors in model.}

  \item{symmetry}{Boolean of whether all combinations
  should be considered (permutation/combination)}

  \item{refs}{List of contributors [[s]] where list-element
  's' is the reference index with genotype vector
  ("A1","A2"). #'}

  \item{condOrder}{Specify conditioning references in refs
  (must be consistent order). For instance
  condOrder=(0,2,1,0) means that we restrict the model such
  that Ref2 and Ref3 are respectively conditioned as 2.
  contributor and 1. contributor in the model.}

  \item{complete}{Boolean for requiring the set of
  combinations to fill the Alist.}
}
\value{
  A matrix where each row represents a genotype
  combination, while each column correspond to a
  contributor
}
\description{
  getContrCombs is a helpfunction to specify all possible
  combined genotypes for a give hypothesis.
}
\details{
  Note: If a reference is unselected, it as numeric() in
  its loci

  The allele-names are coded with a ascii-symbol to have
  threat more than 10-alleles in a loci.
}
\author{
  Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
}

