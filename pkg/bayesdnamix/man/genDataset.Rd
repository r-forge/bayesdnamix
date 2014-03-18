\name{genDataset}
\alias{genDataset}
\title{genDataset}
\usage{
  genDataset(nC, popFreq, mu = 10000, sd = 1000, PrD = 0,
    sorted = FALSE)
}
\arguments{
  \item{nC}{Number of contributors in model.}

  \item{popFreq}{A list of allele frequencies for a given
  population.}

  \item{mu}{Expected peak heights for an allele of a
  contributor.}

  \item{sd}{Standard deviation of peak heights for an
  allele of a contributor.}

  \item{PrD}{Fixed Allele drop-out probability for each
  contributors.}

  \item{sorted}{Boolean for wheter sorting the contributors
  with respect to decreasingly mixture proportions.}
}
\value{
  List with elements Mx,mixData,refData
}
\description{
  genDataset samples a random mixture with Normal[+] peak
  heights.
}
\details{
  Simple function for generating a random peak height
  mixture distributed as normal positive truncated.
}
\author{
  Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
}

