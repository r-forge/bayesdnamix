\name{genDataset}
\alias{genDataset}
\title{genDataset}
\usage{
  genDataset(nC, popFreq, mu = 10000, sd = 1000)
}
\arguments{
  \item{nC}{Number of contributors in model.}

  \item{popFreq}{A list of allele frequencies for a given
  population.}

  \item{mu}{Expected peak heights for an allele of a
  contributor.}

  \item{sd}{Standard deviation of peak heights for an
  allele of a contributor.}
}
\value{
  List with elements trueMx,mixData,refData
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

