% Generated by roxygen2 (4.0.1): do not edit by hand
\name{contLikStutter}
\alias{contLikStutter}
\title{contLikStutter}
\usage{
contLikStutter(nC, mixData, popFreq, refData = NULL, condOrder = NULL,
  xi = NULL, prC = 0, model = 2, pTau = function(x) {     return(1) },
  taumax = 100, maxeval = 5000, threshT = 50, fst = 0, lambda = 0,
  pXi = function(x) {     return(1) })
}
\arguments{
\item{nC}{Number of contributors in model.}

\item{mixData}{Evidence object with list elements adata[[i]] and hdata[[i]]. Each element has a loci-list with list-element 'i' storing qualitative data in 'adata' and quantitative data in 'hdata'.}

\item{popFreq}{A list of allele frequencies for a given population.}

\item{refData}{Reference objects with list element [[s]]$adata[[i]]. The list element has reference-list with list-element 's' having a loci-list adata with list-element 'i storing qualitative data.}

\item{condOrder}{Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. condOrder=-1 means the reference is known-non contributor!}

\item{xi}{A numeric giving stutter-ratio if it is known. Default is NULL, meaning it is integrated out.}

\item{prC}{A numeric for allele drop-in probability. Default is 0.}

\item{model}{A integer for specification of model. See details for more information.}

\item{pTau}{Prior function for tau-parameter. Flat prior is default.}

\item{taumax}{Maximum range of tau-parameter. Default is 1000.}

\item{maxeval}{Maxumum number of evaluations in the interale function. Default is 5000.}

\item{threshT}{The detection threshold given. Used when considering probability of allele drop-outs.}

\item{fst}{is the coancestry coeffecient. Default is 0.}

\item{lambda}{Parameter in modeled peak height shifted exponential model. Default is 0.}

\item{pXi}{Prior function for xi-parameter (stutter). Flat prior on [0,1] is default.}
}
\value{
lik Marginalized likelihood of the hypothesis (model) given observed evidence.
}
\description{
contLikStutter evaluates the marginal likelihood of the STR DNA mixture given some assumed a bayesian model by integrate out parameters.
}
\details{
The procedure are doing numerical integration to approximate the marginal probability by integrate over noisance parameters. Mixture proportions have flat prior.

The user may specify probability of drop-out for each contributors.

Model 1 is gaussian model: yj~N(sum(y)/2*nj*m,sum(y)*tau). Inspired by Tvedebrink.
Model 2 is gamma model: yj~N(sum(y)/(2*tau)*nj*m,tau). Inspired by Cowell.

Function calls procedure in c++ by using the package Armadillo and Boost.
}
\author{
Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
}
\keyword{Bayesian}
\keyword{continuous,}
\keyword{models}

