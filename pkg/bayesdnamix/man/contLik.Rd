\name{contLik}
\alias{contLik}
\title{contLik}
\usage{
  contLik(nC,mixData,popFreq,refData=NULL,condOrder=NULL,model=1,pTau=function(x) { return(1)},taumax=100, maxeval=10000,threshT=50,relaxq=0.2)
}
\arguments{
  \item{nC}{Number of contributors in model.}

  \item{mixData}{Evidence object with list elements
  adata[[i]] and hdata[[i]]. Each element has a loci-list
  with list-element 'i' storing qualitative data in 'adata'
  and quantitative data in 'hdata'.}

  \item{popFreq}{A list of allele frequencies for a given
  population.}

  \item{refData}{Reference objects with list element
  [[s]]$adata[[i]]. The list element has reference-list
  with list-element 's' having a loci-list adata with
  list-element 'i storing qualitative data.}

  \item{condOrder}{Specify conditioning references from
  refData (must be consistent order). For instance
  condOrder=(0,2,1,0) means that we restrict the model such
  that Ref2 and Ref3 are respectively conditioned as 2.
  contributor and 1. contributor in the model.}

  \item{model}{A integer for specification of model. See details for more information.}

  \item{pTau}{Prior function for tau-parameter. Flat prior
  is default.}

  \item{taumax}{Maximum range of tau-parameter. Default is
  10000.}

  \item{maxeval}{Maxumum number of evaluations in the
  interale function.}

  \item{threshT}{The analytical threshold given. Used when
  considering allele drop-outs.}

  \item{relaxq}{A parameter for model 6 which allows for relaxation of variance.}
}
\value{
  lik Marginalized likelihood of the hypothesis (model) given observed evidence.
}
\description{
  contLik evaluates the marginal likelihood of the STR DNA mixture given some assumed model by integrate out parameters.
}
\details{
 The procedure are doing numerical integration to approximate the marginal probability by integrate over noisance parameters. Mixture proportions have flat prior.

 The peak heights (r) are scaled between [0,1] such that the peak heights are not accounted for in the models.

 6 Models implemented:
 Model 1: r ~ Normal(muij,tau)
 Model 2: log(r) ~ Normal(log(muij),tau)
 Model 3: r ~ Normal(muij,tau/muij)
 Model 4: r ~ Dirichlet(tau*mui)
 Model 5: r ~ Normal+(muij,tau)
 Model 6: r ~ Normal(muij,eta/(muij^q)), q is relaxion parameter
 Model 7: r ~ Normal(muij,eta*muij + tau)


 Function calls procedure in c++ by using the package Armadillo
}
\author{
  Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
}
\keyword{continous STR mix data}
\keyword{Bayesian Models}


