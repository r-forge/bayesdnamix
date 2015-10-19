\name{dnamatch-package}
\alias{dnamatch-package}
\alias{dnamatch}
\docType{package}
\title{
A procedure for large scale contamination search
}
\description{
A fast contamination search algorithm to search between trace samples (mixtures) and between trace samples and reference samples where the trace samples are selected within a range of time. The match statistic uses likelihood ratio (LR) calculation based on allele frequencies to create decreasing sorted table of the most likely matches. The user can select threshold for LR to assign matches and change the model to take into account drop-in/out.
}
\details{
\tabular{ll}{
Package: \tab dnamatch\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-10-19\cr
License: \tab   GPL(>=2) \cr
}
dnamatch(fn="testprover",freqfile="ESX17_Norway.csv",reffold="ansattprofiler",TAptrn="TA",SIDpat="-S",BIDpat="_BES")
}
\author{
Oyvind Bleka <Oyvind.Bleka@fhi.no>
}
