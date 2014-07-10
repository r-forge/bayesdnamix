
#' @title contLikMCMC
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description contLikMCMC simulates from the posterior distribution for a bayesian STR DNA mixture model.
#' @details The procedure are doing MCMC to approximate the marginal probability over noisance parameters. Mixture proportions have flat prior.
#' 
#' Function calls procedure in c++ by using the package Armadillo and Boost.
#'
#'
#' @param nC Number of contributors in model.
#' @param mixData Evidence object with list elements adata[[i]] and hdata[[i]]. Each element has a loci-list with list-element 'i' storing qualitative data in 'adata' and quantitative data in 'hdata'.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects with list element [[s]]$adata[[i]]. The list element has reference-list with list-element 's' having a loci-list adata with list-element 'i storing qualitative data.
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. condOrder=-1 means the reference is known-non contributor!
#' @param xi A numeric giving stutter-ratio if it is known. Default is NULL, meaning it is integrated out.
#' @param prC A numeric for allele drop-in probability. Default is 0.
#' @param method Selected MCMC-routine for calculate marginal Likelihood (1=Normal Approximation, 2=Metropolis Hastings).
#' @param pRhoTau Prior function for (rho,tau)-parameters. Flat prior is default.
#' @param rhotaumax Maximum range of rho and tau-parameter. Default is c(1000,1000).
#' @param nDone Maxumum number of random evaluations nlm-optimizing routing. Default is 1.
#' @param threshT The detection threshold given. Used when considering probability of allele drop-outs.
#' @param fst is the coancestry coeffecient. Default is 0.
#' @param lambda Parameter in modeled peak height shifted exponential model. Default is 0.
#' @param pXi Prior function for xi-parameter (stutter). Flat prior on [0,1] is default.
#' @param ximax Maximum range of xi-parameter. Default is 1.
#' @param M Number of evaluations for MCMC-sampling.
#' @param theta0 Start values for parameters in MCMC-sampling. Default is to find these using MLE.
#' @param Sigma Covariance matrix used in kernel-function in MCMC-sampling. Default is to use inverse -hessian from optimization. If theta0 is spesified but not Sigma, Sigma becomes an eye-matrix.
#' @param delta A numerical parameter to scale with Sigma. Useful for better MCMC sampling.
#' @return ret A list(margL,posttheta,postlogL,logpX,nacc) where margL is Marginalized likelihood for hypothesis (model) given observed evidence, posttheta is the posterior samples from a MC routine, postlogL is sampled log-likelihood values, nacc is number of accepted samples.
#' @keywords continuous, Bayesian models
#' @examples
#' \dontrun{
#' #data model:
#' threshT <- 50 #threshold
#' load(system.file("mcData.Rdata", package = "gammadnamix"))
#' popFreq <- data$popFreq
#' mixData <- data$samples$MC15
#' refData <- data$refData
#' #Q-assignation
#' popFreq2 <- popFreq
#' for(i in 1:length(popFreq)) { #make Q-assignation for each loci
#'  tmp <- popFreq[[i]][names(popFreq[[i]])%in%mixData$adata[[i]]]
#'  tmp <- c(tmp,1-sum(tmp))
#'  names(tmp)[length(tmp)] <- "99"
#'  popFreq2[[i]] <- tmp
#' }
#' #Hypothesis:
#' nC <- 2 #number of contributors
#' condhp <- condhd <- rep(0,length(refData)) 
#' condhp [1] <- 1 #condition on first one in refdata
#' #helpfunction for checking samples:
#' validMCMC <- function(X) {
#'  p <- ncol(X)
#'  par(mfrow=c(p,2))
#'  for(i in 1:p) {
#'   plot(density(X[,i],from=0),ty="l")
#'   plot(X[,i],ty="l")
#'  }
#'  par(mfrow=c(1,1))
#' }
#' M <- 1e5 #number of samples from posterior
#' delta = 0.3 #tweak parameter in MCMC
#' #without stutter:
#' hpD4 <- contLikMCMC(nC,mixData,popFreq2,refData=refData,condOrder=condhp,xi=0,threshT=threshT,M=M,delta=delta)
#' hdD4 <- contLikMCMC(nC,mixData,popFreq2,refData=refData,condOrder=condhd,xi=0,threshT=threshT,M=M,delta=delta)
#' validMCMC(hpD4$posttheta)
#' validMCMC(hpD4$posttheta)
#' LRD4 <- hpD4$margL/hdD4$margL #estimated LR
#' #with stutter:
#' hpS4 <- contLikMCMC(nC,mixData,popFreq2,refData=refData,condOrder=condhp,threshT=threshT,M=M,ximax=1,delta=delta)
#' hdS4 <- contLikMCMC(nC,mixData,popFreq2,refData=refData,condOrder=condhd,threshT=threshT,M=M,ximax=1,delta=delta)
#' validMCMC(hpS4$posttheta)
#' validMCMC(hdS4$posttheta)
#' LRS4 <- hpS4$margL/hdS4$margL #estimated LR
#' }

contLikMCMC = function(nC,mixData,popFreq,refData=NULL,condOrder=NULL,xi=NULL,prC=0,method=2,pRhoTau=function(x) {1},rhotaumax=c(1000,1000),nDone=1,threshT=50,fst=0,lambda=0,pXi=function(x) {1},ximax=1,M=1e4,theta0=NULL,Sigma=NULL,delta=0.3){
 #Optimize with MLE if theta0,Sigma not given
 if(is.null(theta0)) {
  mle <- contLikMLE(nC,mixData,popFreq,refData=refData,condOrder=condOrder,xi=xi,prC=prC,threshT=threshT,fst=fst,lambda=lambda,rhotaumax=rhotaumax,ximax=ximax,nDone=1)
  theta0 <- mle$MLE
  Sigma <- mle$Sigma
 }
 if(!is.null(theta0) && is.null(Sigma)) Sigma <- diag( rep(1,length(theta0)) ) #default Sigma inserted
 if(!all(length(theta0)%in%dim(Sigma))) stop("Length of theta0 and dimension of Sigma was not the same!")

 nA = unlist(lapply(mixData$adata,length)) #number of alleles of selected loci
 nL <- length(nA) #number of loci in data
 locnames <- toupper(names(mixData$adata))
 missind <- !locnames%in%toupper(names(popFreq))  #indice of missing loci
 if(any(missind)) {
  stop(paste('Missing locus (',locnames[missind],') in popFreq.',sep=''))
 }

 #convertion of values in popFreq, mixData and Glist$G:
 #loci-order follows as in mixData: "locnames". Rearrange names:
 names(popFreq) <- toupper(names((popFreq))) #make invariant to capital
 popFreq <- popFreq[locnames] #order popFreq to mixData-order

 #get population genotypes:
 getGlist <- function(popFreq) {
  locs <- names(popFreq)
  Glist <- list()
  for(i in 1:length(locs)) {
   G = t(as.matrix(expand.grid(rep(list(as.numeric(names(popFreq[[i]])),as.numeric(names(popFreq[[i]])) ))))) #one genotype per column
   keep = G[2,]>=G[1,] #unique genotypes 
   G <- G[,keep]  #store genotypes
   G <- matrix(as.character(G),nrow=2) #make string names again
   tmpP = t(as.matrix(expand.grid(rep(list(as.numeric(popFreq[[i]]),as.numeric(popFreq[[i]]) )))))
   Gprob = exp(colSums(log(tmpP[,keep]))) #get allele probs
   ishet = G[1,]!=G[2,]
   Gprob[ishet] = 2*Gprob[ishet] #multiply with two to get heterozygote prob
   Glist[[locs[i]]] <- list(G=t(G),Gprob=Gprob)
  }
  return(Glist)
 }

 #Fix references: Assign condition to condM-matrix
 Glist <- getGlist(popFreq) #get population genotype information
 condM <- matrix(-1,nrow=nL,ncol=nC) #default is no references (=-1)
 #assign references to condM-matrix by values of Glist
 if(!is.null(refData) && !is.null(condOrder) && any(condOrder>0)) {
  for(i in 1:nL) {
   subRef <- refData[[locnames[i]]]
   Gset <- Glist[[locnames[i]]]$G #genotype combinations
   if(length(subRef)==0)  stop(paste('Missing locus (',locnames[i],') in refData.',sep=''))
   for(k in 1:length(subRef)) { #for each reference
    if(condOrder[k]>0) {
     Gind1 <- subRef[[k]][1]==Gset[,1] & subRef[[k]][2]==Gset[,2]
     Gind2 <- subRef[[k]][2]==Gset[,1] & subRef[[k]][1]==Gset[,2]
     condM[i,condOrder[k]] = which(Gind1 | Gind2) - 1 #subtract with one since we work from 0-indice
    }
   }
  }
 }

 #fix known sample-information:
 mkvec <- numeric()
 nkval <- rep(0,nL)
 for(i in 1:nL) {
  tmp <- rep(0, length(popFreq[[i]]))
  if(!is.null(condOrder) & !is.null(refData)) {
   for(k in 1:length(condOrder)) {
    if(condOrder[k]!=0) { 
     ind <- which( names(popFreq[[i]])%in%refData[[locnames[i]]][[k]] )
     tmp[ind] = (length(ind)==1) + 1
    }
   }
  }
  nkval[i] <- sum(tmp) #number of sampled (for each loci)
  mkvec <- c(mkvec,tmp) 
 }

 #decode allele-names to index names + Stutter-preparation
 allAbpind <- as.numeric()
 for(i in 1:nL) {
   anames <- names(popFreq[[i]]) #old names
   #find corresponding index of 1-ahead pb-twin
   numAnames <- as.numeric(anames)
   BPind <- rep(0,length(numAnames)) #init 0 if allele doesn't have any bp-twin
   for(j in 1:length(numAnames)) {
    bpind <- which((numAnames[j]-1)==numAnames) #store index of what allele it stutters to
    if(length(bpind)>0) BPind[j] <- bpind #assign index-placement
   }
   allAbpind <- c(allAbpind,BPind)
   anames2 <- 0:(length(popFreq[[i]])-1) #new names
   names(popFreq[[i]]) <- anames2
   for(j in 1:nA[i]) { #for each allele
    mixData$adata[[i]][j] <- anames2[anames==mixData$adata[[i]][j]] #update name
   }
  }

 #fix genotypes:
 Glist <- getGlist(popFreq) #get population genotype information
 Gprob <- lapply(Glist,function(x) return(x$Gprob))
 Gset <- lapply(Glist,function(x) return(x$G))


 #Fix input-variables for Cpp-function
 logPE <- 1
 CnA <- c(0,cumsum(nA))
 allA <- as.integer(unlist(mixData$adata))
 allY <- as.numeric(unlist(mixData$hdata))
 sY <- sapply(mixData$hdata,sum)
 nG <- sapply(Gprob,length) #number of genotype combinations
 CnG <- c(0,cumsum(nG))
 CnG2 <- c(0,cumsum(nG*2)) #note: 2 columns for each genotype!!
 pG <- unlist(Gprob) #vectorize over all loci
 Gvec <- as.integer(rbind(unlist(Gset))) #vectorize a big matrix (loci are put chronologic)
 condRef <- c(condM) #vectorized over all loci
 nAall <- sapply(popFreq,length) #Number of population-alleles on each loci
 CnAall <- c(0,cumsum(nAall)) #cumulative number of alleles
 pA <- unlist(popFreq) #need each allele probability for drop-in probabilities

 #Two cases: Stutter unknown or Stutter known
 loglikYtheta <- function(theta) {   #call c++- function: length(theta)=nC+1
  Cval  <- .C("loglikgammaC",as.numeric(logPE),as.numeric(theta),as.integer(nC),as.integer(nL),as.integer(nA), as.numeric(allY),as.integer(allA),as.integer(CnA),as.integer(allAbpind),as.integer(nAall),as.integer(CnAall),as.integer(Gvec),as.integer(nG),as.integer(CnG),as.integer(CnG2),as.numeric(pG),as.numeric(pA), as.numeric(prC), as.integer(condRef),as.numeric(threshT),as.numeric(fst),as.integer(mkvec),as.integer(nkval),as.numeric(lambda),PACKAGE="gammadnamix")[[1]]
  loglik <- Cval + log(pRhoTau(theta[c(nC,nC+1)])) + log(pXi(theta[nC+2])) #weight with prior of tau and 
  return(loglik) #weight with prior of tau and stutter.
 } 
 loglikYthetaS <- function(theta2) {   #call c++- function: length(theta)=nC
  theta <- c(theta2,xi) #stutter-parameter added as known
  Cval  <- .C("loglikgammaC",as.numeric(logPE),as.numeric(theta),as.integer(nC),as.integer(nL),as.integer(nA), as.numeric(allY),as.integer(allA),as.integer(CnA),as.integer(allAbpind),as.integer(nAall),as.integer(CnAall),as.integer(Gvec),as.integer(nG),as.integer(CnG),as.integer(CnG2),as.numeric(pG),as.numeric(pA), as.numeric(prC), as.integer(condRef),as.numeric(threshT),as.numeric(fst),as.integer(mkvec),as.integer(nkval),as.numeric(lambda),PACKAGE="gammadnamix")[[1]]
  loglik <- Cval + log(pRhoTau(theta[c(nC,nC+1)])) + log(pXi(xi)) #weight with prior of tau and stutter.
  return(loglik)
 }
 logdmvnorm <- function(X,mean,cholC) { #function taken from mvtnorm-package
   p <- nrow(cholC)
   tmp <- backsolve(cholC,t(X)-mean,p,transpose=TRUE)
   rss <- colSums(tmp^2)
   logretval <- -sum(log(diag(cholC))) - 0.5*p*log(2*pi) - 0.5*rss
   return(logretval)
 }
 p <- length(theta0) #dimension
 if(method==1) C <- chol(delta*Sigma) #scale variance with a factor 2: ensures broad posterior
 if(method==2) C <- chol(delta*diag(diag(Sigma))) #scale variance with a factor
 X <- t( t(C)%*%matrix(rnorm(p*M),ncol=M,nrow=p)) #proposal values

 #get boundary of parameters
 if(is.null(xi)) { #if stutter unknown
  lower <- rep(0,nC+2)
  upper = c(rep(1,nC-1),rhotaumax,1)
 } else {
  lower <- rep(0,nC+1)
  upper = c(rep(1,nC-1),rhotaumax)
 }

 #v1: Importance sampling using Normal(theta0,delta*Sigma)
 if(method==1) {
  X <- t( t(X)+ theta0 ) #lazy-bayes simulation
  logpX <- logdmvnorm(X=X,mean=theta0,cholC=C) #get logged distr-values
  indrm <- numeric()
  for(i in 1:p) indrm <- c(indrm, which( X[,i]<lower[i] | X[,i]>upper[i] ) ) #throw away non-acceptable
  if(nC>2) { #check restrictions if nC>2
   indrm <- c(indrm, which( rowSums(X[,1:(nC-1)])>1 ) ) 
  }
  if(length(indrm)>0) {
   X <- X[-indrm,]  #remove proposed
   logpX <- logpX[-indrm] 
  }
  #go through each proposal
  M <- nacc <- nrow(X) #update counter
  postlogL <- rep(NA,M)
  if(is.null(xi)) {
   for(m in 1:M) postlogL[m] <- loglikYtheta(X[m,])
  } else {
   for(m in 1:M) postlogL[m] <- loglikYthetaS(X[m,])
  }
  bigsum <- sum( exp( postlogL - logpX ) )
  margL <- bigsum/M #estimated marginal likelihood
  posttheta <- X #lazy-bayes samples
 #v2: MCMC by Gelfand and Dey (1994), using h() = Normal(theta0,delta*Sigma)
 } else if(method==2) {
  U <- runif(M) #random numbers
  nacc <- 0
  posttheta <- matrix(NA,ncol=p,nrow=M) #accepted theta
  postlogL <- rep(NA,M) #accepted log-likelihoods
  posttheta[1,] <- theta0
  if(is.null(xi)) { #if stutter model
   postlogL[1] <- loglikYtheta(theta0) #get start-likelihood 
   for(m in 2:M) { #for each proposal
    posttheta[m,] <-  X[m,]+posttheta[m-1,] #proposed theta
    if(any(posttheta[m,]<lower | posttheta[m,]>upper)) {
     postlogL[m] <- -Inf
    } else if(nC>2 && sum(posttheta[m,1:(nC-1)])>1)  {
     postlogL[m] <- -Inf
    } else {
     postlogL[m] <- loglikYtheta(posttheta[m,]) #get new-likelihood 
    }
    pr <- exp(postlogL[m]- postlogL[m-1]) #acceptance rate
    if(U[m]>pr) { #if not accepted, i.e. random pr too large
     posttheta[m,] <-  posttheta[m-1,]
     postlogL[m] <- postlogL[m-1 ]
    } else {
     nacc <- nacc + 1
    }
   } #end for each sample
  } else { #if stutter ratio assumed known known
   postlogL[1] <- loglikYthetaS(theta0) #get start-likelihood 
   for(m in 2:M) { #for each proposal
    posttheta[m,] <-  X[m,]+posttheta[m-1,] #proposed theta
    if(any(posttheta[m,]<lower | posttheta[m,]>upper)) {
     postlogL[m] <- -Inf
    } else if(nC>2 && sum(posttheta[m,1:(nC-1)])>1)  {
     postlogL[m] <- -Inf
    } else {
     postlogL[m] <- loglikYthetaS(posttheta[m,]) #get new-likelihood 
    }
    pr <- exp(postlogL[m]- postlogL[m-1]) #acceptance rate
    if(U[m]>pr) { #if not accepted, i.e. random pr too large
     posttheta[m,] <-  posttheta[m-1,]
     postlogL[m] <- postlogL[m-1 ]
    } else {
     nacc <- nacc + 1
    }
   } #end for each sample
  }
  logpX <- logdmvnorm(posttheta,mean=theta0,cholC=chol(Sigma)) #insert with Normal-approx of post-theta
  #plot(postlogL,ty="l")
  #plot(logpX,ty="l")
  margL <- 1/mean(exp(logpX - postlogL)) #estimated marginal likelihood
 }
 return(list(margL=margL,posttheta=posttheta,postlogL=postlogL,logpX=logpX,nacc=nacc))
} #end function

