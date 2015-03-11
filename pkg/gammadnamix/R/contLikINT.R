
#' @title contLikINT
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description contLikINT marginalizes the likelihood of the STR DNA mixture given some assumed a bayesian model by numerical integration.
#' @details The procedure are doing numerical integration to approximate the marginal probability over the noisance parameters. Mixture proportions have flat prior.
#' 
#' Function calls procedure in c++ by using the package Armadillo and Boost.
#'
#' @param nC Number of contributors in model.
#' @param samples A List with samples which for each samples has locus-list elements with list elements adata and hdata. 'adata' is a qualitative (allele) data vector and 'hdata' is a quantitative (peak heights) data vector.
#' @param popFreq A list of allele frequencies for a given population.
#' @param lower Lower bounds of parameters
#' @param upper Upper bounds of parameters
#' @param refData Reference objects has locus-list element [[i]] with a list element 'r' which contains a 2 long vector with alleles for each references.
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. condOrder=-1 means the reference is known-non contributor!
#' @param knownRef Specify known references from refData (index). For instance knownRef=(1,2) means that reference 1 and 2 is known allele samples in the hypothesis. This is affected by fst-correction.
#' @param xi A numeric giving stutter-ratio if it is known. Default is NULL, meaning it is integrated out.
#' @param prC A numeric for allele drop-in probability. Default is 0.
#' @param reltol Required relative tolerance error of evaluations in integration routine. Default is 0.001.
#' @param threshT The detection threshold given. Used when considering probability of allele drop-outs.
#' @param fst is the coancestry coeffecient. Default is 0.
#' @param lambda Parameter in modeled peak height shifted exponential model. Default is 0.
#' @param pXi Prior function for xi-parameter (stutter). Flat prior on [0,1] is default.
#' @param kit shortname of kit: {"ESX17","ESI17","ESI17Fast","ESX17Fast","Y23","Identifiler","NGM","ESSPlex","ESSplexSE","NGMSElect","SGMPlus","ESX16", "Fusion","GlobalFiler"}
#' @return ret A list(margL,deviation,nEvals) where margL is Marginalized likelihood for hypothesis (model) given observed evidence, deviation is the confidence-interval of margL, nEvals is number of evaluations.
#' @export 
#' @references Hahn,T. (2005). CUBA - a library for multidimensional numerical integration. Computer Physics Communications, 168(2),78-95.
#' @keywords continuous, Bayesian models, Marginalized Likelihood estimation


contLikINT = function(nC,samples,popFreq,lower,upper,refData=NULL,condOrder=NULL,knownRef=NULL,xi=NULL,prC=0,reltol=0.01,threshT=50,fst=0,lambda=0,pXi=function(x)1,kit=NULL){
 require(cubature) 
 if(length(lower)!=length(upper)) stop("Length of integral limits differs")
 np2 <- np <- nC + 2 + sum(is.null(xi)) #number of unknown parameters
 if(length(lower)!=np) stop("Length of integral limits differs from number of unknown parameters specified")
 ret <- prepareC(nC,samples,popFreq,refData,condOrder,knownRef,kit)
 nodeg  <- is.null(kit) #boolean whether modeling degradation FALSE=YES, TRUE=NO
 if(nodeg) {
   np2 <- np2 - 1
   if(length(lower)>np2) {
    lower <- lower[-(nC+2)] #remove beta boundary
    upper <- upper[-(nC+2)] #remove beta boundary
   }

 }
 if(length(lower)!=np2) stop("The length integral limits was not the same as number of parameters!")
 liktheta <- function(theta) {   
  theta2 <- theta[1:(nC+1)] #take out mx,mu,sigma
  if(nodeg) {
    theta2 <- c(theta,1) #add beta=1 to parameters
  } else {
    theta2 <- c(theta2,theta[nC+2]) #add beta to parameters
  }
  if(is.null(xi)) {  #if xi unknown
    theta2 <- c(theta2,theta[np2]) #add xi param to parameters
  } else { #if xi known
    theta2 <- c(theta2,xi) #add xi param to parameters
  }
  Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta2),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(prC), ret$condRef,as.numeric(threshT),as.numeric(fst),ret$mkvec,ret$nkval,as.numeric(lambda),ret$bp,as.integer(0),PACKAGE="gammadnamix")[[1]]
  if(is.null(xi)) Cval <- Cval + log(pXi(theta2[np2])) #weight with prior
  return(exp(Cval)) #weight with prior of tau and stutter.
 }
 nU <- nC-ret$nK #number of unknowns
 bisectMx <- (nC==2 && nU==2) #if exact 2 unknown contributors in the hypothesis
 if(bisectMx) lower[1] <- 0.5 #restrict outcome of mixture proportions
 foo <- adaptIntegrate(liktheta, lowerLimit = lower , upperLimit = upper , tol = reltol)
 val <- foo$integral
 dev <- val + c(-1,1)*foo$error
 nEvals <- foo[[3]]
 if(nU>1) { #if more than 1 unknown 
  comb <- factorial(nU)
  val <- comb*val
  dev <- comb*dev
 }
 return(list(margL=val,deviation=dev,nEvals=nEvals))
}

