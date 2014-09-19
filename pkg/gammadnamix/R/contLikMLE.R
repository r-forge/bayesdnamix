
#' @title contLikMLE
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description contLikMLE optimizes the likelihood of the STR DNA mixture given some assumed a bayesian model.
#' @details The procedure are doing numerical optimization to approximate the marginal probabilit over noisance parameters. Mixture proportions have flat prior.
#' 
#' The procedure also does a Laplace Approximation of the marginalized likelihood (theta integrated out) and returns the log-marginal likelihood as logmargL in the fit-list.
#' 
#' Function calls procedure in c++ by using the package Armadillo and Boost.
#'
#' @param nC Number of contributors in model.
#' @param samples A List with samples which for each samples has locus-list elements with list elements adata and hdata. 'adata' is a qualitative (allele) data vector and 'hdata' is a quantitative (peak heights) data vector.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects with list element [[s]]$adata[[i]]. The list element has reference-list with list-element 's' having a loci-list adata with list-element 'i storing qualitative data.
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. condOrder=-1 means the reference is known-non contributor!
#' @param knownRef Specify known references from refData (index). For instance knownRef=(1,2) means that reference 1 and 2 is known allele samples in the hypothesis. This is affected by fst-correction.
#' @param xi A numeric giving stutter-ratio if it is known. Default is NULL, meaning it is integrated out.
#' @param prC A numeric for allele drop-in probability. Default is 0.
#' @param nDone Maxumum number of random evaluations nlm-optimizing routing. Default is 1.
#' @param threshT The detection threshold given. Used when considering probability of allele drop-outs.
#' @param fst is the coancestry coeffecient. Default is 0.
#' @param lambda Parameter in modeled peak height shifted exponential model. Default is 0.
#' @param pXi Prior function for xi-parameter (stutter). Flat prior on [0,1] is default.
#' @param alpha Significant level used for confidence interval for parameter estimation. Default is 0.05.
#' @param delta Standard deviation of normal distribution when drawing random startpoints. Default is 10.
#' @return ret A list(fit,model,nDone,delta) where fit is Maximixed likelihood elements for given model.
#' @export
#' @references Cowell,R.G. et.al. (2014). Analysis of forensic DNA mixtures with artefacts. Applied Statistics, 64(1),1-32.
#' @keywords continuous model, Maximum Likelihood Estimation
contLikMLE = function(nC,samples,popFreq,refData=NULL,condOrder=NULL,knownRef=NULL,xi=NULL,prC=0,nDone=1,threshT=50,fst=0,lambda=0,pXi=function(x)1,delta=10,alpha=0.05){
 ret <- prepareC(nC,samples,popFreq,refData,condOrder,knownRef)

 if(is.null(xi)) {
  negloglikYphi <- function(phi) {   #call c++- function: length(theta)=nC+1
   Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(phi),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(prC), ret$condRef,as.numeric(threshT),as.numeric(fst),ret$mkvec,ret$nkval,as.numeric(lambda),as.integer(1),PACKAGE="gammadnamix")[[1]]
   loglik <- Cval + log(pXi(1/(1+exp(-theta[nC+2])))) #weight with prior of tau and 
   return(-loglik) #weight with prior of tau and stutter.
  }
 } else {  
  negloglikYphi <- function(phi2) {   #call c++- function: length(theta)=nC
   phi<- c(phi2,xi) #stutter-parameter added as known
   Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(phi),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(prC), ret$condRef,as.numeric(threshT),as.numeric(fst),ret$mkvec,ret$nkval,as.numeric(lambda),as.integer(1),PACKAGE="gammadnamix")[[1]]
   return(-Cval)
  }
 }

 np <- nC + 1 + sum(is.null(xi)) #number of unknown parameters
 maxL <- -Inf #
 nOK <- 0 #number of times for reaching optimum
 suppressWarnings({
  while(nOK<nDone) {
    p0 <- rnorm(np,sd=delta) #generate random start value
    if(is.infinite(negloglikYphi(p0))) next #skip to next if still INF
    tryCatch( {
       foo <- nlm(f=negloglikYphi, p=p0,hessian=TRUE)
       Sigma <- solve(foo$hessian)
       if(any(diag(Sigma)<=0)) next; #not a local optimum
       if(foo$code==1 && foo$iterations>2) {
        likval <- -foo$min
        nOK=nOK+1 #it was accepted as an optimum
        if(likval>maxL) {
         maxL <- likval #maximized likelihood
         maxPhi <- foo$est #set as topfoo     
         maxSigma <- Sigma 
        }
       }
    },error=function(e) e) #end trycatch 
  } #end while loop
 })
 #transfer back: 
 mx <- numeric()
 if(nC>1) {
  mx <- 1/(1+exp(-maxPhi[1:(nC-1)]))
  if(nC>2) { #need to transfer back
   for(i in 2:(nC-1)) {
    mx[i] <- mx[i]*(1-sum(mx[1:(i-1)]))
   }
  }
 }
 musigma <- exp(maxPhi[nC:(nC+1)]) #inverse-log
 thetahat <- c(mx,musigma) #last index is removed. This could again be a known contributor
 thetahat2 <- c(mx,1-sum(mx),musigma) #last index is removed. This could again be a known contributor
 if(is.null(xi)) {
  tmp <- 1/(1+exp(-maxPhi[nC+2]))
  thetahat <- c(thetahat,tmp) #add xi to parameters
  thetahat2 <- c(thetahat2,tmp) #add xi to parameters
 }

 #Delta-method to Sigma matrix
 Jacob <- function(phi,mle) { #Jabobian matrix (in value phi)
  J <- matrix(0,length(phi),length(phi))
  if(nC>1) {
   DmDm <- matrix(0,nC-1,nC-1)
   for(i in 1:(nC-1)) {
    mitmp <- 1/(1+exp(-phi[i])) #temporary variable
    for(j in 1:i) {
     if(j==i) {
       DmDm[i,i] <- exp(-phi[i])*mitmp^2
       if(i>1) DmDm[i,i] <- DmDm[i,i]*(1-sum(mle[1:(j-1)])) #note using mle(theta) here!
     } else {
       DmDm[i,j] <- -mitmp*(sum( DmDm[1:i,j] ))
     } #end case
    } #end for each col j (xj)
   } #end for each row i (fi)
  J[1:(nC-1),1:(nC-1)] <- DmDm
  }
  for(i in nC:(nC+1)) J[i,i] <- exp(phi[i])
  if(is.null(xi)) {
   tmp <- exp(-phi[nC+2])
   J[nC+2,nC+2] <- tmp*(1+tmp)^(-2)
  }
  return(J)
 } #end jacobian
 J <- Jacob(maxPhi,mle=thetahat)
 Sigma <- (t(J)%*%maxSigma%*%J) #this is correct covariance of thetahat. Observed hessian is used
 #Sigma <- (J%*%maxSigma%*%t(J))/sqrt(np) #this formula is wrong!

 #get extended Sigma (all parameters)
 Sigma2 <- matrix(NA,nrow=np+1,ncol=np+1) #extended covariance matrix also including mx[nC]
 Sigma2[nC:np+1,nC:np+1] <- Sigma[nC:np,nC:np] 
 Sigma2[nC:np+1,nC:np+1] <- Sigma[nC:np,nC:np] 
 if(nC>1) {
  Sigma2[nC:np+1,1:(nC-1)] <- Sigma[nC:np,1:(nC-1)] 
  Sigma2[1:(nC-1),1:(nC-1)] <- Sigma[1:(nC-1),1:(nC-1)] 
  Sigma2[1:(nC-1),nC:np+1] <- Sigma[1:(nC-1),nC:np] 
  Sigma2[nC,nC] <- sum(Sigma[1:(nC-1),1:(nC-1)])
  for(k in (1:(np+1))[-nC]) {
   Sigma2[nC,k] <- Sigma2[k,nC] <- -sum(Sigma[1:(nC-1),k-sum(k>nC)]) 
  }
 } else {
  Sigma2[1,1:(np+1)] <- Sigma2[1:(np+1),1] <- 0 #no uncertainty
 }
 #CI for theta:
 thetaSD <- sqrt(diag(Sigma2))
 thetaCI <- thetahat2 + cbind(qnorm(alpha/2)*thetaSD,0,qnorm(1-alpha/2)*thetaSD)

 phinames <- c("log(mu)","log(sigma)")
 thetanames <- c("mu","sigma")
 if(nC>1) {
  phinames  <- c(paste0("nu",1:(nC-1)),phinames)
  thetanames <- c(paste0("mx",1:(nC-1)),thetanames)
 }
 thetanames2 <- c(paste0("mx",1:nC),"mu","sigma")
 if(np==(nC+2)) {
  phinames <- c(phinames,"logit(xi)")
  thetanames <- c(thetanames,"xi") 
  thetanames2 <- c(thetanames2,"xi") 
 } 
 colnames(thetaCI) <- c(paste0(alpha/2*100,"%"),"MLE",paste0((1-alpha/2)*100,"%"))
 rownames(thetaCI) <- thetanames2
 colnames(maxSigma) <- rownames(maxSigma) <- phinames 
 colnames(Sigma) <- rownames(Sigma) <- thetanames
 colnames(Sigma2) <- rownames(Sigma2) <- thetanames2
 names(maxPhi) <- phinames
 names(thetahat) <- thetanames
 names(thetahat2) <- thetanames2

 #laplace approx:
 logmargL <- 0.5*(np*log(2*pi)+determinant(Sigma)$mod[1]) + maxL #get log-marginalized likelihood
 nU <- nC-ret$nK #number of unknowns
 if(nU>1) { #if more than 1 unknown 
  logmargL <- log(factorial(nU)) + logmargL #get correct ML adjusting for symmetry
 }
 fit <- list(phihat=maxPhi,thetahat=thetahat,thetahat2=thetahat2,phiSigma=maxSigma,thetaSigma=Sigma,thetaSigma2=Sigma2,loglik=maxL,thetaCI=thetaCI,logmargL=logmargL)
 #store model:
 model <- list(nC=nC,samples=samples,popFreq=popFreq,refData=refData,condOrder=condOrder,knownRef=knownRef,xi=xi,prC=prC,threshT=threshT,fst=fst,lambda=lambda,pXi=pXi)
 ret <- list(fit=fit,model=model,nDone=nDone,delta=delta)
 return(ret)
} #end function

