#setwd("C:/Users/oebl/Dropbox/Forensic/MixtureProj/myArticles/stutterLR/script")
rm(list=ls())
library(gammadnamix)

#load data:
load(system.file("mcData.Rdata", package = "gammadnamix"))
samples <- list(MC15=mcdata$samples$MC15) #use only sample MC15

retlist <- Qassignate(samples,mcdata$popFreq,mcdata$refData) #get Q-assignated allele frequiences
popFreqQ <- retlist$popFreq
refDataQ <- retlist$refData

#Model:
threshT = 50 #peak height threshold
condHp <- c(1,2,3) #Conditioned references under hp
condHd <- c(1,2,0) #Conditioned references under hd
nC <- 4 #number of contributors

#NON Q-assignated data:
set.seed(1)
nDone <- 1 #number of random optimization start points
hptime <- system.time( {  hpSmle <- contLikMLE(nC,samples,mcdata$popFreq,mcdata$refData,condOrder=condHp,threshT=threshT,nDone=nDone) } )[3]
print(hpSmle$fit$thetahat) #MLE under hp
print(sqrt(diag(hpSmle$fit$thetaSigma2))) #standard error
set.seed(1)
hdtime <- system.time( {  hdSmle <- contLikMLE(nC,samples,mcdata$popFreq,mcdata$refData,condOrder=condHd,threshT=threshT,nDone=nDone) } )[3] #slow for replicate version
print(hdSmle$fit$thetahat) #MLE under hd
print(sqrt(diag(hdSmle$fit$thetaSigma2))) #standard error
print(c(hptime,hdtime)) #time usage
LRmle <- exp(hpSmle$fit$loglik-hdSmle$fit$loglik) #MLE optimized 
LRlap <- exp(hpSmle$fit$logmargL-hdSmle$fit$logmargL) #Laplace approximated
print(log10(LRmle))
print(log10(LRlap))


#For Q-assignated data: (much faster than using all data)
set.seed(1)
nDone <- 3 #number of random optimization start points
hptime <- system.time( {  hpSmleQ <- contLikMLE(nC,samples,popFreqQ,refDataQ,condOrder=condHp,threshT=threshT,nDone=nDone) } )[3] #-117.9732
print(hpSmleQ$fit$thetahat) #MLE
print(sqrt(diag(hpSmleQ$fit$thetaSigma2))) #standard error

hdtime <- system.time( {  hdSmleQ <- contLikMLE(nC,samples,popFreqQ,refDataQ,condOrder=condHd,threshT=threshT,nDone=nDone) } )[3] #-129.3092
print(c(hptime,hdtime))
print(hdSmleQ$fit$thetahat) #MLE
print(sqrt(diag(hdSmleQ$fit$thetaSigma2))) #standard error

LRmleQ <- exp(hpSmleQ$fit$loglik-hdSmleQ$fit$loglik) #MLE optimized 
LRlapQ <- exp(hpSmleQ$fit$logmargL-hdSmleQ$fit$logmargL) #Laplace approximated
print(log10(LRmleQ))
print(log10(LRlapQ))



#Integrated likelihood
reltol <- 0.05 #relative error
low <- rep(0,nC+2) #lower boundary of integral
up <- rep(1,nC+2) #upper boundary of integral
up[nC] <- 10000
up[nC+1] <- 1
hptime <- system.time( {  hpSintQ <- contLikINT(nC,samples,popFreqQ,low,up,refDataQ,condOrder=condHp,threshT=threshT,reltol=reltol) } )[3]
hdtime <- system.time( {  hdSintQ <- contLikINT(nC,samples,popFreqQ,low,up,refDataQ,condOrder=condHd,threshT=threshT,reltol=reltol) } )[3]
LRintQ <- hpSintQ$margL/hdSintQ$margL #MLE optimized 
print(log10(LRintQ))




