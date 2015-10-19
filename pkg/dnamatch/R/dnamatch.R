###########
#Changelog#
###########

#' @title dnamatch
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description dnamatch is a function for doing contamination search between trace samples (mixtures) and between trace samples and reference samples where the trace samples are selected within a range of time. 
#' @details dnamatch imports case-stains from genemapper-format automatically and feeds it into a structure for doing comparison matching.
#' 
#' Note: Samplenames must have this format: "SID_BID_CID". The samplename must be unique and consist of SID=Sample ID, BID=Unique ID for all samples, CID=Case ID.
#' 
#' Timestamp="YY-MM-DD-HH-MM-SS" is used as identification key of when the dnamatch function was run.
#' Match results are stored as a table in matchfile.csv with the timestamp. Matches which are earlier found in matchfile.csv will not be stored again. The column "Checked" can be used for comments.
#' More details about a given dnamatch run are stored in the folder "session" with corresponding timestamp.
#'
#' @param fn A folder with stain files. Full directory must be given. 
#' @param freqfile A file containing population frequencies for alleles. Full directory must be given.
#' @param reffold A folder with stored personal references. Default is no references. Full directory must be given.
#' @param TAptrn Filename structure of samples to evaluate. For instance TAptrn="TA-".
#' @param SIDpat Pattern of name which is sample ID (SID_BID_CID)=("-S0001_BES00001-14_2014234231"). Here SIDpat="-S". Can also be a vector.
#' @param BIDpat Pattern of name which is unique for all samples to evaluate. (SID_BID_CID)=("-S0001_BES00001-14_2014234231"). Here BIDpat="_BES". 
#' @param Thist Number of days back in time for file to be imported(stains). Samples older than Tnew are declared as OLD.
#' @param Tnew Number of days back in time which are defined as NEW files. If today is 2014-11-16, Tnew=3 means that files from 2014-11-13 are declared as NEW.
#' @param threshLR Threshold for extracting a match to file.
#' @param threshHeight Acceptable peak height (rfu) in stains.
#' @param threshStutt Acceptable stutter ratio in stains (relative peak heights) .
#' @param threshMaj If second largest allele has ratio (relative to the largest allele) above this threshold, then second allele is part of the major profile (used for extracting major from mixture). If the relative peak height between second and third largest allele has ratio greater than this threshold, no major is assigned.
#' @param pD Assumed drop-out rate per marker (parameter in LR model).
#' @param pC Assumed drop-in rate per marker (parameter in LR model).
#' @param sameKL Boolean whether matches within same case ID number are shown.
#' @param N Database size used for estimating population frequencies. Used to assign new allele frequenceis.
#' @export

#rm(list=ls())
#setwd("C:/Users/oebl/Dropbox/Forensic/MixtureProj/myDev/")
#library(roxygen2);roxygenize("dnamatch")
#library(dnamatch) #source("dnamatch.R")
#dnamatch(fn="testprover",freqfile="ESX17_Norway.csv",reffold="ansattprofiler",TAptrn="TA",SIDpat="-S",BIDpat="_BES")

dnamatch = function(fn,freqfile,reffold=NULL,TAptrn,SIDpat,BIDpat,Thist=1,Tnew=1,threshLR=1e3,threshHeight=10,threshStutt=0.1,threshMaj=0.6,pD=0,pC=0,sameKL=FALSE,N=5000) {
 require(forensim) #used in tippet-evaluations
 newf0 <- 5/(2*N)#new allele-frequence for new alleles

 ###############
 #HELPFUNCTIONS#
 ###############
 tableReader=function(filename) {   #Robust function for reading tables:
  tab <- read.table(filename,header=TRUE,sep="\t",stringsAsFactors=FALSE)
  tryCatch( {  if(ncol(tab)==1) tab <- read.table(filename,header=TRUE,sep=",",stringsAsFactors=FALSE) } ,error=function(e) e) 
  tryCatch( {  if(ncol(tab)==1) tab <- read.table(filename,header=TRUE,sep=";",stringsAsFactors=FALSE) } ,error=function(e) e) 
  if(ncol(tab)==1) tab <- read.table(filename,header=TRUE,sep=";",stringsAsFactors=FALSE)
  return(tab) #need dataframe to keep allele-names correct!!
 }
 strsplit2 <- function(x,spl) {  #function which takes multiple split-seperators
  if(nchar(x)==0) return("")
  txt <- x
  for(j in 1:length(spl)) {
   txt <- unlist(strsplit(txt,split=spl[j]))
  }
  return(txt)
 }
 makematrix <- function(x) {  #function which transposes the vector back to matri
  if(is.null(dim(x))) x<- t(x)
  return(x)
 }
prim = as.integer(c(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997,1009,1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069,1087,1091,1093,1097,1103,1109,1117,1123,1129,1151,1153,1163,1171,1181,1187,1193,1201,1213,1217,1223,1229,1231,1237,1249,1259,1277,1279,1283,1289,1291,1297,1301,1303,1307,1319,1321,1327,1361,1367,1373,1381,1399,1409,1423,1427,1429,1433,1439,1447,1451,1453,1459,1471,1481,1483,1487,1489,1493,1499,1511,1523,1531,1543,1549)) #max number of alleles is 244

#0) Configurate Kit:
tab=tableReader(freqfile)
locs <- toupper(colnames(tab[,-1]))
popFreq <- popFreqP <- list()
for(loc in locs) {
 freqs <- tab[,which(loc==locs)+1]
 popFreq[[loc]] <- popFreqP[[loc]] <- tab[!is.na(freqs ),which(loc==locs)+1]
 names(popFreq[[loc]]) <- tab[!is.na(freqs ),1]
 names(popFreqP[[loc]]) <- prim[1:length(popFreq[[loc]])]
}
names(popFreq) <- names(popFreqP) <- locs

#2.1) IMPORT References:
bigREF <- matrix(nrow=0,ncol=length(locs)+1)
Rfiles <- NULL
if(!is.null(reffold)) Rfiles <- list.files(reffold) #No pattern in ref-files. 

for(Rfile in Rfiles) {
 X=tableReader(paste0(reffold,"/",Rfile))
 #File config:
 cn = colnames(X) #colnames 
 sind = grep("ID",toupper(cn)) #sample col-ind
 if(length(sind)==0) sind = grep("SAMPLE",toupper(cn)) #sample col-ind
 if(length(sind)>1) stop("More than one type of ID was find in the file")
 sn = unique(as.character(X[,sind])) #unique id 
 if(length(sn)==1) { #only one sample was found
  sn <- strsplit(Rfile,"\\.")[[1]][1] #get samplename of filename
  X[,sind] <- sn #update sample name within file!
  # next #this file did not contain reference-samples!
 }
 lind = grep("marker",tolower(cn)) #locus col-ind
 X[,lind] <- toupper(X[,lind]) #make all loci upper-case
 ln = unique(X[,lind]) #locus names: 
 if(length(grep("AM",ln))>0) ln <- ln[-grep("AM",ln)] #remove AMEL
 if(!all(ln%in%locs)) {
  print(paste0("Reference-File '",Rfile,"' had another kit-format than ESX17!"))
  next
 }
 Aind = grep("allele",tolower(cn))[1:2] #allele col-ind. Only consider the two first
 ishom  <- X[,Aind[2]]=="" | is.na(X[,Aind[2]])
 X[ishom ,Aind[2]] <- X[ishom ,Aind[1]]  #add homozygote to both

#NEW CODE: Introduce fast encoding:
 smallREF <- matrix(1,ncol=length(locs),nrow=length(sn))
 for(loc in locs) {
  locind <- which(X[,lind]==loc) #get index of loci
  tmp <- popFreq[[loc]]
  tmpP <- popFreqP[[loc]]
  
  newA1 <- X[locind,Aind[1]][!X[locind,Aind[1]]%in%names(popFreq[[loc]])] #get new Alleles
  newA2 <- X[locind,Aind[2]][!X[locind,Aind[2]]%in%names(popFreq[[loc]])] #get new Alleles
  newA <- unique(c(newA1,newA2)) #get unique alleles
  if( length(newA)>0 ) { 
    tmp <- popFreq[[loc]]
    tmpP <- popFreqP[[loc]]
    popFreqP[[loc]] <- popFreq[[loc]] <- c(popFreq[[loc]],rep(newf0,length(newA))) #insert new freqs.
    names(popFreq[[loc]]) <- c(names(tmp),newA) #insert new names
    names(popFreqP[[loc]]) <- c(names(tmpP),prim[(length(tmp)+1):length(popFreq[[loc]])]) #insert new P-names
  } 
  #convert:
  Pname  <- as.integer(names(popFreqP[[loc]]))  #get primenumbers
  Aname  <- names(popFreq[[loc]])  #get allele names of population
  for(an in Aname) {#for each alleles in population
   for(ai in Aind) {
    aind <- X[locind,ai]==an #get rows which has corresponding allele
    srow <- which(sn%in%X[locind,sind])#get belonging rows in smallREF (some samples may miss a locus!)
    smallREF[srow[aind],which(locs==loc)] <- smallREF[srow[aind],which(locs==loc)]*Pname[which(Aname==an)]
   }
  }#end for each alleles
 } #end for each locus
 smallREF[smallREF==1]=NA #insert missing loci as NA
 bigREF <- rbind(bigREF, cbind(sn,smallREF)) #add to matrix
} #end for each files
if(length(bigREF)>0) colnames(bigREF) <- c("ID",locs)
if(length(bigREF)==0) print("NOTE: No reference samples were imported!")

#Check for duplicated R-ID:
allsn <- bigREF[,1]
if(length(unique(allsn ))!=length(allsn )) stop("Some imported personal references had same R-ID") #is ok

#2.2) Import stains: 
TAfiles <- list.files(fn,pattern=TAptrn) #No pattern in files. SHould check samplenames
bigSS <- numeric() #matrix used to store OLD single source profiles
bigEVID <- numeric() #matrix used to store NEW stains

#Task:
#matching bigREF,bigOLDSS,bigNEWSS -> bigNEWEVID (all)
#matching bigNEWSS -> bigOLDEVID (mix)

for(TAfile in TAfiles) {
 fname <- paste0(fn,"/",TAfile)

 #check time of file:
 Tfile <- file.info(fname)$mtime #time when file was modified!
 Fdate <-  format(Tfile) #get time stample for file
 Tdiff <- difftime(Sys.time(),Tfile,units="days")
 Tdiff <- as.integer(strsplit(format(Tdiff,digits=1)," ")[[1]][1])
 if( Tdiff > Thist ) next #don't import if file was modified more than Thist days ago

 #determine if file is in category NEW
 isNEW <- Tdiff <= Tnew #if not over Tnew days old

 X=tableReader(fname)
 #File config:
 cn = colnames(X) #colnames 
 sind = grep("sample",tolower(cn)) #sample col-ind
 if(length(sind)>1) sind <- sind[grep("name",tolower(cn)[sind])] #if multiple names
 sn = unique(X[,sind]) #sample names
 sn <- sn[grep(BIDpat,sn)] #require Amplification pattern
 if(length(sn)==0) next #this file did not contain samples!

 lind = grep("marker",tolower(cn)) #locus col-ind
 X[,lind] <- toupper(X[,lind]) #make all loci upper-case
 ln = unique(X[,lind]) #locus names:
 ln <- ln[-grep("AMEL",ln)] #remove AMEL
 if(!all(ln%in%locs)) {
  print(paste0("Mixture-File '",TAfile,"' had another kit-format!"))
  next
 }
 Aind = grep("allele",tolower(cn)) #allele col-ind
 Hind = grep("height",tolower(cn)) #height col-ind
 if(length(Aind)==0 | length(Hind)==0) next #if no allele or peak height info

 #extract IDs
 tmpS <- NULL
 for(pat in SIDpat) {
  tmp <- sapply(strsplit(sn,pat ),function(x) x[2]) 
  if(all(!is.na(tmp))) {
   tmpS <- tmp
   break
  }
 }
 if(is.null(tmpS)) next
 tmp2 <- sapply(strsplit(tmpS,BIDpat),function(x) x[2]) 
 if(any(is.na(tmp2)) ) next
 SID <- sapply(strsplit(tmpS,"_"),function(x) x[1])  #extract S-nr
 KID <- sapply(strsplit(tmpS,"_"),function(x) x[3]) #extract KL-nr
 tmp <- strsplit2(TAfile,c(TAptrn,"\\."))[1] #TA-id
 TID <- paste0(TAptrn,tmp)#,substring(TAfile, nchar(TAptrn)+nchar(tmp)+1,nchar(TAptrn)+nchar(tmp)+3)) #get TA-ID

 for(ss  in 1:length(sn)) { #for each combination of SID and BID (unique stains)
    subX <- X[X[,sind ]==sn[ss],] #get submatrix
    MIXvec <- rep(NA,length(locs)) #place to store allele-info into evidence-matrix
    SSvec <- rep(NA,length(locs)) #place to store allele-info into MAJOR single source-matrix
    for(ii in 1:nrow(subX)) { #for each marker-row
      loc <- toupper(subX[ii,lind])
      locind <- grep(loc,locs) #find correct locus
      if(length(locind)==0) next
      Ainfo <- as.character(subX[ii,Aind]) #extract allele infor
      Hinfo <- as.numeric(subX[ii,Hind]) #extract peak height info
      usedInd <- !(is.na(Hinfo) | Ainfo=="" | Ainfo=="NA" | Ainfo=="OL" )
      Ainfo <- Ainfo[usedInd]
      Hinfo <- Hinfo[usedInd]

      #(1) Peak height above threshold
      keep <- Hinfo >= threshHeight  #require minumum peak height:
      Ainfo <- Ainfo[keep] 
      Hinfo <- Hinfo[keep] 
      if(length(Ainfo)==0) next #skip if empty

      #(2) Stutter-filter (again): 
      AllsNum <- as.numeric(Ainfo) #convert to numbers
      stuttindL <- which(as.character(AllsNum)%in%as.character(AllsNum-1)) #stutter-alleles one low bp
      stuttindH <- which(as.character(AllsNum-1)%in%as.character(AllsNum)) #stutter-alleles one high bp
      stuttR <- Hinfo[stuttindL]/Hinfo[stuttindH] #stutter-ratio is comparing observed peak heights
      remove<- stuttindL[ stuttR<threshStutt ] #alleles to remove
      if(length(remove)>0) {
       Ainfo <- Ainfo[-remove] 
       Hinfo <- Hinfo[-remove] 
      }
      #(3) Extract major profiles: Put into an own matrix:

      #decreasingly sort:
      ord <- order(Hinfo,decreasing=TRUE)
      Ainfo <- Ainfo[ord]
      Hinfo <- Hinfo[ord]
      SS <- Ainfo[1] #larges allele always included
      nA <- length(Ainfo)
      if(nA>1 && Hinfo[2]/Hinfo[1]>threshMaj) { #if more than 1 allele and ratios of two largest suffice
        SS <- c(SS,Ainfo[2]) #add second allele
        if(nA>2 && Hinfo[3]/Hinfo[2]>threshMaj) { #if more than 2 major alleles
          SS <- NA #couldn't deside anything for marker
        }
      }  

     #check if allele is missing
      Anew <- unique(Ainfo[!Ainfo%in%names(popFreq[[loc]])])
      if( length(Anew)>0 ) { 
         tmp <- popFreq[[loc]]
         tmpP <- popFreqP[[loc]]
         popFreqP[[loc]] <- popFreq[[loc]] <- c(popFreq[[loc]],rep(newf0,length(Anew))) #insert new freqs.
         names(popFreq[[loc]]) <- c(names(tmp),Anew) #insert new names
         names(popFreqP[[loc]]) <- c(names(tmpP),prim[(length(tmp)+1):length(popFreq[[loc]])]) #insert new P-names
      } 
      Pname  <- as.integer(names(popFreqP[[loc]]))
      MIXvec[locind] <- prod(Pname[names(popFreq[[loc]])%in%Ainfo] ) #use updated primenumbers
      if(!any(is.na(SS))) SSvec[locind] <- prod(Pname[names(popFreq[[loc]])%in%SS] ) #use updated primenumbers
    } #end for each markers
     bigSS <- rbind(bigSS, c(isNEW,Fdate,TID,KID[ss],SID[ss],sn[ss],SSvec) )  #add Extracted Single source profiles 
     bigEVID <- rbind(bigEVID, c(isNEW,Fdate,TID,KID[ss],SID[ss],sn[ss],MIXvec) ) #add to matrix
 } #end for each combination  of SID and BID
 if(which(TAfiles==TAfile)%%10==0) print(paste0(round(which(TAfiles==TAfile)/length(TAfiles)*100),"% import complete"))
} #end for each files
if(length(bigEVID)==0) stop("No evidence was imported! Program stops.")
cn <- c("NEW","Time","TAID","KLID","SID","ID",locs) #also adding KL and TA id
colnames(bigSS) <- colnames(bigEVID) <- cn
#head(bigSS)
startloc <- ncol(bigSS) - length(locs)

#Check for duplicates:
dupind  <- duplicated( bigEVID[,cn=="ID"] )
if(sum(dupind)>0) {
 print(bigEVID[dupind,1:startloc])
 print(paste0("There was ",sum(dupind)," duplicated matches. These are removed..."))
 bigEVID <- bigEVID[!dupind ,]
 bigSS <- bigSS[!dupind ,]
}

#count number of alleles in each marker: (could also be done while importing)
nalleles <- matrix(0,ncol=length(locs),nrow=nrow(bigEVID)) 
for(loc in locs) {
 Pname  <- as.integer(names(popFreqP[[loc]]))
 for(pp in Pname) { #for each primenumber:
  tmp <- round(as.numeric(bigEVID[,which(loc==locs)+startloc]))
  inds <- which(tmp%%pp==0)
  if(length(inds)>0) nalleles[inds,which(loc==locs)] <- nalleles[inds,which(loc==locs)] + 1 #add to counter
 }
}

#Categorize:
########
#Empty:#
########
empty <- which(rowSums(nalleles)==0)
if(length(empty)>0) { #remove empty profiles
 bigSS <- bigSS[-empty,] 
 bigEVID <- bigEVID[-empty,] 
 nalleles <- nalleles[-empty,]
}

########
#SINGLE#
########
#single-source (those having always less than 2 alleles in all markers)
single <- which(apply(nalleles,1, function(x) all(x<=2))) #get indice of single source-profiles

#seperate as partial and full single-source:
singlePart <- single[apply(nalleles[single,],1, function(x) any(x==0))] 
singleFull <-  single[!single%in%singlePart]

#####
#MIX#
#####
#mix-source (with possible stutter): (observation with at least 2)
mixture <- which(apply(nalleles,1, function(x) any(x>2)  )) 

#seperate as partial and full mixtures:
mixturePart <- mixture[apply(nalleles[mixture,],1, function(x) any(x==0)  )]  #those missing any markers
mixtureFull <-  mixture[!mixture%in%mixturePart] #have all markers

#######
#MAJOR#
#######
#Only for mixtures. Used to seperate from single-profiles
#seperate as partial and full single-source:
majorFull <- mixtureFull[which(rowSums(is.na(bigSS[mixtureFull,]))==0)]  #full extracted majors
majorPart <- mixture[which(rowSums(is.na(bigSS[mixture,]))>0 & rowSums(is.na(bigSS[mixture,]))<length(locs))] #partial majors (those with SOME markers)
major <- c(majorFull,majorPart) #majors

############################################################################################################################################
######################################IMPORT + CATEGORIZATION DONE#########################################################################
#############################################################################################################################################

#should be empty:
intersect(singlePart,singleFull) #
intersect(mixturePart,mixtureFull)
intersect(majorFull,majorPart)
intersect(single,mixture)
intersect(single,major)
#indices of those having extracted a major:
#intersect(mixture,major) 

#SEARCH WITH LR:
#Compare aginst Single: Assume C=1, model (pC=0.05,pD=0.1)
#Compare aginst mixture: Estimate C = ceiling(max(nA)/2) and then use model (pC=0.05,pD=0.1 (or estimate? No more consumuing)) 

#Step 1) Calculate LR for each genotype-profiles in population for each loci against each (single,mixture)

#Merge single-profiles and mixtures for each locus:
indNEWREF <- intersect(c(single,major),which(bigEVID[,1]=="TRUE")) #get NEW-Single-sources
indOLDREF <- intersect(c(single,major),which(bigEVID[,1]=="FALSE")) #get OLD-Single-sources
nallREF <- nrow(bigREF)+length(indOLDREF)+length(indNEWREF) #TOTAL NUMBER OF REFERENCES (ref,old and new)

indTAR <- c(which(bigEVID[,1]=="TRUE"), intersect( which(bigEVID[,1]=="FALSE"),mixture) ) #get target samples to search
bigLRtoEVID <- matrix(0,nrow=nallREF ,ncol=length(indTAR)) #storing results in stain-comparisons

print("Calculating LR for all combinations:")
print(paste0("nReferences:",nrow(bigLRtoEVID)))
print(paste0("nSamples:",ncol(bigLRtoEVID)))
print(paste0("nComparisons:",prod(dim(bigLRtoEVID))))

#Step 1: Estimate number of contributors:
Chat <- apply(nalleles[indTAR,],1,function(x) ceiling(max(x)/2)) #estimated number of contr

for(loc in locs) { #for each locus
 print(loc)
 lindS <- which(locs==loc)+startloc #column-index for locus for stains
 lindR <- which(locs==loc)+1 #column-index for locus for references
 Pname  <- as.integer(names(popFreqP[[loc]]))

 #get all C=1 - profiles: (refs,singles,SSmajor)
 allREF <- c(bigREF[,lindR], bigEVID[single,lindS],bigSS[major,lindS]) #Merge ALL types in references 

 NAallREF  <- is.na(allREF) #get NA (empty) references (they are not stored for LR)
 unallREF <- unique(allREF[!NAallREF]) #get unique reference-genotypes

 #for each evidence (mixture evid)

 #Simple model:
 SSevid <- cbind( bigEVID[indTAR,lindS ],Chat ) #consider NEW stains only
 unSSevid <- unique(SSevid)
 unSSevid <- unSSevid[!is.na(unSSevid[,1] ),] #remove NA evidence
 #if(is.null(dim(unSSevid))) unSSevid  <- t(unSSevid )

 for(ss in 1:nrow(unSSevid)) {
     nC <- as.integer(unSSevid[ss,2]) #number of contributors assumed for the evidence
     Ei <- Pname[ round(as.numeric(unSSevid[ss,1]))%%Pname==0] #get P-allele names in evidence
     subfreq <- popFreqP[[loc]][Pname%in%Ei]
     subfreq <- c(subfreq,1-sum(subfreq))
     hd0 <- likEvid(Repliste=Ei,T=NULL,V=NULL,x=nC,theta=0,prDHet=rep(pC,nC),prDHom=rep(pC^2,nC),prC=pC,freq=subfreq ) 
 
     insNEW <- which(SSevid[,1]==unSSevid[ss,1] & SSevid[,2]==unSSevid[ss,2]) #columns where to assign calculation (for all references)
     bigLRtoEVID[!NAallREF,insNEW] <- bigLRtoEVID[!NAallREF,insNEW] - log(hd0)   #NB: DONT insert for NAs
     for(rr in 1:length(unallREF)) {      #for each unique references
      Ri <- Pname[round(as.numeric(unallREF[rr]))%%Pname==0] #get P-allele names in references
      if(length(Ri)==1) Ri  <- rep(Ri,2)
      hp0 <- likEvid(Repliste=Ei,T=Ri,V=NULL,x=nC-1,theta=0,prDHet=rep(pC,nC),prDHom=rep(pC^2,nC),prC=pC,freq=subfreq ) 
      insREF <- which(allREF==unallREF[rr]) #references to insert
      bigLRtoEVID[insREF,insNEW] <- bigLRtoEVID[insREF,insNEW] + log(hp0)
     } #end for each reference
 } #end for each evidence
}
#Finished to calculate LR

#Get ID for all matches:
IDind <- which("ID"==colnames(bigEVID))
SIDind <- which("SID"==colnames(bigEVID))
refID <- c(bigREF[,1],bigEVID[single,IDind],bigSS[major,IDind]) #references to check
refSID <- c(bigREF[,1],bigEVID[single,SIDind],bigSS[major,SIDind]) #STAIN ID (SID) only

tarID <-  bigEVID[indTAR,IDind]
tarSID <-  bigEVID[indTAR,SIDind]

#Get whether ref is (NA,OLD,NEW) :
allref_isNEW <- c(rep("FALSE",nrow(bigREF)),bigEVID[c(single,major),which("NEW"==colnames(bigEVID))]) #boolean whether references isNEW: Personel is OLD
#get type of ref:
allref_type = c(rep("P",nrow(bigREF)),rep("S",length(single)),rep("M",length(major))) #P - personel, S - single source, M - major

#functions to return subset of reference-matrix and target-matrix of type what :
getAllref_what <- function(what,NAref=TRUE) {
 ind <- which(what==colnames(bigEVID))
 vecList <- list()
 if(NAref) vecList$REF <- c(rep(NA,nrow(bigREF)),bigEVID[c(single,major),ind]) #what to references 
 if(!NAref) vecList$REF <- c(rep("0",nrow(bigREF)),bigEVID[c(single,major),ind]) #what to references 
 vecList$TAR <-bigEVID[indTAR,ind] #get what for indTAR  stains
 return(vecList)
} 

vecTime <- getAllref_what("Time")  #Get Time for all matches:
vecTAID <- getAllref_what("TAID",NAref=FALSE)  #Get TA-ID for all matches:
vecKLID <- getAllref_what("KLID",NAref=FALSE)  #Get KL-ID for all matches:


#MATCH LIST: Print a list of all matches (comparisons)
matches  <- which(exp(bigLRtoEVID)>=threshLR,arr.ind = TRUE) #assume always a matrix

#Filter:
#1) Remove because it was the same stain:
matches <- makematrix(matches[!refSID[matches[,1]]==tarSID[matches[,2]],]) #remove same stain

#2) Remove because it was the same KL-number
if(!sameKL) {
 matches <- makematrix(matches[vecKLID$REF[matches[,1]]!=vecKLID$TAR[matches[,2]],]) #keep only those with different KL
}

#get number of days between matches (TA-registed files):
refTimes <- strptime(vecTime$REF,format="%Y-%m-%d %H:%M:%S")
tarTimes <- strptime(vecTime$TAR,format="%Y-%m-%d %H:%M:%S")
diffMatchDays <- difftime(tarTimes[matches[,2]],refTimes[matches[,1]],units="days") #number of days: (target - referanse)
#sort(diffMatchDays)
#matchChron <- matches[diffMatchDays<1,] #Chronologic: require that target is registrated at least -1 days after ref
#matchSameTA <- matches[vecTAID$REF[matches[,1]]==vecTAID$TAR[matches[,2]],] #keep only those with similar TA

#List of matches: 
#TypeMatch: <3 within (W) or >=3 from old to new (B)
#Timediff (days), typeRef (SS,MAJ),ref(ID,TA),tar(ID,TA), (Chat in target)
createList = function(match) {
 tarNEW <- bigEVID[match[,2],which("NEW"==colnames(bigEVID))]
 refNEW <- allref_isNEW[match[,1]]

 timediff <- round(difftime(vecTime$TAR[match[,2]],vecTime$REF[match[,1]],units="days"),1) #number of days back from target-TA to when ref-TA was modified 
 type <- allref_type[match[,1]]
 nCcol <- Chat[match[,2]] #number of estimated contributors in target
 refCol <- cbind(refID[match[,1]],vecTAID$REF[match[,1]])
 tarCol <- cbind(tarID[match[,2]],vecTAID$TAR[match[,2]])
 LRcol <- exp(bigLRtoEVID)[match]
 ord <- order(as.numeric(LRcol),decreasing=TRUE)
 matchlist <- makematrix(cbind(refNEW,tarNEW,timediff,type,refCol,tarCol,nCcol,LRcol)[ord,])
 cn <- c("ID","TA")
 colnames(matchlist) <-  c("refNEW","tarNEW","Timediff","ref-type",paste0( c(rep("ref",length(cn)),rep("target",length(cn))) ,cn),"nC","LR")
 return(matchlist)
}

matchlist <- createList(match=matches)

#FILTER AWAY MATCHES OF TYPE: "REF(P) -> EVID_OLD" and "SS_OLD -> EVID_OLD"
matchlist <- makematrix(matchlist[!(matchlist[,1]=="FALSE" & matchlist[,2]=="FALSE"),]) #do not show results from search within OLD 
#matchSameTAlist <- createList(matchSameTA) #matches with same TA
#matchDifferTAlist <- createList(matchNEW[!(matchNEW[,1]%in%matchNEWsameTA[,1] & matchNEW[,2]%in%matchNEWsameTA[,2]),])

#FILTER AWAY SYMMETRIC MATCHES (i.e. when ref-type="S" and nC=1)
matchlist2 <- matchlist3 <- makematrix(matchlist[,c(5,7)]) #take only unique matches
dupind <- which(duplicated(rbind(matchlist3,matchlist2[,2:1]))) - nrow(matchlist2) #get indices of duplicates
rm <- numeric()
while(1) {
 if(length(dupind)==0) break #stop when all duplicates are considered
 ind <- matchlist3[dupind[1],1]==matchlist3[dupind,2] & matchlist3[dupind[1],2]==matchlist3[dupind,1] 
 rm <- c(rm,dupind[ind])
 dupind <- dupind[-1] #remove corresponding match with lowest LR
}
if(length(rm)>0) matchlist <- makematrix(matchlist[-rm,]) #removing symmetrical matches
stamp <- format(Sys.time(), "%y-%m-%d-%H-%M-%S") #timestamp: Year-month-day-hour-minute-second

###############
#STORE MATCHES#
###############
sfolder <- "sessions" #name of a session folder
if(length( grep(sfolder,list.dirs(full.names=FALSE,recursive=FALSE)) )==0) dir.create(sfolder) #create session folder
save.image(file=paste0("sessions/matching_",stamp,".Rdata")) #store session
#store match-lists:
#store matches not personell:
write.table(makematrix(matchlist[matchlist[,4]!="P",]),file=paste0("sessions/stainmatches_",stamp,".csv"),row.names=FALSE,sep =";") #store session
#store only personell:
write.table(makematrix(matchlist[matchlist[,4]=="P",]),file=paste0("sessions/personelmatches_",stamp,".csv"),row.names=FALSE,sep =";") #store session


################
#STORE Profiles#
################

matchingprofiles = function(matchlist) {
 if(nrow(matchlist)==0) return()
 cn <- colnames(matchlist)
 refind <- which("refID"==cn)
 tarind <- which("targetID"==cn)
 refindTA <- which("refTA"==cn)
 tarindTA <- which("targetTA"==cn)
 EVIDcn <- colnames(bigEVID)

 fulltab <- numeric()
 for(ss in 1:nrow(matchlist)) { #for each matched
  SSrow <- bigSS[bigSS[,IDind]==matchlist[ss,refind],] #get row of reference (bigSS)
  if(length(SSrow)==0) SSrow <- bigREF[bigREF[,1]==matchlist[ss,refind],] #get row of reference (bigREF)
  SScn <- names(SSrow)

  EVIDrow <- bigEVID[bigEVID[,IDind]==matchlist[ss,tarind],] #get row of target stain
  nC <- paste0("nC=",matchlist[ss,"nC"==cn])
  lr <- paste0("LR=",matchlist[ss,"LR"==cn])
  rowtab  <- rbind(c(matchlist[ss,refind], matchlist[ss,tarind]),c(matchlist[ss,refindTA], matchlist[ss,tarindTA]),matchlist[ss,1:2],c(matchlist[ss,4],nC),rep("-------",2))  

  for(loc in locs) { #for each locus
   Pname  <- as.integer(names(popFreqP[[loc]]))
   EPi <- round(as.numeric(EVIDrow[loc==EVIDcn]))
   RPi <- round(as.numeric(SSrow[loc==SScn]))
   Ri  <- Ei <- NA
   if(!is.na(RPi)) Ri <- paste0(names(popFreq[[loc]])[RPi%%Pname==0],collapse="/")
   if(!is.na(EPi)) Ei <- paste0(names(popFreq[[loc]])[EPi%%Pname==0],collapse="/")
   rowtab  <- rbind(rowtab  , c(Ri,Ei))
  }
  fulltab <- rbind(fulltab,rep(paste0("----------",ss,"----------"),2) )
  fulltab <- rbind(fulltab,rowtab )
 }
 return(fulltab)
}

#store details of all matches:
write.table(matchingprofiles(matchlist),file=paste0("sessions/matchinfo_",stamp,".csv"),row.names=FALSE,sep =";") #store session
SID <- makematrix(matchlist[,c(5,7)]) #get SID of ref and tar of current matches
print(paste0("Number of matches=",nrow(matchlist)))

#THIS BLOCK: Loads matches in matchfile. Makes sure that new matches don't duplicate with old
matchfile <- paste0("matchfile.csv")
if(!file.exists(matchfile)) { #create file if not found
 cn <- c("refNEW","tarNEW","Timediff","ref-type","refID","refTA","targetID","targetTA","nC","LR","Time","checked")
 x <- matrix(,ncol=length(cn),nrow=0)
 colnames(x) <- cn
 write.table(x,file=matchfile,row.names=FALSE,col.names = TRUE,sep =";") 
}

if(nrow(matchlist)>0) { #if at least 1 match, 
 #load and store match-file:
 matchlist2 <- read.table(file=matchfile,sep =";",header=TRUE) 
 SID2 <- makematrix(matchlist2[,c(5,7)])   #get SID of ref and tar of previous matches
 isdup <- duplicated(rbind(SID2,SID)) #found (redudance) matches already found before
 matchlist3 <- makematrix(matchlist[!isdup[(nrow(matchlist2)+1):length(isdup)],]) #get only NEW matches!
 if(length(matchlist3)>0) {
   matchlist3 <- cbind(matchlist3,stamp,"") #include Time and checked
   colnames(matchlist3) <-  colnames(matchlist2) 
   write.table(matchlist3,file=matchfile,row.names=FALSE,col.names = FALSE,sep =";",append=TRUE) 
 }
}


} #end dnamatch function