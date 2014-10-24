#' @title euroformixTK
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description euroformixTK is a GUI wrapper for gammadnamix
#' @details The function is a graphical layer for the functions in the package gammadnamix. See ?gammadnamix for more information.
#' @param envirfile A file to a saved environment of a project
#' @export

#library(roxygen2)
#setwd("C:/Users/oebl/Dropbox/Forensic/MixtureProj/myDev/")
#roxygenize("gammadnamix")

#library(gammadnamix); sessionInfo();euroformix()
#setwd("C:/Users/oebl/Dropbox/Forensic/MixtureProj/myDev/quantLR/euroformix0")
#rm(list=ls())
#envirfile=NULL
#source("euroformix.R");euroformix()

euroformix = function(envirfile=NULL) {

 #size of main window
 mwH <- 1000
 mwW <- 1500

 #Required for GUI:
 require(gWidgetstcltk) #requires only gWidgets.

 #type of gwidgets-kit
 options(guiToolkit="tcltk")

 #version:
 version = 1

 #software name:
 softname <- paste0("EuroFormix v",version)


 #####################
 #create environment #
 #####################
 if(is.null(envirfile)) {
  mmTK = new.env() #create new envornment object

  #Toolbar options: can be changed any time by using toolbar
  assign("optFreq",list(freqsize=0),envir=mmTK) #options when new frequencies are found (size of imported database,minFreq)
  assign("optMLE",list(nDone=3,delta=10,dec=4),envir=mmTK) #options when optimizing (nDone,delta)
  assign("optMCMC",list(delta=10,niter=10000),envir=mmTK) #options when running MCMC-simulations (delta, niter)
  assign("optINT",list(reltol=0.005,maxmu=20000,maxsigma=1),envir=mmTK) #options when integrating (reltol and boundaries)
  assign("optDC",list(alphaprob=0.9999,maxlist=1000),envir=mmTK) #options when doing deconvolution (alphaprob, maxlist)
  assign("optDB",list(maxDB=10000,QUALpC=0.05),envir=mmTK) #options when doing LRmix
  assign("optLRMIX",list(range=0.6,nticks=31,nsample=2000,alpha=0.05,ntippets=1e6),envir=mmTK) #options when doing database search (maxDB)


  #initializing environment variables
  assign("workdir",NULL,envir=mmTK) #assign working directory to mmTK-environment
  assign("freqfolder",NULL,envir=mmTK) #assign freqfolder to mmTK-environment
  assign("kits",NULL,envir=mmTK) #assign kitList to mmTK-environment
  assign("selPopKitName",NULL,envir=mmTK) #selected kit and population for popFreq

  #imported data:
  assign("popFreq",NULL,envir=mmTK) #assign to mmTK-environment
  assign("mixData",NULL,envir=mmTK) #assign to mmTK-environment
  assign("refData",NULL,envir=mmTK) #assign to mmTK-environment
  assign("dbData",NULL,envir=mmTK) #assign dbData: matrix referenceDatabase to mmTK-environment (encoded)

  #models: stored setups for model specification
  assign("setEVID",NULL,envir=mmTK) #assign model (evidence weighting)
  assign("setDB",NULL,envir=mmTK) #assign model (database search)
  assign("setDC",NULL,envir=mmTK) #assign model (deconvolution)
  assign("setGEN",NULL,envir=mmTK) #assign model (generation)

  #results: stored results after calculations
  assign("resDB",NULL,envir=mmTK) #assign database search results (i.e. ranked table of results)
  assign("resEVID",NULL,envir=mmTK) #assign evidence weighting results (i.e. calculated LR with MLE estimates)
  assign("resDC",NULL,envir=mmTK) #assign deconvolved result (i.e. ranked table of results)
  assign("resEVIDINT",NULL,envir=mmTK) #assign evidence weighting results - Based on Integration 
  assign("resEVIDLRMIX",NULL,envir=mmTK) #assign evidence weighting results - Based on LRmix
 } else {
  load(envirfile) #loading environment
 }

 ####################################
 #auxiliary functions and variables:#
 ####################################
 prim = c(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113, 127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263, 269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421, 431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593, 599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757, 761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941, 947,953,967,971,977,983,991,997,1009,1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069,1087,1091,1093, 1097,1103,1109,1117,1123,1129,1151,1153,1163,1171,1181,1187,1193,1201,1213,1217,1223,1229,1231,1237,1249, 1259,1277,1279,1283,1289,1291,1297,1301,1303,1307,1319,1321,1327,1361,1367,1373,1381,1399,1409,1423,1427, 1429,1433,1439,1447,1451,1453,1459,1471,1481,1483,1487,1489,1493,1499,1511,1523,1531,1543,1549) 

 #helpfunction to return minimum frequency (used for new alleles)
 getminFreq = function() {
  popFreq <- get("popFreq",envir=mmTK) #get selected popFreq
  freqsize <- get("optFreq",envir=mmTK)$freqsize #get selected size of frequence-database
  if(freqsize>0) {
   return(5/(2*freqsize))
  } else {
   return(min(unlist(popFreq))) #minumum observed frequency was used 
  }
 }

 #Function to get data from environment
 #sel used to select a specific datasubset
 getData = function(type,sel=NULL) {
   Data <- NULL
   if(type=="mix") Data <- get("mixData",envir=mmTK) #assign kit to mmTK-environment
   if(type=="ref") Data <- get("refData",envir=mmTK) #assign kit to mmTK-environment 
   if(type=="db") Data <- get("dbData",envir=mmTK) #assign kit to mmTK-environment 
   if(!is.null(sel)) return(Data[sel]) #returns only selected datasubset
   return(Data)
 }

 #function for inserting sample/ref/db-names into existing gcheckboxgroup
 restoreCheckbox = function(type) {
  subD <- getData(type)
  if(!is.null(subD)) { return( names(subD))
  } else { return(numeric()) }
 }

 #Load kit with allelefreq-info from filefolder:
 loadKitList = function(freqpath) {
  freqfiles = list.files(freqpath)
  kitList <- list()
  for(i in 1:length(freqfiles)) {
   filename = paste(freqpath,"/",freqfiles[i],sep="")
   tab=tableReader(filename)
   Anames = tab[,1]
   tab = tab[,-1]
   freqlist = list()
   for(j in 1:ncol(tab)) { #for each locus
     tmp = tab[,j]
     tmp2 = tmp[!is.na(tmp)]
     names(tmp2) = Anames[!is.na(tmp)]
     freqlist[[j]] = tmp2
   }
   names(freqlist) = toupper(colnames(tab)) #LOCUS-names are assigned as Upper-case
   kit = unlist(strsplit(freqfiles[i],"_"))[1]
   pop = unlist(strsplit(freqfiles[i],"_"))[2]
   pop = unlist(strsplit(pop,"\\."))[1]
   kitind = kit==names(kitList) 
   kitList[[kit]][[pop]] = freqlist
  }
  assign("kits",kitList,envir=mmTK) #assign kit to mmTK-environment
 }

 #Function which takes rownames and adds to first column
 addRownameTable = function(tab) {
  tmp <- colnames(tab)
  tab <- cbind(rownames(tab),tab)
  colnames(tab) <- c("X",tmp)
  return(tab)
 }

 #Robust function for reading tables:
 tableReader=function(filename) {
  tab <- read.table(filename,header=TRUE,sep="\t",stringsAsFactors=FALSE)
  tryCatch( {  if(ncol(tab)==1) tab <- read.table(filename,header=TRUE,sep=",",stringsAsFactors=FALSE) } ,error=function(e) e) 
  tryCatch( {  if(ncol(tab)==1) tab <- read.table(filename,header=TRUE,sep=";",stringsAsFactors=FALSE) } ,error=function(e) e) 
  if(ncol(tab)==1) tab <- read.table(filename,header=TRUE,sep=";",stringsAsFactors=FALSE)
  return(tab) #need dataframe to keep allele-names correct!!
 }

 #save result table to file:
 saveTable = function(tab,sep="txt") {
  tabfile  = gfile(text="Save table",type="save") #csv is correct format!
  if(!is.na(tabfile)) {
   if(length(unlist(strsplit(tabfile,"\\.")))==1) tabfile = paste0(tabfile,".",sep)
   if(sep=="txt" | sep=="tab") write.table(tab,file=tabfile,quote=FALSE,sep="\t",row.names=FALSE) 
   if(sep=="csv") write.table(tab,file=tabfile,quote=FALSE,sep=",",row.names=FALSE) 
   print(paste("Table saved in ",tabfile,sep=""))
  }
 } #end file

 strsplit2 <- function(x,spl) {
  if(nchar(x)==0) return("")
  txt <- x
  for(j in 1:length(spl)) {
   txt <- unlist(strsplit(txt,split=spl[j]))
  }
  return(txt)
 }

 #Helpfunctions for converting profiles from list to table.
 sample_listToTable = function(Y) {
   sn = names(Y) #Y is a list on form Y[[samplename]][[locusname]]$adata,Y[[samplename]][[locusname]]$hdata
   aM = 0   #count number of max allele data:
   hM = 0   #count number of max allele heights:
   for(ss in sn) { #for each sample
    aM = max(unlist( lapply(Y[[ss]],function(x) length(x$adata)) ),aM)
    hM = max(unlist( lapply(Y[[ss]],function(x) length(x$hdata)) ),hM)
   }
   #create tables:
   X=numeric()
   for(ss in sn) { #for each sample
    newsample=numeric() #for allele
    ln = names(Y[[ss]])
    for(loc in ln) {
     newrow = Y[[ss]][[loc]]$adata
     newsample = rbind(newsample, c(newrow,rep("",aM-length(newrow))))
    }
    newsample2=numeric() #for heights
    if(hM>0) {
     for(loc in ln) {
      newrow = Y[[ss]][[loc]]$hdata
      newsample2 = rbind(newsample2, c(newrow,rep("",hM-length(newrow))))
     }      
    }
    X = rbind(X,cbind(ss,ln,newsample,newsample2))
   }
   cn = c("SampleName","Marker", paste("Allele",1:aM,sep=""))
   if(hM>0) cn = c(cn,paste("Height",1:hM,sep=""))
   colnames(X)  = cn
   return(X)
 } #end of functions

 #Helpfunctions for converting profiles from table to list.
 sample_tableToList = function(X) {
  cn = colnames(X) #colnames 
  lind = grep("marker",tolower(cn),fixed=TRUE) #locus col-ind
  if(length(lind)==0) lind = grep("loc",tolower(cn),fixed=TRUE) #try another name
  sind = grep("sample",tolower(cn),fixed=TRUE) #sample col-ind
  if(length(sind)>1)  sind = sind[grep("name",tolower(cn[sind]),fixed=TRUE)] #use only sample name
  A_ind = grep("allele",tolower(cn),fixed=TRUE) #allele col-ind
  H_ind = grep("height",tolower(cn),fixed=TRUE) #height col-ind
  ln = unique(toupper(X[,lind])) #locus names: Convert to upper case
  sn = unique(X[,sind]) #sample names
  I = length(ln)
  Y = list() #insert non-empty characters:
  for(k in 1:length(sn)) { #for each sample in matrix
   Y[[sn[k]]] = list() #one list for each sample
   for(i in 1:I) { #for each locus
     xind = X[,sind]==sn[k] & toupper(X[,lind])==ln[i] #get index in X for given sample and locus
     keep <- which(!is.na(X[xind,A_ind]) & X[xind,A_ind]!="")
     if(length(keep)>0) {
      if(length(A_ind)>0) Y[[sn[k]]][[ln[i]]]$adata = as.character(X[xind,A_ind][keep])
      if(length(H_ind)>0) Y[[sn[k]]][[ln[i]]]$hdata = as.numeric(as.character(X[xind,H_ind][keep]))
     }
   }
  }
  return(Y)
 }

  


###################################################################
###########################GUI#####################################
###################################################################

 #Menu bar file-lists:
 f_setwd = function(h,...) {
  dirfile = gfile(text="Select folder",type="selectdir")
  if(!is.na(dirfile)) {
   setwd(dirfile)
   assign("workdir",dirfile,envir=mmTK) #assign working directory
  }
 }
 f_openproj = function(h,...) {
  projfile = gfile(text="Open project",type="open", filter=list("Project"=list(patterns=list("*.Rdata"))))
  if(!is.na(projfile)) {
   dispose(mainwin)
   euroformix(projfile) #send environment into program
  }
 }
 f_saveproj = function(h,...) {
  projfile = gfile(text="Save project",type="save")
  if(!is.na(projfile)) {
   if(length(unlist(strsplit(projfile,"\\.")))==1) projfile = paste0(projfile,".Rdata")
   print("Size of stored objects (in MB):") #prints size of each stored object
   print(sapply(mmTK,object.size)/1e6) #prints size of each stored object
   save(mmTK,file=projfile) #save environment
   print(paste("Project saved in ",projfile,sep=""))
  }
 }
 f_quitproj = function(h,...) {
  ubool <- gconfirm("Do you want to save project?",title="Quit Program",icon="info")
  if(svalue(ubool)) {
   savefile <- gfile(text="Save file",type="save")
   if(length(unlist(strsplit(savefile ,"\\.")))==1) savefile = paste0(savefile ,".Rdata")
   save(mmTK,file=savefile) 
   print(paste("Project saved in ",savefile,sep=""))
  } else { 
   print("Program terminated without saving")
  }
  dispose(mainwin) #remove window!
 }

 #helpfunction to get value in from user and store
 setValueUser <- function(what1,what2,txt) {
   listopt <- get(what1,envir=mmTK) #get object what 1.
   val <- listopt[[what2]]
   sw <- gwindow(title="User input",visible=FALSE, width=300,height=50)
   grid <- glayout(spacing=0,container=sw )
   grid[1,1] <- glabel(txt, container=grid)
   grid[1,2] <- gedit(text=val,container=grid,width=15)
   grid[2,1] <- gbutton("OK", container=grid,handler = function(h, ...) { 
    listopt[[what2]] <- as.numeric(svalue(grid[1,2])) #insert new value
    assign(what1,listopt,envir=mmTK) #assign user-value to opt-list
    dispose(sw)
   } )
   grid[2,2] <- gbutton("Cancel", container=grid,handler = function(h, ...) { dispose(sw) } )
   visible(sw) <- TRUE
 }

 #helpfunction which checks that at value is in interval of [0,1]
 checkProb = function(x,what) {
  if(x < 0 || x>1) {
   gmessage(message=paste0(what," must be specified in interval [0,1] "),title="Wrong input",icon="error")
   stop("Wrong user-input")
  }
 }
 checkPositive = function(x,what) {
  if(x < 0 ) {
   gmessage(message=paste0(what," cannot be a negative number"),title="Wrong input",icon="error")
   stop("Wrong user-input")
  }
 }
 checkPosInteger = function(x,what) {
  if(x < 1 || round(x)!=x ) {
   gmessage(message=paste0(what," must be a positive integer"),title="Wrong input",icon="error")
   stop("Wrong user-input")
  }
 }



###########GUI WINDOW STARTS#########
 
 ##########
 #Menu bar#
 ##########
 mblst = list( #project saving and so on
  File=list(  
   'Set directory'=list(handler=f_setwd),
   'Open project'=list(handler=f_openproj),
   'Save project'=list(handler=f_saveproj),
   'Quit'=list(handler=f_quitproj,icon="close")
  ),
  Frequencies=list(
   'Set size of frequency database'=list(handler=function(h,...) {  
      setValueUser(what1="optFreq",what2="freqsize",txt="Set size of imported database:") 
    })
  ),
  Optimization=list(
   'Set number of random startpoints'=list(handler=function(h,...) {  
      setValueUser(what1="optMLE",what2="nDone",txt="Set number of random startvalues:") 
    }),
   'Set variance of randomizer'=list(handler=function(h,...) {  
      setValueUser(what1="optMLE",what2="delta",txt="Set variance of randomizer:") 
    })
  ),
  MCMC=list(
   'Set number of samples'=list(handler=function(h,...) {  
      setValueUser(what1="optMCMC",what2="niter",txt="Set number of sample iterations:") 
    }),
   'Set variance of randomizer'=list(handler=function(h,...) {  
      setValueUser(what1="optMCMC",what2="delta",txt="Set variation of randomizer:") 
    })
  ),
  Integration=list(
   'Set relative error requirement'=list(handler=function(h,...) {  
      setValueUser(what1="optINT",what2="reltol",txt="Set relative error:") 
    }),
   'Set maximum of mu-parameter'=list(handler=function(h,...) {  
      setValueUser(what1="optINT",what2="maxmu",txt="Set upper boundary of mu parameter:") 
    }),
   'Set maximum of sigma-parameter'=list(handler=function(h,...) {  
      setValueUser(what1="optINT",what2="maxsigma",txt="Set upper boundary of mu parameter:") 
    })
  ),
  Deconvolution=list(
   'Set required summed probability'=list(handler=function(h,...) {  
      setValueUser(what1="optDC",what2="alphaprob",txt="Set required summed posterior genotype-probability of list:") 
    }),
   'Set max listsize'=list(handler=function(h,...) {  
      setValueUser(what1="optDC",what2="maxlist",txt="Set size of maximum elements in deconvoluted list:") 
    })
  ),
  'Database search'=list(
   'Set maximum view-elements'=list(handler=function(h,...) {  
      setValueUser(what1="optDB",what2="maxDB",txt="Set max size of view elements from database:") 
    }),
   'Set drop-in probability for qualitative model'=list(handler=function(h,...) {  
      setValueUser(what1="optDB",what2="QUALpC",txt="Set allele drop-in probability for qualitative model:") 
    })
  ),
  'Qual LR'=list(
   'Set upper range for sensitivity'=list(handler=function(h,...) {  
      setValueUser(what1="optLRMIX",what2="range",txt="Set upper limit of dropout in sensitivity analysis:") 
    }),
   'Set nticks for sensitivity'=list(handler=function(h,...) {  
      setValueUser(what1="optLRMIX",what2="nticks",txt="Set number of ticks in sensitivity analysis:") 
    }),
   'Set required samples in dropout distr.'=list(handler=function(h,...) {  
      setValueUser(what1="optLRMIX",what2="nsample",txt="Set required number number of samples in dropout distribution:") 
    }),
   'Set significance level in dropout distr.'=list(handler=function(h,...) {  
      setValueUser(what1="optLRMIX",what2="alpha",txt="Set significance level (quantiles) in dropout distribution:") 
    }),
   'Set number of tippets'=list(handler=function(h,...) {  
      setValueUser(what1="optLRMIX",what2="ntippets",txt="Set number of tippet samples in tippet plot:") 
    })
  )
 )

##################################################################################################
########### Program starts #######################################################################
##################################################################################################

 #change working directory to the one stored in mmTK-environment
 wd=get("workdir",envir=mmTK) #assign working directory to mmTK-environment
 if(!is.null(wd)) setwd(wd)
 
 #Main window:
 spc <- 20
 mainwin <- gwindow(softname, visible=FALSE, width=mwW,height=mwH)
 gmenu(mblst,container=mainwin)
 nb = gnotebook(container=mainwin)
 tabGEN = ggroup(horizontal=FALSE,expand=TRUE,spacing=spc,container=nb,label="Generate data") #tab1: Generates data(with peak heights) for a given model (plot EPG in addition)
 tabimport = ggroup(horizontal=FALSE,expand=TRUE,spacing=40,container=nb,label="Import data") #tab2: (imports all files)
 tabmodel = ggroup(horizontal=FALSE,expand=TRUE,spacing=spc,container=nb,label="Model specification") #tab3: specify model used in weight-of-evidence (INT/MLE) or in a Database search 
 tabMLE = ggroup(horizontal=FALSE,expand=TRUE,spacing=spc,container=nb,label="MLE fit") #results from MLE
 tabDC = ggroup(horizontal=FALSE,expand=TRUE,spacing=spc,container=nb,label="Deconvolution") #results from a deconvolution
 tabDB= ggroup(horizontal=FALSE,expand=TRUE,spacing=spc,container=nb,label="Database search") #results from a database search
 tabVAL = ggroup(horizontal=FALSE,expand=TRUE,spacing=spc,container=nb,label="Model validation") #modelvalidation from MLE
 tabLRMIX <- ggroup(horizontal=FALSE,expand=TRUE,spacing=spc,container=nb,label="Qual. LR") #LRmix
 svalue(nb) <- 2 #initial start at second tab


####################################################
###############Tab 1: Create Data:##################
####################################################

 tabGENtmp <- glayout(spacing=0,container=tabGEN)

 refreshTabGEN = function(thlist=list(mu=1000,sigma=0.15,xi=0.1,mx=NULL) ) { #can be repeated
  visible(mainwin) <- FALSE 

  #load/save helpfunctions for generated data
  f_openprof = function(h,...) {
    proffile = gfile(text="Open profile",type="open",filter=list("text"=list(patterns=list("*.txt","*.csv","*.tab")),"all"=list(patterns=list("*"))))
    if(!is.na(proffile)) {
     Data = tableReader(proffile) #load profile
     print(Data)
     setDataToGUI(sample_tableToList(Data)[[1]] ,h$action) #convert from table to list and load into GUI
    }
  }
  f_saveprof = function(h,...) {
    Data <- list(getDataFromGUI(h$action)) #get data from GUI
    if(h$action==0) names(Data) <- "generated"
    if(h$action>0) names(Data) <- paste0("ref",h$action)
    Data = sample_listToTable(Data) #convert from table to list
    saveTable(Data,"csv")
  }

  #helpfunction to set Data to GUI
  setDataToGUI <- function(Data,type) {
    dloc <- names(Data) #locus in Data
    for(i in 1:length(locs)) {
     dind <- grep(toupper(locs[i]),toupper(dloc)) #get dataindex
     if(type==0) { #if evidence
       svalue(tabGENb2[i,1]) <- svalue(tabGENb2[i,2]) <- ""
       if(length(dind)>0) { #if locus found
         if(!is.null(Data[[dind]]$adata)) svalue(tabGENb2[i,1]) <- paste0(Data[[dind]]$adata,collapse=",") #insert alleles
         if(!is.null(Data[[dind]]$hdata)) svalue(tabGENb2[i,2]) <- paste0(Data[[dind]]$hdata,collapse=",") #insert alleles
       }
     } else { #if reference
       svalue(tabGENb3[i,type]) <- ""
       if( length(dind)>0 && !is.null(Data[[dind]]$adata) ) svalue(tabGENb3[i,type]) <- paste0(Data[[dind]]$adata,collapse=",") #insert
     }
   } #end for each locus
  } #end function

  #helpfunction to get data from GUI
  getDataFromGUI <- function(type) {
    outloc <- locs #store locs
    Data <- list()
    for(i in 1:length(locs)) {
     outloc[i] <- svalue(tabGENb1[i,1]) #get new loc-names
     if(type==0) { #store evidence
      Data[[outloc[i]]] <- list( adata=unlist(strsplit(svalue(tabGENb2[i,1]),",")) , hdata=as.numeric(unlist(strsplit(svalue(tabGENb2[i,2]),","))) )
     } else { #store reference
      Data[[outloc[i]]] <- list( adata=unlist(strsplit(svalue(tabGENb3[i,type]),",")) )
     }    
    } #end for each locus
    return(Data)
  }

  #layout: 
  tabGENa = glayout(spacing=0,container=(tabGENtmp[1,1] <-gframe("Parameters",container=tabGENtmp))) 
  tabGENb = glayout(spacing=0,container=(tabGENtmp[2,1] <-gframe("Edit",container=tabGENtmp))) 
  tabGENc = glayout(spacing=3,container=(tabGENtmp[3,1] <-gframe("Import/Export profile",container=tabGENtmp)))  
  tabGENd = glayout(spacing=3,container=(tabGENtmp[4,1] <-gframe("Further action",container=tabGENtmp))) 

  #number of contributors
  set <- get("setGEN",envir=mmTK) #get stored setup
  par <- set$param
  nC <- set$model$nC_hd #number of contributors
 
  #default values
  if(is.null(thlist$mx)) thlist$mx <- format((nC:1)/sum(nC:1),3) #default value


  #user input:
  tabGENa[1,1] <- glabel("mu (amount of dna)",container=tabGENa)
  tabGENa[1,2] <- gedit(thlist$mu,container=tabGENa,width=10)
  tabGENa[2,1] <- glabel("sigma (coeffecient of variation)",container=tabGENa)
  tabGENa[2,2] <- gedit(thlist$sigma,container=tabGENa,width=10)
  tabGENa[3,1] <- glabel("xi (stutter ratio)",container=tabGENa)
  tabGENa[3,2] <- gedit(thlist$xi,container=tabGENa,width=10)
  for(k in 1:nC) { #for each contributors
   tabGENa[k+3,1] <- glabel( paste0("mx",k," (mix-proportion contr. ",k,")") ,container=tabGENa)
   tabGENa[k+3,2] <- gedit(thlist$mx[k],container=tabGENa,width=10)
  }

  #Generate:
  simdata <- genDataset(nC, popFreq=set$popFreq, mu=thlist$mu, sigma=thlist$sigma, sorted=FALSE,threshT=par$threshT, refData=set$refData, mx=thlist$mx/sum(thlist$mx),nrep=1, stutt = thlist$xi, prC=par$prC, lambda=par$lambda)

  #insert data in GUI
  mixData <- simdata$samples[[1]]
  refData <- simdata$refData
  locs <- names(mixData) #get locus names

  #show Loci
  tabGENb1 <-  glayout(spacing=0,container=(tabGENb[1,1] <-gframe("Loci",container=tabGENb))) 
  for(i in 1:length(locs))  tabGENb1[i,1] = gedit(locs[i],container=tabGENb1,width=nchar(locs[i])+3)

  #show allele,heights
  tabGENb2 <-  glayout(spacing=0,container=(tabGENb[1,2] <-gframe("Evidence (allele,heights)",container=tabGENb))) 
  for(i in 1:length(locs)) {
   adata <- mixData[[locs[i]]]$adata
   hdata <- round(mixData[[locs[i]]]$hdata)
   tabGENb2[i,1] = gedit(paste0(adata,collapse=","),container=tabGENb2,width=sum(nchar(adata)) + length(adata))
   tabGENb2[i,2] = gedit(paste0(hdata,collapse=","),container=tabGENb2,width=sum(nchar(hdata)) + length(hdata))
  }

  #show references:
  tabGENb3 <-  glayout(spacing=0,container=(tabGENb[1,3] <-gframe("Reference(s)",container=tabGENb))) 
  for(k in 1:nC) {
   for(i in 1:length(locs)) {
    adata <- refData[[locs[i]]][[k]]
    tabGENb3[i,k] = gedit(paste0(adata,collapse=","),container=tabGENb3,width=sum(nchar(adata)) + length(adata))
   }
  }

  #storage buttons:
  tabGENc[1,1] <- gbutton(text="Store evidence",container=tabGENc,handler=f_saveprof,action=0)
  for(k in 1:nC) tabGENc[1,1+k] <- gbutton(text=paste0("Store ref",k),container=tabGENc,handler=f_saveprof,action=k)
  tabGENc[2,1] <- gbutton(text="Load evidence",container=tabGENc,handler=f_openprof,action=0)
  for(k in 1:nC) tabGENc[2,1+k] <- gbutton(text=paste0("Load ref",k),container=tabGENc,handler=f_openprof,action=k)

  #further action
  tabGENd[1,1] <- gbutton(text="Generate again",container=tabGENd,handler=function(h,...) { 
   mx <- rep(0,nC)
   for(k in 1:nC)  mx[k] <- as.numeric(svalue(tabGENa[3+k,2])) 
   refreshTabGEN(list(mu=as.numeric(svalue(tabGENa[1,2])),sigma=as.numeric(svalue(tabGENa[2,2])),xi=as.numeric(svalue(tabGENa[3,2])),mx=mx/sum(mx)) )
   })

  tabGENd[2,1] <- gbutton(text="Plot EPG",container=tabGENd,handler=function(h,...) {
   kit <- get("selPopKitName",envir=mmTK)[1] #get
   if(is.na(kit)) return()
   Data <- getDataFromGUI(0) #get data from GUI
   plotEPG(Data,kitname=kit,sname="edited") #plot epg's
   focus(mainwin)
  })

  visible(mainwin) <- TRUE #INCREASE SPEED
  focus(mainwin)
 } #end refreshTabGE

####################################################
###############Tab 2: Import Data:##################
####################################################

 #When program starts, import assumed model for EVID.


 #a) button for loading kits from directory:
 f_loadkd = function(h,...) { loadKitList(freqpath=get("freqfolder",envir=mmTK)); }

 #b) load/save profiles/database: Supports any filetype
 f_importprof = function(h,...) {
  type=h$action #get type of profile
#  proffile = gfile(text=paste("Open ",type,"-file",sep=""),type="open",filter=list("text"=list(patterns=list("*.txt","*.csv","*.tab"))))
  proffile = gfile(text=paste("Open ",type,"-file",sep=""),type="open",filter=list("text"=list(patterns=list("*.txt","*.csv","*.tab")),"all"=list(patterns=list("*"))))
  if(!is.na(proffile)) {
   Data = tableReader(proffile) #load profile
   print("Raw fil import:")
   print(Data[1:min(nrow(Data),100),]) #print raw-input data
  ###################
  ##DATABASE IMPORT##
  ###################
   if(type=="db") { #importing database file
    popFreq <- get("popFreq",envir=mmTK) 
    if(is.null(popFreq)) {
     gmessage("Population frequencies needs to be imported for database search",title="Error")
    } else {
     minFreq <- getminFreq()
     #saving MEMORY by convert database here!
     #Data <- dbToGenoTable(Data) #convert database table to genotype matrix

     #Assumption: No reference has same samplename. All samplenames are rowed sequentially
     cn = colnames(Data) #colnames 
     lind = grep("marker",tolower(cn),fixed=TRUE) #locus col-ind
     sind = grep("sample",tolower(cn),fixed=TRUE) #sample col-ind
     aind = grep("allele",tolower(cn),fixed=TRUE) #allele col-ind
     aind <- aind[1:2] #assume only 2 possible alles in reference profiles
     locsDB <- unique(Data[,lind]) #unique markers in DB
     locsPop <- toupper(names(popFreq)) #markers in population
     sn <- unique(Data[,sind]) #unique samples in DB
     newData <- numeric() #encoded reference table
     dbDatalocs <- numeric() #locus in newData
     for(loc in locsDB) { #for each marker in DB:
      locind <- grep(toupper(loc),toupper(names(popFreq)),fixed=TRUE) #get position in popFreq
      if(length(locind)==0) next #if locus not existing in popFreq the locus in not imported
      newCol <- rep(NA,length(sn)) #new column in newData
      subA <- Data[Data[,lind]==loc,aind] #extract allele data with matching marker
      subS <- Data[Data[,lind]==loc,sind] #extract sample names with matching marker
      isHom <- which(subA[,2]=="" | is.na(subA[,2])) #if homozygote assignation:
      if(length(isHom)>0) subA[isHom,2] <- subA[isHom,1] #copy first allele
      okSind <- which(sn%in%subS) #samples to consider for particular marker
      if(length(okSind)==0) next #if no samples exists
      print(paste0("Import data from locus: ",loc)) 
      newCol[okSind] <- 1 #init existing as 1. NA for missing allele-info
      doneEncode <- matrix(FALSE,ncol=2,nrow=length(okSind)) #matrix checking which we are finished with encoding a genotype
      Afreqs <- names(popFreq[[locind]]) #get allele-names. Update for each column
      for(k in 1:length(aind)) { #for each allele-column
       for(j in 1:length(Afreqs)) { #for each unique allele in popFreq:
        okAind <- which(subA[,k]==Afreqs[j]) #find matching alleles in subA[okSind]
        if(length(okAind)==0) next
        doneEncode[okAind,k] = TRUE #check that is finished
        newCol[okSind][okAind] <- newCol[okSind][okAind]*prim[j] #multiply with primenumber
       } #end for each allele j

       #THREAT NEW ALLELES,MISSTYPOS ETC:
       if(any(!doneEncode[,k])) { #if any not encoded
        newA <- unique(subA[!doneEncode[,k],k]) #get new alleles
        newA <- newA[!is.na(newA)] #important to remove NA's
        if(length(newA)==0) next
        tmp <- popFreq[[locind]]
        popFreq[[locind]] <- c(tmp, rep(minFreq,length(newA)))
        names(popFreq[[locind]]) <- c(names(tmp),newA) #add unique
        print(paste0("In locus ",loc,": Allele(s) ",newA," was inserted with min. frequency ",prettyNum(minFreq)))
        for(j in 1:length(newA)) { #for each unique allele in popFreq:
         okAind <- which(subA[,k]==newA[j]) #find matching alleles in subA[okSind]
         if(length(okAind)==0) next
         newCol[okSind][okAind] <- newCol[okSind][okAind] * prim[ which(names(popFreq[[locind]])==newA[j]) ] #multiply with EXTENDED primenumber
        } #end for each allele j
       } #end if not encoded 
       Afreqs <- names(popFreq[[locind]]) #get allele-names again!
      } #end for each column k
      dbDatalocs <- c(dbDatalocs,toupper(names(popFreq)[locind])) #all locus
      newData <- cbind(newData,newCol) #add column
     } #end for each locus loc
   
     #RESCALE popFreq?
     for(i in 1:length(popFreq)) {
      popFreq[[i]] <- popFreq[[i]]/sum(popFreq[[i]])
     }
     print("Frequencies was normalized")
     colnames(newData) <- dbDatalocs
     rownames(newData) <- sn #insert sample names
     print(paste0("Database successfully imported with ",nrow(newData)," samples."))

     tmp <- unlist(strsplit(proffile,"/",fixed=TRUE)) #just label the file
     fname <- tmp[length(tmp)] #get filename
     fname <- unlist(strsplit(fname,"\\."))[1]

     dbData <- get("dbData",envir=mmTK) #get already existing databases
     if(is.null(dbData)) dbData <- list() #create list-object
     dbData[[fname]] <- newData #add database to list

     assign("dbData",dbData,envir=mmTK) #store matrix in environment for later use
     assign("popFreq",popFreq,envir=mmTK) #assign updated popFreq
     tabimportB[2,3][] <- names(dbData)
    } #end if popFreq exist
   } else { 
    Data = sample_tableToList(Data) #convert from table to list 
    #get already stored data:
    if(type=="mix") Data2 <- getData("mix") #get data from mmTK-environment
    if(type=="ref") Data2 <- getData("ref") #get data from mmTK-environment

    if(is.null(Data2)) { #if no previous already there
     Data2 <- Data
    } else {
     for(k in 1:length(Data)) { #for each profile
      Data2[[names(Data)[k]]] <- Data[[k]] #insert dataframe
     }
    }
    if(type=="mix")  assign("mixData",Data2,envir=mmTK) #assign data to mmTK-environment
    if(type=="ref")  assign("refData",Data2,envir=mmTK) #assign data to mmTK-environment
    if(type=="mix")  tabimportB[2,1][] <- names(Data2)
    if(type=="ref")  tabimportB[2,2][] <- names(Data2)
   }
  }
 }

 #c) 
 #prints evidence, references, EPG, databases and population frequencies
 f_viewdata = function(h,...) {
  evidD <- getData("mix") #get selected data
  mixSel  <- numeric()
  if(!is.null(evidD)) mixSel <- svalue(tabimportB[2,1])  #get selected mixtures

  if(h$action=="freq") { #prints EPG
   popFreq <- get("popFreq",envir=mmTK) #get frequencies
   if(is.null(popFreq)) {
    tkmessageBox(message="Please import and select population frequencies!")
    return
   } else {
    locs <- names(popFreq)
    nL <- length(locs)
    unAchr <- unique(unlist(sapply( popFreq,names) )) #get character alleles
    ord <- order(as.numeric(unAchr)) 
    unAchr <- unAchr[ord]  #make increasing order
    outD <- unAchr

    matsi <- matrix(NA,nrow=nL,ncol=length(mixSel))  #vector with summed alleles for each marker
    rownames(matsi) <- locs
    colnames(matsi) <- mixSel
    for(i in 1:nL) {
     freqs <- popFreq[[i]] #get frequencies
     newRow <- rep(NA,length(unAchr))
     for(j in 1:length(freqs)) {
      rowind <- which(unAchr==names(freqs[j] ))
      newRow[rowind] <- freqs[j]
     }
     outD <- cbind(outD,newRow)
     for(msel in mixSel) {
      adata <- evidD[[msel]][[locs[i]]]$adata #selected samples
      if(length(adata)>0) matsi[i,which(msel==mixSel)] <- sum(freqs[names(freqs)%in%adata])
     } #end for each samples
    } #end for each locus

    #plot table with frequncies
    colnames(outD) = c("Allele",locs) 
    dbwin <- gwindow("Population frequencies", visible=FALSE,height=mwH)
    gtable(outD ,container=dbwin) #create table
    visible(dbwin) <- TRUE

    #Calculate random man probability of matching for each samples
    for(msel in mixSel) { 
     cind <- which(msel==mixSel)
     print(paste0("Calculating exact MAC random man distribution for sample ",msel,". This could take a while..."))
     si <- matsi[!is.na(matsi[,cind]),cind]
     macdist <- exactMACdistr(si,2*length(si)-4) #allow down to 4 missmatches
     cumprob <- rev(cumsum(rev(macdist)))
     ymax <- max(cumprob)
     delta <- 0.03
     if(msel!=mixSel[1]) dev.new() #create new plot for next EPG
     barplot(cumprob,main="Random man probability having number of allele matches>=k",xlab="k number of allele matches",space=0,ylim=c(0,ymax+ymax*2*delta))
     text(1:length(cumprob)-0.5,cumprob+ymax*delta,labels=format(cumprob,digits=2))
     mtext(paste0("Sample: ",msel))
    }
   } #end if popFreq
  } #end freq

  if(h$action=="mix") { #prints EPG
     kitname <- svalue(tabimportA[3,1]) #get kit name. Must be same name as in generateEPG
     for(msel in mixSel) {
      subD <- evidD[[msel]] #selected samples
      print("------------------------------------")
      print(paste("Samplename: ",msel,sep=""))
      print(subD)
      if(which(msel==mixSel)>1) dev.new() #create new plot for next EPG
      plotEPG(Data=subD,kitname,sname=msel) #plot epg's
     } 
   focus(mainwin)
  }

  if(h$action=="ref") { #print tables only
     refD <- getData("ref")
     refSel <- numeric()
     if(!is.null(refD))  refSel <- svalue(tabimportB[2,2])  #get selected references
     nR <- length(refSel) #number of selected references

     #first: Print references:
     for(rsel in refSel) {
      print("------------------------------------")
      print(paste("Referencename: ",rsel,sep=""))
      tab <- matrix(unlist(refD[[rsel]]),ncol=2,byrow=TRUE)
      ref <- paste0(tab[,1],"/",tab[,2])
      names(ref) <- names(refD[[rsel]])
      print(ref)
     }
   
     #second: Print #matches to selected samples
     for(msel in mixSel) { #for each selected 
      print("------------------------------------")
      print(paste0("Number of matching alleles with samplename ",msel,":"))
      subD <- evidD[[msel]] #selected samples
      locs <- names(subD)
      tab <- matrix(NA,ncol=nR,nrow=length(locs))
      rownames(tab) <- locs
      colnames(tab) <- refSel 
      for(loc in  locs) { #for each locus
        for(rsel in refSel) {
          refA <- refD[[rsel]][[loc]]$adata
          if(!is.null(refA)) tab[which(loc==locs),which(rsel==refSel)] <- sum(refA%in%subD[[loc]]$adata)
        }
      } 
      print(t(tab))
     } #end for each mix    
   } #end if references

   if(h$action=="db") {  #view imported databases
    popFreq <- get("popFreq",envir=mmTK)
    if(length(tabimportB[2,3][])==0) return();
    dbSel <- svalue(tabimportB[2,3])  #get selected database
    for(dsel in dbSel) { #for each selected databases
     subD <- getData("db",dsel)[[dsel]] #get selected data
     dblocs <- toupper(colnames(subD)) #get database locs
     outD <- rownames(subD) #will be character
     macD <- matrix(0,nrow=length(outD),ncol=length(mixSel)) #matching allele counter for each reference
     nlocs <- rep(0,length(outD))
     for(i in 1:length(dblocs)) { #for each locus in db     
      Ainfo <- names(unlist(popFreq[[dblocs[i]]])) #extract allele-info
      #translate to genotypes
      Pinfo <- prim[1:length(Ainfo)]
      G = t(as.matrix(expand.grid(rep(list(Ainfo,Ainfo )))))
      GP = t(as.matrix(expand.grid(rep(list(Pinfo,Pinfo )))))
      keep = GP[2,]>=GP[1,] #unique genotypes
      G0 <- G[,keep]  #store genotypes

      G <- paste0(G0[1,],paste0("/",G0[2,]))
      GP <- GP[,keep]  #store genotypes
      GP <- GP[1,]*GP[2,]
      newRow <- rep(NA,nrow(subD))  
      tmpmac <- matrix(0,nrow=length(GP),ncol=length(mixSel)) #genotype match for each samples
      for(j in 1:length(GP)) { #for each genotype
       #Always: Find corresponding genotype name by looking up stored primenumbers (same order as in popFreq!)
       rowind <- which(subD[,i]==GP[j]) #samples with this genotype
       if(length(rowind)==0) next
       newRow[rowind] <- G[j]
       #Optional (if mixtures selected): Count matching alleles:
       hasloc <- FALSE #total number locus checked over evidence
       for(msel in mixSel) { #for each selected mixture
        evid0 <- evidD[[msel]][[dblocs[i]]]$adata
        if(!is.null(evid0)) hasloc <- TRUE
        tmpmac[j,which(msel==mixSel)] <- sum(G0[,j]%in%evid0)
       } 
       if(ncol(macD)>0) { #if compared with evidence
        macD[rowind,] = macD[rowind,] + tmpmac[j,] #add match counter 
        nlocs[rowind] <- nlocs[rowind] + as.integer(hasloc)
       } #end if evidence had locus
      } #end for each genotype
      outD <- cbind(outD,newRow) #extend
     } #end for each locus
     colnames(outD) <- c("Reference",dblocs)
     nmax <- min(nrow(outD), get("optDB",envir=mmTK)$maxDB) #max limit of number in showing dropdown
     if(nrow(outD)>nmax ) print(paste0("Database contained more than ",nmax," samples. Showing first ",nmax," samples only!"))


     #if selected Mixtures: show ranked matches in database-reference:
     if(ncol(macD)>0) { #
      ord <- order(rowSums(macD),decreasing=TRUE)
      outD2 <- cbind(outD[ord,1],macD[ord,],nlocs)
      colnames(outD2) <- c("Reference",mixSel,"nLocs") 
      dbwin2 <- gwindow(paste0("Number of sample matching alleles in references in database ",dsel), visible=FALSE,height=mwH)
      if(nmax<=1) gtable(outD2,container=dbwin2) #create table
      if(nmax>1) gtable(outD2[1:nmax,],container=dbwin2) #create table
      visible(dbwin2) <- TRUE
     }
     dbwin <- gwindow(paste0("References in imported database ",dsel), visible=FALSE,height=mwH)
     if(nmax<=1) gtable(outD,container=dbwin) #create table
     if(nmax>1) gtable(outD[1:nmax,],container=dbwin) #create table
     visible(dbwin) <- TRUE        

    } #end for each databases
   } #end if db -case
 }  #end viewdata

 ###############
 #start layout:#
 ###############
 tabimportA = glayout(spacing=5,container=gframe("Step 1) Import and select Population frequencies",container=tabimport)) #kit and population selecter
 tabimportA2 = glabel("",container=tabimport) #evidence,ref dataframe
 tabimportB = glayout(spacing=5,container=gframe("Step 2) Import and select Evidence, Reference, Database",container=tabimport)) #evidence,ref dataframe
 tabimportB2 = glabel("",container=tabimport) #evidence,ref dataframe
 tabimportC = glayout(spacing=20,container=gframe("Step 3) Select Interpretation",container=tabimport)) #Tasks button

 #Choose box and import button
 tabimportA[1,1] = gbutton(text="1) Select directory\n(with frequency files)",container=tabimportA,handler=
  function(h,...) {
   dirfile = gfile(text="Select folder",type="selectdir")
   if(!is.na(dirfile)) assign("freqfolder",dirfile,envir=mmTK) #assign freqfolder
  }
 )
 tabimportA[1,2] = gbutton(text="2) Import from directory\n(with frequency files)",container=tabimportA,handler=
  function(h,...) {
   loadKitList(freqpath=get("freqfolder",envir=mmTK))
   kitList <- get("kits",envir=mmTK)
   tabimportA[3,1][] <- names(kitList)
  }
 ) 
 tabimportA[2,1] <- glabel(text="Select kit:",container=tabimportA)
 tabimportA[2,2] <- glabel(text="Select population:",container=tabimportA)

 tabimportA[3,1] <- gcombobox(items="", width=100, selected =0 , editable = FALSE, container = tabimportA, handler=
    function(h,...) {
     if(!is.null(get("dbData",envir=mmTK))) {
      print("You can not change selected kit after loading a database!") 
     } else {
      kitList <- get("kits",envir=mmTK)
      if(!is.null(kitList)) tabimportA[3,2][] <- names(kitList[[svalue(tabimportA[3,1])]])
     }
    })
 #population-selection
 tabimportA[3,2] <- gcombobox(items="", width=100, selected = 0 , editable = FALSE , container = tabimportA, handler=
    function(h,...) {
     if(!is.null(get("dbData",envir=mmTK))) {
      print("You can not change selected population after loading a database!") 
     } else {
      kitList <- get("kits",envir=mmTK)
      if(!is.null(kitList)) {
       popList <- kitList[[svalue(tabimportA[3,1])]][[svalue(tabimportA[3,2])]] #get selected frequencies
       assign("popFreq",popList,envir=mmTK) #assign popFreq get("popFreq",envir=mmTK)
       popkitname <- c(svalue(tabimportA[3,1]),svalue(tabimportA[3,2])) #store selected kit and popnames
       assign("selPopKitName",popkitname,envir=mmTK) 
      }
     }
    })
#previous stored kit/pop-names:
 popkitname <- get("selPopKitName",envir=mmTK) 
 if(!is.null(popkitname)) {
  tabimportA[3,1][] <- popkitname[1]
  tabimportA[3,2][] <- popkitname[2]
 }
 tabimportA[2,3] <-  gbutton(text="View frequencies",container=tabimportA,handler=f_viewdata,action="freq")  #view popFreq-data

 #Choose box and import button
 tabimportB[1,1] = gbutton(text="Import evidence",container=tabimportB,handler=f_importprof,action="mix")
 tabimportB[2,1] = gcheckboxgroup(items="", container = tabimportB)
 tabimportB[2,1][] = restoreCheckbox("mix")

 #Choose box and import button
 tabimportB[1,2] = gbutton(text="Import reference",container=tabimportB,handler=f_importprof,action="ref")
 tabimportB[2,2] = gcheckboxgroup(items="", container = tabimportB)
 tabimportB[2,2][] = restoreCheckbox("ref")

 #Choose box and import button
 tabimportB[1,3] = gbutton(text="Import database",container=tabimportB,handler=f_importprof,action="db")
 tabimportB[2,3] = gcheckboxgroup(items="", container = tabimportB)
 tabimportB[2,3][] = restoreCheckbox("db")

 #view data:
 tabimportB[3,1] = gbutton(text="View evidence",container=tabimportB,handler=f_viewdata,action="mix")
 tabimportB[3,2] = gbutton(text="View references",container=tabimportB,handler=f_viewdata,action="ref")
 tabimportB[3,3] = gbutton(text="View database",container=tabimportB,handler=f_viewdata,action="db")

 #helpfunction used to extract selected importdata-elements to further model-setup
 selectDataToModel <- function(h,....) {
   #All: popFreq must be imported!
   #GEN: No other requirement
   #EVID: must have both mixture and reference profiles
   #DB: Database search require both mixture and database, reference profiles is optional
   #DC: Deconvolution requires only mixtures. Reference profiles is optional
   popFreq <- get("popFreq",envir=mmTK)
   mixSel <- refSel <- dbSel <- numeric()
   if(length(tabimportB[2,1][])>0) mixSel <- svalue(tabimportB[2,1])  #get selected mixtures
   if(length(tabimportB[2,2][])>0) refSel <- svalue(tabimportB[2,2])  #get selected references
   if(length(tabimportB[2,3][])>0) dbSel <- svalue(tabimportB[2,3])  #get selected databases

   if(is.null(popFreq)) {
    gmessage("No frequencies was specified!\n Please import table with population frequencies.")
   } else if(h$action!="GEN" && length(mixSel)==0) {
    gmessage(message="Please import and select mixture-profile!")
   } else if(h$action=="EVID" && length(refSel)==0) {
    gmessage(message="Please import and select reference-profiles for weight of evidence!")
   } else if(h$action=="DB" && length(dbSel)==0) {
    gmessage(message="Please import and select reference-database for database search!")
   } else {
#    delete(tabmodeltmp,tabmodel)  
#    dispose(tabmodel)  
    refreshTabmodel(mixSel,refSel,dbSel,h$action) #refresh table with selected data
    svalue(nb) <- 3 #change tab of notebook
   }
 }

 #Button-choices further:
 tabimportC[1,1] = gbutton(text="Generate sample",container=tabimportC,handler=selectDataToModel,action="GEN")
 tabimportC[1,2] = gbutton(text="Deconvolution",container=tabimportC,handler=selectDataToModel,action="DC")
 tabimportC[1,3] = gbutton(text="Weight-of-Evidence",container=tabimportC,handler=selectDataToModel,action="EVID")
 tabimportC[1,4] = gbutton(text="Database search",container=tabimportC,handler=selectDataToModel,action="DB")


####################################################################################################################
#######################################Tab 3: Model setup:##########################################################
#####################################################################################################################

tabmodeltmp <- glayout(spacing=spc,container=tabmodel)


  #helpfunction to get boundary from Toolbar:
  getboundary = function(nC,xi=NULL) {
    optint <- get("optINT",envir=mmTK) #options when integrating (reltol and boundaries)
    lower <- rep(0,nC+2)
    upper <- rep(1,nC+2)
    upper[nC] <- optint$maxmu
    upper[nC+1] <- optint$maxsigma
    if(!is.null(xi)) { #must remove stutter ratio if known
     lower <- lower[-nC+2]
     upper <- lower[-nC+2]
    }
    return(bound=list(lower=lower,upper=upper))
  }


  ##EVID INTEGRATION (can be done anywhere after model setup)##
  doINT <- function(type) { #Used either with EVID or DB
     #sig = number of decimals
     optint <- get("optINT",envir=mmTK) #options when integrating (reltol and boundaries)
     if(type=="EVID") {
       set <- get("setEVID",envir=mmTK)
       par <- set$param
       mod <- set$model
       print(paste0("Calculating integrals with relative error ",optint$reltol))
       print("This may take a while...")   
       bhp <- getboundary(mod$nC_hp,par$xi) #get boundaries under hp
       bhd <- getboundary(mod$nC_hd,par$xi) #get boundaries under hd

       #integrate:
       int_hd <- contLikINT(mod$nC_hd, set$samples, set$popFreqQ, bhd$lower, bhd$upper, set$refDataQ, mod$condOrder_hd, mod$knownRef_hd, par$xi, par$prC, optint$reltol, par$threshT, par$fst, par$lambda, par$pXi) 
       print("Calculation under Hd done...")

       int_hp <- contLikINT(mod$nC_hp, set$samples, set$popFreqQ, bhp$lower, bhp$upper, set$refDataQ, mod$condOrder_hp, mod$knownRef_hp, par$xi, par$prC, optint$reltol, par$threshT, par$fst, par$lambda, par$pXi) 
       LR <- int_hp$margL/int_hd$margL
       dev <- range(c(int_hp$deviation/int_hd$deviation,int_hp$deviation/rev(int_hd$deviation))) #get deviation interval of LR
       res <- list(LR=LR,dev=dev)
       assign("resEVIDINT",res,envir=mmTK) #assign deconvolved result to environment
       #Print a GUI message:
       gmessage(message=paste0("The LR (integrated Likelihood based)\nwas calculated as \n",format(LR,digits=4)," [",format(dev[1],digits=4)," , ",format(dev[2],digits=4),"]"),title="Continuous LR (Integrated Likelihood based)",icon="info")
     } 
     if(type=="DB") { #Case of DB-search
       doDB("INT") #do database search with integration
     }
  } #end Integration


  refreshTabmodel = function(mixSel,refSel,dbSel,type) { 
   #type={"GEN","EVID","DB","DC"}
   visible(mainwin) <- FALSE
   popFreq <- get("popFreq",envir=mmTK)
   locs <- names(popFreq)  #get names of loci for imported population frequencies. used as default in list
   mixD = getData("mix")
   refD = getData("ref") 
#mixSel = names(mixD)
#refSel = names(refD)
#type="EVID"
   nM = length(mixSel) #number of mix-profiles
   nR = length(refSel) #number of ref-profiles

   tabmodelA = glayout(spacing=5,container=(tabmodeltmp[1,1] <-gframe("Model specification",container=tabmodeltmp))) 
   tabmodelB = glayout(spacing=0,container=(tabmodeltmp[1,2] <-gframe("Data selection",container=tabmodeltmp))) 
   tabmodelC = glayout(spacing=0,container=(tabmodeltmp[2,1] <-gframe("Show selected data",container=tabmodeltmp)))  
   tabmodelD = glayout(spacing=5,container=(tabmodeltmp[2,2] <-gframe("Calculations",container=tabmodeltmp)))  

   #Hypothesis selection: subframe of A
   txt <- "Contributor(s) under Hp:"
   if(type=="DB") txt <- paste0(txt, "\n(DB-reference already included)")
   if(type!="GEN") { #if not generating samples
    tabmodelA1 = glayout(spacing=0,container=(tabmodelA[1,1] <-gframe("Evidence(s)",container=tabmodelA))) 
    if(type!="DC") tabmodelA2 = glayout(spacing=0,container=(tabmodelA[2,1] <-gframe(txt,container=tabmodelA))) 
   } 
   tabmodelA3 = glayout(spacing=0,container=(tabmodelA[3,1] <-gframe("Contributor(s) under Hd:",container=tabmodelA)))
   tabmodelA4a = glayout(spacing=0,container=(tabmodelA[4,1] <-gframe("Continuous Model Parameters",container=tabmodelA))) 
   if(type=="EVID") tabmodelA4b = glayout(spacing=0,container=(tabmodelA[5,1] <-gframe("Qualitative Model Parameters",container=tabmodelA))) 
   tabmodelA5 = glayout(spacing=0,container=(tabmodelA[6,1] <-gframe("Advanced Parameters",container=tabmodelA))) 

   if(type=="DB") { #add database-names if included:
    tabmodelA6 = glayout(spacing=0,container=(tabmodelA[7,1] <-gframe("Database(s) to search",container=tabmodelA))) 
    for(dsel in dbSel) tabmodelA6[which(dsel==dbSel),1] =  glabel(text=dsel,container=tabmodelA6)
   }
   if(type!="GEN") { #don't show evidence if GEN
    for(msel in mixSel) {
     tabmodelA1[which(msel==mixSel),1] <- gcheckbox(text=msel,container=tabmodelA1,checked=TRUE)
     enabled(tabmodelA1[which(msel==mixSel),1]) <- FALSE
    }
   }
   for(rsel in refSel) {
    if(!type%in%c("DC","GEN")) tabmodelA2[which(rsel==refSel),1]  <- gcheckbox(text=rsel,container=tabmodelA2,checked=TRUE) #Check as default under Hp
    tabmodelA3[which(rsel==refSel),1]  <- gcheckbox(text=rsel,container=tabmodelA3,checked=!(type=="EVID")) #references under Hd (not checked if evidnece)
   }
   #specify number of unknowns
   if(!type%in%c("DC","GEN")) {
    tabmodelA2[nR+1,1] <- glabel(text="#unknowns (Hp): ",container=tabmodelA2)
    tabmodelA2[nR+1,2] <- gedit(text="1",container=tabmodelA2,width=4)
   }
   tabmodelA3[nR+1,1] <- glabel(text="#unknowns (Hd): ",container=tabmodelA3)
   tabmodelA3[nR+1,2] <- gedit(text="2",container=tabmodelA3,width=4)

   #Model parameters for continuous model:
   tabmodelA4a[1,1] <- glabel(text="Probability of Dropin: ",container=tabmodelA4a)
   tabmodelA4a[1,2] <- gedit(text="0",container=tabmodelA4a,width=4)
   tabmodelA4a[2,1] <- glabel(text="fst-correction: ",container=tabmodelA4a)
   tabmodelA4a[2,2] <- gedit(text="0",container=tabmodelA4a,width=4)

   #Model parameters for qualitative model:
   if(type=="EVID") {
    tabmodelA4b[1,1] <- glabel(text="Probability of Dropin: ",container=tabmodelA4b)
    tabmodelA4b[1,2] <- gedit(text="0.05",container=tabmodelA4b,width=4)
    tabmodelA4b[2,1] <- glabel(text="fst-correction: ",container=tabmodelA4b)
    tabmodelA4b[2,2] <- gedit(text="0.02",container=tabmodelA4b,width=4)
   }

   #Advanced parameters:
   tabmodelA5[1,1] <- gcheckbox(text="Q-assignation",container=tabmodelA5,checked=TRUE,horisontal=TRUE) #checked only if not generating
   tabmodelA5[2,1] <- glabel(text="Detection threshold: ",container=tabmodelA5)
   tabmodelA5[2,2] <- gedit(text="150",container=tabmodelA5,width=4)
   tabmodelA5[3,1] <- glabel(text="Stutter ratio (xi): ",container=tabmodelA5)
   tabmodelA5[3,2] <- gedit(text="",container=tabmodelA5,width=4)
   tabmodelA5[4,1] <- glabel(text="Dropin peak height \n hyperparam (lambda):",container=tabmodelA5)
   tabmodelA5[4,2] <- gedit(text="0",container=tabmodelA5,width=4)

   tabmodelA5[5,1] <- glabel(text="Stutter prior function p(x): ",container=tabmodelA5)
   tabmodelA5[6,1] <- gedit(text="1",container=tabmodelA5,width=20)
   enabled(tabmodelA5[6,1]) <- FALSE #NOT IMPLEMENTED

   if(type=="GEN") { #deactivate options for generation:
    enabled(tabmodelA4a[2,2]) <- FALSE #deactivate fst-correction
    enabled(tabmodelA5[1,1]) <- FALSE #deactivate Q-assignation fst-correction
   }

   #Data selection
   tabmodelB[1,1] <- glabel(text="Loci:",container=tabmodelB)
   for(loc in locs) { #insert locus names from popFreq
    tabmodelB[1+which(loc==locs),1] <- loc  #insert loc-name
   }
   for(msel in mixSel) { #for each selected mixture
    tabmodelB[1,1 + which(msel==mixSel)] <- glabel(text=msel,container=tabmodelB) #get selected mixturenames
    for(loc in locs) { #for each locus
     exist <- !is.null(mixD[[msel]][[loc]]$adata) ##&& !is.null(mixD[[msel]][[loc]]$hdata) #check if exist alleles!
     tabmodelB[1+which(loc==locs),1 + which(msel==mixSel)]  <- gcheckbox(text="",container=tabmodelB,checked=exist)
     if(!exist) enabled(tabmodelB[1+which(loc==locs),1 + which(msel==mixSel)]) <- FALSE #deactivate non-existing locus
    }
   }  
   for(rsel in refSel) { #for each selected reference
    tabmodelB[1,1 + nM + which(rsel==refSel)] <- glabel(text=rsel,container=tabmodelB) #name of reference
    for(loc in locs) { #for each locus
     exist <- !is.null(refD[[rsel]][[loc]]$adata) #check if allele exists
     tabmodelB[1+which(loc==locs),1 + nM + which(rsel==refSel)]  <- gcheckbox(text="",container=tabmodelB,checked=exist)
     if(!exist) enabled(tabmodelB[1+which(loc==locs),1 + nM + which(rsel==refSel)]) <- FALSE #deactivate non-existing locus
    }
   }

   #helpfunction which takes GUI settings and stores them in "set'type'"
   storeSettings <- function(lrtype="PLOT") {
     #lrtype={"CONT"},{"QUAL"},{"PLOT"}
      sellocs <- numeric() #Selected loci (which all mixtures, references has)
      for(loc in locs) { #for each locus in popFreq
       isOK <- TRUE
       for(msel in mixSel) isOK <- isOK && svalue(tabmodelB[1+which(loc==locs),1 + which(msel==mixSel)])  #check if locus checked for samples
       for(rsel in refSel) isOK <- isOK && svalue(tabmodelB[1+which(loc==locs),1 + nM + which(rsel==refSel)]) #check if locus checked for references
       if(isOK) sellocs <- c(sellocs,loc) #locus can be evaluated
      }
      if(length(sellocs)==0) { #don't do anything if no loci will be evaluated
       gmessage(message="No loci are evaluated! Be sure that all selected data have valid data in their loci.",title="No loci found!",icon="error")
       stop("No evaluation done.")
      }
      popFreq <- popFreq[sellocs] #consider only relevant loci in popFreq
      print(c("Locs to be evaluated:",paste0(sellocs,collapse=",")))

      #Check if samples have peak heights if cont. model is used:
      if(lrtype=="CONT") {
        hdatas <- sapply( mixD[mixSel], function(x) sum(sapply(x,function(y) length(y$hdata)) ) ) #number of alleles for each samples
        if(any(hdatas==0)) {
          gmessage(paste0("The sample(s) ", paste0(mixSel[hdatas==0],collapse=",")," did not have peak heights! \nEvaluate with qualitative LR model"),title="Wrong data input!",icon="error")
          stop("No evaluation done.")
        }
      }

      #Insert missing alleles to the popFreq-object:
      minFreq <- getminFreq() #get frequency used 
      for(loc in sellocs) {
       tmp <- popFreq[[loc]] #get allele-names in popFreq
       newA <- numeric()
       for(ss in mixSel) {
        adata <- mixD[[ss]][[loc]]$adata
        newA <- c(newA, adata[!adata%in%names(tmp)])
       }
       for(rr in refSel) {
        adata <- refD[[rr]][[loc]]$adata
        newA <- c(newA, adata[!adata%in%names(tmp)])
       }
       newA <- unique(newA)
       if(length(newA)>0) { #if new allele found
        popFreq[[loc]] <- c(tmp, rep(minFreq,length(newA)))
        popFreq[[loc]] <- popFreq[[loc]]/sum(popFreq[[loc]]) #normalize
        names(popFreq[[loc]]) <- c(names(tmp),newA) #add unique
        print(paste0("In locus ",loc,": Allele(s) ",paste0(newA,collapse=",")," was inserted with min. frequency ",prettyNum(minFreq)))
       }
      } #end for each loci

      #READ FROM GUI:
      threshT <- as.numeric(svalue(tabmodelA5[2,2]))
      xi <- svalue(tabmodelA5[3,2]) #stutter ratio
      lambda <-   as.numeric(svalue(tabmodelA5[4,2])) #lambda is hyperparameter to dropin-peak height model
      pXi = function(x)1  #default
      if(lrtype!="QUAL") prC <-  as.numeric(svalue(tabmodelA4a[1,2])) #dropin probability
      if(lrtype=="QUAL") prC <-  as.numeric(svalue(tabmodelA4b[1,2])) #dropin probability
      if(lrtype!="QUAL") fst <-  as.numeric(svalue(tabmodelA4a[2,2])) #correction
      if(lrtype=="QUAL") fst <-  as.numeric(svalue(tabmodelA4b[2,2])) #correction

      #CHECK VARIABLES:
      checkPositive(threshT,"Threshold")
      checkPositive(lambda,"Dropin-peak height hyperparameter")
      checkProb(prC,"Allele dropin-probability")
      checkProb(fst,"fst-correction")
      if(xi=="") {
       xi <- NULL #assume unknown stutter ratio
      } else {
       checkProb(xi,"Stutter ratio (xi)")
      }

      #prepare data for function in gammadnamix! Data-object must have correct format!
      #go through selected data from GUI:
      samples <- list()
      refData <- NULL
      if(nR>0) refData <- list()
      for(loc in sellocs) { #for each selected locus in popFreq
       for(msel in mixSel) { #for each mixture samples
         subD <- mixD[[msel]][[loc]] #get selected data
         if(lrtype=="CONT" && is.null(subD$hdata)) { #peak height not found!
            gmessage(message=paste0("The evidence ",msel," didn't contain peak heights for locus ",loc,". Please unselect locus"),title="Wrong input",icon="error")
            stop("Unselect loci which does not have peak heights.")
         }
         if(!is.null(subD$hdata)) {
          keep <- subD$hdata>=threshT #allele to keep (above threshold)
          subD$hdata <- subD$hdata[keep]
          subD$adata <- subD$adata[keep]
         }
         samples[[msel]][[loc]] <- subD #insert samples
       }
       if(nR>0) refData[[loc]] <- list()
       for(rsel in refSel) refData[[loc]][[rsel]] <- refD[[rsel]][[loc]]$adata #insert references
      } #end for each locus

      #get specified preposition 
      condOrder_hp <- condOrder_hd <- rep(0,nR)
      if(type=="DC") condOrder_hp <- NULL
      for(rsel in refSel) { #for each reference under hp and hd
        if(!type%in%c("DC","GEN")) {
         valhp <- as.integer(svalue(tabmodelA2[which(rsel==refSel),1])) 
         condOrder_hp[which(rsel==refSel)] <- valhp +  valhp*max(condOrder_hp)
        }
        valhd <- as.integer(svalue(tabmodelA3[which(rsel==refSel),1])) 
        condOrder_hd[which(rsel==refSel)] <- valhd + valhd*max(condOrder_hd)
      }
     
      #number of contributors in model:
      nC_hp <- NULL
      if(!type%in%c("DC","GEN")) nC_hp <-  as.integer(svalue(tabmodelA2[nR+1,2])) + sum(condOrder_hp>0)
      nC_hd <-  as.integer(svalue(tabmodelA3[nR+1,2])) + sum(condOrder_hd>0)

      #get model parameters:
      popFreqQ <- popFreq
      refDataQ <- refData
      if(svalue(tabmodelA5[1,1])) { #if Q-assignation
       ret <- Qassignate(samples, popFreq, refData) #call function in gammadnamix
       popFreqQ <- ret$popFreq
       refDataQ <- ret$refData
      }
  
      #get specified preposition 
      knownref_hp <- knownref_hd <- NULL #known non-contributors under Hp always NULL (since they exist under condOrder)
      if(type=="EVID") { #only for Evidence
       knownref_hd <- which(condOrder_hp>0 & condOrder_hd==0) #those references conditioned under hp but not hd
       if(length(knownref_hd)==0) knownref_hd <- NULL
      }

      #get input to list: note: "fit_hp" and "fit_hd" are list-object from fitted model
      model <- list(nC_hp=nC_hp,nC_hd=nC_hd,condOrder_hp=condOrder_hp,condOrder_hd=condOrder_hd,knownref_hp=knownref_hp,knownref_hd=knownref_hd) #proposition
      param <- list(xi=xi,prC=prC,threshT=threshT,fst=fst,lambda=lambda,pXi=pXi) 
      set <- list(samples=samples,refData=refData,popFreq=popFreq,model=model,param=param,refDataQ=refDataQ,popFreqQ=popFreqQ)     
      if(type=="DB") set$dbSel <- dbSel #store name of selected databases
      assign(paste0("set",type),set,envir=mmTK) #store setup for relevant type
   } #end store settings from GUI to environment


   #Calculation button:  
   #Plot EPG of selected data (calls storeSettings first and then plotEPG by importing selected data)
   if(type!="GEN") {
     tabmodelC[1,1] = gbutton(text="Plot EPG",container=tabmodelC,handler= function(h,...) { 
      storeSettings("PLOT") #store settings
      #loads each selections and plot epg:
      set = get(paste0("set",type),set,envir=mmTK) #store setup for relevant type
      print("Assumed population frequencies:")
      print(set$popFreqQ)
      print("Considered references:")
      print(set$refDataQ)
      print("Considered Evidence samples:")
      print(set$samples)

      #plot EPG:
      kit <- get("selPopKitName",envir=mmTK)[1] #get selected kit
      if(is.na(kit)) return()
      samples <- set$samples
      for(sn in names(samples)) {
       if(which(sn==names(samples))>1) dev.new() #create new plot for next EPG
       plotEPG(samples[[sn]],kitname=kit,sname=sn) #plot epg's
      }
     })
   }

   #Calculation button:  
   if(type=="GEN") {
    tabmodelD[1,1] = gbutton(text="Generate sample",container=tabmodelD,handler= function(h,...) { 
     storeSettings("CONT") #store settings
     refreshTabGEN() #generate a dataset based on stored settings
     svalue(nb) <- 1 #go to data generation-window
    })
   } else {
    tabmodelD[1,1] = gbutton(text="Continuous LR\n(Maximum Likelihood based)",container=tabmodelD,handler=
	function(h,...) {
      storeSettings("CONT") #store settings
      refreshTabMLE(type) #refresh MLE fit tab (i.e. it fits the specified model)
      svalue(nb) <- 4 #go to mle-fit window (for all cases) when finished
    }) #end cont. calculation button
    if(type%in%c("EVID","DB")) tabmodelD[2,1] = gbutton(text="Continuous LR \n(Integrated Likelihood based)",container=tabmodelD,handler=
	function(h,...) {
      storeSettings("CONT") #store model-settings
      doINT(type) #Integrate either for EVID or DB search
    }) #end cont. calculation button
    if(type=="EVID") tabmodelD[3,1] = gbutton(text="Qualitative LR\n(semi-continuous)",container=tabmodelD,handler=
	function(h,...) {
      storeSettings("QUAL") #store model-settings (use other input)
      refreshTabLRMIX() #refresh LRmix calculation tab (i.e. it fits the specified model)
      svalue(nb) <- 8 #go to mle-fit window (for all cases) when finished
    }) #end cont. calculation button
   } #end type-cases
   visible(mainwin) <- TRUE
   focus(mainwin)
  } #end refresh setup tab-frame


#################################################
##########CONT DB-SEARCHING (MLE or INT)#########
#################################################
   doDB <- function(ITYPE="MLE") {  #take a given inference-type {"MLE","INT","QUAL"} #qual means that only qualitative LR is done
     require(forensim)  
     set <- get("setDB",envir=mmTK) #get setup for DB
     opt <- get("optDB",envir=mmTK) #options when doing LRmix
     popFreq <- set$popFreq #get original saved popfreq 
     popFreqQ <- set$popFreqQ #get popFreq saved by setup
     locs_hd <- names(popFreqQ) #get locus analysed under hd
     mixSel <- names(set$samples) #get name of selected mixtures
     mod <- set$model #take out model specifications
     par <- set$param #take out param specifications
     refData <- set$refDataQ #take out selected references         

     if(ITYPE=="MLE") {
        mleobj <- set$mlefit_hd #get object under hd
        mleopt <- get("optMLE",envir=mmTK)
        logLi_hd <- logLiki(mleobj) #get log-probabilities for each locus (under Hd)
     }
     if(ITYPE=="INT") { #Calculate with INT
       optint <- get("optINT",envir=mmTK)
       bhp <- getboundary(mod$nC_hp+1,par$xi) #get boundaries under hp
       bhd <- getboundary(mod$nC_hd,par$xi) #get boundaries under hd
       hd0stored <- list() #A list to store previous calculations
     }

     #LRmix settings:
     nsample <- get("optLRMIX",envir=mmTK)$nsample
     totA <-  sapply(  set$samples, function(x) sum(sapply(x,function(y) length(y$adata)) ) ) #number of alleles for each samples
     refHd <- NULL
     if( any(mod$condOrder_hd>0) ) refHd <- lapply(set$refData ,function(x) x[mod$condOrder_hd]) #must have the original refData!
     pDhat <- rep(NA,length(totA))
     print("Estimating allele dropout probability based on MC method...")
     for(ss in 1:length(totA)) { #for each sample (do dropout estimation) under Hd
       DOdist <- simDOdistr(totA=totA[ss],nC=mod$nC_hd,popFreq=popFreq,refData=refHd,minS=nsample,prC=opt$QUALpC)
       pDhat[ss] <- quantile(DOdist ,0.5) #get quantiles
     }
     print(paste0("Median(s) is given as pDhat=",pDhat))
     pDhat <- round(mean(pDhat),2) #used pDhat in database search
     print(paste0("Mean of the medians over each sample was given as pDhat=",pDhat))
     pDvec <- rep(pDhat,max(mod$nC_hp+1,mod$nC_hd))
     #Dropout estimation finished 

     DBtab <- numeric()  #used to store table when searched
     for(dsel in set$dbSel) { #for each selected databases
       subD <- getData("db",dsel)[[dsel]] #get selected database
       dblocs <- toupper(colnames(subD)) #get database locs
       indD <- rownames(subD) #get individual names in database
       macD <- rep(0,length(indD)) #matching allele counter for each reference
       LRD <- rep(0,length(indD)) #matching allele counter for each reference
       nlocs <- rep(0,length(indD)) #Number of loci which are used for calculating score - Note: Require that reference in DB has a genotype

       #convert allele-names of elements in database to one in popFreq
       for(loc in dblocs) { #for each locus in db     
         if(is.null(popFreq[[loc]])) next #skip to next locus
         Ainfo <- names(unlist(popFreq[[loc]])) #extract allele-info of frequncies
         #translate database to original genotypes
         Pinfo <- prim[1:length(Ainfo)]
         G = t(as.matrix(expand.grid(rep(list(Ainfo,Ainfo )))))
         GP = t(as.matrix(expand.grid(rep(list(Pinfo,Pinfo )))))
         keep = GP[2,]>=GP[1,] #unique genotypes
         G <- G[,keep]  #store genotypes
         GP <- GP[1,keep]*GP[2,keep] #get prime product
         G[!G%in%names(popFreqQ[[loc]])] <- "99"         #Rename missing alleles in popFreqQ to "99":
         G0 <- paste0(G[1,],paste0("/",G[2,])) #make db-ref in one vector only
  
         newRow <- rep(NA,length(indD))  
         for(j in 1:length(GP)) { #for each genotype in population: Check in database
           #Always: Find corresponding genotype name by looking up stored primenumbers (same order as in popFreq!)
           rowind <- which(subD[,which(loc==dblocs)]==GP[j]) #samples with this genotype
           if(length(rowind)==0) next
           newRow[rowind] <- G0[j]
           #Get MAC of references:
           tmpmac <- 0
           #Count matching alleles over all mixtures:
           for(msel in mixSel) { #for each selected mixture
            evid0 <- set$samples[[msel]][[loc]]$adata
            tmpmac <- tmpmac + sum(G[,j]%in%evid0)
           } 
           macD[rowind] = macD[rowind] + tmpmac #add match counter  
           nlocs[rowind] <- nlocs[rowind] + 1 #counted only once!
          } #end for each genotype
          subD[,which(loc==dblocs)] <- newRow #force insertion of genotype-names
        } #end for each locus

        #determine individuals which will with LR=0 when pC=0 (i.e. no peak explained by unknowns or stutter)
        #can combine it to calculate qualitative LR
        LR0bool <- rep(FALSE,nrow(subD)) #boolean for reference which is not necessary to calculate for contLR
        LR1 <- rep(1,nrow(subD)) #LRmix vec
        nU_hp <- mod$nC_hp - sum(mod$condOrder_hp>0) #number of unknowns under Hp                    
        nU_hd <- mod$nC_hd - sum(mod$condOrder_hd>0) #number of unknowns under Hd                    
        for(loc in dblocs ) { #for each locus in db 
          if(is.null(popFreq[[loc]])) next #skip to next locus
              evidlist <- lapply( set$samples, function(x) x[[loc]]$adata ) #take out sample data:
              condR <- unlist(refData[[loc]][mod$condOrder_hp] ) #take out known refs under Hp 
              dbR <- subD[,which(loc==dblocs)] #take out DB-refs
              isNA <- is.na(dbR) #take out missing references
              if(all(isNA)) next #skipt locus if none to calculate
              dbR2 <- matrix(NA,nrow=nrow(subD),ncol=2) #create a matrix with NA
              dbR2[!isNA,] <- t(matrix(unlist(strsplit(dbR[!isNA] , "/")) ,nrow=2)) #store into new matrix
              move = as.numeric(dbR2[,2])<as.numeric(dbR2[,1])
              dbR2[move[!isNA],] <- dbR2[move[!isNA],c(2,1)]   #sort they are same genotype           
              dbR2 <- dbR2[!isNA,]
              if(sum(!isNA)==1) dbR2 <- rbind(dbR2) #require conversion if one possible combination
              undbR <- unique(dbR2) #get unique genotypes
              for(j in 1:nrow(undbR)) {
                dbind <-  which(dbR2[,1]==undbR[j,1] & dbR2[,2]==undbR[j,2]) #get index of matching genotypes
                ref0 <- c(undbR[j,],condR ) #conditional references
                Evid <- NULL
                for(ss in length(evidlist)) { #for each evidence
                  Ei <- evidlist[[ss]]	
                  if(par$prC==0) {
                   Ei2 <- Ei[!Ei%in%ref0] #set of unknown alleles (not explained by ref0) 
	             Ei3 <- Ei2[!as.numeric(Ei2)%in%(as.numeric(ref0)-1)] #set of unknown al leles (after explained by being stutter from ref0)
                   Ei4 <- Ei3[(as.numeric(Ei3)-1)%in%as.numeric(Ei3)] #set of unknown alleles (after explaining possible stutter from unknown)
                   Ei5 <- unique(c(Ei4,Ei3[!as.numeric(Ei3)%in%(as.numeric(Ei3)-1)])) #Alleles to explain removed stutter 
                   if(length(Ei5)>(2*nU_hp)) LR0bool[dbind] <- TRUE #if more alleles than unknown
                  } #end case of no drop-in
                  if(ss>1) Ei <- c(Ei,"0") #insert zero
                  Evid <- c(Evid,Ei)
                } #end for each evidence
                hp0 <- likEvid( Evid,T=ref0,V=NULL,x=nU_hp,theta=par$fst, prDHet=pDvec, prDHom=pDvec^2, prC=opt$QUALpC, freq=popFreqQ[[loc]])
                if(par$fst>0 | j==1) hd0 <- likEvid( Evid,T=condR,V=undbR[j,],x=nU_hd,theta=par$fst, prDHet=pDvec, prDHom=pDvec^2, prC=opt$QUALpC, freq=popFreqQ[[loc]])
                LR1[dbind] <- LR1[dbind]*hp0/hd0   #if more alleles than unknown
              }#end for each genotypes
        } #end for each locus

        print(paste0("Calculating continuous LR for ",sum(!LR0bool)," individual(s) in database ",dsel,"..."))
        #CONT LR calculation for each reference in table: 
        #unsubD <- unique( subD ) #get unique values. Not in use
        for(rind in 1:length(indD)) { #for each individual in database
         if(LR0bool[rind]) next #skip individual in database (which will have LR=0)
         Dind <- subD[rind,] #take out individual
         dblocs2 <- dblocs[!is.na(Dind)] #take out loci which the reference in database have
         locevalind <- locs_hd%in%dblocs2
         loceval <- locs_hd[locevalind] #locus to evaluate

         #setup for hp:
         #insert Dind to refData         
         if(is.null(refData)) { #
          refData2 <- list()
          condOrder_hp <- 1 #put conditional-index to model  
          condOrder_hd <- 0 #put conditional-index to model 
         } else { 
          refData2 <- refData[loceval] #take out only relevant loci to analyse
          condOrder_hp <- c(mod$condOrder_hp,max(mod$condOrder_hp)+1) #put conditional-index to model 
          condOrder_hd <- c(mod$condOrder_hd,0) #put conditional-index to model 
         }
         nR <- length(condOrder_hp) #number of references in refData2
         for(loc in loceval) refData2[[loc]]$ref99 <- unlist(strsplit(Dind[ which(loc==dblocs) ], "/"))  #insert data into a new ref
         samples <- lapply( set$samples, function(x) x[loceval] ) #take only relevant mixture data:
         

         if(ITYPE=="MLE") { #calculate with MLE
           logLi_hdeval <- logLi_hd[locevalind] #take out relevant values
           mlefit_hp <- contLikMLE(mod$nC_hp+1,samples,popFreqQ[loceval],refData2,condOrder_hp,mod$knownRef_hp,par$xi,par$prC,mleopt$nDone,par$threshT,par$fst,par$lambda,delta=mleopt$delta,pXi=par$pXi)
           if(par$fst>0) { #must calculate Hd once again (assume Rj is known)
            mlefit_hdj <- contLikMLE(mod$nC_hd,samples,popFreqQ[loceval],refData2,condOrder_hd,nR,par$xi,par$prC,mleopt$nDone,par$threshT,par$fst,par$lambda,delta=mleopt$delta,pXi=par$pXi)
            LRD[rind] <- exp(mlefit_hp$fit$loglik - mlefit_hdj$fit$loglik) #insert calculated LR adjusted by fst-correction
           } else {
            LRD[rind] <- exp(mlefit_hp$fit$loglik - sum(logLi_hdeval)) #insert calculated LR:
           }  
         } #END DB WITH TYPE MLE
         if(ITYPE=="INT") { #Calculate with INT
           int_hp <- contLikINT(mod$nC_hp+1, samples, popFreqQ[loceval], bhp$lower, bhp$upper, refData2, condOrder_hp, mod$knownRef_hp, par$xi, par$prC, optint$reltol, par$threshT, par$fst, par$lambda, par$pXi)     
           hd0INT <- numeric()
           if(par$fst==0 && length(hd0stored)>0) { #If any previous calculated values
             if(par$fst==0) { #can use stored value if fst=0
              for(l in 1:length(hd0stored)) { #for each element in list
               hd0locs <- hd0stored[[l]][-1]
               if(all(hd0locs%in%loceval) && all(loceval%in%hd0locs)) hd0INT <-  as.numeric(hd0stored[[l]][1]) #get stored hd-value
              }
             }
           }
           if(length(hd0INT)==0) { #calculate and store
             int_hd <- contLikINT(mod$nC_hd, samples, popFreqQ[loceval], bhp$lower, bhp$upper, refData2, condOrder_hd,nR, par$xi, par$prC, optint$reltol, par$threshT, par$fst, par$lambda, par$pXi) 
             hd0INT <- int_hd$margL
             hd0stored[[length(hd0stored) + 1]] <- c(hd0INT,loceval)
           }
           LRD[rind] <- int_hp$margL/hd0INT
         } #END DB WITH TYPE INT
         if(rind%%50==0) print(paste0(round(rind/length(indD)*100),"% finished"))
       }
       print(paste0(100,"% finished for database ",dsel))
      DBtab <- rbind(DBtab , cbind(indD,LRD,LR1,macD,nlocs) ) #add to DBtab
    } #end for each databases
    colnames(DBtab) <- c("Referencename","cont LR","qual LR","MAC","nLocs")
    assign("resDB",DBtab ,envir=mmTK) #assign deconvolved result to environment
    refreshTabDB("LRC") #LR is order to sort with
    svalue(nb) <- 6 #go to database search results window when finished     
   } #end doDB

##########################################################################################
############################Tab 4: MLE estimation:########################################
##########################################################################################
#HERE: WE DO BASICLY MLE-FITTING HERE: 
#BUT DECONVOLUTION-function AND DATABASE SEARCHING is implemented here (saves memory usage)!

tabMLEtmp <- glayout(spacing=30,container=tabMLE)

  f_savetableEVID = function(h,...) { #function for storing LR
   resEVID <- get("resEVID",envir=mmTK) #get EVID calculations when GUI starts
   if(is.null(resEVID)) {
    tkmessageBox(message="There is no Weight-of-Evidence results available!")
   } else {
    #tab <- c(resEVID$LRi,resEVID$LRmle,resEVID$LRlap)
    #tab <- cbind(c(names(resEVID$LRi),"JointMLE","JointLaplace"),tab,log10(tab))
    tab <- c(resEVID$LRi,resEVID$LRmle)
    tab <- cbind(c(names(resEVID$LRi),"JointMLE"),format(tab,digits=4),format(log10(tab),digits=4))
    colnames(tab) <- c("Marker","LR","log10LR")
    saveTable(tab, "txt") 
   }
  }

  f_savetableALL = function(h,...) { #function for storing MLE estimates of fitted models
   set <- get(paste0("set",h$action),envir=mmTK) #get all setup-object 

   printMLE <- function(mlefit,hyp,sig=4,colps="-") {
    mle <- cbind(mlefit$thetahat2,sqrt(diag(mlefit$thetaSigma2))) #standard deviation
    log10LR <- mlefit$loglik/log(10) #=log10(exp(mlefit$loglik))
    txt0 <- paste0("-------Estimates under ",hyp,"---------\n\n")
    txt1 <- paste0(c("param","MLE","Std.Err."),collapse=colps)
    for(i in 1:nrow(mle)) txt1 <- paste0(txt1,"\n",paste0( c(rownames(mle)[i],format(mle[i,],digits=sig)),collapse=colps) )
    txt2 <- "\n\n"
    txt2 <- paste0(txt2, "log10Lik=",format(log10LR,digits=sig))
    txt2 <- paste0(txt2, "\nLik=",format(10^log10LR,digits=sig),"\n\n")
    txt <- paste0(txt0,txt1,txt2)
    return(txt)
   }
   txt <- ""
   if(!is.null(set$mlefit_hd)) txt <- paste0(txt,printMLE(set$mlefit_hd$fit,"Hd"))
   if(!is.null(set$mlefit_hp)) txt <- paste0(txt,printMLE(set$mlefit_hp$fit,"Hp"))
   saveTable(txt, "txt") 
  } #end savetableALL


  #helpfunction ran when call deconvolution
  doDC <- function(mleobj) {
     dcopt <- get("optDC",envir=mmTK) #options when Deconvolution
     DCtable <- deconvolve(mleobj,dcopt$alphaprob,dcopt$maxlist)$table1 #be sure of having ri
     rownames(DCtable) <- 1:nrow(DCtable)
     DCtable<-addRownameTable(DCtable)
     colnames(DCtable)[1] <- "rank"
     assign("resDC",DCtable,envir=mmTK) #assign deconvolved result to environment
     refreshTabDC() #update table with deconvolved results
     svalue(nb) <- 5 #go to deconvolution results window (for all cases) when finished     
  }

  #helpfunction ran when call deconvolution
  doMCMC <- function(mleobj) { 
     optlist <- get("optMCMC",envir=mmTK)  #options for MCMC 
     if(any(is.na(mleobj$fit$phiSigma))) return();
     print(paste0("Sampling ",optlist$niter," samples with variation ",optlist$delta,". This could take a while... "))
     print("Note: You can change default number of iterations in toolbar menu.")
     mcmcfit <- contLikMCMC(mleobj,optlist$niter,optlist$delta)
     print(paste0("Sampling acceptance ratio=",round(mcmcfit$accrat,2),". This should be around 0.2"))
     validMCMC(mcmcfit,acf=FALSE) #don't plot acf
  }

  refreshTabMLE = function(type) { 
    #type={"EVID","DB","DC"}
   
    #optimizing options
    opt<- get("optMLE",envir=mmTK) #options when optimizing (nDone,delta)
    dec <- opt$dec #number of significant desimals to have in MLE print

    checkPositive(opt$delta,"Variance parameter of randomizer")
    checkPosInteger(opt$nDone,"Number of random startpoints")

    if(type!="START") {
     print(paste0(opt$nDone," random startpoints with variation ",opt$delta," are applied in the optimizer.")) 
     print("This could take a while...")
    }

    #helpfunction used to show MLE fit 
    tableMLE <- function(mlefit,tabmleX) {
     tabmleX1 = glayout(spacing=0,container=(tabmleX[1,1] <-gframe("Parameter estimates:",container=tabmleX))) 
     tabmleX2 = glayout(spacing=0,container=(tabmleX[2,1] <-gframe("Maximum Likelihood value",container=tabmleX))) 
     mle <- cbind(mlefit$thetahat2,sqrt(diag(mlefit$thetaSigma2)))
     log10LR <- mlefit$loglik/log(10) #=log10(exp(mlefit$loglik))
     tab <- cbind(rownames(mle),format(mle,digits=dec))
     colnames(tab) <- c("param","MLE","Std.Err.")
     tabmleX1[1,1] <- gtable(tab,container=tabmleX1,width=200,height=nrow(tab)*26)#,noRowsVisible=TRUE) #add to frame
     tabmleX2[1,1] =  glabel(text="log10lik=",container=tabmleX2)
     tabmleX2[1,2] =  glabel(text=format(log10LR,digits=dec),container=tabmleX2)
     tabmleX2[2,1] =  glabel(text="Lik=",container=tabmleX2)
     tabmleX2[2,2] =  glabel(text=format(10^log10LR,digits=dec),container=tabmleX2)
#     tabmleX2[3,1] =  glabel(text="Laplace P(E)=",container=tabmleX2)
#     tabmleX2[3,2] =  glabel(text=format(exp(fit$logmargL),digits=sig),container=tabmleX2)
    }

    if(type=="START") { #loads already calculated results if program starts
     set <- get("setEVID",envir=mmTK) #get setup for EVID
     mlefit_hd <- set$mlefit_hd
     mlefit_hp <- set$mlefit_hp
     if(is.null(mlefit_hd)) return(); #LR has not been calculated, return out of function!
    } else { #otherwise, function was called to make new calculations
     if(type=="EVID") set <- get("setEVID",envir=mmTK) #get setup for EVID
     if(type=="DB") set <- get("setDB",envir=mmTK) #get setup for DB
     if(type=="DC") set <- get("setDC",envir=mmTK) #get setup for DC

     #take out relevant parameters from stored list
     mod <- set$model
     par <- set$param     

     #fit under hd: (does it for all methods)
     mlefit_hd <- contLikMLE(mod$nC_hd,set$samples,set$popFreqQ,set$refDataQ,mod$condOrder_hd,mod$knownRef_hd,par$xi,par$prC,opt$nDone,par$threshT,par$fst,par$lambda,delta=opt$delta,pXi=par$pXi)
     if(!is.null(set$mlefit_hd) && set$mlefit_hd$fit$loglik>mlefit_hd$fit$loglik )  mlefit_hd <- set$mlefit_hd #the old model was better

     #fit under hp: (only for evidence)
     if(type=="EVID") {
      mlefit_hp <- contLikMLE(mod$nC_hp,set$samples,set$popFreqQ,set$refDataQ,mod$condOrder_hp,mod$knownRef_hp,par$xi,par$prC,opt$nDone,par$threshT,par$fst,par$lambda,delta=opt$delta,pXi=par$pXi)
      if(!is.null(set$mlefit_hp) && set$mlefit_hp$fit$loglik>mlefit_hp$fit$loglik )  mlefit_hp <- set$mlefit_hp #the old model was better
     } else {
      mlefit_hp <- NULL #not used otherwise
     }

     #store MLE result:
     #store best mle-values once again
     set$mlefit_hp=mlefit_hp #store fitted mle-fit
     set$mlefit_hd=mlefit_hd #store fitted mle-fit
     if(type=="EVID") assign("setEVID",set,envir=mmTK) #store setup for EVID
     if(type=="DB") assign("setDB",set,envir=mmTK) #store setup for DB
     if(type=="DC") assign("setDC",set,envir=mmTK) #store setup for DC
    }

    #GUI:
    tabmleA = glayout(spacing=0,container=(tabMLEtmp[1,1] <- gframe("Estimates under Hd",container=tabMLEtmp))) 
    tableMLE(mlefit_hd$fit,tabmleA)
    tabmleA3 = glayout(spacing=0,container=(tabmleA[3,1] <-gframe("Further Action",container=tabmleA))) 
    tabmleA3[1,1] <- gbutton(text="MCMC simulation",container=tabmleA3,handler=function(h,...) { doMCMC(mlefit_hd) } )
    tabmleA3[2,1] <- gbutton(text="Deconvolution",container=tabmleA3,handler=function(h,...) { doDC(mlefit_hd) }  )
    tabmleA3[3,1] <- gbutton(text="Model validation",container=tabmleA3,handler=function(h,...) { } )

    if(type=="EVID") { #used only for weight-of-evidence
     tabmleB = glayout(spacing=0,container=(tabMLEtmp[1,2] <-gframe("Estimates under Hp",container=tabMLEtmp))) 
     tableMLE(mlefit_hp$fit,tabmleB)
     tabmleB3 = glayout(spacing=0,container=(tabmleB[3,1] <-gframe("Further Action",container=tabmleB))) 
     tabmleB3[1,1] <- gbutton(text="MCMC simulation",container=tabmleB3,handler=function(h,...) { doMCMC(mlefit_hp) } )
     tabmleB3[2,1] <- gbutton(text="Deconvolution",container=tabmleB3,handler=function(h,...) {  doDC(mlefit_hp) }  )
     tabmleB3[3,1] <- gbutton(text="Model validation",container=tabmleB3,handler=function(h,...) { } )
    }

    #We show weight-of-evidence
    tabmleD = glayout(spacing=5,container=(tabMLEtmp[2,1] <-gframe("Further evaluation",container=tabMLEtmp))) 
    tabmleD[1,1] <- gbutton(text="Optimize model more",container=tabmleD,handler=function(h,...) { refreshTabMLE(type)  } )
    if(type=="EVID") {
     logLRmle <- mlefit_hp$fit$loglik - mlefit_hd$fit$loglik
     log10LRmle <- logLRmle/log(10) 
     LRmle <- exp(logLRmle)
     LRlap <- exp(mlefit_hp$fit$logmargL - mlefit_hd$fit$logmargL)
     LRi <- exp(logLiki(mlefit_hp)-logLiki(mlefit_hd))
     resEVID <- list(LRmle=LRmle,LRlap=LRlap,LRi=LRi) 
     assign("resEVID",resEVID,envir=mmTK) #store EVID calculations
    } 
    if(type=="START") {
     resEVID <- get("resEVID",envir=mmTK) #get EVID calculations when GUI starts
     if(!is.null(resEVID)) { #put variables in environment
       LRmle <- resEVID$LRmle
       LRlap <- resEVID$LRlap
       LRi <- resEVID$LRi
     }
    } #end if start
    if(type=="EVID" || type=="START") {
     tabmleC = glayout(spacing=0,container=(tabMLEtmp[1,3] <-gframe("Weight-of-evidence\n(MLE based)",container=tabMLEtmp))) 
     tabmleC1 = glayout(spacing=0,container=(tabmleC[1,1] <-gframe("Joint LR",container=tabmleC))) 
     tabmleC1[1,1] =  glabel(text="LR=",container=tabmleC1)
     tabmleC1[1,2] =  glabel(text=format(LRmle,digits=dec),container=tabmleC1)
     tabmleC1[2,1] =  glabel(text="log10LR=",container=tabmleC1)
     tabmleC1[2,2] =  glabel(text=format(log10LRmle,digits=dec),container=tabmleC1)
     tabmleC3 = glayout(spacing=0,container=(tabmleC[2,1] <-gframe("LR for each loci",container=tabmleC))) 
     for(i in 1:length(LRi)) {
      tabmleC3[i,1] =  glabel(text=names(LRi)[i],container=tabmleC3)
      tabmleC3[i,2] =  glabel(text=format(LRi[i],digits=dec),container=tabmleC3)
     }
#     tabmleC2 = glayout(spacing=0,container=(tabmleC[3,1] <-gframe("Laplace approximation based:",container=tabmleC))) 
#     tabmleC2[1,1] =  glabel(text="LR=",container=tabmleC2)
#     tabmleC2[1,2] =  glabel(text=LRlap,container=tabmleC2)
#     tabmleC2[2,1] =  glabel(text="log10LR=",container=tabmleC2)
#     tabmleC2[2,2] =  glabel(text=log10(LRlap),container=tabmleC2)
     tabmleD[2,1] <- gbutton(text="Continuous LR\n(Integrated Likelihood based)",container=tabmleD,handler=function(h,...) { doINT("EVID") } ) 
    } #end if EVID or START
    if(type=="DB") tabmleD[2,1] <- gbutton(text="Database search",container=tabmleD,handler=function(h,...) { doDB("MLE")} )
    tabmleE = glayout(spacing=0,container=(tabMLEtmp[2,2] <-gframe("Save results to file",container=tabMLEtmp))) 
    tabmleE[1,1] <- gbutton(text="All results",container=tabmleE,handler=f_savetableALL,action=type)
    if(type=="EVID") tabmleE[2,1] <- gbutton(text="Only LR results",container=tabmleE,handler=f_savetableEVID)

    focus(mainwin) #focus window after calculations are done
  } #end refresh tab-frame of MLE-fit

  refreshTabMLE(type="START") #Show already calculted evidence-results when program starts


##############################################################
###############Tab 5: Deconvolution results:##################
##############################################################
 tabDCtmp <- glayout(spacing=30,container=tabDC)

 f_savetableDC = function(h,...) {
   DCtable <- get("resDC",envir=mmTK) #get deconvolved result 
   if(is.null(DCtable)) {
    tkmessageBox(message="There is no deconvolution results available!")
   } else {
    saveTable(DCtable, "txt") #save 
   }
 }
 refreshTabDC = function() {
   DCtable <- get("resDC",envir=mmTK) #get deconvolved result 
   if(!is.null(DCtable)) {
    tabDCa = glayout(spacing=1,container=(tabDCtmp[1,1] <-glayout(spacing=0,container=tabDCtmp)),expand=TRUE) #table layout
    tabDCb = glayout(spacing=1,container=(tabDCtmp[2,1] <-glayout(spacing=0,container=tabDCtmp)),expand=TRUE) #table layout
    tabDCa[1,1] <- gtable(DCtable,container=tabDCa,height=mwH*0.5,width=mwW,height=mwH-2*mwH/3,do.autoscroll=TRUE,noRowsVisible=TRUE) #add to frame
    tabDCb[2,1] <- gbutton(text="Save table",container=tabDCb,handler=f_savetableDC)  
   }
 }
 refreshTabDC() #open results when program starts


##############################################################
###############Tab 6: Database search:########################
##############################################################

 tabDBtmp <- glayout(spacing=30,container=tabDB) #grid-layout

 f_savetableDB = function(h,...) {
   DBsearch <-get("resDB",envir=mmTK) #load results from environment
   if(is.null(DBsearch)) {
    tkmessageBox(message="There is no database search results available.")
    return()
   }
   sep="txt"
   tabfile = gfile(text="Save table",type="save")
   if(!is.na(tabfile)) {
    if(length(unlist(strsplit(tabfile,"\\.")))==1) tabfile = paste0(tabfile,".",sep)
    ord <- order(as.numeric(DBsearch[,2]),decreasing=TRUE) 
    if(h$action=="LRC") ord <- order(as.numeric(DBsearch[,3]),decreasing=TRUE) 
    if(h$action=="MAC") ord <- order(as.numeric(DBsearch[,4]),decreasing=TRUE) 
    if(length(ord)<=1) write.table(DBsearch,file=tabfile,quote=FALSE,sep="\t",row.names=FALSE) #load environment
    if(length(ord)>1) write.table(DBsearch[ord,],file=tabfile,quote=FALSE,sep="\t",row.names=FALSE) #load environment
    print(paste("Full ranked table saved in ",tabfile,sep=""))
   }
  }

 refreshTabDB = function(rank) {
   DBtable <- get("resDB",envir=mmTK) #get deconvolved result 
   if(!is.null(DBtable)) {
    if(rank=="LRC") rankind <- 2
    if(rank=="LRQ") rankind <- 3
    if(rank=="MAC") rankind <- 4
    ord <- order(as.numeric(DBtable[,rankind]),decreasing=TRUE) #need to convert to numeric!
    tabDBa = glayout(spacing=1,container=(tabDBtmp[1,1] <-glayout(spacing=0,container=tabDBtmp)),expand=TRUE) #table layout
    tabDBb = glayout(spacing=1,container=(tabDBtmp[2,1] <-glayout(spacing=0,container=tabDBtmp)),expand=TRUE) #table layout
    if(length(ord)<=1) tabDBa[1,1] <- gtable(DBtable,container=tabDBa,width=mwW,height=mwH*0.5,do.autoscroll=TRUE,noRowsVisible=TRUE) #add to frame
    if(length(ord)>1) tabDBa[1,1] <- gtable(DBtable[ord[1:min(get("optDB",envir=mmTK)$maxDB,length(ord))],] ,container=tabDBa,width=mwW,height=mwH*0.5,do.autoscroll=TRUE,noRowsVisible=TRUE) #add to frame
    tabDBb[2,1] <- gbutton(text="Save table",container=tabDBb,handler=f_savetableDB,action=rank)  
   }
 }
 refreshTabDB("LRC") #when program starts
 visible(mainwin) <- TRUE
 focus(mainwin)


###############################################################
###############Tab 7: Model validation:########################
###############################################################



###############################################################
###############Tab 8: LRmix module:############################
###############################################################
 #uses only qualitative information

 tabLRMIXtmp <- glayout(spacing=30,container=tabLRMIX)

 f_savetableEVIDLRMIX = function(h,...) { #function for storing LR
   LRi <- get("resEVIDLRMIX",envir=mmTK) #get EVID calculations when GUI starts
   if(is.null(LRi)) {
    tkmessageBox(message="There was no Weight-of-Evidence results available!")
   } else {
    tab <- c(LRi,prod(LRi))
    tab <- cbind(c(names(LRi),"Joint"),tab,log10(tab))
    colnames(tab) <- c("Locus","LR","log10LR")
    saveTable(tab, "txt") 
   }
  }
  noSamples = function(hyp) { #helpfunction for tell user that wrong model assumption was used.
    gmessage(message=paste0("No samples was accepted out of the first 1000.\nPlease change hypothesis under ",hyp),title="Wrong model specification",icon="error")
    stop()
  }  

 refreshTabLRMIX = function() {
  require(forensim)
 
  #helpfunction to make GUI-table with LR calculations
  tableLR = function(LRi) { 
   sig <- 4 #number of signif to show
   tabLRmixB1[1,1] = glabel(text="Locus",container=tabLRmixB1)
   tabLRmixB1[1,2] = glabel(text="LR",container=tabLRmixB1)
   tabLRmixB1[1,3] = glabel(text="log10LR",container=tabLRmixB1)

   for(loc in locs) {
    i = which(loc==locs)
    tabLRmixB1[i+1,1] = glabel(text=loc,container=tabLRmixB1)
    tabLRmixB1[i+1,2] = glabel(text=format(LRi[i],digits=sig),container=tabLRmixB1)
    tabLRmixB1[i+1,3] = glabel(text=format(log10(LRi[i]),digits=sig),container=tabLRmixB1)
   }

   #show jointly:
   tabLRmixB2[1,1] = glabel(text="LR",container=tabLRmixB2)
   tabLRmixB2[2,1] = glabel(text="log10LR",container=tabLRmixB2)
   tabLRmixB2[1,2] = glabel(text=format(prod(LRi),digits=sig),container=tabLRmixB2)
   tabLRmixB2[2,2] = glabel(text=format(sum(log10(LRi),digits=sig)),container=tabLRmixB2)
  }

  #helpfunction for calculating LR for each given dropout pD (takes a numeric)
  doLR = function(pD) {
    pDhp <- rep(pD,mod$nC_hp)
    pDhd <- rep(pD,mod$nC_hd)
    hpvec <- hdvec <- rep(1,length(locs))
    for(loc in locs) {
      hpvec[which(loc==locs)] <- likEvid( Evidlist[[loc]],T=refList_hp[[loc]]$Ri,V=refList_hp[[loc]]$Ki,x=mod$nC_hp-refList_hp[[loc]]$nR,theta=par$fst, prDHet=pDhp, prDHom=pDhp^2, prC=par$prC, freq=set$popFreqQ[[loc]])
      hdvec[which(loc==locs)] <- likEvid( Evidlist[[loc]],T=refList_hd[[loc]]$Ri,V=refList_hd[[loc]]$Ki,x=mod$nC_hd-refList_hd[[loc]]$nR,theta=par$fst, prDHet=pDhd, prDHom=pDhd^2, prC=par$prC, freq=set$popFreqQ[[loc]])
    }
    LRi <- hpvec/hdvec
    names(LRi) <- locs
    return(LRi)
  }

 #helpfunction to get conditional refs under a hypothesis
  getConds <- function(condOrder,knownref) {
    cond <- which(condOrder>0) #ind of conditional refs (they are increasingly sorted)
    Ri <- Ki <- NULL
    for(rr in cond ) Ri <- c(Ri,set$refDataQ[[loc]][[rr]])
    for(rr in knownref) Ki <- c(Ki,set$refDataQ[[loc]][[rr]])
    return(list(Ri=Ri,Ki=Ki,nR=length(cond) ))
  }

  #take out relevant parameters from stored list
  set <- get("setEVID",envir=mmTK) #get setup for EVID
  mod <- set$model
  par <- set$param     

  #Data:
  locs <- names(set$popFreqQ) #get analysing loci
  nS <- length(set$samples) #number of samples
 
  #Prepare Evidence and refs under each hypothesis:
  Evidlist <- list()
  refList_hp <- list()
  refList_hd <- list()
  for(loc in locs) {
    Ei <- NULL #get evidence
    for(ss in 1:nS) {
     if(ss>1) Ei <- c(Ei,0) #seperate with 0  
     adata <- set$samples[[ss]][[loc]]$adata
     if(length(adata)==0) adata=0 #is empty
     Ei <- c(Ei,adata)
    } 
    Evidlist[[loc]] <- Ei
    refList_hp[[loc]] <- getConds(mod$condOrder_hp,mod$knownref_hp) #under hp
    refList_hd[[loc]] <- getConds(mod$condOrder_hd,mod$knownref_hd) #under hd
  }


  #GUI:
  tabLRmixA = glayout(spacing=0,container=(tabLRMIXtmp[1,1] <-gframe("Analysis of qualitative LR",container=tabLRMIXtmp))) 
  tabLRmixB = glayout(spacing=0,container=(tabLRMIXtmp[1,2] <-gframe("Weight-of-Evidence",container=tabLRMIXtmp))) 

  tabLRmixA1 = glayout(spacing=0,container=(tabLRmixA[1,1] <-gframe("Preanalysis",container=tabLRmixA)))  
  tabLRmixA2 = glayout(spacing=0,container=(tabLRmixA[2,1] <-gframe("Calculation",container=tabLRmixA))) 
  tabLRmixA3 = glayout(spacing=0,container=(tabLRmixA[3,1] <-gframe("Postanalysis",container=tabLRmixA))) 

  tabLRmixB1 = glayout(spacing=0,container=(tabLRmixB[1,1] <-gframe("Loci",container=tabLRmixB)))  
  tabLRmixB2 = glayout(spacing=0,container=(tabLRmixB[2,1] <-gframe("Jointly",container=tabLRmixB)))  


  #Preanalysis (sensistivity and dropout plots)
  tabLRmixA1[1,1] <- gbutton(text="Sensitivity",container=tabLRmixA1,handler=function(x) {
    #get range from options under toolbar:
    optLRMIX <- get("optLRMIX",envir=mmTK) #options when integrating (reltol and boundaries)
    range <- optLRMIX$range
    nTicks <- optLRMIX$nticks
    checkProb(range,"The range")
    checkPosInteger(nTicks, "The number of ticks")  
    pDvec <- seq(1e-6,range,l=nTicks) #set very small dropout probability
    LRivec <- Vectorize(doLR) (pDvec) #return LR(pD,loc)
    logLR <<- colSums(log10(LRivec)) #joint LR
    ylims <- range(logLR[!(is.nan(logLR) | is.infinite(logLR))]) 
    if(ylims[2]>-1) ylims[1] <- -1  #lower limit
    plot(pDvec ,logLR ,xlab="Pr(Dropout)",ylab="log_10(LR)",ylim=ylims + c(0,0.5),main="Sensitivity plot")
    lines(spline(pDvec ,logLR,n=nTicks )) #use spline to interpolate points
    abline(h=0,col="gray")
   })
  tabLRmixA1[2,1] <- gbutton(text="Conservative LR",container=tabLRmixA1,handler=function(x) {
    optLRMIX <- get("optLRMIX",envir=mmTK) 
    nsample <- optLRMIX$nsample
    alpha <- optLRMIX$alpha
    totA <-  sapply(  set$samples, function(x) sum(sapply(x,function(y) length(y$adata)) ) ) #number of alleles for each samples
    print("Total number of observed alleles for sample(s):")
    print(totA)
    refHp <- refHd <- NULL
    if( any(mod$condOrder_hp>0) ) refHp <- lapply(set$refData ,function(x) x[mod$condOrder_hp]) #must have the original refData!
    if( any(mod$condOrder_hd>0) ) refHd <- lapply(set$refData ,function(x) x[mod$condOrder_hd]) #must have the original refData!

  
    for(ss in 1:nS) { #for each sample (do dropout calculation)
     print(paste0("For evidence ",names(set$samples)[[ss]],":"))
     print("Estimating quantiles from allele dropout distribution under Hp...")
     DOdist <- simDOdistr(totA=totA[ss],nC=mod$nC_hp,popFreq=set$popFreq,refData=refHp,minS=nsample,prC=par$prC)
     if(length(DOdist)==0) noSamples("Hp")
     qqhp <- quantile(DOdist ,c(alpha,1-alpha)) #get quantiles
     print(qqhp)
     print("Estimating quantiles from allele dropout distribution under Hd...")
     DOdist <- simDOdistr(totA=totA[ss],nC=mod$nC_hd,popFreq=set$popFreq,refData=refHd,minS=nsample,prC=par$prC)
     if(length(DOdist)==0) noSamples("Hd")
     qqhd <- quantile(DOdist ,c(alpha,1-alpha)) #get quantiles
     print(qqhd)
     if(ss==1) consDO <- rbind(qqhp,qqhd)
     if(ss>1) {
      if(qqhp[1]<consDO[1,1]) consDO[1,1] <- qqhp[1]
      if(qqhp[2]>consDO[1,2]) consDO[1,2] <- qqhp[2]
      if(qqhd[1]<consDO[2,1]) consDO[2,1] <- qqhd[1]
      if(qqhd[2]>consDO[2,2]) consDO[2,2] <- qqhd[2]
     }
    }
    consDO <- signif(consDO,2) #2 significant numbers to use
    print(consDO) #print out

    #calculate corresponding LR's
    consLRi <- Vectorize(doLR) (c(consDO))
    conslog10LR <- colSums(log10(consLRi)) #get joint log10LR
    selInd <- which.min(conslog10LR) #select reporting LR (which gives min LR)
    selDO <- consDO[selInd] #take out selected dropout
    LRi <- consLRi[,selInd] #take out selected LR-values
    svalue(tabLRmixA2[1,2]) <- selDO  #update selected dropout
    assign("resEVIDLRMIX",LRi,envir=mmTK) #assign evidence weighting results - Based on LRmix
    tableLR(LRi) #update table with calculated LR
  })
  #Calculation
  tabLRmixA2[1,1] <-  glabel(text="Dropout prob:",container=tabLRmixA2)
  tabLRmixA2[1,2] <-  gedit(text="0.05",container=tabLRmixA2,width=8) #this is updated after dropout distr is ran
  tabLRmixA2[2,1] <-  gbutton(text="Calculate LR",container=tabLRmixA2,handler=function(x) {
    pD <- as.numeric(svalue(tabLRmixA2[1,2]))
    checkProb(pD,"The allele dropout probability")
    LRi <- doLR(pD) 
    assign("resEVIDLRMIX",LRi,envir=mmTK) #assign evidence weighting results - Based on LRmix
    tableLR(LRi) #update table with calculated LR
  }) 
  tabLRmixA2[2,2] <-  gbutton(text="Save table",container=tabLRmixA2,handler=f_savetableEVIDLRMIX)

  #Tippet-analysis frame:
  tippets <- mod$knownref_hd #known non-contributors under Hd
  if(!is.null(tippets)) {
   tN <- names(set$refData[[1]][tippets]) #tippet names
   tabLRmixA3[1,1] <- glabel( "Select reference to tippet:",container=tabLRmixA3)
   tabLRmixA3[2,1] <- gcombobox( items=tN ,container=tabLRmixA3)
   tabLRmixA3[3,1] <- gbutton(text="Do Tippet",container=tabLRmixA3,handler=function(x) {

     #calculate LR for all genotypes in original popFreq.
     pD <- as.numeric(svalue(tabLRmixA2[1,2])) #take dropout-value as given in GUI
     tipref <- svalue(tabLRmixA3[2,1]) #get name of reference to tippet
     Glist <- getGlist(set$popFreqQ) #get random man-Glist 

     print("Precalculating for tippet plot...")
     #calculate LRs directly here: 
     tipsel <- which(tN==tipref) #index of tippet to select
     tipind <- mod$knownref_hd[tipsel] #get tip-ind in refData
     modtipind <- mod$condOrder_hp[tipind] #get position in system of tippet
     pDhp <- rep(pD,mod$nC_hp)
     pDhd <- rep(pD,mod$nC_hd)
     for(loc in locs) { #Calcualte for each locus:
       nG <- length(Glist[[loc]]$Gprob) #number of genotypes
       Glist[[loc]]$LR <- rep(NA,nG) #init space for LR
       refhptmp <- refList_hp[[loc]]$Ri  #take out contributing replicates under Hp
       nrefhdtmp <- refList_hd[[loc]]$Ki  #take out non-contributing replicates under Hd
       for(j in 1:nG) { #for each genotypes
        refhptmp[ 2*modtipind -c(1,0) ] <- Glist[[loc]]$G[j,] #insert genotype to reference
        nrefhdtmp[ 2*tipsel-c(1,0) ] <- Glist[[loc]]$G[j,] #insert genotype to reference
        hp0 <- likEvid( Evidlist[[loc]],T=refhptmp,V=refList_hp[[loc]]$Ki,x=mod$nC_hp-refList_hp[[loc]]$nR,theta=par$fst, prDHet=pDhp, prDHom=pDhp^2, prC=par$prC, freq=set$popFreqQ[[loc]])
        hd0 <- likEvid( Evidlist[[loc]],T=refList_hd[[loc]]$Ri,V=nrefhdtmp,x=mod$nC_hd-refList_hd[[loc]]$nR,theta=par$fst, prDHet=pDhd, prDHom=pDhd^2, prC=par$prC, freq=set$popFreqQ[[loc]])
        Glist[[loc]]$LR[j] <- hp0/hd0 #store LR
       }
     } #end for each locus
     ntippet <- get("optLRMIX",envir=mmTK)$ntippets #get number of tippets to simulate
     print(paste0("Simulating ",ntippet," tippets..."))
     RMLR <- rep(0,ntippet) #vector of tippets
     for(loc in locs) {
      RMLR <- RMLR + sample(log10(Glist[[loc]]$LR),ntippet,replace=TRUE,prob=Glist[[loc]]$Gprob)
     }
     minmax = range(RMLR[!is.infinite(RMLR)])
     xbar <- sum(10^RMLR)/ntippet #mean LR
     print( paste0("Mean of samples = ", xbar ))
     empvarLR <- (sum((10^(2*RMLR))) - ntippet*xbar^2)/(ntippet-1)  #emperical variance
     varxbar <- empvarLR/ntippet
     print( paste0("Standard Error of samples = ",sqrt(varxbar)) )
     quantiles<-quantile(RMLR,c(0.01,0.05,0.5,0.95,0.99))
     print(quantiles)
     if(ntippet>5e6) stop("Number of tippets was above 5mill. This is too large to plot!")
     lr0 <- NULL #log10 LR
     LRi <- get("resEVIDLRMIX",envir=mmTK) #get EVID calculations when GUI starts
     if(!is.null(LRi))  lr0 <- sum(log10(LRi))
     mtxt = paste0("Tippet calculation for ",tipref," with ",ntippet," samples.")
     plot(ecdf(RMLR),xlim=c(minmax[1],max(minmax[2],lr0)) , main=mtxt,xlab="log10(LR)")
     if(!is.null(lr0)) {
      points(lr0,1,pch=10,col="blue")
      lines(rep(lr0,2),c(1,0),lty=1,col="blue",lwd=0.5)
      print(paste0("Discriminatory metric (log10(LR) - q99) = ",lr0 - quantiles[length(quantiles)] ))
     }
   }) #end TIPPET BUTTON
  } #end if not tippet possible
 } #end refresh funtion

} #end funcions

