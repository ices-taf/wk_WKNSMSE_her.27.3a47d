#-------------------------------------------------------------------------------
# WKNSMSE
#
# Author: Benoit Berges
#         WMR, The Netherland
# email: benoit.berges@wur.nl
#
#  MSE of North Sea Herring
#
# Date: 2018/11/18
#
# Build for R3.5.1, 64bits
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 1) load packages 
#    setup paths 
#    load functions
#-------------------------------------------------------------------------------

#rm(list=ls())


args=(commandArgs(TRUE))
#args <- 'ftar=0.24_btrig=1.4e6_HCR=1_IAV=0_BB=0'
args    <- strsplit(args,"_")
ftarget <- as.numeric(substr(args[[1]][1],6,9))
btrigger<- as.numeric(substr(args[[1]][2],7,11))

# HCR
HCR     <- ifelse(as.numeric(substr(args[[1]][3],5,nchar(args[[1]][3])))==1,"A","B")

# IAV
if(as.numeric(substr(args[[1]][4],5,nchar(args[[1]][4]))) == 0){
  IAV <- 0
}
if(as.numeric(substr(args[[1]][4],5,nchar(args[[1]][4]))) == 1){
  IAV <- 'A'
}
if(as.numeric(substr(args[[1]][4],5,nchar(args[[1]][4]))) == 2){
  IAV <- 'B'
}

# BB
if(as.numeric(substr(args[[1]][5],4,nchar(args[[1]][4]))) == 0){
  BB <- 0
}
if(as.numeric(substr(args[[1]][5],4,nchar(args[[1]][5]))) == 1){
  BB <- 'A'
}
if(as.numeric(substr(args[[1]][5],4,nchar(args[[1]][5]))) == 2){
  BB <- 'B'
}

if(IAV == 0){
  IAV <- NULL
}else if(IAV == "B"){
  IAV <- c("A","B")
}

if(BB == 0){
  BB  <- NULL
}else if(BB == "B"){
  BB  <- c("A","B")
}

cat(ftarget,"\t",btrigger,"\t",HCR,"\t",IAV,"\t",BB,"\n")


library(FLSAM)
library(FLEDA)
library(FLFleet)
library(minpack.lm)  # install.packages("minpack.lm")
library(stats)

# define path to directory
#path          <- "D:/Work/Herring MSE/NSAS/"
#path              <- "D:/git/wk_WKNSMSE_her.27.3a47d/R/"
#path              <- "F:/WKNSMSE/wk_WKNSMSE_her.27.3a47d/R"
#path <- 'E:/wk_WKNSMSE_her.27.3a47d/R'
#path <- 'D:/Repository/NSAS_MSE/wk_WKNSMSE_her.27.3a47d/R/'
path <- "/home/hintz001/wk_WKNSMSE_her.27.3a47d/R"
assessment_name   <- "NSAS_WKNSMSE2018"
try(setwd(path),silent=TRUE)

# paths to different subfolders
dataPath      <- file.path(".","data/")
outPath       <- file.path(".","results/")
scriptPath    <- file.path(".","side_scripts/")
functionPath  <- file.path(".","functions/")

source(file.path(functionPath,"MSE_assessment.R"))
source(file.path(path,"04_forecastScenarios.r"))
source(file.path(path,"forecastFunctions.r"))

#-------------------------------------------------------------------------------
# 2) Initialize
#
# Define MSE parameters, 
# load objects initialized previously
#   * Biology
#     - biol
#     - surveys
#   * Fisheries
#     - FCProp
#     - catchD
#     - F sel: FAsel, FCsel, FBDsel
#-------------------------------------------------------------------------------

nits                <- 200
# load object
load(file.path(outPath,paste0(assessment_name,'_init_MSE_',ac(nits),'.RData')))
stkAssessment.ctrl <- NSH.ctrl
if(nits == 1000)
  load(file.path(outPath,"stkAssessment2018.init1000.RData"))
if(nits == 200)
  load(file.path(outPath,"stkAssessment2018.init.RData"))

# load MSE parameters
load(file.path(outPath,paste0(assessment_name,'_parameters_MSE_',ac(nits),'.RData')))
fishery@landings.sel[,projPeriod] <- sweep(fishery@landings.sel[,projPeriod],c(2:6),quantMeans(fishery@landings.sel[,projPeriod]),"/")
biol@harvest.spwn[] <- 0.67

strFleet    <- c('A','B','C','D')
nFleets     <- length(strFleet)
nAges       <- dim(biol@stock.n)[1]
surveyNames <- names(surveys)
escapeRuns  <- numeric()


#-------------------------------------------------------------------------------
# 2) ref points
# Ftarget = FMSY and Btrigger = MSYBtrigger for now
#
#-------------------------------------------------------------------------------

referencePoints <- list(Fmsy = 0.26,
                        Fsq  = NA,
                        Flim = 0.34,
                        Fpa  = 0.30,
                        Blim = 800000,
                        Bpa  = 900000,
                        MSYBtrigger = 1400000,
                        Ftarget = ftarget,
                        F01 = 0.05,
                        Btrigger = btrigger)

managementRule  <- list(HCR = HCR,
                        TACIAV=IAV, #"A","A+B","NULL"
                        BB = BB)
                        
runName         <- paste0("NSAS_Ftar_",referencePoints$Ftarget,
                              "_Btrig_",referencePoints$Btrigger,
                              "_HCR_",managementRule$HCR,
                              "_TACIAV_",paste(managementRule$TACIAV,collapse=""),
                              "_BB_",paste(managementRule$BB,collapse=""),
                              "_",nits,"iters.RData")

newUptakes      <- F
#------------------------------------------------------------------------------#
# 3) Define TACs for A, B and D fleets. 
# TACs for A and B fleets are taken out of HAWG2018. This needs updating
#
# Note 1: The Cfleet is defined as a proportion of F.
# Note 2: TAC for C and D fleets are for the WB
# Note 3: TACs for the D fleet is kept constant for future years
#------------------------------------------------------------------------------#

TAC                       <- FLQuant(NA,dimnames=list(age='all',
                                      year=histMinYr:(futureMaxYr+3),
                                      unit=c('TAC'),
                                      season='all',
                                      area=c('A','B','C','D'),
                                      iter=1:nits))

#- TAC for A fleet in NS
TAC_A                     <- read.table(file.path(dataPath,'TAC_A.csv'),sep = ",")
TAC[,ac(TAC_A[,1]),,,"A"] <- TAC_A[,2]
#- TAC for B fleet in NS
TAC_B                     <- read.table(file.path(dataPath,'TAC_B.csv'),sep = ",")
TAC[,ac(TAC_B[,1]),,,"B"] <- TAC_B[,2]
#- TAC for D fleet in WB
TAC_D                     <- read.table(file.path(dataPath,'TAC_D.csv'),sep = ",")
TAC[,ac(TAC_D[,1]),,,"D"] <- TAC_D[,2]
TAC[,ac((max(TAC_D[,1])+1):(futureMaxYr+3)),,,"D"] <- TAC[,ac(max(TAC_D[,1])),,,"D"] # TAC is fixed for the D fleet

#- TAC for C fleet in WB
TAC_C                     <- read.table(file.path(dataPath,'TAC_C.csv'),sep = ",")
TAC[,ac(TAC_C[,1]),,,"C"] <- TAC_C[,2]
#- Fixed TAC C in WB (used for transfer to the A fleet)
TAC[,ac((max(TAC_C[,1])+1):(futureMaxYr+3)),,,"C"] <- TAC_C[dim(TAC_C)[1],2]

if(newUptakes){
  #- setup variables for transfer, uptake and split
  uptakeFleets              <- read.table(file.path(dataPath,'over_underfishing2017.csv'),sep = ",")

  #- Transfer from C fleet TAC to fleet A
  Ctransfer                 <- matrix(runif((length(projPeriod)+3)*nits,min=0.4, max=0.5),nrow=nits,ncol=length(projPeriod)+3)    # Transfer of TAC from IIIa to IVa for C fleet in assessment year. Set between 0.4 and 0.5
  # update for the D fleets
  Duptake                   <- matrix(1,nrow=nits,ncol=length(projPeriod)+3)    # assume full uptake for the D fleet
  DSplitHist                <- read.table(file.path(dataPath,'D_split.csv'),sep = ",") # get mean and sd from historical data for NSAS/WBSS split for the D fleet
  Dsplit                    <- matrix(rnorm((length(projPeriod)+3)*nits,
                                      mean=mean(DSplitHist$V2),
                                      sd=sd(DSplitHist$V2)/2),
                                      nrow=nits,ncol=length(projPeriod)+3)
  #Dsplit                    <- matrix(rnorm((length(projPeriod)+3)*nits,mean=0.6,sd=0.1),nrow=nits,ncol=length(projPeriod)+3)
  # update for the B fleet
  Buptake                   <- matrix(rnorm ((length(projPeriod)+3)*nits,
                                      mean(an(as.vector(uptakeFleets[2:16,3])),na.rm=TRUE), # mean over available historical values
                                      sd(an(as.vector(uptakeFleets[2:16,3])),na.rm=TRUE)/2),   # sd over available historical values
                                      nrow=nits,ncol=length(projPeriod)+3)
  TAC_var                   <- array(NA,
                                     dim=c(length(projPeriod)+3,nits,4),
                                     dimnames=list('years' = ac(an(projPeriod)[1]:(an(projPeriod)[length(projPeriod)]+3)),
                                                   'iter' = 1:nits,
                                                   'var' = c('Ctransfer','Duptake','Dsplit','Buptake')))
  TAC_var[,,'Ctransfer']    <- t(Ctransfer)
  TAC_var[,,'Duptake']      <- t(Duptake)
  TAC_var[,,'Dsplit']       <- t(Dsplit)
  TAC_var[,,'Buptake']      <- t(Buptake)

  save(Ctransfer,Duptake,DSplitHist,Dsplit,Buptake,TAC_var,file=file.path(outPath,paste0("SplitUptakes",nits,".RData")))
} else {
  load(file.path(outPath,paste0("SplitUptakes",nits,".RData")))
}
CATCH                     <- TAC
FHCR                      <- FLQuant(NA, dimnames=list(age=0:8,year=projPeriod,unit="unique",season="all",area="unique",iter=1:nits))
SSBHCR                    <- FLQuant(NA, dimnames=list(age="all",year=projPeriod,unit="unique",season=c("FcY","CtY"),area="unique",iter=1:nits))


#------------------------------------------------------------------------------#
# 2) Start running the MSE
#------------------------------------------------------------------------------#

start.time          <- Sys.time()
for (iYr in an(projPeriod)){
  cat(iYr,"\n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))

  #----------------------------------------
  # define year names
  #----------------------------------------
  TaY <- ac(iYr-1)  # terminal year in the assessment
  ImY <- ac(iYr)    # intermediate year in the short term forecast, ie, current year
  FcY <- ac(iYr+1)  # year for which the advice is given, ie, forecast two years ahead the last year in the assessment
  FuY <- c(ImY,FcY) # combination of future years

  #----------------------------------------
  # update the biol number at age in ImY
  #----------------------------------------

  #- Define mortality rates for iYr-1 to calculate survivors to iYr
  m           <- m(biol)[,ac(iYr-1),,]
  z           <- areaSums(landings.sel(fishery)[,ac(iYr-1),,,,]) + m

  #- Update biological model to iYr
  survivors   <- stock.n(biol)[,ac(iYr-1)] * exp(-z)
  stock.n(biol)[ac((range(biol,"min")+1):range(biol,"max")),ac(iYr),,] <- survivors[-dim(survivors)[1],,,,,]@.Data
  biol@harvest[,ac(iYr-1)] <- areaSums(landings.sel(fishery)[,ac(iYr-1),,,,])

  #- Update recruitment
  recruitBio <- array( 0, dim=c(1,nits)) # initialize array
  for(idxIter in 1:nits){
    if(idxIter %in% itersSR)
      recruitBio[idxIter] <- ifelse(  c(ssb(iter(biol[,ac(iYr-1)],idxIter)))<=params(iter(biol.sr,idxIter))["b"],
                                      params(iter(biol.sr,idxIter))["a"] * c(ssb(iter(biol[,ac(iYr-1)],idxIter))),
                                      params(iter(biol.sr,idxIter))["a"] * params(iter(biol.sr,idxIter))["b"])
    if(idxIter %in% itersRI)
      recruitBio[idxIter] <- params(iter(biol.sr,idxIter))["a"] * c(ssb(iter(biol[,ac(iYr-1)],idxIter))) * exp(-params(iter(biol.sr,idxIter))["b"] * c(ssb(iter(biol[,ac(iYr-1)],idxIter))))
  }
  recruitBio     <- recruitBio * exp(sr.res[,ac(iYr),drop=T])
  stock.n(biol)[1,ac(iYr)] <- recruitBio

  #- Plusgroup
  if (!is.na(range(biol,"plusgroup"))){
    stock.n(biol)[ac(range(biol,"max")),ac(iYr),] <- stock.n(biol)[ac(range(biol,"max")),ac(iYr),] + survivors[ac(range(biol,"max"))]
  }

  #- Apply process error per cohort age 1 to 8
  for(idxAge in 2:nAges){
    biol@stock.n[idxAge,ac(iYr)] <- biol@stock.n[idxAge,ac(iYr)]*varProccError[idxAge-1,ac(an(ac(iYr))-(idxAge-1))]
  }

  cat("\n Finished biology \n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))

  #- Update fishery to year iYr-1
  landings.n(fishery)[,ac(iYr-1)]     <- sweep(sweep(landings.sel(fishery)[,ac(iYr-1),,,,],c(1:4,6),z,"/"),c(1:4,6),stock.n(biol)[,ac(iYr-1)]*(1-exp(-z)),"*")
  catch.n(biol)[,ac(iYr-1)]           <- areaSums(landings.n(fishery[,ac(iYr-1)]))
  print(computeLandings(fishery[,ac(iYr-1)])/TAC[,ac(iYr-1)])

  #-------------------------------------------------------------------------------
  # Assessment
  #-------------------------------------------------------------------------------

  # filter stock object up to intermediate year to include biological variables
  stkAssessment            <- window(biol,end=an(TaY))
  stkAssessment@catch.n    <- stkAssessment@catch.n * catchVar[,ac(histMinYr:TaY),,1]
  stkAssessment@landings.n <- stkAssessment@catch.n
  stkAssessment@landings   <- computeLandings(stkAssessment)

  # smooth M prior to running the assessment, median filter of order 5
  for(idxIter in 1:nits){
    for(idxAge in 1:nAges){
      stkAssessment@m[idxAge,,,,,idxIter] <- fitted(loess(data ~ year,data=data.frame(data=stkAssessment@m[idxAge,,,,,idxIter,drop=T],year=histMinYr:TaY),span=0.5))
    }
  }

  # update surveys
  for(idxSurvey in surveyNames){
    agesSurvey  <- an(rownames(surveys[[idxSurvey]]@index))
    yearSurvey  <- an(colnames(surveys[[idxSurvey]]@index))
    surveyProp  <- mean(c(surveys[[idxSurvey]]@range[6],surveys[[idxSurvey]]@range[7]))
    surveys[[idxSurvey]]@index[,TaY] <- surveyVars[ac(agesSurvey),TaY,idxSurvey,'catchabilities']*
                                        exp(-z[ac(agesSurvey),TaY]*surveyProp)*
                                        biol@stock.n[ac(agesSurvey),TaY]*surveyVars[ac(agesSurvey),TaY,idxSurvey,'residuals']
  }

  stkAssessment.tun                 <- window(surveys,end=an(TaY)+1)
  stkAssessment.tun[["HERAS"]]      <- window(surveys[["HERAS"]],end=an(TaY))
  stkAssessment.tun[["IBTS-Q3"]]    <- window(surveys[["IBTS-Q3"]],end=an(TaY))

  stkAssessment.ctrl@range[5]       <- an(TaY)+1
  stkAssessment.ctrl@residuals      <- F

  # Run stock assessment
  ret <- MSE_assessment(  stkAssessment,
                          stkAssessment.tun,
                          stkAssessment.ctrl,
                          escapeRuns,
                          stkAssessment.init)
  stkAssessment.init      <- ret$resInit
  escapeRuns              <- ret$escapeRuns
  stkAssessment           <- ret$stk

  print(escapeRuns)
  print(stkAssessment@stock.n[,ac(TaY)] / biol@stock.n[,ac(TaY)])

  cat("\n Finished stock assessment \n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))

  #-------------------------------------------------------------------------------
  # Forecast
  #-------------------------------------------------------------------------------
  #(iStocks,iFishery,iYr,iTAC,iHistMaxYr,mpPoints,managementRule)
  projNSAS                  <- projectNSH(stkAssessment,fishery,iYr,TAC,histMaxYr,referencePoints,managementRule)
  TAC[,FcY,,,c("A","B")]    <- projNSAS$TAC[,,,,c("A","B")]
  FHCR[,FcY]                <- projNSAS$Fbar
  SSBHCR[,FcY,,"FcY"]       <- projNSAS$SSB$FcY #store HCR SSB in the forecast year
  SSBHCR[,FcY,,"CtY"]       <- projNSAS$SSB$CtY #store HCR SSB in the continuation year

  cat("\n Finished forecast \n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))

  #-------------------------------------------------------------------------------
  # TAC to real F again
  #-------------------------------------------------------------------------------

  #-Calculate effort accordingly (assuming constant catchability)
  
  #- Get a good starting condition
  #mults <- matrix(NA,nrow=nits,ncol=4)
  #for(idxIter in 1:nits)
  #  mults[idxIter,] <- optim(par=runif(4),fn=TAC2sel_V2,iYr=ImY,iBiol=biol[,ImY],iFishery=fishery[,ImY],iTAC=TAC[,ImY],catchVar=catchVar,TAC_var=TAC_var,iTer=idxIter,control=list(maxit=1000),lower=rep(1e-8,4),method="L-BFGS-B")$par
  CATCH[,ImY,,,"A"]        <- TAC[,ImY,,,"A"] + TAC_var[ImY,,'Ctransfer'] * TAC[,ImY,,,"C"]
  CATCH[,ImY,,,"B"]        <- TAC[,ImY,,,"B",drop=T] * TAC_var[ImY,,'Buptake',drop=T]

  multso <- matrix(NA,nrow=nits,ncol=4)
  for(idxIter in 1:nits)
    multso[idxIter,] <- nls.lm(par=runif(4),fn=TAC2sel,iYr=ImY,iBiol=biol[,ImY],iFishery=fishery[,ImY],iTAC=CATCH[,ImY],catchVar=catchVar,TAC_var=TAC_var,iTer=idxIter,control=nls.lm.control(maxiter=1000),lower=rep(1e-8,4))$par

  #Check for very high F
  idx <- which(quantMeans(sweep(landings.sel(fishery[,ac(ImY)]),3:6,t(multso),"*")[ac(2:6),,,,"A"]) > 5)
  if(length(idx)>0){
    print(idx)
    fishery@landings.sel[,ac(ImY),,,,-idx] <- sweep(landings.sel(fishery[,ac(ImY),,,,-idx]),3:6,t(multso)[,-idx],"*")
    fishery@landings.sel[,ac(ImY),,,,idx]  <- landings.sel(fishery[,ac(an(ImY)-1),,,,idx])
  } else {
    #When setting a very high F this indicates that the optimiser doesn't work well, so replace with last year F
    fishery@landings.sel[,ac(ImY)] <- sweep(landings.sel(fishery[,ac(ImY)]),3:6,t(multso),"*")
  }

  cat("\n Finished effort calc \n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
}
save.image(file=paste(outPath,runName,".RData",sep=""))


