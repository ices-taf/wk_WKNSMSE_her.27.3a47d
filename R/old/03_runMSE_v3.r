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

rm(list=ls())

library(FLSAM)
library(FLEDA)
library(minpack.lm)  # install.packages("minpack.lm")
library(stats)

# define path to directory
#path          <- "D:/Work/Herring MSE/NSAS/"
#path              <- "D:/git/wk_WKNSMSE_her.27.3a47d/R/"
#path              <- "F:/WKNSMSE/wk_WKNSMSE_her.27.3a47d/R"
path <- 'E:/wk_WKNSMSE_her.27.3a47d/R'
assessment_name   <- "NSAS_WKNSMSE2018"
try(setwd(path),silent=TRUE)

# paths to different subfolders
dataPath      <- file.path(".","data/")
outPath       <- file.path(".","results/")
scriptPath    <- file.path(".","side_scripts/")
functionPath  <- file.path(".","functions/")

source(file.path(functionPath,"stf_ImY.R"))
source(file.path(functionPath,"optF_ImY.R"))
source(file.path(functionPath,"stf_FcY.R"))
source(file.path(functionPath,"optF_FcY.R"))
source(file.path(functionPath,"MSE_assessment.R"))

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

nits                <- 10
# load object
load(file.path(outPath,paste0(assessment_name,'_init_MSE_',ac(nits),'.RData')))
stkAssessement.ctrl <- NSH.ctrl

# load MSE parameters
load(file.path(outPath,paste0(assessment_name,'_parameters_MSE_',ac(nits),'.RData')))

strFleet    <- c('A','B','C','D')
nFleets     <- length(strFleet)
nAges       <- dim(biol)[1]
surveyNames <- names(surveys)


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
                        Ftarget = 0.26,
                        F01 = 0.05,
                        Btrigger = 1400000)

managementRule  <- list(HCR = "A",
                        TACIAV="A", #"A","A+B","NA"
                        BB = NA)

#------------------------------------------------------------------------------#
# 3) Define TACs for A, B and D fleets. 
# TACs for A and B fleets are taken out of HAWG2018. This needs updating
#
# Note 1: The Cfleet is defined as a proportion of F.
# Note 2: TAC for C and D fleets are for the WB
# Note 3: TACs for the D fleet is kept constant for future years
#------------------------------------------------------------------------------#

TAC                                   <- FLQuant(NA,dimnames=list(age='all',
                                                                  year=histMinYr:(futureMaxYr+3),
                                                                  unit=c('TAC'),
                                                                  season='all',
                                                                  area=c('A','B','C','D'),
                                                                  iter=1:nits))

# TAC for A fleet in NS
TAC_A                     <- read.table(file.path(dataPath,'TAC_A.csv'),sep = ",")
TAC[,ac(TAC_A[,1]),,,"A"] <- TAC_A[,2]
# TAC for B fleet in NS
TAC_B                     <- read.table(file.path(dataPath,'TAC_B.csv'),sep = ",")
TAC[,ac(TAC_B[,1]),,,"B"] <- TAC_B[,2]
# TAC for D fleet in WB
TAC_D                     <- read.table(file.path(dataPath,'TAC_D.csv'),sep = ",")
TAC[,ac(TAC_D[,1]),,,"D"] <- TAC_D[,2]
TAC[,ac((max(TAC_D[,1])+1):(futureMaxYr+3)),,,"D"] <- TAC[,ac(max(TAC_D[,1])),,,"D"] # TAC is fixed for the D fleet

# TAC for C fleet in WB
TAC_C                     <- read.table(file.path(dataPath,'TAC_C.csv'),sep = ",")
TAC[,ac(TAC_C[,1]),,,"C"] <- TAC_C[,2]
# fixed TAC C in WB (used for transfer to the A fleet)
TAC[,ac((max(TAC_C[,1])+1):(futureMaxYr+3)),,,"C"] <- TAC_C[dim(TAC_C)[1],2]
# randomize TAC C in WB (used for transfer to the A fleet)
#TAC[,ac((max(TAC_C[,1])+1):(futureMaxYr+3)),,,"C"] <- rnorm(length((max(TAC_C[,1])+1):(futureMaxYr+3)),
#                                                            mean(TAC_C[,2]),
#                                                           sd(TAC_C[,2])) # set random TAC for the C fleet in WB

# setup variables for transfer, uptake and split
uptakeFleets <- read.table(file.path(dataPath,'over_underfishing2017.csv'),sep = ",")
# need to input the split for the D fleet - ask Henrik

# transfer from C fleet TAC to fleet A
Ctransfer   <- matrix(runif((length(projPeriod)+3)*nits,min=0.4, max=0.5),nrow=nits,ncol=length(projPeriod)+3)    # Transfer of TAC from IIIa to IVa for C fleet in assessment year. Set between 0.4 and 0.5
# update for the D fleeta
Duptake     <- matrix(1,nrow=nits,ncol=length(projPeriod)+3)    # assume full uptake for the D fleet
Dsplit      <- matrix(rnorm((length(projPeriod)+3)*nits,mean=0.6,sd=0.1),nrow=nits,ncol=length(projPeriod)+3)
# update for the B fleet
Buptake     <- matrix(rnorm ((length(projPeriod)+3)*nits,
                      mean(an(as.vector(uptakeFleets[2:16,3])),na.rm=TRUE), # mean over available historical values
                      sd(an(as.vector(uptakeFleets[2:16,3])),na.rm=TRUE)/2),   # sd over available historical values
                      nrow=nits,ncol=length(projPeriod)+3)
TAC_var     <- array(NA,
                     dim=c(length(projPeriod)+3,nits,
                           4),
                     dimnames=list('years' = ac(an(projPeriod)[1]:(an(projPeriod)[length(projPeriod)]+3)),
                                   'iter' = 1:nits,
                                   'var' = c('Ctransfer','Duptake','Dsplit','Buptake')))
TAC_var[,,'Ctransfer'] <- t(Ctransfer)
TAC_var[,,'Duptake']   <- t(Duptake)
TAC_var[,,'Dsplit']    <- t(Dsplit)
TAC_var[,,'Buptake']   <- t(Buptake)

#------------------------------------------------------------------------------#
# 5) run stf to get TAC in 2019. This is because the MP is different than what
# was calculated during HAWG2018
#
# This is used to update biol which is the bio object.
#
# Note: the recruitment is inferred from SAM. This needs to be changed with a 
# stock recruitment relationship.
#------------------------------------------------------------------------------#

# create stf object and calculate stock in ImY
##############!!!!!!!!!!!!Modify path to function folder to make this work!!!!!!!!!!!!#########
stf <- stf_ImY( stkAssessement,
                fishery,
                TAC,
                TAC_var,
                catchVar[1,projPeriod,,'FCprop'],
                recFuture, # recruitment is know in 2018 for the initial assessment
                c('2018','2019','2020'))

fishery@landings.sel[,'2018'] <- stf@harvest[,'2018']
fishery@landings.n[,'2018']   <- stf@landings.n[,'2018']
fishery@landings[,'2018']     <- stf@landings[,'2018']

#------------------------------------------------------------------------------#
# 6) Start running the MSE
# MSE starts in 2018
#------------------------------------------------------------------------------#

F2plusIter  <-  array(NA,dim=c(length(projPeriod),nits))
F01Iter     <-  array(NA,dim=c(length(projPeriod),nits))

escapeRuns <- numeric()

start.time          <- Sys.time()
for (iYr in an(projPeriod)){
  print(iYr)
  cat(iYr,"\n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
  
  # define years for fisheries
  DtY <- ac(iYr)  # terminal year in the assessment
  ImY <- ac(an(iYr)+1)    # intermediate year in the short term forecast, ie, current year
  FcY <- ac(an(iYr)+2)  # year for which the advice is given, ie, forecast two years ahead the last year in the assessment
  CtY <- ac(an(DtY)+3)             #Continuation year
  FuY <- c(ImY,FcY,CtY) # combination of future years
  
  #-------------------------------------------------------------------------------
  # update Biology object.
  #
  # Need to include S-R function. For now, this is using the recruitment estimation
  # from SAM.
  #
  #-------------------------------------------------------------------------------
  ###################################################################################
  ############################ Recruitment ##########################################
  ###################################################################################
  #- in the OM predict recruitment following either ricker or segreg and multiply that prediction with sr.res
  # e.g. yr = 2018, then ssb of 2017 produced rec in 2018
  
  recruitBio <- array( 0, dim=c(1,nits)) # initialize array
  
  for(idxIter in 1:nits){
    if(idxIter %in% itersSR)
      recruitBio[idxIter] <- ifelse(  c(ssb(iter(biol[,ac(iYr-1)],idxIter)))<=params(iter(biol.sr,idxIter))["b"],
                                      params(iter(biol.sr,idxIter))["a"] * c(ssb(iter(biol[,ac(iYr-1)],idxIter))),
                                      params(iter(biol.sr,idxIter))["a"] * params(iter(biol.sr,idxIter))["b"])
    if(idxIter %in% itersRI)  
      recruitBio[idxIter] <- params(iter(biol.sr,idxIter))["a"] * c(ssb(iter(biol[,ac(iYr-1)],idxIter))) * exp(-params(iter(biol.sr,idxIter))["b"] * c(ssb(iter(biol[,ac(iYr-1)],idxIter))))
  }
  
  recruitBio     <- recruitBio * exp(drop(sr.res[,DtY]))
  
  ###################################################################################
  ############################ Update F ##########################################
  ###################################################################################
  biol@harvest[,DtY]  <- apply(fishery@landings.sel[,DtY],c(1,2,4,5,6),'sum')
  
  Z <- biol@harvest[,ac(c(an(DtY)-1,DtY))] + biol@m[,ac(c(an(DtY)-1,DtY))]
  
  # compute stock
  # propagate stock number with Z, only fill first slot
  survivors                           <- drop(biol@stock.n[,ac(an(DtY)-1),1])*exp(-drop(Z[,ac(an(DtY)-1)])) # stock.n is the same for all fleets in the stf object, taking first element
  biol@stock.n[2:nAges,DtY]           <- survivors[1:(nAges-1),]
  biol@stock.n[nAges,DtY]             <- biol@stock.n[nAges,DtY,1] + survivors[nAges,]
  biol@stock.n[1,DtY]                 <- recruitBio
  # apply process error per cohort age 1 to 8
  for(idxAge in 2:nAges){
    biol@stock.n[idxAge,DtY] <- biol@stock.n[idxAge,DtY]*varProccError[idxAge-1,ac(an(DtY)-(idxAge-1))]
  }
  
  biol@stock[,DtY]                    <- computeStock(biol[,DtY])
  
  # compute catch and landings
  biol@catch.n[,DtY]      <- drop(biol@harvest[,DtY])*drop(biol@stock.n[,DtY])*drop((1-exp(-Z[,DtY]))/Z[,DtY])*drop(catchVar[,DtY,,'residuals'])
  biol@landings.n[,ImY]   <- drop(biol@harvest[,DtY])*drop(biol@stock.n[,DtY])*drop((1-exp(-Z[,DtY]))/Z[,DtY])
  
  biol@catch       <- computeCatch(biol)
  biol@landings    <- computeLandings(biol)
  
  # update surveys
  for(idxSurvey in surveyNames){
    agesSurvey  <- an(rownames(surveys[[idxSurvey]]@index))
    yearSurvey  <- an(colnames(surveys[[idxSurvey]]@index))
    
    surveyProp  <- mean(c(surveys[[idxSurvey]]@range[6],surveys[[idxSurvey]]@range[7]))
    
    surveys[[idxSurvey]]@index[,DtY] <- drop(surveyVars[ac(agesSurvey),DtY,idxSurvey,'catchabilities'])*
                                        drop(exp(-Z[ac(agesSurvey),DtY]*surveyProp)*
                                        biol@stock.n[ac(agesSurvey),DtY])*drop(surveyVars[ac(agesSurvey),DtY,idxSurvey,'residuals'])
  }
  
  cat("\n Finished biology \n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
  
  #-------------------------------------------------------------------------------
  # Assessment
  #-------------------------------------------------------------------------------
  
  # filter stock object up to intermediate year to include biological variables
  stkAssessement <- window(  biol,
                             start=an(fullPeriod[1]),
                             end=an(DtY))
  
  stkAssessement@landings.n <- stkAssessement@catch.n
  stkAssessement@landings   <- computeLandings(stkAssessement)

  # smooth M prior to running the assessment, median filter of order 5
  for(idxIter in 1:nits){
    for(idxAge in 1:nAges){
      stkAssessement@m[idxAge,,,,,idxIter] <- runmed(stkAssessement@m[idxAge,,,,,idxIter],k=5)
    }
  }
  
  # select tuning object for assessment. filter up to terminal year
  stkAssessement.tun <- surveys
  
  for(idxSurvey in surveyNames){
    minYearSurvey     <- min(as.numeric(colnames(stkAssessement.tun[[idxSurvey]]@index)))
    stkAssessement.tun[[idxSurvey]] <- window( stkAssessement.tun[[idxSurvey]],
                                               start=minYearSurvey,
                                               end=an(DtY))
    stkAssessement.tun[[idxSurvey]]@effort[,DtY] <- 1
    stkAssessement.tun[[idxSurvey]]@catch.n <- stkAssessement.tun[[idxSurvey]]@index
  }
  
  stkAssessement.ctrl@range[5] <- an(DtY)
  stkAssessement.ctrl@residuals <- F
  
  # stock assessment
  if(DtY == projPeriod[1]){
    ret <- MSE_assessment(  stkAssessement,
                            stkAssessement.tun,
                            stkAssessement.ctrl,
                            escapeRuns)
    stkAssessement.init <- ret$resInit
    ret <- MSE_assessment(  stkAssessement,
                            stkAssessement.tun,
                            stkAssessement.ctrl,
                            escapeRuns,
                            stkAssessement.init)
   }else{
    ret <- MSE_assessment(  stkAssessement,
                            stkAssessement.tun,
                            stkAssessement.ctrl,
                            escapeRuns,
                            stkAssessement.init)
  }
  
  # update stock assessment with results
  escapeRuns      <- c(escapeRuns,ret$escapeRuns)
  
  stkAssessement  <- ret$stk
  stkAssessement@stock   <- computeStock(stkAssessement)
  
  cat("\n Finished stock assessment \n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
  
  #-------------------------------------------------------------------------------
  # Fishery, i.e. Short term forecast
  #-------------------------------------------------------------------------------
  
  ############# Recruitment #############
  # Recruitment for future years are computed trough a 10 year weighted average
  # recruitment is the same for Intermediate, forecast and continuation years
  
  # finding standard deviation on recruitment for the weighted average.
  # This is using the iterations (as no residuals)
  recSd <- array(NA,dim=c(1,dim(stkAssessement)[2]),
                 dimnames = list('sdlogR','year' = an(fullPeriod[1]):an(DtY)))
  for(idxYear in 1:(dim(stkAssessement)[2])){
    recSd[idxYear] <- sd(log(rec(stkAssessement[,idxYear])))
  }
  
  # recruitment for stf as weighted mean over 10 years
  recFuture <- exp( apply( log(rec(stkAssessement)[,ac((an(DtY)-10):DtY)]),
                           6,
                           weighted.mean,
                           w=1/recSd[,ac((an(DtY)-10):DtY)],
                           na.rm=T))
  recFuture <- as.array(drop(recFuture))
  
  ############# setting up stf and compute Intermediate year #############
  ##############!!!!!!!!!!!!Modify path to function folder to make this work!!!!!!!!!!!!#########
  cat(paste("\n Time running Start - sft ImY",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))

  stf <- stf_ImY( stkAssessement,
                  fishery,
                  TAC,
                  TAC_var,
                  catchVar[1,projPeriod,,'FCprop'],
                  recFuture,
                  FuY)
  cat(paste("\n Time running End - sft ImY",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
  
  # update fishery object
  fishery@landings.sel[,ImY] <- stf@harvest[,ImY]
  fishery@landings.n[,ImY]   <- stf@landings.n[,ImY]
  fishery@landings[,ImY]     <- stf@landings[,ImY]

  
  ############# Compute Forecast year #############
  # use HCR and TAC C and D fleets to define TAcs for FcY
  ##############!!!!!!!!!!!!Modify path to function folder to make this work!!!!!!!!!!!!#########
  cat(paste("\n Time running Start - sft FcY",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))

  stf <- stf_FcY( stf,
                  fishery,
                  TAC,
                  TAC_var,
                  catchVar[1,projPeriod,,'FCprop'],
                  recFuture,
                  HCR,
                  referencePoints,
                  FuY)
  
  cat(paste("\n Time running End - sft FcY",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
  
  # update TAC
  TAC[,FcY,,,'A'] <- stf@catch[,FcY,'A']
  TAC[,FcY,,,'B'] <- stf@catch[,FcY,'B']

  cat("\n Finished forecast \n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
}




plot(apply(drop(iter(biol@stock.n[,ac(2016:2019)],8)),2,'sum'))

plot(iter(biol@stock[,ac(2016:2019)],8))
ssb(iter(biol,8))
N <- iter(biol@stock.n,8)
wts <- iter(biol@stock.wt,8)
Ftot <- iter(biol@harvest,8)
F.spwn <- iter(biol@harvest.spwn,8)
M <- iter(biol@m,8)
M.spwn <- iter(biol@m.spwn,8)
mats <- iter(biol@mat,8)
ssb(iter(biol,8))
ssb_FcY               <- sum( drop(N[,ac(iYr)]*wts[,ac(iYr)])*
                                exp(-drop(Ftot[,ac(iYr)]*F.spwn[,ac(iYr)]-M[,ac(iYr)]*M.spwn[,ac(iYr)]))*drop(mats[,ac(iYr)]))

#save.image(file=paste(outPath,runName,"_",settings$RecRegime,".RData",sep=""))


