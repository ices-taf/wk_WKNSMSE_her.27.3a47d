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

source(file.path(functionPath,"optF_TACdiff.R"))
source(file.path(functionPath,"rescaleF_Ftargets.R"))
source(file.path(functionPath,"stf_ImY.R"))
source(file.path(functionPath,"stf_FcY.R"))
source(file.path(functionPath,"MSE_assessment.R"))

#-------------------------------------------------------------------------------
# 2) Initialize
#
# Define MSE parameters, 
# load objects initialized previously
#   * Biology
#     - stocks
#     - surveys
#   * Fisheries
#     - FCProp
#     - catchD
#     - F sel: FAsel, FCsel, FBDsel
#-------------------------------------------------------------------------------

# load object
load(file.path(outPath,paste0(assessment_name,'_init_MSE.RData')))
stkAssessement.ctrl <- NSH.ctrl

# load MSE parameters
load(file.path(outPath,paste0(assessment_name,'_parameters_MSE.RData')))

strFleet    <- c('A','B','C','D')
nFleets     <- length(strFleet)
nAges       <- dim(stocks)[1]
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

HCR <- 'A'

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
Ctransfer   <- runif(length(projPeriod),min=0.4, max=0.5)    # Transfer of TAC from IIIa to IVa for C fleet in assessment year. Set between 0.4 and 0.5
# update for the D fleeta
Duptake     <- rep(1, length(projPeriod)) # assume full uptake for the D fleet
Dsplit      <- rep(0.60, length(projPeriod)) # NSAS/WBSS split randomization based on historical records. Fixed for now, need to contact henrik
# update for the B fleet
Buptake     <- rnorm (length(projPeriod), 
                      mean(an(as.vector(uptakeFleets[2:16,3])),na.rm=TRUE), # mean over available historical values
                      sd(an(as.vector(uptakeFleets[2:16,3])),na.rm=TRUE))   # sd over available historical values

TAC_var     <- array(NA,
                     dim=c(length(projPeriod),
                           4),
                     dimnames=list('years' = projPeriod,
                                   'var' = c('Ctransfer','Duptake','Dsplit','Buptake')))
TAC_var[,'Ctransfer'] <- Ctransfer
TAC_var[,'Duptake']   <- Duptake
TAC_var[,'Dsplit']    <- Dsplit
TAC_var[,'Buptake']   <- Buptake

#------------------------------------------------------------------------------#
# 5) run stf to get TAC in 2019. This is because the MP is different than what
# was calculated during HAWG2018
#
# This is used to update stocks which is the bio object.
#
# Note: the recruitment is inferred from SAM. This needs to be changed with a 
# stock recruitment relationship.
#------------------------------------------------------------------------------#

recFuture <- as.array(subset(rec(stkAssessement),year==2018)$data)

# create stf object and calculate stock in ImY
stf <- stf_ImY( stocks,
                fisheryFuture,
                TAC,
                TAC_var,
                FCPropIts,
                recFuture, # recruitment is know in 2018 for the initial assessment
                c('2018','2019','2020'))

#------------------------------------------------------------------------------#
# 6) Start running the MSE
# MSE starts in 2019
#------------------------------------------------------------------------------#

escapeRuns <- numeric()

start.time          <- Sys.time()
bunit               <- dimnames(biol@n)$unit
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
  
  # calculate residuals to add to the catches
  resiCatch <- array(NA,dim=c(nAges,nits),
                     dimnames = list('ages' = 0:(nAges-1),'iter' = 1:nits))
  for(idxIter in 1:nits){
    resiCatch[,idxIter] <- rnorm( length(varCatchMat[,2,idxIter]), 
                                  (1:length(varCatchMat[,2,idxIter]))*0, 
                                  varCatchMat[,2,idxIter]) # residual to add to the catch
  }
  
  # catch numbers with added residuals
  stocks@catch.n[,DtY]  <- drop(apply(stf@catch.n[,DtY],c(1,2,4,5,6),'sum'))*exp(resiCatch)
  stocks@catch          <- computeCatch(stocks)
  
  # landing numbers
  stocks@landings.n[,DtY]  <- apply(stf@landings.n[,DtY],c(1,2,4,5,6),'sum')
  stocks@landings          <- computeLandings(stocks)
  
  # stock numbers, stock number are the same for all flets, we use 1st field here
  stocks@stock.n[,DtY]            <- stf@stock.n[,DtY,1]
  stocks@stock[,DtY]              <- computeStock(stf[,DtY,1])
  
  # update F
  stocks@harvest[,DtY]           <- apply(stf@harvest[,DtY],c(1,2,4,5,6),'sum')
  
  # update surveys
  
  for(idxSurvey in surveyNames){
    qSelect     <- qMat[rownames(qMat) == idxSurvey,,,drop=FALSE]
    varSurv     <- varSurvMat[rownames(qMat) == idxSurvey,,,drop=FALSE]
    agesSurvey  <- an(rownames(surveys[[idxSurvey]]@index))
    nAgesSurvey <- length(agesSurvey)
    yearSurvey  <- an(colnames(surveys[[idxSurvey]]@index))
    
    resiSurv    <- array(NA,
                         dim=c(nAgesSurvey,nits),
                         dimnames = list('ages' = agesSurvey,'iter' = 1:nits))
    # residuals at age for the selected survey
    for(idxIter in 1:nits){
      resiSurv[,idxIter] <- rnorm( length(varSurv[,2,idxIter]), 
                                    (1:length(varSurv[,2,idxIter]))*0, 
                                    varSurv[,2,idxIter]) # residual to add to the catch
    }
    
    # total mortality
    Z <- stocks@m[ac(agesSurvey),DtY] + stocks@harvest[ac(agesSurvey),DtY]
    surveyProp  <- mean(c(surveys[[idxSurvey]]@range[6],surveys[[idxSurvey]]@range[7]))
    
    surveys[[idxSurvey]]@index[,DtY] <- qSelect[,2,]*drop(exp(-Z*surveyProp)*stocks@stock.n[ac(agesSurvey),DtY])*exp(resiSurv)
  }
  
  cat("\n Finished biology \n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
  
  #-------------------------------------------------------------------------------
  # Assessment
  #-------------------------------------------------------------------------------
  
  # filter stock object up to intermediate year to include biological variables
  stkAssessement <- window(  stocks,
                             start=an(fullPeriod[1]),
                             end=an(DtY))
  
  stkAssessement@catch.n  <- stkAssessement@landings.n
  stkAssessement@catch    <- computeCatch(stkAssessement)
  
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
  
  stkAssessement <- window(  stocks,
                             start=an(fullPeriod[1]),
                             end=an(CtY))
  
  stkAssessement[,ac(an(fullPeriod[1]):an(DtY))] <- ret$stk
  
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
  recSd <- array(NA,dim=c(1,dim(stkAssessement)[2]-3),
                 dimnames = list('sdlogR','year' = an(fullPeriod[1]):an(DtY)))
  for(idxYear in 1:(dim(stkAssessement)[2]-3)){
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
  stf <- stf_ImY( stkAssessement,
                  fisheryFuture,
                  TAC,
                  TAC_var,
                  FCPropIts,
                  recFuture,
                  FuY)

  
  ############# Compute Forecast year #############
  # use HCR and TAC C and D fleets to define TAcs for FcY
  stf <- stf_FcY( stf,
                  fisheryFuture,
                  TAC,
                  TAC_var,
                  FCPropIts,
                  recFuture,
                  HCR,
                  referencePoints,
                  FuY)
  
  TAC[,FcY,,,'A'] <- stf@catch[,FcY,'A']
  TAC[,FcY,,,'B'] <- stf@catch[,FcY,'B']

  cat("\n Finished forecast \n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
}
save.image(file=paste(outPath,runName,"_",settings$RecRegime,".RData",sep=""))


