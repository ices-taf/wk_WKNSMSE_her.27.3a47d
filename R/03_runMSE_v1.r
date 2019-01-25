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
path              <- "D:/git/wk_WKNSMSE_her.27.3a47d/R/"
#path              <- "F:/WKNSMSE/wk_WKNSMSE_her.27.3a47d/R"
assessment_name   <- "NSAS_WKNSMSE2018"
try(setwd(path),silent=TRUE)

# paths to different subfolders
dataPath      <- file.path(".","data/")
outPath       <- file.path(".","results/")
scriptPath    <- file.path(".","side_scripts/")
functionPath  <- file.path(".","functions/")

source(file.path(functionPath,"optF_TACdiff.R"))
source(file.path(functionPath,"stf_ImY.R"))

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

# load MSE parameters
load(file.path(outPath,paste0(assessment_name,'_parameters_MSE.RData')))

strFleet  <- c('A','B','C','D')
nFleets   <- length(strFleet)
nAges     <- dim(stocks[[1]])[1]


#-------------------------------------------------------------------------------
# 2) Organize objects
# here, 2 main objects are defined:
#   * biology
#   * fishery object for each fleet
#-------------------------------------------------------------------------------

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
# 4) run stf to get TAC in 2019. This is because the MP is different than what
# was calculated during HAWG2018
#
# This is used to update stocks which is the bio object.
#
# Note: the recruitment is inferred from SAM. This needs to be changed with a 
# stock recruitment relationship.
#------------------------------------------------------------------------------#

# create stf object and calculate stock in ImY
stf <- stf_ImY( NSH.sim,
                stocks,
                fisheryFuture,
                TAC,
                TAC_var,
                FCPropIts,
                c('2018','2019','2020'))

#------------------------------------------------------------------------------#
# 5) Start running the MSE
# MSE starts in 2019
#------------------------------------------------------------------------------#

escapeRuns <- numeric()

start.time          <- Sys.time()
bunit               <- dimnames(biol@n)$unit
for (iYr in an(projPeriod)){
  cat(iYr,"\n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
  
  #-------------------------------------------------------------------------------
  # update Biology object.
  #
  # Need to include S-R function. For now, this is using the recruitment estimation
  # from SAM.
  #
  #-------------------------------------------------------------------------------
  
  # Catches are updated using the stf object for ImY, simply adding residuals to 
  # the catches. One need to include the S-R function
  for(idxIter in 1:nits){
    resiCatch <- rnorm( length(varCatchMat[,2,idxIter]), 
                        (1:length(varCatchMat[,2,idxIter]))*0, 
                        varCatchMat[,2,idxIter]) # residual to add to the catch
    
    # catches = landings + deviation
    stocks[[idxIter]]@catch.n[,iYr]  <- rowSums(stf@catch.n[,iYr,,,,idxIter])*exp(resiCatch)
    stocks[[idxIter]]@catch          <- computeCatch(stocks[[idxIter]])
    #stocks[[idxIter]]@catch[,iYr]    <- sum(stocks[[idxIter]]@catch.n[,iYr]*
    #                                        apply(drop(fisheryFuture[,iYr,,'catch.wt',,1]),1,mean))
    #stocks[[idxIter]]@catch[,iYr]    <- sum( rowSums(stf@catch.n[,iYr,,,,idxIter]* # catch.n per fleet
    #                                                      stf[,iYr,,,,idxIter]@catch.wt)* # catch.wt per fleet
    #                                              exp(resiCatch)) # sum catches accross fleets
    #stocks[[idxIter]]@catch.wt[,iYr] <- stocks[[idxIter]]@catch.n[,iYr]/stocks[[idxIter]]@catch[,iYr]
    #stocks[[idxIter]]@catch.wt[,iYr] <- apply(drop(fisheryFuture[,iYr,,'catch.wt',,1]),1,mean) # weight at age as mean accross the fleets
    
    
    stocks[[idxIter]]@landings.n[,iYr]        <- rowSums(stf@catch.n[,iYr,,,,idxIter])
    stocks[[idxIter]]@landings                <- computeLandings(stocks[[idxIter]])
    
    stocks[[idxIter]]@stock.n[,iYr]           <- stf@stock.n[,iYr,1,,,idxIter]
    
    stocks[[idxIter]]@harvest[,iYr]           <- rowSums(stf@harvest[,iYr,,,,idxIter])
    
    surveyNames <- unique(rownames(qMat)) # get all the survey names
    
    for(idxSurvey in 1:length(surveyNames)){
      qSelect     <- subset(qMat[,,idxIter],rownames(qMat) == surveyNames[idxSurvey])
      varSurv     <- subset(varSurvMat[,,idxIter],rownames(varSurvMat) == surveyNames[idxSurvey])
      
      resiSurv <- rnorm(  length(varSurv[,2]), 
                          (1:length(varSurv[,2]))*0, 
                          varSurv[,2]) # residual to add to the catch
      
      # select number of age for the corresponding ages to the survey
      NSelect     <- stocks[[idxIter]]@stock.n[match(as.character(qSelect[,1]),rownames(stocks[[idxIter]]@stock.n)), # filter ages
                                               iYr]  # current year
      NSelect     <- drop(NSelect) # drop dimensions with 1 level
      FSelect     <-  stocks[[idxIter]]@harvest[match(as.character(qSelect[,1]),rownames(stocks[[idxIter]]@harvest)), # filter ages for F
                                                iYr] # filter year
      FSelect     <- drop(FSelect) # drop dimensions with 1 level
      Z           <-  stocks[[idxIter]]@m[match(as.character(qSelect[,1]),rownames(stocks[[idxIter]]@m)), # filter ages for M
                                          iYr]  # filter years for M
      Z           <- Z + FSelect
      surveyProp  <- mean(c(surveys[[surveyNames[idxSurvey]]]@range[6],surveys[[surveyNames[idxSurvey]]]@range[7]))

      # filling survey index for current year
      surveys[[surveyNames[idxSurvey]]]@index[,iYr,,,,idxIter] <- qSelect[,2]*exp(-Z*surveyProp)*NSelect*exp(resiSurv)
    }
  }
  
  cat("\n Finished biology \n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
  
  #-------------------------------------------------------------------------------
  # Assessment
  #-------------------------------------------------------------------------------
  
  C1 <- dim(stocks[[1]]@catch.wt)
  C1[6]<- 10
  C2 <- C1
  C2[1]<- 1
  
  dmns        <- dimnames(stocks[[1]]@harvest)
  dmns$iter   <- 1:nits
  dmns2        <- dimnames(stocks[[1]]@catch)
  dmns2$iter   <- 1:nits
  
  NSH <- FLStock( catch.wt=FLQuant(dim=C1,dimnames=dmns),
                  catch.n=FLQuant(dim=C1,dimnames=dmns),
                  catch=FLQuant(dim=C2,dimnames=dmns2),
                  landings.wt=FLQuant(dim=C1,dimnames=dmns),
                  landings.n=FLQuant(dim=C1,dimnames=dmns),
                  landings=FLQuant(dim=C2,dimnames=dmns2),
                  discards.wt=FLQuant(dim=C1,dimnames=dmns),
                  discards.n=FLQuant(dim=C1,dimnames=dmns),
                  discards=FLQuant(dim=C2,dimnames=dmns2),
                  stock.wt=FLQuant(dim=C1,dimnames=dmns),
                  stock.n=FLQuant(dim=C1,dimnames=dmns),
                  stock=FLQuant(dim=C2,dimnames=dmns2),
                  m=FLQuant(dim=C1,dimnames=dmns),
                  mat=FLQuant(dim=C1,dimnames=dmns),
                  harvest=FLQuant(dim=C1,dimnames=dmns),
                  harvest.spwn=FLQuant(dim=C1,dimnames=dmns),
                  m.spwn=FLQuant(dim=C1,dimnames=dmns))
  NSH@range <- stocks[[1]]@range
  
  for(idxIter in 1:nits){
    iter(NSH, idxIter) <- stocks[[idxIter]]
  }
  
  NSH <- window(  NSH,
                   start=an(fullPeriod[1]),
                   end=an(iYr))
  
  NSH.tun <- surveys
  
  for(idxSurvey in 1:length(surveys)){
    minYearSurvey     <- min(as.numeric(colnames(NSH.tun[[surveyNames[idxSurvey]]]@index)))
    NSH.tun[[idxSurvey]] <- window( NSH.tun[[idxSurvey]][,,,,,idxIter],
                                     start=minYearSurvey,
                                     end=an(iYr))
    NSH.tun[[idxSurvey]]@effort[,iYr] <- 1
    NSH.tun[[idxSurvey]]@catch.n <- NSH.tun[[idxSurvey]]@index
  }
  
  NSH.ctrl@range[5] <- an(iYr)
  NSH.ctrl@residuals <- F
  
  res <- FLSAM.MSE(NSH,NSH.tun,NSH.ctrl,return.sam=T)
  
  continueRuns        <- which(!(1:dims(NSH)$iter) %in% escapeRuns)
  if (iay == iy)  res <- FLSAM.MSE(NSH,idx0,sam0.ctrl,return.sam=T) #first year only, we can create an FLSAMs object for all iterations similar to the final NSAS 2018 assessment too if that is needed
  if (iay > iy)   res <- FLSAM.MSE(stk0,idx0,sam0.ctrl,starting.sam=res,return.sam=T) #here we add the starting values for each iteration
  
  trouble <- data.frame(iter = 1:it , failure =  unlist(lapply(res , function(x) is.na(x))))
  trouble <- trouble[trouble$failure == T,]
  
  #if we only have one iteration with problems, we try to run that one again, but now without starting values
  if (dim(trouble)[1] == 1){
    tres    <-  try(FLSAM(iter(stk0,trouble$iter),
                          FLIndices(lapply(idx0 , function(x) iter(x,trouble$iter))),
                          sam0.ctrl,silent=T))
    #if that didn't work, we'll just save it as an empty FLSAM object
    if(class(tres)=="try-error"){
      res[[trouble$iter]] <- new("FLSAM")
    } else {
      res[[trouble$iter]] <- tres
    }
  }
  #if we have more iterations with trouble, we can make use of the FLSAM.MSE routine
  if (dim(trouble)[1] > 1){
    resTrouble          <- FLSAM.MSE(iter(stk0,trouble$iter),
                                     FLIndices(lapply(idx0,function(x) iter(x,trouble$iter))),
                                     sam0.ctrl,return.sam=T)
    counter <- 1
    for(ii in trouble$iter){
      if(is.na(resTrouble[[counter]])){
        res[[ii]] <- new("FLSAM")
      } else {
        res[[ii]] <-  resTrouble[[counter]]
      }
      counter <- counter + 1
    }
  }
  
  
  #Update the stock object and store which runs to escape
  for(ii in 1:it){
    if(!is.na(res[[ii]]@harvest[1,1,drop=T])){
      iter(stk0@harvest,ii) <- res[[ii]]@harvest
      iter(stk0@stock.n,ii) <- res[[ii]]@stock.n
    } else {
      escapeRuns <- sort(unique(c(escapeRuns,ii)))
    }
  }
  
  
  
  
  
  
  
  
  
  
  iYr <- '2017'
  NSH.ctrl@range[5] <- an(iYr)
  NSH.ctrl@residuals <- F
  idxIter <- 2
  
  for(idxIter in 1:nits){
    NSH2 <- window(  stocks[[idxIter]],
                    start=an(fullPeriod[1]),
                    end=an(iYr))
    NSH2@catch.n  <- NSH2@landings.n
    NSH2@catch    <- computeCatch(NSH2)
    NSH2@landings <- computeLandings(NSH2)
    
    # run a 7 years median filter to prevent any big discontinuity in M
    #for(idxAge in dimnames(NSH2@m)$age){
    #  NSH2@m[idxAge,] <- runmed(NSH2@m[idxAge,],k=11)
    #  NSH2@stock.wt[idxAge,] <- runmed(NSH2@stock.wt[idxAge,],k=11)
    #  NSH2@mat[idxAge,] <- runmed(NSH2@mat[idxAge,],k=11)
    #}
    
    NSH.tun2 <- surveys
    
    for(idxSurvey in 1:length(surveys)){
      minYearSurvey     <- min(as.numeric(colnames(NSH.tun2[[surveyNames[idxSurvey]]]@index)))
      NSH.tun2[[idxSurvey]] <- window( NSH.tun2[[idxSurvey]][,,,,,idxIter],
                                      start=minYearSurvey,
                                      end=an(iYr))
      NSH.tun2[[idxSurvey]]@effort[,iYr] <- 1
      NSH.tun2[[idxSurvey]]@catch.n <- NSH.tun2[[idxSurvey]]@index
    }

    NSH.sim[[idxIter]] <- FLSAM(NSH2,NSH.tun2,NSH.ctrl)
  }
  
  cat("\n Finished stock assessment \n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
  
  #-------------------------------------------------------------------------------
  # Fishery, i.e. Short term forecast
  #-------------------------------------------------------------------------------
  
  # define years for fisheries
  DtY <- ac(iYr-1)  # terminal year in the assessment
  ImY <- ac(iYr)    # intermediate year in the short term forecast, ie, current year
  FcY <- ac(iYr+1)  # year for which the advice is given, ie, forecast two years ahead the last year in the assessment
  CtY <- ac(an(DtY)+3)             #Continuation year
  FuY <- c(ImY,FcY,CtY) # combination of future years
  
  ############# Recruitment #############
  # Recruitment for intermediate year, forecast year and continuation year
  # intermediate year: taken from SAM
  # forecast and continuation years: 10 years weighted average
 
  
  # apply HCR to obtain F2-6 and F0-1
  
  cat("\n Finished forecast \n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
}
save.image(file=paste(outPath,runName,"_",settings$RecRegime,".RData",sep=""))


