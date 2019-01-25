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

FuY <- c('2018','2019','2020')


ImY <- FuY[1]
FcY <- FuY[2]
CtY <- FuY[3]

############# initialize stf FLStock object #############
dsc         <- "North Sea Herring"
nam         <- "NSAS"
dms         <- dimnames(stocks[[1]]@m)
dms$year    <- ac((an(ImY)-3):an(CtY))
dms$unit    <- c("A","B","C","D")
dms$iter    <- 1:nits

nAges       <- length(dms$age)
nits        <- length(dms$iter)

# Create the stf object 
stf         <- FLStock(name=nam,desc=dsc,FLQuant(NA,dimnames=dms))

for(idxIter in 1:nits){
  for(idxFleet in dms$unit){
    # update stf object with current stock object
    stf[,,idxFleet,,,idxIter]          <- window( stocks[[idxIter]],
                                                  start=an(dms$year)[1],
                                                  end=rev(an(dms$year))[1])
    # update stf object with stock numbers from SAM object
    stf[,,idxFleet,,,idxIter]@stock.n  <- window( NSH.sim[[idxIter]]@stock.n,
                                                  start=an(dms$year)[1],
                                                  end=rev(an(dms$year))[1])
    
    # update catch.wt
    if(idxFleet == 'B' || idxFleet == 'D'){
      stf[,,idxFleet,,,idxIter]@catch.wt <- fisheryFuture[,FuY,'BD','catch.wt',,idxIter]
    }else{
      stf[,,idxFleet,,,idxIter]@catch.wt <- fisheryFuture[,FuY,idxFleet,'catch.wt',,idxIter]
    }
  }
}

# Fill slots that have no meaning for NSAS
stf@discards.n[]          <- 0
stf@discards[]            <- 0
stf@discards.wt[]         <- 0


############################ END TEST

############# Compute F in intermediate year #############

# update intermediate year with results from assessment

source(file.path(functionPath,"optF_TACdiff.R"))

Fscalor <- array( 0, dim=c(nFleets,nits)) # initialize array
for(idxIter in 1:nits){
  # find F for the different fleets
  Fscalor[,idxIter] <- nls.lm(  par=runif(4), # starting point
                                lower=rep(1e-8,4),
                                upper=NULL,
                                optF_TACdiff,                                  # function to optimize
                                harvest_sf  = NSH.sim[[idxIter]]@harvest,       # single fleet FLStock object
                                catch.wt_mf = fisheryFuture[,,,'catch.wt',,idxIter],       # catch weight at age single fleet. Using stock weight at age for now. How to get catch weight for 2018?
                                stock.n_sf  = NSH.sim[[idxIter]]@stock.n,       # stock at age
                                M           = stocks[[idxIter]]@m,              # natural mortality
                                Fsel        = fisheryFuture[,,,'sel',,idxIter], # selectivity stored as FLQuant object. Normalized between 0 and 1.
                                iYr         = ImY,                              # year of interest
                                TACs        = TAC[,,,,,idxIter],             # TAC FLQuant object for fleets A, B and D
                                FCProp      = FCPropIts[,idxIter],
                                TAC_var     = TAC_var,
                                recruit     = recFuture[idxIter],
                                nls.lm.control(ftol = (.Machine$double.eps),maxiter = 1000), # optimizer control object
                                jac=NULL)$par
}


stf@harvest[,ImY,'A'] <- t(apply(fisheryFuture[,ImY,'A','sel'],1,'*',Fscalor[1,]))
stf@harvest[,ImY,'B'] <- t(apply(fisheryFuture[,ImY,'BD','sel'],1,'*',Fscalor[2,]))
stf@harvest[,ImY,'C'] <- t(apply(fisheryFuture[,ImY,'C','sel'],1,'*',Fscalor[3,]))
stf@harvest[,ImY,'D'] <- t(apply(fisheryFuture[,ImY,'BD','sel'],1,'*',Fscalor[4,]))































recFuture <- as.array(subset(rec(NSH.sim),year==2018)$value)

# create stf object and calculate stock in ImY
stf <- stf_ImY( NSH.sim,
                stocks,
                fisheryFuture,
                TAC,
                TAC_var,
                FCPropIts,
                recFuture, # recruitment is know in 2018 for the initial assessment
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
  
  # Catches are updated using the stf object for ImY, simply adding residuals to 
  # the catches. One need to include the S-R function
  for(idxIter in 1:nits){
    resiCatch <- rnorm( length(varCatchMat[,2,idxIter]), 
                        (1:length(varCatchMat[,2,idxIter]))*0, 
                        varCatchMat[,2,idxIter]) # residual to add to the catch
    
    # catches = landings + deviation
    stocks[[idxIter]]@catch.n[,DtY]  <- rowSums(stf@catch.n[,DtY,,,,idxIter])*exp(resiCatch)
    stocks[[idxIter]]@catch          <- computeCatch(stocks[[idxIter]])
    
    stocks[[idxIter]]@landings.n[,DtY]        <- rowSums(stf@catch.n[,DtY,,,,idxIter])
    stocks[[idxIter]]@landings                <- computeLandings(stocks[[idxIter]])
    
    stocks[[idxIter]]@stock.n[,DtY]           <- stf@stock.n[,DtY,1,,,idxIter]
    
    stocks[[idxIter]]@harvest[,DtY]           <- rowSums(stf@harvest[,DtY,,,,idxIter])
    
    surveyNames <- unique(rownames(qMat)) # get all the survey names
    
    for(idxSurvey in 1:length(surveyNames)){
      qSelect     <- subset(qMat[,,idxIter],rownames(qMat) == surveyNames[idxSurvey])
      varSurv     <- subset(varSurvMat[,,idxIter],rownames(varSurvMat) == surveyNames[idxSurvey])
      
      resiSurv <- rnorm(  length(varSurv[,2]), 
                          (1:length(varSurv[,2]))*0, 
                          varSurv[,2]) # residual to add to the catch
      
      # select number of age for the corresponding ages to the survey
      NSelect     <- stocks[[idxIter]]@stock.n[match(as.character(qSelect[,1]),rownames(stocks[[idxIter]]@stock.n)), # filter ages
                                               DtY]  # current year
      NSelect     <- drop(NSelect) # drop dimensions with 1 level
      FSelect     <-  stocks[[idxIter]]@harvest[match(as.character(qSelect[,1]),rownames(stocks[[idxIter]]@harvest)), # filter ages for F
                                                DtY] # filter year
      FSelect     <- drop(FSelect) # drop dimensions with 1 level
      Z           <-  stocks[[idxIter]]@m[match(as.character(qSelect[,1]),rownames(stocks[[idxIter]]@m)), # filter ages for M
                                          DtY]  # filter years for M
      Z           <- Z + FSelect
      surveyProp  <- mean(c(surveys[[surveyNames[idxSurvey]]]@range[6],surveys[[surveyNames[idxSurvey]]]@range[7]))

      # filling survey index for current year
      surveys[[surveyNames[idxSurvey]]]@index[,DtY,,,,idxIter] <- qSelect[,2]*exp(-Z*surveyProp)*NSelect*exp(resiSurv)
    }
  }
  
  cat("\n Finished biology \n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
  
  #-------------------------------------------------------------------------------
  # Assessment
  #-------------------------------------------------------------------------------
  
  # prepare objects to enter stock assessment
  C1 <- dim(stocks[[1]]@catch.wt)
  C1[6]<- 10
  C2 <- C1
  C2[1]<- 1
  
  dmns        <- dimnames(stocks[[1]]@harvest)
  dmns$iter   <- 1:nits
  dmns2        <- dimnames(stocks[[1]]@catch)
  dmns2$iter   <- 1:nits
  
  stkAssessement <- FLStock(catch.wt=FLQuant(dim=C1,dimnames=dmns),
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
  
  stkAssessement@range <- stocks[[1]]@range
  
  # little hack to go from list of FLStock objects to using the $iter. This will need fixing later
  for(idxIter in 1:nits){
    iter(stkAssessement, idxIter) <- stocks[[idxIter]]
  }
  
  stkAssessement <- window(  stkAssessement,
                             start=an(fullPeriod[1]),
                             end=an(DtY))
  
  stkAssessement.tun <- surveys
  
  for(idxSurvey in 1:length(surveys)){
    minYearSurvey     <- min(as.numeric(colnames(stkAssessement.tun[[surveyNames[idxSurvey]]]@index)))
    stkAssessement.tun[[idxSurvey]] <- window( stkAssessement.tun[[idxSurvey]][,,,,,idxIter],
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
  
  stkAssessement@stock   <- computeStock(stkAssessement)
  
  cat("\n Finished stock assessment \n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
  
  #-------------------------------------------------------------------------------
  # Fishery, i.e. Short term forecast
  #-------------------------------------------------------------------------------
  
  ############# Recruitment #############
  # Recruitment for future years are computed trough a 10 year weighted average
  # recruitment is the same for Intermediate, forecast and continuation years
  
  # finding standard deviation on recruitment using the iterations (no residuals)
  recSd <- array(NA,dim=c(1,dim(stkAssessement)[2]),
                 dimnames = list('sdlogR','year' = an(fullPeriod[1]):an(DtY)))
  for(idxYear in 1:dim(stkAssessement)[2]){
    recSd[idxYear] <- sd(log(rec(stkAssessement[,idxYear])))
  }
  
  # recruitment for stf as weighted mean over 10 years
  recFuture <- exp( apply( log(rec(stkAssessement)[,ac((an(DtY)-10):DtY)]),
                           6,
                           weighted.mean,
                           w=1/recSd[,ac((an(DtY)-10):DtY)],
                           na.rm=T))
  recFuture <- as.array(drop(recFuture))
  
  ############# Intermediate year #############
  # hack to be fixed here. Going from $iter to list of FLStock objects
  for(idxIter in 1:nits){
    NSH.sim[[idxIter]] <- iter(stkAssessement,idxIter)
  }
  
  
  ############# initialize stf FLStock object #############
  dsc         <- "North Sea Herring"
  nam         <- "NSAS"
  dms         <- dimnames(stocks[[1]]@m)
  dms$year    <- ac((an(ImY)-3):an(CtY))
  dms$unit    <- c("A","B","C","D")
  dms$iter    <- 1:nits
  
  nAges       <- length(dms$age)
  nits        <- length(dms$iter)
  
  # Create the stf object 
  stf         <- FLStock(name=nam,desc=dsc,FLQuant(NA,dimnames=dms))
  
  for(idxIter in 1:nits){
    for(idxFleet in dms$unit){
      # update stf object with current stock object
      stf[,,idxFleet,,,idxIter]          <- window( stocks[[idxIter]],
                                                    start=an(dms$year)[1],
                                                    end=rev(an(dms$year))[1])
      # update stf object with stock numbers from SAM object
      stf[,,idxFleet,,,idxIter]@stock.n  <- window( NSH.sim[[idxIter]]@stock.n,
                                                    start=an(dms$year)[1],
                                                    end=rev(an(dms$year))[1])
      
      # update catch.wt
      if(idxFleet == 'B' || idxFleet == 'D'){
        stf[,,idxFleet,,,idxIter]@catch.wt <- fisheryFuture[,FuY,'BD','catch.wt',,idxIter]
      }else{
        stf[,,idxFleet,,,idxIter]@catch.wt <- fisheryFuture[,FuY,idxFleet,'catch.wt',,idxIter]
      }
    }
  }
  
  # Fill slots that have no meaning for NSAS
  stf@discards.n[]          <- 0
  stf@discards[]            <- 0
  stf@discards.wt[]         <- 0
  
  
  ############################ END TEST
  
  ############# Compute F in intermediate year #############
  
  # update intermediate year with results from assessment
  
  source(file.path(functionPath,"optF_TACdiff.R"))
  
  optF_TACdiff( c(1,1,1,1),         # scalor 4x1
                NSH.sim[[idxIter]]@harvest,   # single fleet FLStock object
                fisheryFuture[,,,'catch.wt',,idxIter],   # stock number single fleet
                NSH.sim[[idxIter]]@stock.n,  # catch weight at age single fleet
                stocks[[idxIter]]@m,            # natural mortality
                fisheryFuture[,,,'sel',,idxIter],         # selectivity stored as FLQuant object. Normalized between 0 and 1.
                ImY,          # year of interest
                TAC[,,,,,idxIter],         # TAC FLQuant object for fleets A, B and D
                FCPropIts[,idxIter],
                TAC_var,
                recFuture[idxIter])
  
  Fscalor <- array( 0, dim=c(nFleets,nits)) # initialize array
  for(idxIter in 1:nits){
    # find F for the different fleets
    Fscalor[,idxIter] <- nls.lm(  par=runif(4), # starting point
                                  lower=rep(1e-8,4),
                                  upper=NULL,
                                  optF_TACdiff,                                  # function to optimize
                                  harvest_sf  = NSH.sim[[idxIter]]@harvest,       # single fleet FLStock object
                                  catch.wt_mf = fisheryFuture[,,,'catch.wt',,idxIter],       # catch weight at age single fleet. Using stock weight at age for now. How to get catch weight for 2018?
                                  stock.n_sf  = NSH.sim[[idxIter]]@stock.n,       # stock at age
                                  M           = stocks[[idxIter]]@m,              # natural mortality
                                  Fsel        = fisheryFuture[,,,'sel',,idxIter], # selectivity stored as FLQuant object. Normalized between 0 and 1.
                                  iYr         = ImY,                              # year of interest
                                  TACs        = TAC[,,,,,idxIter],             # TAC FLQuant object for fleets A, B and D
                                  FCProp      = FCPropIts[,idxIter],
                                  TAC_var     = TAC_var,
                                  recruit     = recFuture[idxIter],
                                  nls.lm.control(ftol = (.Machine$double.eps),maxiter = 1000), # optimizer control object
                                  jac=NULL)$par
  }
  
  ############# update stf object #############
  
  # update F
  #stf@harvest[,ImY,'A'] <- t(apply(fisheryFuture[,ImY,'A','sel'],1,'*',Fscalor[1,]))
  #stf@harvest[,ImY,'B'] <- t(apply(fisheryFuture[,ImY,'BD','sel'],1,'*',Fscalor[2,]))
  #stf@harvest[,ImY,'C'] <- t(apply(fisheryFuture[,ImY,'C','sel'],1,'*',Fscalor[3,]))
  #stf@harvest[,ImY,'D'] <- t(apply(fisheryFuture[,ImY,'BD','sel'],1,'*',Fscalor[4,]))
  
  
  # compute stock.n, catch.n and landing.n
  for(idxIter in 1:nits){
    stf@harvest[,ImY,'A',,,idxIter] <- fisheryFuture[,ImY,'A','sel',,idxIter]*Fscalor[1,idxIter]
    stf@harvest[,ImY,'B',,,idxIter] <- fisheryFuture[,ImY,'BD','sel',,idxIter]*Fscalor[2,idxIter]
    stf@harvest[,ImY,'C',,,idxIter] <- fisheryFuture[,ImY,'C','sel',,idxIter]*Fscalor[3,idxIter]
    stf@harvest[,ImY,'D',,,idxIter] <- fisheryFuture[,ImY,'BD','sel',,idxIter]*Fscalor[4,idxIter]
    
    Z <-  rowSums(drop(stf@harvest[,ImY,,,,idxIter])) + # sum accross the fleets
      drop(stf[,ImY,1,,,idxIter]@m) # M is the same for all fleets in the stf object
    
    # propagate stock number with Z
    survivors                             <- stf[,ac(an(ImY)-1),1,,,idxIter]@stock.n*exp(-Z) # stock.n is the same for all fleets in the stf object
    stf@stock.n[2:nAges,ImY,,,,idxIter]   <- survivors[1:(nAges-1)]
    stf@stock.n[nAges,ImY,,,,idxIter]     <- stf[nAges,ImY,1,,,idxIter]@stock.n + survivors[nAges]
    stf@stock.n[1,ImY,,,,idxIter]         <- subset(recruit,year==an(ImY))$value[idxIter]
    
    for(idxFleet in 1:nFleets){  
      stf@catch.n[,ImY,idxFleet,,,idxIter]     <-  stf@stock.n[,ImY,idxFleet,,,idxIter]*
        (1-exp(-Z))*
        stf@harvest[,ImY,idxFleet,,,idxIter]/Z
      stf@catch[,ImY,idxFleet,,,idxIter]       <- computeCatch(stf[,ImY,idxFleet,,,idxIter])
      stf@landings.n[,ImY,idxFleet,,,idxIter]  <- stf@catch.n[,ImY,idxFleet,,,idxIter]
      stf@landings[,ImY,idxFleet,,,idxIter]    <- computeLandings(stf[,ImY,idxFleet,,,idxIter])
    }
  }
  
  
  
  
  
  # initialize stf object and fill in intermediate year
  stf <- stf_ImY( NSH.sim,
                  stocks,
                  fisheryFuture,
                  TAC,
                  TAC_var,
                  FCPropIts,
                  recFuture, # recruitment is know in 2018 for the initial assessment
                  FuY)
 
  
  # apply HCR to obtain F2-6 and F0-1
  
  cat("\n Finished forecast \n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
}
save.image(file=paste(outPath,runName,"_",settings$RecRegime,".RData",sep=""))


