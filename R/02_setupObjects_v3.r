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
library(FLFleet)

# define path to directory
#path          <- "D:/Work/Herring MSE/NSAS/"
#path              <- "D:/git/wk_WKNSMSE_her.27.3a47d/R/"
#path              <- "F:/WKNSMSE/wk_WKNSMSE_her.27.3a47d/R"
path <- 'E:/wk_WKNSMSE_her.27.3a47d/R'
#path <- '/home/berge057/ICES/wk_WKNSMSE_her.27.3a47d/R/'
assessment_name   <- "NSAS_WKNSMSE2018"
try(setwd(path),silent=TRUE)

# paths to different subfolders
dataPath      <- file.path(".","data/")
outPath       <- file.path(".","results/")
scriptPath    <- file.path(".","side_scripts/")
functionPath  <- file.path(".","functions/")

flag_saveAll <- FALSE # flag to indicate whether we want to save every step of the initialization

# loading function
source(file.path(functionPath,"randBlocks.R"))
source(file.path(functionPath,"randNums.R"))
source(file.path(functionPath,"make.arma.resid.R"))
source(file.path(functionPath,"make.arma.resid.lst.R"))


#-------------------------------------------------------------------------------
# 2) load assessment objects (single and multi fleet)
#    define MSE parameters
#    load raw M
#
# Note 1: the assessments we use is without the LAI index. The assessments that
# are ran during HAWG is using the LAI so results are slightly different. See 
# 00_test_no_LAI.R for a comparison of the assessments. This is for convenience
# as the LAI is a component index and is weekly structured, therefore 
# complicated to implement
#-------------------------------------------------------------------------------

#- Load single fleet and multi fleet assessment objects - using assessments without the LAI
#load(file.path(outPath,paste0(assessment_name,'_mf.Rdata')))
#load(file.path(outPath,paste0(assessment_name,'_sf.Rdata')))
load(file.path(outPath,paste0(assessment_name,'_mf_noLAI.Rdata')))
load(file.path(outPath,paste0(assessment_name,'_sf_noLAI.Rdata')))

# parameters
n.retro.years       <-  7                                       # Number of years for which to run the retrospective
nFutureyrs          <- 20 + 3
histMinYr           <- dims(NSH)$minyear
histMaxYr           <- dims(NSH)$maxyear
yearCurrent         <- histMinYr:histMaxYr # vector the years
futureMaxYr         <- histMaxYr + nFutureyrs
histPeriod          <- ac(histMinYr:histMaxYr)
projPeriod          <- ac((histMaxYr+1):futureMaxYr)
fullPeriod          <- c(histPeriod,projPeriod)
recrPeriod          <- ac(2007:2017)
selPeriod           <- ac(2007:2017)
fecYears            <- ac(2007:2017)
nits                <- 1000 # number of random samples

# reading the raw M and applying plus group
#raw_M             <- read.csv(file.path(dataPath,"Smoothed_span50_M_NotExtrapolated_NSASSMS2016.csv"),header=TRUE)
raw_M             <- read.csv(file.path(dataPath,"Raw_NotExtrapolated_NSAS_SMS2016.csv"),header=TRUE)
colnames(raw_M)   <- sub("X","",colnames(raw_M))
rownames(raw_M)   <- raw_M[,1]
raw_M             <- raw_M[,-1]# Trim off first column as it contains 'ages'
raw_M             <- cbind(replicate(as.numeric(colnames(raw_M)[1])-histMinYr,raw_M[,1]), raw_M)
raw_M             <- cbind(raw_M,raw_M[,dim(raw_M)[2]])
colnames(raw_M)   <- histMinYr:histMaxYr
raw_M             <- raw_M[1:9,] + 0.11 # trim age 9 and add 0.11 (assessment profiling from WKPELA)

#-------------------------------------------------------------------------------
# 3)  create random samples using variance/covariance matrix
#     initialize biol object
# 
# Note: this takes forever, therefore saving the object after this process
#-------------------------------------------------------------------------------

NSH.sim         <- simulate(NSH,NSH.tun,NSH.ctrl,n=nits)
names(NSH.sim)  <- paste0('iter',1:nits)

C1              <- dim(NSH@catch.wt)
C1[6]           <- nits
C2              <- C1
C2[1]           <- 1

dmns          <- dimnames(NSH@catch.wt)
dmns$iter     <- 1:nits
dmns2         <- dimnames(NSH@catch)
dmns2$iter    <- 1:nits

biol  <- FLStock( catch.wt=FLQuant(dim=C1,dimnames=dmns),
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

units(biol) <- units(NSH)

for(idxIter in 1:nits){
  print(paste('init step 3 random samples - iter=',idxIter))
  iter(biol,idxIter) <- NSH + NSH.sim[[idxIter]]
}

biol          <- window(window(biol,end=histMaxYr+1),start=histMinYr,end=futureMaxYr) # extend the FLStock object to the full projection period
biol@m.spwn[,ac(2018:2040)] <- 0.67

if(flag_saveAll)
  save( biol,
        NSH.sim,
        file=file.path(outPath,paste0(assessment_name,'init_MSE_3.RData')))

#-------------------------------------------------------------------------------
# 4) create FLStocks object using random samples (with future years as NA)
#-------------------------------------------------------------------------------

biol@harvest.spwn[,projPeriod]  <- biol@harvest.spwn[,ac(histMaxYr)] # propagate Fprop before spawning
biol@m.spwn[,projPeriod]        <- biol@m.spwn[,ac(histMaxYr)] # propagate Fprop before spawning
biol@stock                      <- computeStock(biol)

#-------------------------------------------------------------------------------
# 5) allocating future maturity, stock weight and M at age
#
# w@a and mat are on the same randomization
# M gets an independent randomization
# we randomize the following:
# number of years in the chain
# start year in the chain
# reversing of the chain
#
# Note 1: M is from raw_M, not the M (smoothed) from the assessment
# Note 2: this part takes a while
#-------------------------------------------------------------------------------

# define fishery object
dmns      <- dimnames(biol@m)
dmns$area <- c("A","B","C","D")

fishery   <- FLCatch(price=FLQuant(NA,dimnames=dmns))

# generate random blocks for weight at age and maturity
yrChain   <- randBlocks(an(fecYears),an(projPeriod),nits)
yrChainM  <- randBlocks(an(fecYears),an(projPeriod),nits)

for(idxIter in 1:nits){
  print(paste('init step 5 bio variables - iter=',idxIter))
  
  # future maturity at age
  biol@mat[,projPeriod,,,,idxIter]          <- array( iter(biol@mat[,ac(yrChain[[idxIter]])],idxIter),
                                                      dim=dim(iter(biol@mat[,projPeriod],idxIter)))
  
  # future stock weight at age
  biol@stock.wt[,projPeriod,,,,idxIter]     <- array( iter(biol@stock.wt[,ac(yrChain[[idxIter]])],idxIter),
                                                      dim=dim(iter(biol@stock.wt[,projPeriod],idxIter)))
  
  # future catch weight at age
  biol@catch.wt[,projPeriod,,,,idxIter]     <- array( iter(biol@catch.wt[,ac(yrChain[[idxIter]])],idxIter),
                                                      dim=dim(iter(biol@catch.wt[,projPeriod],idxIter)))
  
  # future natural mortality at age, on a different block chain
  biol@m[,projPeriod,,,,idxIter]            <- as.matrix( raw_M[,ac(yrChainM[[idxIter]])])
  
  # # multi fleet landing weight at age
  # A fleet
  fishery@landings.wt[,colnames(NSHs3$residual@catch.wt[,,,,'A']),
                      ,,'A',idxIter]                      <- NSHs3$residual@catch.wt[,,,,'A']
  fishery@landings.wt[,projPeriod,,,'A',idxIter]          <- NSHs3$residual@catch.wt[,ac(yrChain[[idxIter]]),,,'A']
  
  # B fleet
  fishery@landings.wt[,colnames(NSHs3$residual@catch.wt[,,,,'BD']),
                      ,,'B',idxIter]                      <- NSHs3$residual@catch.wt[,,,,'BD']
  fishery@landings.wt[,projPeriod,,,'B',idxIter]          <- NSHs3$residual@catch.wt[,ac(yrChain[[idxIter]]),,,'BD']
  
  # C fleet
  fishery@landings.wt[,colnames(NSHs3$residual@catch.wt[,,,,'C']),
                      ,,'C',idxIter]                      <- NSHs3$residual@catch.wt[,,,,'C']
  fishery@landings.wt[,projPeriod,,,'C',idxIter]          <- NSHs3$residual@catch.wt[,ac(yrChain[[idxIter]]),,,'C']
  
  # D fleet
  fishery@landings.wt[,colnames(NSHs3$residual@catch.wt[,,,,'BD']),
                      ,,'D',idxIter]                      <- NSHs3$residual@catch.wt[,,,,'BD']
  fishery@landings.wt[,projPeriod,,,'D',idxIter]          <- NSHs3$residual@catch.wt[,ac(yrChain[[idxIter]]),,,'BD']
  
  # loop to delete zero weights. One uses the mean over the projected years to fill in the gaps
  for(idxFleet in 1:dim(fishery@landings.wt)[3]){
    for(idxAges in 1:dim(fishery@landings.wt)[1]){
      # find indices where weights are 0
      idxZeros <- which(drop(fishery@landings.wt[idxAges,
                                                 projPeriod,
                                                 ,,idxFleet,idxIter])==0,arr.ind = T)
      # find indices where weights are not 0
      idxNonZeros <- which(drop(fishery@landings.wt[idxAges,
                                                    projPeriod,
                                                    ,,idxFleet,idxIter])!=0,arr.ind = T)
      
      # put mean to the years where catch weight is zero
      fishery@landings.wt[idxAges,
                          projPeriod[idxZeros], # subset years that are zero
                          ,,idxFleet,idxIter] <- mean(drop(fishery@landings.wt[ idxAges,
                                                                                projPeriod[idxNonZeros], # subset years that are non zero
                                                                                ,,idxFleet,idxIter]))
      
    }
  }
}

# filling in catch.wt in previous years
yearsMulti <- colnames(NSHs3$residual@catch.wt)

fishery@landings.wt[,yearsMulti,,,'A']           <- NSHs3$residual@catch.wt[,yearsMulti,,,'A']
fishery@landings.wt[,yearsMulti,,,'B']           <- NSHs3$residual@catch.wt[,yearsMulti,,,'BD']
fishery@landings.wt[,yearsMulti,,,'C']           <- NSHs3$residual@catch.wt[,yearsMulti,,,'C']
fishery@landings.wt[,yearsMulti,,,'D']           <- NSHs3$residual@catch.wt[,yearsMulti,,,'BD']

fishery@landings.sel[,yearsMulti,,,'A']          <- NSH3f.sam@harvest[,yearsMulti,,,'A']
fishery@landings.sel[,yearsMulti,,,'B']          <- NSH3f.sam@harvest[,yearsMulti,,,'BD']
fishery@landings.sel[,yearsMulti,,,'C']          <- NSH3f.sam@harvest[,yearsMulti,,,'C']
fishery@landings.sel[,yearsMulti,,,'D']          <- NSH3f.sam@harvest[,yearsMulti,,,'BD']

# landing.wt = catch.wt
biol@landings.wt                                 <- biol@catch.wt

if(flag_saveAll)
  save( biol,
        NSH.sim,
        fishery,
        file=file.path(outPath,paste0(assessment_name,'init_MSE_5.RData')))

#-------------------------------------------------------------------------------
# 6) creating survey indices
# 
# survey indices are created for each random sample using Q, N, F and M at age
# residuals are added using a normal distribution with sigma as sd of observations
# I(a,y) = Q(a)*exp(-Z(A,y)*surveyProp)*N(a,y)*res(a,y)
# Z = M + F
#
# The following survey indices are to be generatred
# HERAS age 1 to 8
# IBTS-Q1 age 1
# IBTS0 age 0
# IBTS-Q3 age 0 to 5
#
# Note 1: the raw M is used here, as opposed to the smoothed one used in the 
# assessment
# Note 2: we don't use the LAI index here, see assessment part
# Note 3: filling up matrices of residuals and catchabilities takes a while
#-------------------------------------------------------------------------------

#load(file.path(outPath,paste0(assessment_name,'init_MSE_5.RData')))

surveys     <- lapply(NSH.tun,propagate,iter=nits)
surveys     <- window(window(surveys,end=histMaxYr+1),start=histMinYr,end=futureMaxYr) # extend the FLStock object to the full projection period


############# initialize FLQuant object containing catchabilities and residuals #############
dmns        <- dimnames(NSH@harvest)
dmns$unit   <- names(NSH.tun)
dmns$season <- c('catchabilities','residuals')
dmns$year   <- fullPeriod
dmns$iter   <- 1:nits

surveyVars    <- FLQuant(array( NA, # covariance matrix using a period of 10 years for all the ages
                                dim=c(length(dmns$age), # ages
                                      length(fullPeriod),   # years
                                      length(NSH.tun),            # number of surveys
                                      2,            # quantity stored (catchability and residuals)
                                      1,
                                      nits)), # iterations
                          dimnames=dmns)

# get catchabilities and residuals for current and future years
for(idxIter in 1:nits){
  print(paste('init step 6 survey indices residuals & catchabilities - iter=',idxIter))
  
  sdAll <- obs.var(NSH.sim[[idxIter]]) # getting observation variance from SAM object for the current iteration
  qAll  <- catchabilities(NSH.sim[[idxIter]]) # getting catchabilities from SAM object for the current iteration
  
  surveyNames <- as.character(unique(qAll$fleet)) # get all the survey names
  
  # creating indices for surveys
  # loop on all available surveys
  for(surveyCurrent in surveyNames){
    sdSelect    <-  subset(sdAll, sdAll$fleet == surveyCurrent) # subset observation variance
    qSelect     <- subset(qAll,qAll$fleet == surveyCurrent)
    
    # building residuals resi(a,y)
    # initialize residual array
    maxYearSurvey     <- futureMaxYr
    minYearSurvey     <- min(as.numeric(colnames(NSH.tun[[surveyCurrent]]@index)))
    yearCurrentSurvey <- minYearSurvey:maxYearSurvey # vector the years in the survey
    resi              <- array(0, dim=c(dim(sdSelect)[1],length(yearCurrentSurvey))) # initialize array nAges x nYears
    colnames(resi)    <- yearCurrentSurvey
    rownames(resi)    <- qSelect$age
    # generate the residuals using a normal distribution - are residuals log or linear??? Probably log
    for(idxResi in 1:dim(sdSelect)[1]){
      resi[idxResi,]  <- rlnorm(length(yearCurrentSurvey), 
                                0, 
                                sdSelect$value[idxResi])
    }
    
  
    # fill in object
    surveyVars[ac(qSelect[,2]),,surveyCurrent,'catchabilities',,idxIter]  <- replicate(length(fullPeriod),qSelect[,3])
    surveyVars[ac(qSelect[,2]),ac(yearCurrentSurvey),surveyCurrent,'residuals',,idxIter]       <- resi
  }
}

surveyNames <- names(NSH.tun)

# calculate new survey indices for each survey
for(surveyCurrent in surveyNames){
  maxYearSurvey     <- max(as.numeric(colnames(NSH.tun[[surveyCurrent]]@index)))
  if(maxYearSurvey > histMaxYr) 
    maxYearSurvey <- histMaxYr
  minYearSurvey     <- min(as.numeric(colnames(NSH.tun[[surveyCurrent]]@index)))
  yearCurrentSurvey <- minYearSurvey:maxYearSurvey # vector the years in the survey
  ageSurvey <- NSH.tun[[surveyCurrent]]@range[1]:NSH.tun[[surveyCurrent]]@range[2]
  
  
  NSelect     <- biol@stock.n[ac(ageSurvey),ac(yearCurrentSurvey)]  # filter years
  NSelect     <- drop(NSelect) # drop dimensions with 1 level
  FSelect     <- biol@harvest[ac(ageSurvey),ac(yearCurrentSurvey)] # filter years for F
  FSelect     <- drop(FSelect) # drop dimensions with 1 level
  Z           <- raw_M[ac(ageSurvey),ac(yearCurrentSurvey)]  # filter years for M
  Z           <- drop(replicate(nits,as.matrix(Z))) + FSelect
  surveyProp  <- mean(c(NSH.tun[[surveyCurrent]]@range[6],NSH.tun[[surveyCurrent]]@range[7]))
  
  
  # calculate survey indices
  surveys[[surveyCurrent]]@index[,ac(yearCurrentSurvey)] <- drop(surveyVars[ac(ageSurvey),ac(yearCurrentSurvey),surveyCurrent,'catchabilities',,])*
                                                            exp(-Z*surveyProp)*NSelect*
                                                            drop(surveyVars[ac(ageSurvey),ac(yearCurrentSurvey),surveyCurrent,'residuals',,])
  if(surveyCurrent == 'HERAS') 
    surveys[[surveyCurrent]]@index[1,ac(1989:1996)] <- -1 # trim age 1 from 1989 to 1996 specifically for HERAS
  
}

if(flag_saveAll)
  save( biol,
        NSH.sim,
        fishery,
        surveys,
        surveyVars,
        file=file.path(outPath,paste0(assessment_name,'init_MSE_6.RData')))

# plotting survey indices
#surveyName <- 'IBTS-Q1'

#par(mfrow=c(4,2))
#for(ageIdx in surveys[[surveyName]]@range[1]:surveys[[surveyName]]@range[2]){
  #a     <-drop(NSH.tun[[surveyName]]@index)
#a     <- NSH.tun[[surveyName]]@index # in case of only 1 age in survey object
#  years <- as.numeric(colnames(NSH.tun[[surveyName]]@index))
#  plot(years,a[ac(ageIdx),],type='l',ylab=c('age',ac(ageIdx)))
#  
#  for(idxIter in 1:100){
#    b     <-surveys[[surveyName]]@index[ac(ageIdx),,,,,idxIter]
#    b     <- b[,match(years,colnames(b))]
#    lines(years,b,col='green')
#  }
#}


#-------------------------------------------------------------------------------
# 7) creating catches from each random samples
# 
# survey indices are created for each random sample using F and N and M at age
# residuals are added using a normal distribution with sigma as sd of observations
# C(a,y) = F(a,y)/Z(a,y)*(1-exp(-Z(a,y)))*N(a,y)*res(a,y)
#
# Note 1: the raw M is used here, as opposed to the smoothed one used in the 
# assessment
#-------------------------------------------------------------------------------

#load(file.path(outPath,paste0(assessment_name,'init_MSE_6.RData')))

############# initialize FLQuant object containing catch residuals #############
dmns        <- dimnames(NSH@harvest)
dmns$unit   <- 'catch unique'
dmns$season <- c('residuals','FCprop') # residuals and proportion of F
dmns$year   <- fullPeriod
dmns$iter   <- 1:nits

catchVar    <- FLQuant(array( NA, # covariance matrix using a period of 10 years for all the ages
                              dim=c(length(dmns$age), # ages
                                    length(fullPeriod),   # years
                                    1,            # fleet (4)
                                    2,            # quantity stored (residuals, Fcprop)
                                    1,
                                    nits)), # iterations
                         dimnames=dmns)

# fill in residuals for the catches for current and future years
for(idxIter in 1:nits){
  print(paste('init step 7 catch residuals - iter=',idxIter))
  
  sdAll     <- obs.var(NSH.sim[[idxIter]]) # getting observation variance from SAM object for the current iteration
  sdSelect  <- subset(sdAll, sdAll$fleet == 'catch unique') # subset observation variance
  
  # fill the residual array
  resi              <- array(0, dim=c(dim(sdSelect)[1],length(fullPeriod))) # initialize array nAges x nYears
  colnames(resi)    <- fullPeriod
  rownames(resi)    <- sdSelect$age
  for(idxResi in 1:dim(sdSelect)[1]){
    resi[idxResi,]  <- rlnorm(length(fullPeriod), 
                              0, 
                              sdSelect$value[idxResi]) # lognormal distribution
  }
  
  catchVar[,,,'residuals',,idxIter]       <- resi
}

# compute catches
#C(a,y) = F(a,y)/Z(a,y)*(1-exp(-Z(a,y)))*N(a,y)*res(a,y)

NSelect     <- biol@stock.n[,histPeriod]  # filter years
NSelect     <- drop(NSelect) # drop dimensions with 1 level
FSelect     <- biol@harvest[,histPeriod] # filter years for F
FSelect     <- drop(FSelect) # drop dimensions with 1 level
Z           <- raw_M[,histPeriod]  # filter years for M
Z           <- drop(replicate(nits,as.matrix(Z))) + FSelect

# catches, with added error based on observation variance
biol@catch.n[,histPeriod] <-  FSelect/Z*(1-exp(-Z))*NSelect*drop(catchVar[,histPeriod,,'residuals'])
biol@catch                <- computeCatch(biol)
#We don't believe the closure catch data, so put it to NA
biol@catch.n[,ac(1978:1979)]           <- NA

# landings, here modelled as catches without error
biol@landings.n[,histPeriod] <-  FSelect/Z*(1-exp(-Z))*NSelect
biol@landings <- computeLandings(biol)

if(flag_saveAll)
  save( biol,
        catchVar,
        NSH.sim,
        fishery,
        surveys,
        surveyVars,
        file=file.path(outPath,paste0(assessment_name,'init_MSE_7.RData')))

# plotting the catches
#par(mfrow=c(3,3))
#for(ageIdx in 1:dim(NSH@catch.n)[1]){
#  a     <- drop(NSH@catch.n)
#  years <- histMinYr:histMaxYr # vector the years
#  plot(years,a[ageIdx,],type='l',ylab=c('age',ac(ageIdx-1)))
#  for(idxIter in 1:100){
#    b     <- iter(biol@catch.n[ageIdx],idxIter)
#    b     <- b[,match(years,colnames(b))]
#    lines(years,b,col='green')
#  }
#}

#-------------------------------------------------------------------------------
# 8) C fleet: proportion of F of the C fleet for the future years
# The proportion (all ages combined) is obtained from the multi-fleet assessment 
# using a random draw similar to M and weight at age
# For the given quantity, we randomize the following:
# number of years in the chain
# start year in the chain
# reversing of the chain
#
# Note 1: multi fleet assessment is performed as a multi-fleet from 1997 onward 
# and as a single fleet assessment prior to 1997, therefore an FLStock object 
# with 2 fields
# Note 2: NSAS/WBSS split for the C fleet will be kept constant to 30 NSAS/70 WBSS
# Note 3: we don't have random samples here as we don't simulate different time
# series for the multi-fleet assessment
# Note 4: for now, we add F for all the ages but I don't this is correct,
# depending on whether F accross the ages is additive. Should not this be based 
# on the catch.
#-------------------------------------------------------------------------------

#load(file.path(outPath,paste0(assessment_name,'init_MSE_7.RData')))

yrChainFC   <- randBlocks(an(fecYears),an(projPeriod),nits)
ages        <- NSH3f.sam@range[1]:NSH3f.sam@range[2]

# initializing object
FA            <- array(0, dim=c(length(ages),length(yearCurrent))) # initialize array nAges x nYears
colnames(FA)  <- yearCurrent
rownames(FA)  <- ages

FBD            <- array(0, dim=c(length(ages),length(yearCurrent))) # initialize array nAges x nYears
colnames(FBD)  <- yearCurrent
rownames(FBD)  <- ages

FC            <- array(0, dim=c(length(ages),length(yearCurrent))) # initialize array nAges x nYears
colnames(FC)  <- yearCurrent
rownames(FC)  <- ages

for(idxAge in 1:length(ages)){
  Ftot              <-  NSH3f.sam@harvest[idxAge,1:length(yearCurrent),,,1] + # fleet A
    NSH3f.sam@harvest[idxAge,1:length(yearCurrent),,,2] + # fleet BD
    NSH3f.sam@harvest[idxAge,1:length(yearCurrent),,,3]   # fleet C
  Ftot              <-  drop(Ftot)
  FA[idxAge,]   <-  drop(NSH3f.sam@harvest[idxAge,1:length(yearCurrent),,,1])
  FBD[idxAge,]  <-  drop(NSH3f.sam@harvest[idxAge,1:length(yearCurrent),,,2])
  FC[idxAge,]   <-  drop(NSH3f.sam@harvest[idxAge,1:length(yearCurrent),,,3])
}

FCProp            <- colSums(FC)/(colSums(FA)+colSums(FBD)+colSums(FC))

# fill in object for proportion of F in the NS from the C fleet
# fill in historical period, obviously no randominzation in these.
catchVar[1,histPeriod,,'FCprop'] <- FCProp

# fill in proportion for future years, only using age 0 in the FLQuant object
for(idxIter in 1:nits){
  print(paste('init step 8 catch FcProp - iter=',idxIter))
  # fill in projection period with randomization (see above comments)
  catchVar[1,projPeriod,,'FCprop',,idxIter] <-   as.matrix(FCProp[ac(yrChainFC[[idxIter]])])
}

#-------------------------------------------------------------------------------
# 9) create selection patterns for the different fleets
# Selectivity of fleet projected forward using a random walk, using results from
# the multi-fleet assessment
#-------------------------------------------------------------------------------

fleets <- c('A','C','BD')


############# compute selectivity with random walk for each fleet #############
dmns        <- dimnames(NSH@harvest)
dmns$year   <- projPeriod
dmns$iter   <- 1:nits

# loop through the fleets to compute selectivity in futyre years for each fleet
for(idxFleet in 1:length(fleets)){

  currentHarvest <- NSH3f.sam@harvest[,ac((max(yearCurrent)-10):max(yearCurrent)),,,fleets[idxFleet]] # get F for the current fleet
  #currentHarvest <- NSH@harvest[,ac((max(yearCurrent)-10):max(yearCurrent))] # get F for the current fleet

  #- Create random walk over Fs (as multiplier from last years selection pattern)
  # covariance in F (log) over 10 years for age 0 to 2
  covmat1                   <- cov(apply(log(drop(currentHarvest)), # covariance for the last 10 years
                                         1,
                                         diff))
  covmat10Y                 <- covmat1 # covariance computed for all the ages over the 10 years
  covmat1[ac(0:2),ac(3:8)]  <- 0
  covmat1[ac(3:8),ac(0:2)]  <- covmat1[ac(0:2),ac(3:8)]
  covmat1[ac(3:8),ac(3:8)]  <- covmat1[ac(3:8),ac(0:2)]
  # covariance over 20 years for age 3 and up
  covmat2                   <- cov(apply(log(drop(currentHarvest)), # covariance over the last 20 years
                                         1,
                                         diff))
  covmat2[ac(0:2),ac(0:2)]  <- 0;
  covmatMix                 <- covmat1 + covmat2 # mix of year period for different ages
  covmatMix[is.na(covmatMix)] <- 0

  # create random samples using a multivariate normal distribution for the log covariance between ages for the future years for each iteration.
  # using 10 year period for the log covariance for age 0 to 2 and 20 year period for ages 3+
  wF                                                            <- FLQuant(array(t(mvrnorm(nits*nFutureyrs,
                                                                                           rep(0,length(ages)),
                                                                                           covmatMix)), # covariance matrix using 10 years for age 0 to 2 and 20 years for ages 3+
                                                                                 dim=c(length(ages),nFutureyrs,1,1,1,nits)),
                                                                           dimnames=dmns)

  # handle outliers for the random log covariances samples
  qtil                                                          <- quantile(c(wF),probs=c(0.05,0.95)) # 5/95 percentiles for all ages.
  qtilold                                                       <- quantile(c(wF[ac(4:8),]),probs=c(0.25,0.75)) # 25/75 percentiles for ages 4 to 8
  # set outliers to 0 based on 5/95 percentiles across all ages
  wF@.Data[which(wF<qtil[1] | wF>qtil[2])][]                    <- 0
  # set outliers to 0 based on 5/95 percentiles across all ages
  wF[ac(4:8),]@.Data[which(wF[ac(4:8),]<qtilold[1] | wF[ac(4:8),]>qtilold[2])][] <- 0

  # mimicing random walk through cumsum of the variances (i.e. F residuals) through the year for each age
  Ftemp <- apply(log(drop(currentHarvest)),1,mean) # mean of F at age over the selected number of years selPeriod
  rwF   <- wF
  for(idxIter in 1:nits){
    print(paste('init step 9 sel pattern - iter=',idxIter,' - fleet=', fleets[idxFleet]))

    for(idxAge in 1:length(ages)){
      # compute F at age with residuals estimated through a random walk (cumsum across the years)
      rwF[idxAge,,,,,idxIter] <-  Ftemp[idxAge] + cumsum(drop(wF[idxAge,,,,,idxIter]))
    }
    # define ages for Fbar, 2-6 for the A fleet, 0-1 for B C and D fleets
    if(fleets[idxFleet]=='A') Fbarages <- ac(2:6) else Fbarages <- ac(0:1)
    # compute Fbar for each year in the current iteration
    Fbar <- apply(exp(rwF[Fbarages,,,,,idxIter]),2,mean)
    # compute selectivity as S(a,y) = F(a,y)/Fbar
    for(idxYear in 1:dim(rwF)[2]){
      rwF[,idxYear,,,,idxIter] <- exp(rwF[,idxYear,,,,idxIter])/drop(Fbar[,idxYear])
      rwF[,idxYear,,,,idxIter] <- rwF[,idxYear,,,,idxIter]/max(rwF[,idxYear,,,,idxIter])
    }
    # fill in final FLQuant object
    #fisheryFuture[,,fleets[idxFleet],'sel',,idxIter] <- rwF[,,,,,idxIter]
    if(fleets[idxFleet] == 'BD'){
      fishery@landings.sel[,projPeriod,,,'B',idxIter] <- rwF[,,,,,idxIter]
      fishery@landings.sel[,projPeriod,,,'D',idxIter] <- rwF[,,,,,idxIter]
    }else{
      fishery@landings.sel[,projPeriod,,,fleets[idxFleet],idxIter] <- rwF[,,,,,idxIter]
    }
  }
}

#-------------------------------------------------------------------------------
# 10) Future recruitment
#-------------------------------------------------------------------------------

recPeriod <- ac(2002:2016)

biol.sr <- fmle(as.FLSR(biol,model='segreg')) # just to populate the structure

itersSR <- sample(1:dims(biol)$iter,round(0.15*dims(biol)$iter),replace=F)
itersRI <- which(!(1:dims(biol)$iter) %in% itersSR)

for (its in itersSR){
  iter(params(biol.sr),its)  <- params(fmle(FLSR( rec = rec(iter(biol,its))[,ac(an(recPeriod)[2]:max(an(recPeriod)))], # rec from 1948 to 2017
                                                  ssb = ssb(iter(biol,its))[,ac(an(recPeriod)[1]:(max(an(recPeriod))-1))], # ssb from 1947 to 2016
                                                  model='segreg')))
}  
for (its in itersRI){
  iter(params(biol.sr),its)  <- params(fmle(FLSR(rec = rec(iter(biol,its))[,ac(an(recPeriod)[2]:max(an(recPeriod)))], # rec from 1948 to 2017
                                                  ssb = ssb(iter(biol,its))[,ac(an(recPeriod)[1]:(max(an(recPeriod))-1))], # ssb from 1947 to 2016
                                                  model='ricker')))
}

# THIS IS MODIFIED SO THAT AN ARIMA MODEL IS FITTED FOR THE RESIDUALS OF EACH ITERATION
# AND USE TO PRODUCE THE FUTURE DEVIATIONS FOR THE CORRESPONDING ITERATION
# I want to add something that takes autocorrelation in SR relationship into account
# to do this I use an arima model
### S/R residuals - with autocorrelation
rec.res <- residuals(biol.sr)[,ac(an(recPeriod)[2]:(max(an(recPeriod)))-1)]

# autoregressive model order 1
set.seed(108)

# a list with one model per iteration

arima.fit.lst <- list()
for(its in 1:dims(biol)$iter)
  arima.fit.lst[[its]] <- try(arima(an(iter(rec.res,its)),order=c(1,0,0)))
idx <- which(unlist(lapply(arima.fit.lst,function(x){class(x)=="try-error"}))==T)
for(its in idx)
  arima.fit.lst[[its]] <- try(arima(an(iter(rec.res,its))))

#ny <- 20        # number of years to project - Usually 20
#dy <- range(stkMC)["maxyear"]       # data year
#ay <- dy                            # assessment year
#iy <- ay+1                          # initial projections year (also intermediate)
#fy <- iy + ny -1                    # final year

# create autocorrelation in residuals and propagate throughout stock into the future
# from initial year of projections (iy) to last of projections (ny-1)
sr.res  <- make.arma.resid.lst(arima.fit.lst, age = 0, years = an(projPeriod[1]):max(an(projPeriod)), rec.res)

#-------------------------------------------------------------------------------
# 11) process error
#-------------------------------------------------------------------------------

projYearsCohort <- (an(projPeriod)[1]-8):(max(an(projPeriod)))

############# initialize FLQuant object containing catch residuals #############
dmns        <- dimnames(NSH@harvest)
dmns$age    <- ac(1:8)
dmns$season <- c('procError') # residuals and proportion of F
dmns$year   <- ac(projYearsCohort)
dmns$iter   <- 1:nits

varProccError   <- FLQuant( array(  NA, # covariance matrix using a period of 10 years for all the ages
                                    dim=c(length(dmns$age), # ages
                                    length(dmns$year),   # years
                                    1,            # fleet (4)
                                    1,            # quantity stored (residuals, Fcprop)
                                    1,
                                    nits)), # iterations
                            dimnames=dmns)

# commpute survivors
surv <- biol@stock.n[,histPeriod]*exp(-biol@harvest[,histPeriod]-biol@m[,histPeriod]) # effectively, this is age 1 to 8 in year + 1
surv[dim(surv)[1]-1] <- surv[dim(surv)[1]-1] + surv[dim(surv)[1]]
dimnames(surv)$age <- ac(1:9)

# process error
procError <-  surv[ac(1:8),histPeriod[1:(length(histPeriod)-1)]]/ # survivors age 1 to 8 (0 to 7 in surv object) year 1948 to 2017
              biol@stock.n[ac(1:8),histPeriod[2:length(histPeriod)]] # numbers at age, age 1 to 8


for(idxIter in 1:nits){
  print(paste('init step 11 process error - iter=',idxIter))
  
  # covariance accross the ages using a 10 year period of full cohorts
  covMat  <- cov(t(FLCohort(log(procError))[,ac(1999:2008),,,,idxIter,drop=T]))
  # draw covariates accross the ages for each cohort
  res     <- exp(mvrnorm(length(projYearsCohort),rep(0,dim(covMat)[1]),covMat))
  
  varProccError[,,,,,idxIter] <- t(res)
}


#-------------------------------------------------------------------------------
# 12) tidying up and saving objects for next step
#-------------------------------------------------------------------------------

# prepare stock object
stkAssessement <- biol[,c(histPeriod)]
units(biol) <- units(NSH)
units(stkAssessement) <- units(NSH)

# compute recruitment for 2018 from SAM object (containing information from IBTS0 and Q1 in 2018)
recFuture <- array(NA,dim=c(nits,1))
for(idxIter in 1:nits){
  # recruitment out of SAM for 2018
  recFuture[idxIter] <- as.array(subset(rec(NSH.sim[[idxIter]]),year==2018)$value)
}


# saving object to workspace for future use in MSE
save( biol,            # biology object
      stkAssessement,
      varProccError,   # future process error (across cohorts)
      catchVar,        # future observation variance for the catches
      recFuture,       # recruitment for 2018 for stf from SAM object
      sr.res,          # residuals for stock-recruitment
      biol.sr,         # stock-recruitment fits
      itersSR,         # iterations for segmented regression
      itersRI,         # iterations for Ricket 
      fishery,         # fishery object (fleet wise), contains selection patterns + catch.wt
      surveys,         # survey object
      surveyVars,      # future catchabilities and residuals for the surveys
      NSH.ctrl,        # SAM control object
      file=file.path(outPath,paste0(assessment_name,'_init_MSE_',ac(nits),'.RData')))

# resetting parameters
nFutureyrs          <- 20
yearCurrent         <- histMinYr:histMaxYr # vector the years
futureMaxYr         <- histMaxYr + nFutureyrs
projPeriod          <- ac((histMaxYr+1):futureMaxYr)
fullPeriod          <- c(histPeriod,projPeriod)

# save parameters
save(n.retro.years, 
     nFutureyrs,
     histMinYr,
     histMaxYr,
     yearCurrent,
     futureMaxYr, 
     histPeriod,
     projPeriod,
     fullPeriod,
     recrPeriod,
     selPeriod,
     fecYears,
     nits,
     file=file.path(outPath,paste0(assessment_name,'_parameters_MSE_',ac(nits),'.RData')))
