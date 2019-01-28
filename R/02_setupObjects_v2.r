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

# loading function
source(file.path(functionPath,"randBlocks.R"))
source(file.path(functionPath,"randNums.R"))


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
nits                <- 10 # number of random samples

# reading the raw M and applying plus group
raw_M             <- read.csv(file.path(dataPath,"Smoothed_span50_M_NotExtrapolated_NSASSMS2016.csv"),header=TRUE)
colnames(raw_M)   <- sub("X","",colnames(raw_M))
rownames(raw_M)   <- raw_M[,1]
raw_M             <- raw_M[,-1]# Trim off first column as it contains 'ages'
raw_M             <- cbind(replicate(as.numeric(colnames(raw_M)[1])-histMinYr,raw_M[,1]), raw_M)
raw_M             <- cbind(raw_M,raw_M[,dim(raw_M)[2]])
colnames(raw_M)   <- histMinYr:histMaxYr
# hack to set plus group, converting M array into an FLStock object, using the setPlusGroup, then back to array
# !!!!!! to be updated. Right now one uses an empty FLStock object. This is wrong as I think the setting of the plus group 
# needs the catches as input
NSHM2             <- readFLStock(file.path(dataPath,"index.txt"),no.discards=TRUE,quiet=FALSE)
NSHM2@m[]         <- as.matrix(raw_M)
pg                <- NSH@range['max']
NSHM2             <- setPlusGroup(NSHM2,pg) # really wonder if the setPlusGroup does anything... Needs clarifying.
raw_M             <- drop(NSHM2@m)
raw_M             <- raw_M + 0.11

#-------------------------------------------------------------------------------
# 3) create random samples using variance/covariance matrix
#-------------------------------------------------------------------------------

NSH.sim         <- simulate(NSH,NSH.tun,NSH.ctrl,n=nits)
names(NSH.sim)  <- paste0('iter',1:nits)

#-------------------------------------------------------------------------------
# 4) create FLStocks object using random samples (with future years as NA)
#-------------------------------------------------------------------------------

stocks          <- NSH + NSH.sim # from FLSAMs to FLStocks
stocks          <- window(window(stocks,end=histMaxYr+1),start=histMinYr,end=futureMaxYr) # extend the FLStock object to the full projection period

# update FLStocks object with random samples infered from variance/co-variance matrix
for(idxIter in 1:nits){
  stocks[[idxIter]]@catch.n                    <- stocks[[idxIter]]@stock.n * stocks[[idxIter]]@harvest / 
                                                  (stocks[[idxIter]]@harvest + stocks[[idxIter]]@m) * 
                                                  (1 - exp(-stocks[[idxIter]]@harvest - stocks[[idxIter]]@m)) # compute catch numbers 
  stocks[[idxIter]]@catch.n[,ac(1978:1979)]    <- NA # fill in Na for the the closure catch data
  stocks[[idxIter]]@landings.n                 <- stocks[[idxIter]]@catch.n
  
  stocks[[idxIter]]@harvest.spwn[,projPeriod]  <- stocks[[idxIter]]@harvest.spwn[,ac(histMaxYr)] # propagate Fprop before spawning
  stocks[[idxIter]]@m.spwn[,projPeriod]        <- stocks[[idxIter]]@m.spwn[,ac(histMaxYr)] # propagate Mprop before spawning
  
  stocks[[idxIter]]@stock <- computeStock(stocks[[idxIter]])
}

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
# Note: M is from raw_M, not the M (smoothed) from the assessment
#-------------------------------------------------------------------------------

# generate random blocks for weight at age and maturity
yrChain   <- randBlocks(an(fecYears),an(projPeriod),nits)
yrChainM  <- randBlocks(an(fecYears),an(projPeriod),nits)
# catch weight
multiFleet_catch.wt         <- array( 0, # initialize array, nAges x nYears x nFleets x nits. 3 fleets in assessment (A,BD,C)
                                      dim=c(dim(NSHs3$residual@catch.wt)[1], # nAges
                                            length(projPeriod),              # nYears
                                            3,                               # nFleets
                                            nits),                           # nits.
                                      dimnames = list(0:(dim(NSHs3$residual@catch.wt)[1]-1),
                                                      as.numeric(projPeriod),
                                                      1:3,
                                                      1:nits))

# update FLStocks object
for(idxIter in 1:nits){
  # future maturity at age
  stocks[[idxIter]]@mat      [,projPeriod][]               <- array(stocks[[idxIter]]@mat[,ac(yrChain[[idxIter]])],
                                                                    dim=dim(stocks[[idxIter]]@mat[,projPeriod]))
  # future weight at age
  stocks[[idxIter]]@stock.wt [,projPeriod][]               <- array(stocks[[idxIter]]@stock.wt[,ac(yrChain[[idxIter]])],
                                                                    dim=dim(stocks[[idxIter]]@stock.wt[,projPeriod]))
  
  # future catch weight at age
  stocks[[idxIter]]@catch.wt [,projPeriod][]               <- array(stocks[[idxIter]]@catch.wt[,ac(yrChain[[idxIter]])],
                                                                    dim=dim(stocks[[idxIter]]@catch.wt[,projPeriod]))
  
  # future catch weight at age
  stocks[[idxIter]]@landings.wt [,projPeriod][]             <- array(stocks[[idxIter]]@landings.wt[,ac(yrChain[[idxIter]])],
                                                                     dim=dim(stocks[[idxIter]]@landings.wt[,projPeriod]))
  
  # multi fleet catch weight at age
  multiFleet_catch.wt[,projPeriod,,idxIter] <- drop(NSHs3$residual@catch.wt[,ac(yrChain[[idxIter]])])[,1:length(projPeriod),]
  
  # future natural mortality at age, on a different block chain
  stocks[[idxIter]]@m        [,projPeriod][]               <- array(raw_M[,ac(yrChainM[[idxIter]])],
                                                                    dim=dim(stocks[[idxIter]]@m[,projPeriod]))
}

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
#-------------------------------------------------------------------------------

surveys     <- lapply(NSH.tun,propagate,iter=nits)
surveys     <- window(window(surveys,end=histMaxYr+1),start=histMinYr,end=futureMaxYr) # extend the FLStock object to the full projection period

# initialize array to store catchabilities
temp        <- catchabilities(NSH.sim[[idxIter]])
qMat        <- array(NA,
                     dim=c(dim(temp)[1],
                           2,
                           nits),
                     dimnames=list('fleet' = temp$fleet,
                                   'var' = c('ages','value')))
# initialize array to store observation variance for surveys
temp        <- obs.var(NSH.sim[[idxIter]]) # getting observation variance from SAM object for the current iteration
temp        <- subset(temp, temp$fleet != 'catch unique') # subset observation variance
varSurvMat  <- array(NA,
                     dim=c(dim(temp)[1],
                           2,
                           nits),
                     dimnames=list('fleet' = temp$fleet,
                                   'var' = c('ages','value')))

for(idxIter in 1:nits){
  
  sdAll <- obs.var(NSH.sim[[idxIter]]) # getting observation variance from SAM object for the current iteration
  qAll  <- catchabilities(NSH.sim[[idxIter]]) # getting catchabilities from SAM object for the current iteration
  qMat[,1,idxIter]  <- qAll$age
  qMat[,2,idxIter]  <- qAll$value
  
  varSurvMat[,1,idxIter] <- subset(sdAll, sdAll$fleet != 'catch unique')$age
  varSurvMat[,2,idxIter] <- subset(sdAll, sdAll$fleet != 'catch unique')$value
  
  surveyNames <- as.character(unique(qAll$fleet)) # get all the survey names
  
  # creating indices for surveys
  # loop on all available surveys
  for(idxSurvey in 1:length(surveyNames)){
    sdSelect    <-  subset(sdAll, sdAll$fleet == surveyNames[idxSurvey]) # subset observation variance
    qSelect     <- subset(qAll,qAll$fleet == surveyNames[idxSurvey])
    
    # building residuals r(a,y)
    # initialize residual array
    maxYearSurvey     <- max(as.numeric(colnames(NSH.tun[[surveyNames[idxSurvey]]]@index)))
    if(maxYearSurvey > histMaxYr) 
      maxYearSurvey <- histMaxYr
    minYearSurvey     <- min(as.numeric(colnames(NSH.tun[[surveyNames[idxSurvey]]]@index)))
    yearCurrentSurvey <- minYearSurvey:maxYearSurvey # vector the years in the survey
    resi              <- array(0, dim=c(dim(sdSelect)[1],length(yearCurrentSurvey))) # initialize array nAges x nYears
    colnames(resi)    <- yearCurrentSurvey
    rownames(resi)    <- qSelect$age
    # generate the residuals using a normal distribution - are residuals log or linear??? Probably log
    for(idxResi in 1:dim(sdSelect)[1]){
      resi[idxResi,]  <- rnorm(length(yearCurrentSurvey), 
                               0, 
                               sdSelect$value[idxResi])
    }
    
    NSelect     <- stocks[[idxIter]]@stock.n[match(as.character(sdSelect$age),rownames(stocks[[idxIter]]@stock.n)), # filter ages
                                             match(as.character(yearCurrentSurvey),colnames(stocks[[idxIter]]@stock.n))]  # filter years
    NSelect     <- drop(NSelect) # drop dimensions with 1 level
    FSelect     <-  stocks[[idxIter]]@harvest[match(as.character(sdSelect$age),rownames(stocks[[idxIter]]@stock.n)), # filter ages for F
                                              match(as.character(yearCurrentSurvey),colnames(stocks[[idxIter]]@stock.n))] # filter years for F
    FSelect     <- drop(FSelect) # drop dimensions with 1 level
    Z           <-  raw_M[match(as.character(sdSelect$age),rownames(raw_M)), # filter ages for M
                          match(as.character(yearCurrentSurvey),colnames(raw_M))]  # filter years for M
    Z           <- Z + FSelect
    surveyProp  <- mean(c(NSH.tun[[surveyNames[idxSurvey]]]@range[6],NSH.tun[[surveyNames[idxSurvey]]]@range[7]))
    
    # update survey object
    surveys[[surveyNames[idxSurvey]]]@index[,match(as.character(yearCurrentSurvey),colnames(surveys[[surveyNames[idxSurvey]]]@index)) # filling only the years available for the survey
                                            ,,,,idxIter] <- as.matrix(replicate(length(yearCurrentSurvey), qSelect$value)*exp(-Z*surveyProp)*NSelect*exp(resi))
    surveys[[surveyNames[idxSurvey]]] <- surveys[[surveyNames[idxSurvey]]][,ac(minYearSurvey:futureMaxYr)]
  }
}

# plotting survey indices
surveyName <- 'IBTS-Q3'

par(mfrow=c(3,2))
for(ageIdx in 1:dim(NSH.tun[[surveyName]]@index)[1]){
  #ageIdx  <- 4
  a     <-drop(NSH.tun[[surveyName]]@index)
  #a     <- NSH.tun[[surveyName]]@index # in case of only 1 age in survey object
  years <- as.numeric(colnames(NSH.tun[[surveyName]]@index))
  plot(years,a[ageIdx,],type='l',ylab=c('age',ac(ageIdx-1)))
  
  for(idxIter in 1:nits){
    b     <-surveys[[surveyName]]@index[ageIdx,,,,,idxIter]
    b     <- b[,match(years,colnames(b))]
    lines(years,b,col='green')
  }
}


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

# initialize array to store observation variance
temp        <- obs.var(NSH.sim[[idxIter]]) # getting observation variance from SAM object for the current iteration
temp        <- subset(temp, temp$fleet == 'catch unique') # subset observation variance
varCatchMat <- array(NA,
                     dim=c(dim(temp)[1],
                           2,
                           nits),
                     dimnames=list('fleet' = temp$fleet,
                                   'var' = c('ages','value')))

# looping through all random samples
for(idxIter in 1:nits){
  
  sdAll     <- obs.var(NSH.sim[[idxIter]]) # getting observation variance from SAM object for the current iteration
  sdSelect  <- subset(sdAll, sdAll$fleet == 'catch unique') # subset observation variance
  
  varCatchMat[,1,idxIter] <- sdSelect$age
  varCatchMat[,2,idxIter] <- sdSelect$value
  
  resi              <- array(0, dim=c(dim(sdSelect)[1],length(yearCurrent))) # initialize array nAges x nYears
  colnames(resi)    <- yearCurrent
  rownames(resi)    <- sdSelect$age
  # fill the residual array
  for(idxResi in 1:dim(sdSelect)[1]){
    resi[idxResi,]  <- rnorm(length(yearCurrent), 
                             0, 
                             sdSelect$value[idxResi])
  }
  
  #C(a,y) = F(a,y)/Z(a,y)*(1-exp(-Z(a,y)))*N(a,y)*res(a,y)
  NSelect   <- stocks[[idxIter]]@stock.n[,match(as.character(yearCurrent),colnames(stocks[[idxIter]]@stock.n))] # number at age filtered to current years
  NSelect   <- drop(NSelect) # drop dimensions with 1 level
  FSelect   <- stocks[[idxIter]]@harvest[,match(as.character(yearCurrent),colnames(stocks[[idxIter]]@stock.n))] # F for current years
  FSelect   <- drop(FSelect) # drop dimensions with 1 level
  Z         <-  raw_M[,match(as.character(yearCurrent),colnames(stocks[[idxIter]]@stock.n))]  # current years for M
  Z         <- Z + FSelect # adding M and F
  
  # catches, with added error based on observation variance
  stocks[[idxIter]]@catch.n[,match(as.character(yearCurrent),colnames(stocks[[idxIter]]@stock.n))] <-  # filter only the current years
    FSelect/Z*(1-exp(-Z))*NSelect*exp(resi)
  
  stocks[[idxIter]]@catch <- computeCatch(stocks[[idxIter]])
  #We don't believe the closure catch data, so put it to NA
  stocks[[idxIter]]@catch.n[,ac(1978:1979)]           <- NA
  
  # landings, i.e. here modelled as catches without error
  stocks[[idxIter]]@landings.n[,match(as.character(yearCurrent),colnames(stocks[[idxIter]]@stock.n))] <-  # filter only the current years
    FSelect/Z*(1-exp(-Z))*NSelect
  stocks[[idxIter]]@landings <- computeLandings(stocks[[idxIter]])
  #We don't believe the closure catch data, so put it to NA
  #stocks[[idxIter]]@landings.n[,ac(1978:1979)]           <- NA
}

# plotting the catches
ageIdx  <- 3

par(mfrow=c(3,3))
for(ageIdx in 1:dim(NSH@catch.n)[1]){
  a     <- drop(NSH@catch.n)
  years <- histMinYr:histMaxYr # vector the years
  plot(years,a[ageIdx,],type='l',ylab=c('age',ac(ageIdx-1)))
  for(idxIter in 1:nits){
    b     <- stocks[[idxIter]]@catch.n[ageIdx,,,,]
    b     <- b[,match(years,colnames(b))]
    lines(years,b,col='green')
  }
}

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

FCPropIts         <- array( 0, # initialize array 
                            dim=c(length(fullPeriod),nits), # nAges x nYears x nits
                            dimnames = list(as.numeric(fullPeriod),
                                            1:nits))

# fill in proportion for future years
for(idxIter in 1:nits){
  # fill in historical period, no randominzation in these.
  FCPropIts[histPeriod,idxIter] <- FCProp
  
  # fill in projection period with randomization (see above comments)
  FCPropIts[projPeriod,idxIter] <-   array(FCProp[ac(yrChainFC[[idxIter]])],
                                           dim=length(FCPropIts[projPeriod,idxIter]))
}

#-------------------------------------------------------------------------------
# 9) D fleet: fixed catch per year (following TAC). One will vary the split 
# between NSAS and WBSS
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# 10) create selection patterns for the different fleets
# Selectivity of fleet projected forward using a random walk, using results from
# the multi-fleet assessment
#-------------------------------------------------------------------------------

fleets <- c('A','C','BD')

############# initialize FLQuant object containing F and selec 193tivities #############
dmns        <- dimnames(NSH@harvest)
dmns$unit   <- fleets
dmns$season <- c('F','sel','catch.wt')
dmns$year   <- projPeriod
dmns$iter   <- 1:nits

fisheryFuture    <- FLQuant(array( NA, # covariance matrix using a period of 10 years for all the ages
                                   dim=c(length(ages), # ages
                                         nFutureyrs,   # years
                                         3,            # fleet (4)
                                         3,            # quantity stored (F, sel)
                                         1,
                                         nits)), # iterations
                            dimnames=dmns)


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
    fisheryFuture[,,fleets[idxFleet],'sel',,idxIter] <- rwF[,,,,,idxIter]
  }
}

############# plot selectivity #############


############# fill in FLQuant object with Fs #############
for(idxIter in 1:nits){
  fisheryFuture[,1,'A','F',,idxIter]   <- NSH3f.sam@harvest[,'2018',,,'A']
  fisheryFuture[,1,'C','F',,idxIter]   <- NSH3f.sam@harvest[,'2018',,,'C']
  fisheryFuture[,1,'BD','F',,idxIter]  <- NSH3f.sam@harvest[,'2018',,,'BD']
  
  fisheryFuture[,,'A','catch.wt',,idxIter]   <- multiFleet_catch.wt[,,1,idxIter]
  fisheryFuture[,,'BD','catch.wt',,idxIter]   <- multiFleet_catch.wt[,,2,idxIter]
  fisheryFuture[,,'C','catch.wt',,idxIter]  <- multiFleet_catch.wt[,,3,idxIter]
}

#-------------------------------------------------------------------------------
# 10) tidying up and saving objects for next step
#-------------------------------------------------------------------------------

# prepare objects for MSE: 
# stkAssessement with the current assessment, including 2018
# newstk up to 2017
C1 <- dim(stocks[[1]]@catch.wt)
C1[6]<- 10
C2 <- C1
C2[1]<- 1

dmns        <- dimnames(stocks[[1]]@harvest)
dmns$iter   <- 1:nits
dmns2        <- dimnames(stocks[[1]]@catch)
dmns2$iter   <- 1:nits

stkAssessement  <- FLStock( catch.wt=FLQuant(dim=C1,dimnames=dmns),
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

newstk          <- FLStock( catch.wt=FLQuant(dim=C1,dimnames=dmns),
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

units(newstk) <- units(NSH)

stkAssessement <- window(  stkAssessement,
                           start=an(histPeriod[1]),
                           end=max(an(histPeriod))+1)

for(idxIter in 1:nits){
  iter(stkAssessement,idxIter)          <- stocks[[idxIter]][,c(histPeriod,projPeriod[1])]
  iter(stkAssessement@harvest,idxIter)  <- NSH.sim[[idxIter]]@harvest
  iter(stkAssessement@stock.n,idxIter)  <- NSH.sim[[idxIter]]@stock.n
  
  iter(newstk,idxIter)                  <- stocks[[idxIter]]
}

stkAssessement@stock <- computeStock(stkAssessement)

stocks <- newstk

units(stkAssessement) <- units(NSH)
units(stocks)         <- units(NSH)

# save objects for MSE
save(stocks, 
     surveys,
     FCPropIts,
     fisheryFuture,
     stkAssessement,
     NSH.ctrl,
     qMat,
     varCatchMat,
     varSurvMat,
     file=file.path(outPath,paste0(assessment_name,'_init_MSE.RData')))

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
     file=file.path(outPath,paste0(assessment_name,'_parameters_MSE.RData')))
