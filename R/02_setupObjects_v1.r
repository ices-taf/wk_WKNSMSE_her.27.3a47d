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
#library(R.utils)

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

#sourceDirectory(functionPath) # loading function

source(file.path(functionPath,"randBlocks.R"))
source(file.path(functionPath,"randNums.R"))


#-------------------------------------------------------------------------------
# 2) load assessment objects (single and multi fleet)
#    define MSE parameters
#    load raw M
#-------------------------------------------------------------------------------

#- Load single fleet and multi fleet assessment objects - using assessment without the LAI
#load(file.path(outPath,paste0(assessment_name,'_mf.Rdata')))
#load(file.path(outPath,paste0(assessment_name,'_sf.Rdata')))
load(file.path(outPath,paste0(assessment_name,'_mf_noLAI.Rdata')))
load(file.path(outPath,paste0(assessment_name,'_sf_noLAI.Rdata')))

NSH.ctrlnoLAI <- drop.from.control(NSH.ctrl,fleet="LAI-CNS") # test for no LAI

# parameters
n.retro.years       <-  7                                       # Number of years for which to run the retrospective
nFutureyrs          <- 20
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
# hack to set plus group
NSHM2             <- readFLStock(file.path(dataPath,"index.txt"),no.discards=TRUE,quiet=FALSE)
NSHM2@m[]         <- as.matrix(raw_M)
pg                <- NSH@range['max']
NSHM2             <- setPlusGroup(NSHM2,pg)
raw_M             <- drop(NSHM2@m)
  
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

# update FLStocks object
for(idxIter in 1:nits){
  stocks[[idxIter]]@catch.n                    <- stocks[[idxIter]]@stock.n * stocks[[idxIter]]@harvest / 
                                                  (stocks[[idxIter]]@harvest + stocks[[idxIter]]@m) * 
                                                  (1 - exp(-stocks[[idxIter]]@harvest - stocks[[idxIter]]@m)) # compute catch numbers 
  stocks[[idxIter]]@catch.n[,ac(1978:1979)]    <- NA # fill in Na for the the closure catch data
  stocks[[idxIter]]@landings.n                 <- stocks[[idxIter]]@catch.n

  stocks[[idxIter]]@harvest.spwn[,projPeriod]  <- stocks[[idxIter]]@harvest.spwn[,ac(histMaxYr)] # propagate Fprop before spawning
  stocks[[idxIter]]@m.spwn[,projPeriod]        <- stocks[[idxIter]]@m.spwn[,ac(histMaxYr)] # propagate Mprop before spawning
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
#-------------------------------------------------------------------------------

# generate random blocks for weight at age and maturity
yrChain   <- randBlocks(an(fecYears),an(projPeriod),nits)
yrChainM  <- randBlocks(an(fecYears),an(projPeriod),nits)

# update FLStocks object
for(idxIter in 1:nits){
  # future maturity at age
  stocks[[idxIter]]@mat      [,projPeriod][]               <- array(stocks[[idxIter]]@mat[,ac(yrChain[[idxIter]])],
                                                                    dim=dim(stocks[[idxIter]]@mat[,projPeriod]))
  # future weight at age
  stocks[[idxIter]]@stock.wt [,projPeriod][]               <- array(stocks[[idxIter]]@stock.wt[,ac(yrChain[[idxIter]])],
                                                                    dim=dim(stocks[[idxIter]]@stock.wt[,projPeriod]))
  # future natural mortality at age
  stocks[[idxIter]]@m        [,projPeriod][]               <- array(stocks[[idxIter]]@m[,ac(yrChainM[[idxIter]])],
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
# Note: the raw M is used here, as opposed to the smoothed one used in the 
# assessment
# 
# The following survey indices are to be generatred
# HERAS age 1 to 8
# IBTS-Q1 age 1
# IBTS0 age 0
# IBTS-Q3 age 0 to 5
# LAI-ORSH age 0, 2 weeks
# LAI-BUN age 0, 2 weeks
# LAI- CNS age 0, 4 weeks
# LAI-SNS age 0, 3 weeks
#-------------------------------------------------------------------------------
surveys     <- lapply(NSH.tun,propagate,iter=nits)
surveys     <- window(window(surveys,end=histMaxYr+1),start=histMinYr,end=futureMaxYr) # extend the FLStock object to the full projection period

for(idxIter in 1:nits){
  
  sdAll <- obs.var(NSH.sim[[idxIter]]) # getting observation variance from SAM object for the current iteration
  qAll  <- catchabilities(NSH.sim[[idxIter]]) # getting catchabilities from SAM object for the current iteration
  
  surveyNames <- as.character(unique(qAll$fleet))
  
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
    surveys[[surveyNames[idxSurvey]]]@index[,match(as.character(yearCurrentSurvey),colnames(stocks[[idxIter]]@stock.n)) # filling only the years available for the survey
                                            ,1,1,1,idxIter] <- as.matrix(replicate(length(yearCurrentSurvey), qSelect$value)*exp(-Z*surveyProp)*NSelect*exp(resi))
  }
}

# plotting survey indices
ageIdx  <- 4
surveyName <- 'HERAS'
a     <-drop(NSH.tun[[surveyName]]@index)
#a     <- NSH.tun[[surveyName]]@index # in case of only 1 age in survey object
years <- as.numeric(colnames(NSH.tun[[surveyName]]@index))
plot(years,a[ageIdx,],type='l')

for(idxIter in 1:nits){
  b     <-surveys[[surveyName]]@index[ageIdx,,,,,idxIter]
  b     <- b[,match(years,colnames(b))]
  lines(years,b,col='green')
}


#-------------------------------------------------------------------------------
# 7) creating catches from each random samples
# 
# survey indices are created for each random sample using F and N and M at age
# residuals are added using a normal distribution with sigma as sd of observations
# C(a,y) = F(a,y)/Z(a,y)*(1-exp(-Z(a,y)))*N(a,y)*res(a,y)
#-------------------------------------------------------------------------------

# looping through all random samples
for(idxIter in 1:nits){
  
  sdAll     <- obs.var(NSH.sim[[idxIter]]) # getting observation variance from SAM object for the current iteration
  sdSelect  <- subset(sdAll, sdAll$fleet == 'catch unique') # subset observation variance
  
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
  FSelect   <-  stocks[[idxIter]]@harvest[,match(as.character(yearCurrent),colnames(stocks[[idxIter]]@stock.n))] # F for current years
  FSelect   <- drop(FSelect) # drop dimensions with 1 level
  Z         <-  raw_M[,match(as.character(yearCurrent),colnames(stocks[[idxIter]]@stock.n))]  # current years for M
  Z         <- Z + FSelect # adding M and F
  
  # catches, with added error based on observation variance
  stocks[[idxIter]]@catch.n[,match(as.character(yearCurrent),colnames(stocks[[idxIter]]@stock.n))] <-  # filter only the current years
                                                                                                      FSelect/Z*(1-exp(-Z))*NSelect*exp(resi)
  #We don't believe the closure catch data, so put it to NA
  stocks[[idxIter]]@catch.n[,ac(1978:1979)]           <- NA
  # landings
  stocks[[idxIter]]@landings.n[,match(as.character(yearCurrent),colnames(stocks[[idxIter]]@stock.n))] <-  # filter only the current years
                                                                                                      FSelect/Z*(1-exp(-Z))*NSelect
  #We don't believe the closure catch data, so put it to NA
  stocks[[idxIter]]@landings.n[,ac(1978:1979)]           <- NA
}

# plotting the catches
ageIdx  <- 3
a     <- drop(NSH@catch.n)
years <- histMinYr:histMaxYr # vector the years
plot(years,a[ageIdx,],type='l')
for(idxIter in 1:nits){
  b     <- stocks[[idxIter]]@catch.n[ageIdx,,,,]
  b     <- b[,match(years,colnames(b))]
  lines(years,b,col='green')
}

#-------------------------------------------------------------------------------
# 8) C fleet: proportion of F of the C fleet for the future years
# The proportion is obtained from the multi-fleet assessment using a random
# draw similar to M and weight at age
# For the given quantity, we randomize the following:
# number of years in the chain
# start year in the chain
# reversing of the chain
#
# Note 1: multi fleet assessment is performed as a multi-fleet from 1997 onward 
# and as a single fleet assessment prior to 1997 (therefore an FLStock object 
# with) 2 fields
# Note 2: NSAS/WBSS split for the C fleet will be kept constant to 
# 30 NSAS/70 WBSS
# Note 3: we don't randomize the current years as we don't simulate different 
# time series for the multi-fleet assessment
#-------------------------------------------------------------------------------

yrChainFC   <- randBlocks(an(fecYears),an(projPeriod),nits)
ages        <- NSH3f.sam@range[1]:NSH3f.sam@range[2]

# initializing object
FAProp            <- array(0, dim=c(length(ages),length(yearCurrent))) # initialize array nAges x nYears
colnames(FAProp)  <- yearCurrent
rownames(FAProp)  <- ages

FBDProp            <- array(0, dim=c(length(ages),length(yearCurrent))) # initialize array nAges x nYears
colnames(FBDProp)  <- yearCurrent
rownames(FBDProp)  <- ages

FCProp            <- array(0, dim=c(length(ages),length(yearCurrent))) # initialize array nAges x nYears
colnames(FCProp)  <- yearCurrent
rownames(FCProp)  <- ages

for(idxAge in 1:length(ages)){
  Ftot              <-  NSH3f.sam@harvest[idxAge,1:length(yearCurrent),,,1] + # fleet A
                        NSH3f.sam@harvest[idxAge,1:length(yearCurrent),,,2] + # fleet BD
                        NSH3f.sam@harvest[idxAge,1:length(yearCurrent),,,3]   # fleet C
  Ftot              <-  drop(Ftot)
  FAProp[idxAge,]   <-  drop(NSH3f.sam@harvest[idxAge,1:length(yearCurrent),,,1])/Ftot
  FBDProp[idxAge,]  <-  drop(NSH3f.sam@harvest[idxAge,1:length(yearCurrent),,,2])/Ftot
  FCProp[idxAge,]   <-  drop(NSH3f.sam@harvest[idxAge,1:length(yearCurrent),,,3])/Ftot
}

FCProp            <-  cbind(FCProp,
                            array(0,dim=c(length(ages),nFutureyrs))) # adding future years

FCPropIts         <- array( 0, # initialize array 
                            dim=c(length(ages),length(fullPeriod),nits), # nAges x nYears x nits
                                  dimnames = list(ages,
                                                  as.numeric(fullPeriod),
                                                  1:nits))

# fill in proportion for future years
for(idxIter in 1:nits){
  # fill in historical period, no randominzation in these.
  FCPropIts[,histPeriod,idxIter] <- FCProp
  
  # fill in projection period with randomization (see above comments)
  FCPropIts[,projPeriod,idxIter] <- array(FCProp[,ac(yrChainFC[[idxIter]])],
                                          dim=dim(FCPropIts[,projPeriod,idxIter]))
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


dmns                                                          <- dimnames(NSH@harvest)
dmns$year                                                     <- projPeriod
dmns$iter                                                     <- 1:nits
ages                                                          <- dimnames(stocks@stock.n)$age

#- Create random walk over Fs (as multiplier from last years selection pattern)
# covariance over 10 years for age 0 to 2
covmat1                                                        <- cov(apply(log(NSH@harvest[,ac((max(yearCurrent)-10):max(yearCurrent)),drop=T]), # covariance for the last 10 years
                                                                            1,
                                                                            diff))
covmat10Y                 <- covmat1 # covariance computed for all the ages over the 10 years
covmat1[ac(0:2),ac(3:8)]  <- 0
covmat1[ac(3:8),ac(0:2)]  <- covmat1[ac(0:2),ac(3:8)]
covmat1[ac(3:8),ac(3:8)]  <- covmat1[ac(3:8),ac(0:2)]
# covariance over 20 years for age 3 and up
covmat2                                                       <- cov(apply(log(NSH@harvest[,ac((max(yearCurrent)-20):max(yearCurrent)),drop=T]), # covariance over the last 20 years
                                                                           1,
                                                                           diff))
covmat2[ac(0:2),ac(0:2)]  <- 0;
covmatMix                 <- covmat1 + covmat2 # mix of year period for different ages

# multivariate normal distribution for the covariance
mvrnorm(nits*nFutureyrs,rep(0,length(ages)),cov(apply(log(NSH@harvest[,selPeriod,drop=T]),1,diff)))

# using 10 year period for the covariance for all the ages
wF                                                            <- FLQuant(array(t(mvrnorm(nits*nFutureyrs,
                                                                                         rep(0,length(ages)),
                                                                                         covmat10Y)), # covariance matrix using a period of 10 years for all the ages
                                                                               dim=c(length(ages),nFutureyrs,1,1,1,nits)),
                                                                         dimnames=dmns)

# using 10 year period for the covariance for age 0 to 2 and 20 year period for ages 3+
wF                                                            <- FLQuant(array(t(mvrnorm(nits*nFutureyrs,
                                                                                         rep(0,length(ages)),
                                                                                         covmatMix)), # covariance matrix using 10 years for age 0 to 2 and 20 years for ages 3+
                                                                               dim=c(length(ages),nFutureyrs,1,1,1,nits)),
                                                                         dimnames=dmns)

qtil                                                          <- quantile(c(wF),probs=c(0.05,0.95))
qtilold                                                       <- quantile(c(wF[ac(4:8),]),probs=c(0.25,0.75))
wF@.Data[which(wF<qtil[1] | wF>qtil[2])][]                    <- 0
wF[ac(4:8),]@.Data[which(wF[ac(4:8),]<qtilold[1] | wF[ac(4:8),]>qtilold[2])][] <- 0
rwF                                                           <- FLQuant(aperm(apply(wF,c(1,3:6),cumsum),c(2,1,3:6)),dimnames=dmns)
rwF                                                           <- sweep(rwF,c(1,3:5),apply(log(NSH@harvest[,selPeriod]),1,mean),"+")
fbarages                                                      <- ac(range(NSH)["minfbar"]:range(NSH)["maxfbar"])
landsel                                                       <- sweep(exp(rwF),c(2:6),apply(exp(rwF)[fbarages,],2:6,mean),"/")

for(iFsh in 1:dims(fishery)$unit)
  landings.sel(fishery)[,projPeriod,iFsh]                     <- FLQuant(sweep(landsel,c(1),propN[,,iFsh],"*"),
                                                                         dimnames=dimnames(stocks[,ac(projPeriod)]@stock.n))
plot(landings.sel(fishery[,projPeriod,1]))


#-------------------------------------------------------------------------------
# 10) tidying up and saving objects for next step
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#
#
#
#
#
#
#
#
#
#
#
#
# Old code from there (2015)
#
#
#
#
#
#
#
#
#
#
#
#
#-------------------------------------------------------------------------------

  #-------------------------------------------------------------------------------
  # 1): Create biological population object
  #-------------------------------------------------------------------------------

#HAWG TEST
#load(file.path(inPath,"scaleNsNSAS.RData"))
#load(file.path(inPath,"scaleFsNSAS.RData"))
#scaleNsNew                <- FLQuant(NA,dimnames=dimnames(iter_retro[[1]][[ac(2013)]]@stock.n))
#scaleNsNew[,ac(1947:2013)]<- scaleNs
#scaleNsNew[,ac(2014)]     <- scaleNs[,ac(2013)]
#scaleFsNew                <- FLQuant(NA,dimnames=dimnames(iter_retro[[1]][[ac(2013)]]@harvest))
#scaleFsNew[,ac(1947:2013)]<- scaleFs
#scaleFsNew[,ac(2014)]     <- scaleFs[,ac(2013)]
#stocks2                   <- stocks
#for(iTer in 1:nits){
#  stocks2@harvest[,1:68,,,,iTer]      <- iter_retro[[iTer]][[ac(2013)]]@harvest * scaleFsNew
#  stocks2@stock.n[,1:68,,,,iTer]      <- iter_retro[[iTer]][[ac(2013)]]@stock.n * scaleNsNew
#}
#END HAWG TEST
stocks2                   <- stocks
for(iTer in 1:nits){
  stocks2@harvest[,1:68,,,,iTer]      <- iter_retro[[iTer]][[ac(2013)]]@harvest
  stocks2@stock.n[,1:68,,,,iTer]      <- iter_retro[[iTer]][[ac(2013)]]@stock.n
}

save(stocks,stocks2,file=file.path(outPath,"stocks_stocks2.RData"))
stocks                    <- stocks2

biol                      <- as.FLBiol(stocks)
  #- Random draw from lognormal distribution for new recruitment, estimate lognormal parameters first
recrAge                   <- dimnames(rec(stocks))$age
pars                      <- optim(par=c(17.1,0.20),fn=optimRecDistri,recs=sort(c(rec(NSH[,ac(recrPeriod)]))),
                                   method="Nelder-Mead")$par

biol@n[1,projPeriod]      <- rtlnorm(length(projPeriod)*nits,mean=pars[1],sd=pars[2],lower=0.01*min(biol@n[recrAge,],na.rm=T))


  #-------------------------------------------------------------------------------
  # 2): Create fisheries object
  #-------------------------------------------------------------------------------

dmns                      <- dimnames(m(biol))
dmns$unit                 <- c("A","B","C","D")
fishery                   <- FLCatch(price=FLQuant(NA,dimnames=dmns))
name(fishery)             <- "catches"
desc(fishery)             <- "North Sea Herring"
fishery@range             <- range(biol)

  #-------------------------------------------------------------------------------
  #- Get the proportional contribution of each fishery to the landings and weight
  #-------------------------------------------------------------------------------

    #-------------------------------------------------------------------------------
    #- Partial Ns per fleet and plusgroup setting
    #-------------------------------------------------------------------------------
partialN                  <- read.csv(paste(inPath,"partial_ns.csv",sep=""),header=T)
  dimnames(partialN)[[1]] <- ac(0:9) #Biological sampling is up to age 9
  Ns                      <- partialN
  partialN[partialN==0]   <- NA
  pg                      <- range(stocks)["plusgroup"]
  pgplus                  <- dimnames(partialN)[[1]][which(dimnames(partialN)[[1]] > pg)]
  partialN[ac(pg),]       <- apply(Ns[ac(c(pg:pgplus)),],2,sum,na.rm=T); partialN <- partialN[ac(dimnames(partialN)[[1]][1]:pg),]
  idx                     <- lapply(as.list((histMaxYr-2):histMaxYr),function(x){grep(x,colnames(partialN))})
  dmns$year               <- ac(histMaxYr+1); dmns$iter <- 1;
  propN                   <- FLQuant(NA,dimnames=dmns); propWt <- FLQuant(NA,dimnames=dmns)
  N                       <- array(unlist(lapply(idx,function(x){sweep(partialN[,x],1,rowSums(partialN[,x],na.rm=T),"/")})),
                                   dim=c(length(dimnames(partialN)[[1]]),length(dmns$unit),3))
  N[is.na(N)]             <- 0
                          #check: apply(N,c(1,3),sum,na.rm=T) #must equal 1
propN[]                   <- apply(N,1:2,mean,na.rm=T)

    #-------------------------------------------------------------------------------
    #- Partial Wts per fleet and plusgroup setting
    #-------------------------------------------------------------------------------
partialWt                 <- read.csv(paste(inPath,"partial_ws.csv",sep=""),header=T)
  dimnames(partialWt)[[1]]<- ac(0:9) #Biological sampling is up to age 9
  partialWt[partialWt==0] <- NA
  #- Define the average (weighted) weight-at-age and calculate the deviation of each fleet from that average weight
  #  We assume that the combination of the A,B,C and D fleet together make up the average weight at age in the fishery

  #- Plusgroup correction
  res                     <- partialWt * Ns
  res[ac(pg),]            <- colSums(res[ac(c(pg:pgplus)),],na.rm=T); res <- res[ac(dimnames(res)[[1]][1]:pg),]
  res                     <- res / partialN
  partialWt               <- res
  
  Wt                      <- array(unlist(lapply(idx,function(x){sweep(partialWt[,x],1,(rowSums(partialWt[,x]*partialN[,x],na.rm=T)/rowSums(partialN[,x],na.rm=T)),"/")})),
                                   dim=dim(N))
                          #check:  apply(Wt * N,c(1,3),sum,na.rm=T) #must equal 1
propWt[]                  <- apply(Wt,1:2,mean,na.rm=T)

  #- Put all proportions equal to NA to zero, so that they don't get any weight
propN@.Data[is.na(propN)==T][]  <- 0; propWt@.Data[is.na(propWt)==T][] <- 0

  #-Take single fleet weights and numbers and multiply by the proportions
for(iFsh in dimnames(fishery@landings.n)$unit){
  fishery@landings.n[,  ac(histMinYr:histMaxYr),iFsh]         <- sweep(NSH@landings.n,1,propN[,,iFsh],"*")*NSH@landings.wt / sweep(NSH@landings.wt,1,propWt[,,iFsh],"*")
  fishery@landings.wt[, ac(histMinYr:histMaxYr),iFsh]         <- sweep(NSH@landings.wt,1,propWt[,,iFsh],"*")
  fishery@discards.n[,  ac(histMinYr:histMaxYr),iFsh]         <- 0
  fishery@discards.wt[, ac(histMinYr:histMaxYr),iFsh]         <- 0
}
fishery@landings.n@.Data[is.infinite(fishery@landings.n)==T]  <- 0
fishery@landings.n@.Data[is.na(fishery@landings.n)==T]        <- 0
fishery@landings[,    ac(histMinYr:histMaxYr)]                <- computeLandings(fishery[,ac(histMinYr:histMaxYr)])
fishery@discards[,    ac(histMinYr:histMaxYr)]                <- computeDiscards(fishery[,ac(histMinYr:histMaxYr)])
                          #check: computeLandings(NSH) / window(unitSums(fishery@landings),1947,2014) #must equal 1

  #-Calculate deterministic landing.sel
for(iFsh in dimnames(fishery@landings.sel)$unit)
  landings.sel(fishery)[,ac(histMinYr:histMaxYr),iFsh]        <- FLQuant(sweep(sweep(harvest(stocks[,ac(histMinYr:histMaxYr)]),c(1),propN[,,iFsh],"*"),2:6,
                                                                               fbar(stocks[,ac(histMinYr:histMaxYr)]),"/"),
                                                                         dimnames=dimnames(stocks[,ac(histMinYr:histMaxYr)]@stock.n))

catch.q(     fishery)[]                                       <- 1
discards.sel(fishery)[]                                       <- 0
fishery@discards.wt[]                                         <- 0
fishery@discards.n[]                                          <- 0

  #-------------------------------------------------------------------------------
  #- Vary selectivity of fleet (add random walk to last year selection pattern)
  #-------------------------------------------------------------------------------

dmns                                                          <- dimnames(NSH@harvest)
dmns$year                                                     <- projPeriod
dmns$iter                                                     <- 1:nits
ages                                                          <- dimnames(stocks@stock.n)$age

#- Create random walk over Fs (as multiplier from last years selection pattern)
covmat1                                                        <- cov(apply(log(NSH@harvest[,ac(2003:2013),drop=T]),1,diff))
covmat1[ac(3:8),ac(3:8)] <- covmat1[ac(3:8),ac(0:2)] <- covmat1[ac(0:2),ac(3:8)] <- 0
covmat2                                                       <- cov(apply(log(NSH@harvest[,ac(1997:2013),drop=T]),1,diff))
covmat2[ac(0:2),ac(0:2)] <- 0; covmat <- covmat1 + covmat2
wF                                                            <- FLQuant(array(t(mvrnorm(nits*nFutureyrs,rep(0,length(ages)),cov(apply(log(NSH@harvest[,selPeriod,drop=T]),1,diff)))),
                                                                               dim=c(length(ages),nFutureyrs,1,1,1,nits)),
                                                                         dimnames=dmns)
wF                                                            <- FLQuant(array(t(mvrnorm(nits*nFutureyrs,rep(0,length(ages)),covmat)),
                                                                               dim=c(length(ages),nFutureyrs,1,1,1,nits)),
                                                                         dimnames=dmns)
qtil                                                          <- quantile(c(wF),probs=c(0.05,0.95))
qtilold                                                       <- quantile(c(wF[ac(4:8),]),probs=c(0.25,0.75))
wF@.Data[which(wF<qtil[1] | wF>qtil[2])][]                    <- 0
wF[ac(4:8),]@.Data[which(wF[ac(4:8),]<qtilold[1] | wF[ac(4:8),]>qtilold[2])][] <- 0
rwF                                                           <- FLQuant(aperm(apply(wF,c(1,3:6),cumsum),c(2,1,3:6)),dimnames=dmns)
rwF                                                           <- sweep(rwF,c(1,3:5),apply(log(NSH@harvest[,selPeriod]),1,mean),"+")
fbarages                                                      <- ac(range(NSH)["minfbar"]:range(NSH)["maxfbar"])
landsel                                                       <- sweep(exp(rwF),c(2:6),apply(exp(rwF)[fbarages,],2:6,mean),"/")

for(iFsh in 1:dims(fishery)$unit)
  landings.sel(fishery)[,projPeriod,iFsh]                     <- FLQuant(sweep(landsel,c(1),propN[,,iFsh],"*"),
                                                                         dimnames=dimnames(stocks[,ac(projPeriod)]@stock.n))
plot(landings.sel(fishery[,projPeriod,1]))
  #-------------------------------------------------------------------------------
  #- Vary landing weights (sample from observed landing weights but with some sort
  #  of autorcorrelation to biol)
  #-------------------------------------------------------------------------------
projFishLandwt                                                <- array(iter(NSH@landings.wt,1)[,ac(yrStrngsC)],dim=dim(fishery@landings.wt[,projPeriod,1]))
for(iFsh in dimnames(fishery@landings.wt)$unit)
  fishery@landings.wt[,projPeriod,iFsh]                       <- sweep(projFishLandwt,1,propWt[,,iFsh],"*")


  #-------------------------------------------------------------------------------
  # 4): Save the objects
  #-------------------------------------------------------------------------------
save(biol          ,file=paste(outPath,"biol.RData",          sep=""))
  save(pars        ,file=paste(outPath,"recPars.RData",       sep=""))
save(fishery       ,file=paste(outPath,"fishery.RData",       sep=""))
  save(propN       ,file=paste(outPath,"propN.RData",         sep=""))
  save(propWt      ,file=paste(outPath,"propWt.RData",        sep=""))
  save(ctch        ,file=paste(outPath,"ctch.RData",          sep=""))
  save(landsel     ,file=paste(outPath,"landsel.RData",       sep=""))
save(surveys       ,file=paste(outPath,"surveys.RData",       sep=""))
  save(surv        ,file=paste(outPath,"surv.RData",          sep=""))
  save(surveyQ     ,file=paste(outPath,"surveyQ.RData",       sep=""))
  save(surveyK     ,file=paste(outPath,"surveyK.RData",       sep=""))
save(stocks        ,file=paste(outPath,"stocks.RData",        sep=""))
save(settings      ,file=paste(outPath,"settings.RData",      sep=""))
save.image(         file=paste(outPath,"setup06052015.RData", sep=""))

