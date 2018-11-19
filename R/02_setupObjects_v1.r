#-------------------------------------------------------------------------------
# WKNSMSE2018
#
# Author: Niels Hintzen
#         IMARES, The Netherland
#
#  MSE of North Sea Herring
#
# Date: 2018/11/18
#
# Build for R3.4.3, 64bits
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# WKNSMSE
#
# Author: Niels Hintzen
#         IMARES, The Netherland
#
#  MSE of North Sea Herring
#
# Date: 2018/11/18
#
# Build for R3.4.3, 64bits
#-------------------------------------------------------------------------------

rm(list=ls())

library(FLSAM)
library(FLEDA)

#path          <- "D:/Work/Herring MSE/NSAS/"
path              <- "D:/git/wk_WKNSMSE_her.27.3a47d/R/"
assessment_name   <- "NSAS_WKNSMSE2018"
try(setwd(path),silent=TRUE)

# paths to different subfolders
dataPath      <- file.path(".","data/")
outPath       <- file.path(".","results/")
scriptPath    <- file.path(".","side_scripts/")

#- Load objects
load(file.path(outPath,paste0(assessment_name,'_mf.Rdata')))
load(file.path(outPath,paste0(assessment_name,'_sf.Rdata')))

#- Settings
nyrs        <- 27
histMinYr   <- dims(NSH)$minyear
histMaxYr   <- dims(NSH)$maxyear
futureMaxYr <- histMaxYr + nyrs
histPeriod  <- ac(histMinYr:histMaxYr)
projPeriod  <- ac((histMaxYr+1):futureMaxYr)
recrPeriod  <- ac(2007:2017)
selPeriod   <- ac(2007:2017)
fecYears    <- ac(2007:2017)
niter        <- 1000
settings    <- list(histMinYr=histMinYr,
                    histMaxYr=histMaxYr,
                    futureMaxYr=futureMaxYr,
                    histPeriod=histPeriod,
                    projPeriod=projPeriod,
                    recrPeriod=recrPeriod,
                    nyrs=nyrs,
                    nits=nits,
                    fecYears=fecYears)

source(paste(codePath,"functions.r",sep=""))

  #-------------------------------------------------------------------------------
  # 0): Create stock object & use vcov for new realisations
  #-------------------------------------------------------------------------------
stocks                            <- monteCarloStock(NSH,NSH.sam,niter) # function not working

stocks                            <- window(window(stocks,end=histMaxYr+1),start=histMinYr,end=futureMaxYr)
stocks@catch.n                    <- stocks@stock.n * stocks@harvest / (stocks@harvest + stocks@m) * (1 - exp(-stocks@harvest - stocks@m))
stocks@catch.n[,ac(1978:1979)]    <- NA
stocks@landings.n                 <- stocks@catch.n

stocks@harvest.spwn[,projPeriod]  <- stocks@harvest.spwn[,ac(histMaxYr)]
stocks@m.spwn[,projPeriod]        <- stocks@m.spwn[,ac(histMaxYr)]

#- Sample from selected year range to fill future time period of m, weight, fec and landings.wt
#  Retain a certain degree of auto-correlation and take fec, weight, m and landings
#  weight from the same years to retain correlation between fec, wt, landings.wt
#  and m (if there is any with m...)
#  Take blocks of years instead of randomly years glued together. In some instances,
#  reverse the blocks of year for variation

  #- Take blocks of years, sample blocks from fecYears and add up till length of timeseries
yrs       <- projPeriod
saveBlcks <- matrix(NA,nrow=nits,ncol=length(yrs))

full      <- F
rw        <- 1
# replaced the version below because of occassional errors
while(full==F){
  tmpx <- 0
  while(tmpx<1) {
  set.seed(as.integer(Sys.time()))
  sam   <- sample(1:length(fecYears),nits*1000,replace=T)
  samM  <- matrix(sam,nrow=nits,ncol=1000)
  idx   <- which(t(apply(samM,1,cumsum)) == length(yrs),arr.ind=T)
  tmpx <- length(idx)
  }
  if((rw+nrow(idx) - 1) > nrow(saveBlcks)){ stprow <- nrow(saveBlcks); full <- T} else { stprow <- rw+nrow(idx) -1}
  for(iRow in 1:(stprow-rw+1))
    if ((rw+iRow-1) <= nrow(saveBlcks) && iRow>0) saveBlcks[(rw+iRow-1),1:idx[iRow,2]] <- samM[idx[iRow,1],1:idx[iRow,2]]
  rw    <- stprow + 1
}

  #- Take the sampled blocks and assign years to it
saveBlcksYrs <- array(NA,dim=c(nits,length(yrs),2),dimnames=list(nits=1:nits,blocks=1:length(yrs),strtstp=c("start","stop")))
for(iCol in 1:ncol(saveBlcksYrs)){
  #-saveBlcks - 1 because fecYears[1] is the first possible year
  #-saveBlcks + 2 because rev(fecYears)[1] is the last possible year and correction for as.integer
  #-If blck is reversed is decided afterwards, hence the tight bounds in runif
  strstp  <- as.integer(runif(nits,an(fecYears[1])+saveBlcks[,iCol]-1,an(rev(fecYears)[1])-saveBlcks[,iCol]+2))
  rv      <- sample(c(T,F),nits,replace=T)
  strt    <- ifelse(rv==F,strstp,strstp-saveBlcks[,iCol]+1)
  stp     <- strt + saveBlcks[,iCol] - 1
  saveBlcksYrs[,iCol,"start"]       <- ifelse(rv==F,strt,stp)
  saveBlcksYrs[,iCol,"stop"]        <- ifelse(rv==F,stp,strt)

  #-Correct for those records where only one option in year is possible
  idx     <- which((an(fecYears[1])+saveBlcks[,iCol]-1) == an(rev(fecYears)[1])-saveBlcks[,iCol]+2)
  saveBlcksYrs[idx,iCol,"start"]    <- (an(fecYears[1])+saveBlcks[,iCol]-1)[idx]
  saveBlcksYrs[idx,iCol,"stop"]     <- ifelse((saveBlcksYrs[idx,iCol,"start"] + saveBlcks[idx,iCol] - 1) > an(rev(fecYears)[1]),
                                              (saveBlcksYrs[idx,iCol,"start"] - saveBlcks[idx,iCol] + 1),
                                              (saveBlcksYrs[idx,iCol,"start"] + saveBlcks[idx,iCol] - 1))

  #-add start-stop for those blocks equal to length(fecYears)
  idx     <- which(saveBlcks[,iCol] == length(fecYears))
  if(length(idx) > 0){
    saveBlcksYrs[idx,iCol,"start"]  <- ifelse(rv[idx]==F,an(fecYears[1]),an(rev(fecYears)[1]))
    saveBlcksYrs[idx,iCol,"stop"]   <- ifelse(rv[idx]==F,an(rev(fecYears)[1]),an(fecYears[1]))
  }

  #-If block is large, sampled size might exceed bounds, so take only possible options
  if(length(idx) > 0) idx           <- which(is.na(strstp)==T & is.na(saveBlcks[,iCol])==F)[which(!which(is.na(strstp)==T &
                                                                                                         is.na(saveBlcks[,iCol])==F) %in% idx)]
  if(length(idx) == 0)idx           <- which(is.na(strstp)==T & is.na(saveBlcks[,iCol])==F)
  if(length(idx) > 0){
    #-Calculate options that are possible
    opts                      <- mapply(function(strt,stp){
                                        return(fecYears[which(!fecYears %in% strt:stp)])},
                                        strt=an(rev(fecYears)[1])-saveBlcks[idx,iCol]+2,
                                        stp =an(fecYears[1])     +saveBlcks[idx,iCol]-2,SIMPLIFY=F)
    #-Sample from those options and define start and stop points
    res     <- an(unlist(lapply(opts,sample,1)))
    saveBlcksYrs[idx,iCol,"start"]  <- res
    saveBlcksYrs[idx,iCol,"stop"]   <- ifelse((res + saveBlcks[idx,iCol] - 1) > an(rev(fecYears)[1]),
                                               res - saveBlcks[idx,iCol] + 1,
                                               res + saveBlcks[idx,iCol] - 1)
  }
}

  #-Create year strings and bind them together
yrStrngsIter<- apply(saveBlcksYrs,1,function(x){mapply(seq,from=na.omit(x[,1]),to=na.omit(x[,2]))})
# if each sampled block was the same length, it creates a matrix instead of a list, gives an error below.
# so added fix to change any matrices into lists first
for (i in 1:nits) {
  if (!is.list(yrStrngsIter[[i]]))  {
    tmp <- list()
    for (cC in 1:ncol(yrStrngsIter[[i]])) tmp[[cC]]<- yrStrngsIter[[i]][,cC]
    yrStrngsIter[[i]] <- tmp   
    rm(tmp)
  } 
}

yrStrngs    <- lapply(yrStrngsIter,function(x){do.call(c,x)})
yrStrngsC   <- do.call(c,yrStrngs)

stocks@mat      [,projPeriod][]               <- array(iter(stocks@mat,1)     [,ac(yrStrngsC)],dim=dim(stocks@mat     [,projPeriod]))
stocks@stock.wt [,projPeriod][]               <- array(iter(stocks@stock.wt,1)[,ac(yrStrngsC)],dim=dim(stocks@stock.wt[,projPeriod]))
stocks@m        [,projPeriod][]               <- array(iter(stocks@m,1)       [,ac(yrStrngsC)],dim=dim(stocks@m       [,projPeriod]))


  #-------------------------------------------------------------------------------
  # 1): Create survey object & use vcov for new realisations + error on realisations
  #-------------------------------------------------------------------------------
surveys           <- lapply(NSH.tun,propagate,iter=nits)
for(iSurv in names(surveys)) surveys[[iSurv]] <- window(surveys[[iSurv]],start=range(surveys[[iSurv]])["minyear"],end=futureMaxYr)
dmns              <- dimnames(surveys[["HERAS"]]@index)
dmns$year         <- dimnames(surveys[["SCAI"]]@index)$year
dmns$unit         <- names(surveys)
surv              <- FLQuant(NA,dimnames=dmns)

  #- Get redrawn survey Qs and Ks
load(file=paste(outPath,"random.param.RData",sep=""))

  #- Get the index of each parameter in the random.param object
Qidx              <- unlist(apply(NSH.ctrl@catchabilities,1,function(x)c(na.omit(x))))

  #- Create objects for surveyQ and surveyK's
surveyQ           <- FLQuants("SCAI"=   FLQuant(NA,dimnames=dimnames(surveys[["SCAI"]]@index)),
                              "IBTS0"=  FLQuant(NA,dimnames=dimnames(surveys[["IBTS0"]]@index)),
                              "IBTS-Q1"=FLQuant(NA,dimnames=dimnames(surveys[["IBTS-Q1"]]@index)),
                              "HERAS"=  FLQuant(NA,dimnames=dimnames(surveys[["HERAS"]]@index)))


surveyK           <- surveyQ
for(iSurv in names(surveys)) surveyK[[iSurv]][]         <- 1


  #- Fill the Qs by survey
for(iYr in dimnames(surveyQ[["SCAI"]])$year)
  surveyQ[["SCAI"]][,iYr]       <- exp(random.param[,which(colnames(random.param) == "logScaleSSB")])
for(iYr in dimnames(surveyQ[["IBTS0"]])$year)
  surveyQ[["IBTS0"]][,iYr]      <- exp(random.param[,which(colnames(random.param) %in% "logFpar")[Qidx[grep("IBTS0",names(Qidx))]]])
for(iYr in dimnames(surveyQ[["IBTS-Q1"]])$year)
  surveyQ[["IBTS-Q1"]][,iYr]    <- exp(random.param[,which(colnames(random.param) %in% "logFpar")[Qidx[grep("IBTS-Q1",names(Qidx))]]])
for(iYr in dimnames(surveyQ[["HERAS"]])$year)
  surveyQ[["HERAS"]][,iYr]      <- t(exp(random.param[,which(colnames(random.param) %in% "logFpar")[Qidx[grep("HERAS",names(Qidx))]]]))

  #- Index var no longer used but filled anyway
for(iSurv in names(surveys)) surveys[[iSurv]]@index.var[,projPeriod][] <- surveys[[iSurv]]@index.var[,ac(histMaxYr)]

  #- Vary selectivity of survey
for(iSurv in names(surveys)){
  # iSurv <- names(surveys)[4]
  cat("Vary selection of:",iSurv,"\n")
  #- Get residuals
  Resids                                            <- subset(residuals(NSH.sam),fleet==iSurv & year <= range(NSH.tun[[iSurv]])["maxyear"])
  iResids                                           <- FLQuant(1,dimnames=c(dimnames(NSH.tun[[iSurv]]@index)[1:5],iter="1"))

  #- Substract residuals (log scale)
  for(i in 1:nrow(Resids)){
    if(iSurv == "SCAI") iResids[1,ac(Resids$year[i]),]                    <- exp(Resids$log.obs[i] - Resids$log.mdl[i])
    if(iSurv != "SCAI") iResids[ac(Resids$age[i]),ac(Resids$year[i]),]    <- exp(Resids$log.obs[i] - Resids$log.mdl[i])
  }

  #- Take blocks of residuals, sample blocks from 1-10 and add up till length of timeseries
  yrs       <- range(surveys[[iSurv]])["minyear"]:futureMaxYr
  saveBlcks <- matrix(NA,nrow=nits,ncol=length(yrs))

  full      <- F
  rw        <- 1
  cat("Sampling iters:",iSurv,"\n")
  while(full==F){
    tmpx <- 0
    while(tmpx<1) {
      set.seed(as.integer(Sys.time()))
      sam   <- sample(1:10,nits*1000,replace=T)
      samM  <- matrix(sam,nrow=nits,ncol=1000)
      idx   <- which(t(apply(samM,1,cumsum)) == length(yrs),arr.ind=T)
      tmpx <- length(idx)
    }
    if((rw+nrow(idx) - 1) > nrow(saveBlcks)){ stprow <- nrow(saveBlcks); full <- T} else { stprow <- rw+nrow(idx) -1}
    for(iRow in 1:(stprow-rw+1))
      if ((rw+iRow-1) <= nrow(saveBlcks) && iRow>0) saveBlcks[(rw+iRow-1),1:idx[iRow,2]] <- samM[idx[iRow,1],1:idx[iRow,2]]
    rw    <- stprow + 1
  }

  #- Take the sampled blocks and assign years to it
  saveBlcksYrs <- array(NA,dim=c(nits,length(yrs),2),dimnames=list(nits=1:nits,yrs=yrs,strtstp=c("start","stop")))
  for(iCol in 1:ncol(saveBlcksYrs)){
    strstp  <- as.integer(runif(nits,range(surveys[[iSurv]])["minyear"]+saveBlcks[,iCol]-1,range(NSH.tun[[iSurv]])["maxyear"]-saveBlcks[,iCol]+2))
    rv      <- sample(c(T,F),nits,replace=T)
    strt    <- ifelse(rv==F,strstp,strstp-saveBlcks[,iCol]+1)
    stp     <- strt + saveBlcks[,iCol] - 1
    saveBlcksYrs[,iCol,"start"]  <- ifelse(rv==F,strt,stp)
    saveBlcksYrs[,iCol,"stop"]   <- ifelse(rv==F,stp,strt)
  }

  cat("Populating surv:",iSurv,"\n")
  for(iTer in 1:nits){
    blk <- which(is.na(saveBlcksYrs[iTer,,1])==F)
    idx <- ac(unlist(mapply(seq,from=saveBlcksYrs[iTer,blk,"start"],to=saveBlcksYrs[iTer,blk,"stop"])))
    #- Fill survey pattern with random draws of historic years
    if(iSurv == "SCAI")   iter(surv[1,      ac(yrs),iSurv],iTer)   <- iResids[,idx]
    if(iSurv == "IBTS0")  iter(surv[1,      ac(yrs),iSurv],iTer)   <- iResids[,idx]
    if(iSurv == "IBTS-Q1")iter(surv[1,      ac(yrs),iSurv],iTer)   <- iResids[,idx]
    if(iSurv == "HERAS")  iter(surv[ac(1:8),ac(yrs),iSurv],iTer)   <- iResids[,idx]
  }
}


  #- Recalculate survey index values (also historic)
SCAIfactor                                                                            <- c(exp(yearMeans(
                                                                                               log(window(quantSums(stock.n(NSH)* exp(-NSH@harvest*harvest.spwn(NSH)-m(NSH)*m.spwn(NSH)) * stock.wt(NSH) *mat(NSH)),1972,histMaxYr) *
                                                                                                   subset(catchabilities(NSH.sam),fleet=="SCAI")$value) -
                                                                                               subset(residuals(NSH.sam),fleet=="SCAI")$log.mdl)))

surveys[["SCAI"]]@index[,ac(range(surveys[["SCAI"]])["minyear"]:histMaxYr)]           <- sweep(quantSums(stock.n(stocks)    * exp(-stocks@harvest*harvest.spwn(stocks)-m(stocks)*m.spwn(stocks)) * stock.wt(stocks) *
                                                                                               mat(stocks))[,ac(range(surveys[["SCAI"]])["minyear"]:histMaxYr)],
                                                                                               c(2,6),  surveyK[["SCAI"]][,ac(range(surveys[["SCAI"]])["minyear"]:histMaxYr)],"^")            *
                                                                                               surveyQ[["SCAI"]][,    ac(range(surveys[["SCAI"]])   ["minyear"]:histMaxYr)] * 1/SCAIfactor    * surv[1,      ac(range(surveys[["SCAI"]])    ["minyear"]:histMaxYr),"SCAI"]

surveys[["IBTS0"]]@index[,ac(range(surveys[["IBTS0"]])["minyear"]:(histMaxYr+1))]     <- sweep((stock.n(stocks)["0",]       * exp(-stocks@harvest["0",] * mean(range(surveys[["IBTS0"]])  [c("startf","endf")]) -
                                                                                               m(stocks)["0",]  * mean(range(surveys[["IBTS0"]])   [c("startf","endf")])))[,ac(range(surveys[["IBTS0"]])["minyear"]:(histMaxYr+1))],
                                                                                               c(2,6),  surveyK[["IBTS0"]][,ac(range(surveys[["IBTS0"]])["minyear"]:(histMaxYr+1))],"^")      *
                                                                                               surveyQ[["IBTS0"]][,   ac(range(surveys[["IBTS0"]])  ["minyear"]:(histMaxYr+1))]               * surv[1,      ac(range(surveys[["IBTS0"]])   ["minyear"]:(histMaxYr+1)),"IBTS0"]

surveys[["IBTS-Q1"]]@index[,ac(range(surveys[["IBTS-Q1"]])["minyear"]:(histMaxYr+1))] <- sweep((stock.n(stocks)             * exp(-stocks@harvest * mean(range(surveys[["IBTS-Q1"]])[c("startf","endf")]) -
                                                                                               m(stocks)        * mean(range(surveys[["IBTS-Q1"]]) [c("startf","endf")])))[ac(1),ac(range(surveys[["IBTS-Q1"]])["minyear"]:(histMaxYr+1))],
                                                                                               c(1,2,6),surveyK[["IBTS-Q1"]][,ac(range(surveys[["IBTS-Q1"]])["minyear"]:(histMaxYr+1))],"^")  *
                                                                                               surveyQ[["IBTS-Q1"]][, ac(range(surveys[["IBTS-Q1"]])["minyear"]:(histMaxYr+1))]               * surv[ac(1),  ac(range(surveys[["IBTS-Q1"]]) ["minyear"]:(histMaxYr+1)),"IBTS-Q1"]

surveys[["HERAS"]]@index[,ac(range(surveys[["HERAS"]])["minyear"]:histMaxYr)]         <- sweep(setPlusGroup(stock.n(stocks) * exp(-stocks@harvest * mean(range(surveys[["HERAS"]])  [c("startf","endf")]) -
                                                                                               m(stocks)        * mean(range(surveys[["HERAS"]])   [c("startf","endf")])),8)[ac(1:8),ac(range(surveys[["HERAS"]])  ["minyear"]:histMaxYr)],
                                                                                               c(1,2,6),surveyK[["HERAS"]][,ac(range(surveys[["HERAS"]])    ["minyear"]:histMaxYr)],"^")      *
                                                                                               surveyQ[["HERAS"]][,   ac(range(surveys[["HERAS"]])  ["minyear"]:histMaxYr)]                   * surv[ac(1:8),ac(range(surveys[["HERAS"]])   ["minyear"]:histMaxYr),"HERAS"]
surveys[["HERAS"]]@index[ac(1),ac(range(surveys[["HERAS"]])["minyear"]:1996)]         <- -1


  #-------------------------------------------------------------------------------
  #- 2): Use stock object + error on realisations
  #-------------------------------------------------------------------------------

dmns              <- dimnames(NSH@catch.n)
dmns$year         <- dmns$year[1]:futureMaxYr
dmns$iter         <- 1:nits
ctch              <- FLQuant(NA,dimnames=dmns)

  #- Take blocks of residuals, sample blocks from 1-10 and add up till length of timeseries
yrs       <- range(NSH)["minyear"]:futureMaxYr
saveBlcks <- matrix(NA,nrow=nits,ncol=length(yrs))

full      <- F
rw        <- 1
cat("Sampling iters:","catch.n","\n")
while(full==F){
  tmpx <- 0
  while(tmpx<1) {
    set.seed(as.integer(Sys.time()))
    sam   <- sample(1:10,nits*1000,replace=T)
    samM  <- matrix(sam,nrow=nits,ncol=1000)
    idx   <- which(t(apply(samM,1,cumsum)) == length(yrs),arr.ind=T)
    tmpx <- length(idx)
  }
  if((rw+nrow(idx) - 1) > nrow(saveBlcks)){ stprow <- nrow(saveBlcks); full <- T} else { stprow <- rw+nrow(idx) -1}
  for(iRow in 1:(stprow-rw+1))
    if ((rw+iRow-1) <= nrow(saveBlcks) && iRow>0) saveBlcks[(rw+iRow-1),1:idx[iRow,2]] <- samM[idx[iRow,1],1:idx[iRow,2]]
  rw    <- stprow + 1
}

  #- Take the sampled blocks and assign years to it
saveBlcksYrs <- array(NA,dim=c(nits,length(yrs),2),dimnames=list(nits=1:nits,yrs=yrs,strtstp=c("start","stop")))
for(iCol in 1:ncol(saveBlcksYrs)){
  strstp  <- as.integer(runif(nits,range(NSH)["minyear"]+saveBlcks[,iCol]-1,range(NSH)["maxyear"]-saveBlcks[,iCol]+2))
  rv      <- sample(c(T,F),nits,replace=T)
  strt    <- ifelse(rv==F,strstp,strstp-saveBlcks[,iCol]+1)
  stp     <- strt + saveBlcks[,iCol] - 1
  saveBlcksYrs[,iCol,"start"]  <- ifelse(rv==F,strt,stp)
  saveBlcksYrs[,iCol,"stop"]   <- ifelse(rv==F,stp,strt)
}

  #- Substract and calculate residuals (non-standardized)
Resids                                                        <- subset(residuals(NSH.sam),fleet=="catch")
iResids                                                       <- FLQuant(NA,dimnames=c(dimnames(stocks@stock.n)[1:5],iter="1"))

for(i in 1:nrow(Resids))
  iResids[ac(Resids$age[i]),ac(Resids$year[i]),]              <- exp(Resids$log.obs[i] - Resids$log.mdl[i])

  #- Fill the object with the residuals
for(iTer in 1:nits){
  blk <- which(is.na(saveBlcksYrs[iTer,,1])==F)
  idx <- ac(unlist(mapply(seq,from=saveBlcksYrs[iTer,blk,"start"],to=saveBlcksYrs[iTer,blk,"stop"])))
    #- Fill survey pattern with random draws of historic years
  iter(ctch[,      ac(yrs),],iTer)   <- iResids[,idx]
}
  #- Because residuals are not estimated everywhere, some are NA, replace with 1
ctch@.Data[which(is.na(ctch))] <- 1

#- Add error to catch.n observations
stocks@catch.n    <- stocks@catch.n * ctch
stocks@landings.n <- stocks@catch.n

save.image(file.path(outPath,"catchSurveys.RData"))
  #-------------------------------------------------------------------------------
  #- 3): Perform starting point assessment + retrospective to measure assessment
  #      error
  #-------------------------------------------------------------------------------

dms   <- dims(stocks)
dmns  <- dimnames(stocks@stock.n)
resN  <- array(NA,dim=c(dms$age,dms$year,1,1,nyrs,nits),dimnames=list(age=dmns$age,year=dmns$year,unit="unique",season="all",area=projPeriod,iter=1:nits))
resF  <- array(NA,dim=c(dms$age,dms$year,1,1,nyrs,nits),dimnames=list(age=dmns$age,year=dmns$year,unit="unique",season="all",area=projPeriod,iter=1:nits))

#===============================================================================
# DO THIS ON NEMO
#===============================================================================

  #- Prepare objects for retro, no hess estimation
iter_retro                  <- FLSAMs()
base.assess                 <- NSH.sam
base.assess@control@nohess  <- T
ctrl                        <- NSH.ctrl
ctrl@nohess                 <- T

  #- Start max number of cores (12 - 1)
cl                          <- startCluster(11)
exportCluster(cl,lst=list("stocks","base.assess","surveys","ctrl","retroLB","retro","histMaxYr","NSH.tun"))
  #- Call retros
start.time                  <- Sys.time()
retros                      <- clusterApplyLB(cl,1:nits,fun=retroLB)
end.time                    <- Sys.time()

stopCluster(cl)
#===============================================================================
# END DO THIS ON NEMO
#===============================================================================

  #- Take blocks of 20 years, sample blocks from retro years and add up till length of timeseries
  #  So from 1957 - 2013 (see explanation below

#DDDMMMM change 11 to 27...

saveBlcksYrs <- array(NA,dim=c(nits,2),dimnames=list(nits=1:nits,strtstp=c("start","stop")))
strstp  <- as.integer(runif(nits,1957+20-1,histMaxYr-20+2))
rv      <- sample(c(T,F),nits,replace=T)
strt    <- ifelse(rv==F,strstp,strstp-20+1)
stp     <- strt + 20 - 1
saveBlcksYrs[,"start"]  <- ifelse(rv==F,strt,stp)
saveBlcksYrs[,"stop"]   <- ifelse(rv==F,stp,strt)
saveBlcksYrs  <- as.data.frame(saveBlcksYrs)
saveBlcksYrs  <- cbind(saveBlcksYrs,retroYr = as.integer(runif(nits,1,11)))
saveBlcksYrs  <- cbind(saveBlcksYrs,mult=rbinom(nits,1,0.5))
saveBlcksYrs$mult[which(saveBlcksYrs$mult==0)] <- -1

retros[[65]]  <- retros[[sample(c(1:nits)[-c(65,708,883)],1)]]
retros[[708]] <- retros[[sample(c(1:nits)[-c(65,708,883)],1)]]
retros[[883]] <- retros[[sample(c(1:nits)[-c(65,708,883)],1)]]
iter_retro    <- as(retros,"FLSAMs")

for(iRun in 1:nits){
  print(iRun)
  iter_ns    <- lapply(iter_retro[[iRun]],function(x){return(window(slot(x,"stock.n"),start=histMinYr,end=histMaxYr))})
  iter_fs    <- lapply(iter_retro[[iRun]],function(x){return(window(slot(x,"harvest"),start=histMinYr,end=histMaxYr))})
  iter_nsFLQ <- FLQuant(array(unlist(iter_ns),dim=c(dms$age,length(histPeriod),1,1,1,nyrs)),dimnames=list(age=dms$min:dms$max,year=histPeriod,unit="unique",season="all",area="unique",iter=1:nyrs))
  iter_fsFLQ <- FLQuant(array(unlist(iter_fs),dim=c(dms$age,length(histPeriod),1,1,1,nyrs)),dimnames=list(age=dms$min:dms$max,year=histPeriod,unit="unique",season="all",area="unique",iter=1:nyrs))
  for(i in 2:11){
    iter_nsFLQ[,ac(histMaxYr-i+2),,,,i] <- NA
    iter_fsFLQ[,ac(histMaxYr-i+2),,,,i] <- NA
  }

  #- Drop histMaxYr assessment because it is treated as 'the truth' so no error there
  iter_errorN <- exp(sweep(log(iter_nsFLQ[,,,,,2:11]),1:5,log(stocks@stock.n[,histPeriod,,,,iRun]),"-"))
  iter_errorF <- exp(sweep(log(iter_fsFLQ[,,,,,2:11]),1:5,log(stocks@harvest[,histPeriod,,,,iRun]),"-"))
  #- Align error terms to what they are in relation to terminal year in their own retro aspect
  #  So the 2001 assessment has terminal year 2010 but is lined up as histMaxYr now to make calcs easier
  for(i in 1:10){
    iter_errorN[,(1+i):67,,,,i] <- iter_errorN[,1:(67-i),,,,i]
    iter_errorF[,(1+i):67,,,,i] <- iter_errorF[,1:(67-i),,,,i]
  }

  #- Now fill resN and resF: resN only needs 20 years of data, for 11 years ahead
  #  So sample 11 times blocks of length 20 years from a yearspan of 1957:histMaxYr and
  #  allow these errors to be mirrored (multiply with -1) to assure retrospective bias
  #  can go in two ways.
  for(i in 1:27){
    refyr <- histMaxYr-1+i
    resN[,ac((refyr-19):(refyr)),,,i,iRun] <- exp(saveBlcksYrs[iRun,"mult"] * log(iter_errorN[,ac(saveBlcksYrs[iRun,"start"]:saveBlcksYrs[iRun,"stop"]),,,,saveBlcksYrs[iRun,"retroYr"]]))
    resF[,ac((refyr-19):(refyr)),,,i,iRun] <- exp(saveBlcksYrs[iRun,"mult"] * log(iter_errorF[,ac(saveBlcksYrs[iRun,"start"]:saveBlcksYrs[iRun,"stop"]),,,,saveBlcksYrs[iRun,"retroYr"]]))
  }
}
resN[is.na(resN)][] <- 1
resF[is.na(resF)][] <- 1
resN                <- FLQuant(resN,dimnames=list(age=dmns$age,year=dmns$year,unit="unique",season="all",area=projPeriod,iter=1:nits))
resF                <- FLQuant(resF,dimnames=list(age=dmns$age,year=dmns$year,unit="unique",season="all",area=projPeriod,iter=1:nits))
save(resN,      file=paste(outPath,"resNFinal.RData",sep=""))
save(resF,      file=paste(outPath,"resFFinal.RData",sep=""))
save(iter_retro,file=paste(outPath,"retros.RData"   ,sep=""))

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
wF                                                            <- FLQuant(array(t(mvrnorm(nits*nyrs,rep(0,length(ages)),cov(apply(log(NSH@harvest[,selPeriod,drop=T]),1,diff)))),
                                                                               dim=c(length(ages),nyrs,1,1,1,nits)),
                                                                         dimnames=dmns)
wF                                                            <- FLQuant(array(t(mvrnorm(nits*nyrs,rep(0,length(ages)),covmat)),
                                                                               dim=c(length(ages),nyrs,1,1,1,nits)),
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

