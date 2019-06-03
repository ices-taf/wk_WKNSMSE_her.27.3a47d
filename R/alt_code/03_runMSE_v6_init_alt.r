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

IAV       <- NULL
BB        <- NULL
HCR       <- 'A'
ftarget   <- 0.26
btrigger  <- 1.4e06

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
#path <- 'E:/git/wk_WKNSMSE_her.27.3a47d/R'
#path <- 'D:/Repository/NSAS_MSE/wk_WKNSMSE_her.27.3a47d/R/'
#path <- "/home/hintz001/wk_WKNSMSE_her.27.3a47d/R"
path <- '/home/berge057/ICES/wk_WKNSMSE_her.27.3a47d/R/'
#path <- 'E:/git/wk_WKNSMSE_her.27.3a47d/R'
assessment_name   <- "NSAS_WKNSMSE2018"
try(setwd(path),silent=TRUE)

# paths to different subfolders
dataPath      <- file.path(".","data/")
outPath       <- file.path(".","results/")
scriptPath    <- file.path(".","side_scripts/")
functionPath  <- file.path(".","functions/")

source(file.path(functionPath,"MSE_assessment.R"))
source(file.path(functionPath,"forecastScenarios.r"))
source(file.path(functionPath,"forecastFunctions.r"))


load('/home/berge057/ICES/wk_WKNSMSE_her.27.3a47d/R/results/stkAssessment2018_2000_20y_corrected.RData')

startYear <- iYr

start.time          <- Sys.time()
for (iYr in (startYear:max(projPeriod))){
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
  
  paramRec  <- params(biol.sr)
  ssbRec  <- drop(ssb(biol[,ac(iYr-1)]))
  
  # Ricker
  recruitBio[itersRI] <- paramRec['a',itersRI]*ssbRec[itersRI]*exp(-paramRec['b',itersRI]*ssbRec[itersRI])
  # Segmented regression
  idxSSR1 <- which(ssbRec[itersSR] <= paramRec['b',itersSR])
  # if loop to cover the case of idxSSR1 being empty
  if(length(idxSSR1)!=0){
    recruitBio[itersSR[idxSSR1]]    <- paramRec['a',itersSR[idxSSR1]]*ssbRec[itersSR[idxSSR1]] # SSB < b (slope)
    recruitBio[itersSR[-idxSSR1]]   <- paramRec['a',itersSR[-idxSSR1]]*paramRec['b',itersSR[-idxSSR1]] # SSB > b (plateau)
  }else{
    recruitBio[itersSR[idxSSR1]]    <- paramRec['a',itersSR]*paramRec['b',itersSR]
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
  require(doParallel); ncores <- detectCores()-1; ncores <- ifelse(nits<ncores,nits,ncores);cl <- makeCluster(ncores); registerDoParallel(cl)
  dat     <- as.data.frame(stkAssessment@m)
  datS    <- split(dat,as.factor(paste(dat$age,dat$iter)))
  res     <- foreach(i = 1:length(datS)) %dopar% fitted(loess(data ~ year,data=datS[[i]],span=0.5))
  stkAssessment@m <- FLQuant(c(aperm(array(unlist(res),dim=c(length(1947:TaY),nits,nAges)),c(3,1,2))),dimnames=dimnames(stkAssessment@m))
  stopCluster(cl); detach("package:doParallel",unload=TRUE); detach("package:foreach",unload=TRUE); detach("package:iterators",unload=TRUE)
  
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
                          escapeRuns)
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
  
  projNSAS                  <- projectNSH(stkAssessment,
                                          fishery,
                                          iYr,
                                          TAC,
                                          histMaxYr,
                                          referencePoints,
                                          managementRule)
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
  
  require(doParallel); ncores <- detectCores()-1; ncores <- ifelse(nits<ncores,nits,ncores);cl <- makeCluster(ncores); clusterEvalQ(cl,library(FLCore)); clusterEvalQ(cl,library(minpack.lm)); registerDoParallel(cl)
  multso <- do.call(rbind,foreach(idxIter = 1:nits) %dopar% nls.lm(par=runif(4),
                                                                   fn=TAC2sel,
                                                                   iYr=ImY,
                                                                   iBiol=biol[,ImY],
                                                                   iFishery=fishery[,ImY],
                                                                   iTAC=CATCH[,ImY],
                                                                   catchVar=catchVar,
                                                                   TAC_var=TAC_var,
                                                                   iTer=idxIter,
                                                                   control=nls.lm.control(maxiter=1000),
                                                                   lower=rep(1e-8,4))$par)
  stopCluster(cl); detach("package:doParallel",unload=TRUE); detach("package:foreach",unload=TRUE); detach("package:iterators",unload=TRUE)
  
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


save.image(file=paste(outPath,'stkAssessment2018_2000_40y',".RData",sep=""))

#save.image(file=paste(outPath,'stkAssessment2018_2000_20y_corrected',".RData",sep=""))

#save(stkAssessment.init,file=paste(outPath,runName,".RData",sep=""))


