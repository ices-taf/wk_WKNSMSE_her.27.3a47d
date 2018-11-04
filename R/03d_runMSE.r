#-------------------------------------------------------------------------------
#
# Code for HERTAC to simulate different scenarios of population
#  structure and see if it can be managed sustainably
#
# Code by: Niels Hintzen
#
# To be used with R3.0.x on both Windows as Linux
#-------------------------------------------------------------------------------

  rm(list=ls())
  library(FLCore)
  library(FLFleet)
  library(minpack.lm)

  #- Set paths
  codePath        <- "I:/WKHerTAC/NSAS/R/"
  dataPath        <- "I:/WKHerTAC/NSAS/Data/"
  outPath         <- "I:/WKHerTAC/NSAS/Results/"
  NSASPath        <- "I:/WKHerTAC/NSAS/Results/"
  WBSSPath        <- "I:/WKHerTAC/WBSS/Results/"
  combPath        <- "I:/WKHerTAC/NSASWBSS/Results/"
  FLMETApath      <- "I:/WKHerTAC/FLMeta/"

  #- Paths for Nemo
  if((substr(R.Version()$os,1,3)== "lin" & length(dir("/media"))>0) | (substr(R.Version()$os,1,3)=="min" & length(dir("/media"))>0)){
    codePath      <- sub("I:/","/media/w/",codePath)
    dataPath      <- sub("I:/","/media/w/",dataPath)
    outPath       <- sub("I:/","/media/w/",outPath)
    NSASPath      <- sub("I:/","/media/w/",NSASPath)
    WBSSPath      <- sub("I:/","/media/w/",WBSSPath)
    combPath      <- sub("I:/","/media/w/",combPath)
    FLMETApath    <- sub("I:/","/media/w/",FLMETApath)
  }

  #- Source the FLMeta package code
  source(file.path(FLMETApath,"FLMETA.r"))
  source(paste(codePath,"04_forecastScenarios.r", sep="")) #NSAS

  #- Settings
  histMinYr       <- 1991
  histMaxYr       <- 2013
  nyrs            <- 27
  futureMaxYr     <- histMaxYr + nyrs
  histPeriod      <- ac(histMinYr:histMaxYr)
  projPeriod      <- ac((histMaxYr+1):futureMaxYr)
  nits            <- 10
  settings        <- list(histMinYr=histMinYr,histMaxYr=histMaxYr,futureMaxYr=futureMaxYr,
                          histPeriod=histPeriod,projPeriod=projPeriod,nyrs=nyrs,nits=nits)
  datum           <- date()

#--------------------------------------------------------------
#- 0) Meta match case scenario: 2 biological unit, 5 fisheries, 2 stock assessments
#
# 5 fishing fleets
# 2 populations, one in North Sea, one in Western Baltic
#--------------------------------------------------------------

  #- Load data from north-stock assessments (WP4north)
  biolsL          <- list();
  load(file=file.path(NSASPath,"biol.RData"));         biolsL[["NS"]] <- biol
  load(file=file.path(WBSSPath,"biol.RData"));         biolsL[["WB"]] <- biol
  biols           <- biolsL; rm(biolsL)
  stocksL         <- list();
  load(file=file.path(NSASPath,"stocks.RData"));       stocksL[["NS"]] <- stocks
  load(file=file.path(WBSSPath,"stocks.RData"));       stocksL[["WB"]] <- stocks
  stocks          <- stocksL; rm(stocksL)
  fisheriesL      <- list();
  load(file=file.path(NSASPath,"fishery.RData"));      fisheriesL[["NS"]] <- fishery
  load(file=file.path(WBSSPath,"fishery.RData"));      fisheriesL[["WB"]] <- fishery
  fisheries       <- fisheriesL; rm(fisheriesL)
  devNs           <- list();
  load(file=file.path(NSASPath,"resNFinal.RData"));    devNs[["NS"]] <- resN
  load(file=file.path(WBSSPath,"resNFinal.RData"));    devNs[["WB"]] <- resN
  devFs           <- list();
  load(file=file.path(NSASPath,"resFFinal.RData"));    devFs[["NS"]] <- resF
  load(file=file.path(WBSSPath,"resNFinal.RData"));    devFs[["WB"]] <- resF
  propNs          <- list();
  load(file=file.path(NSASPath,"propN.RData"));        propNs[["NS"]] <- propN
  load(file=file.path(WBSSPath,"propN.RData"));        propNs[["WB"]] <- propN
  catchResids     <- list();
  load(file=file.path(NSASPath,"ctch.RData"));         catchResids[["NS"]] <- ctch
  load(file=file.path(WBSSPath,"ctch.RData"));         catchResids[["WB"]] <- ctch

  #- Fill the fisheries object with the right data and merge the C & D fleet together
  fisheries[[1]]                <- window(fisheries[[1]][,,,,,1:nits],start=histMinYr,end=futureMaxYr)
  fisheries[[2]]                <- window(fisheries[[2]][,,,,,1:nits],start=histMinYr,end=futureMaxYr)
  dmns                          <- dimnames(fisheries[[1]]@landings.n)
  dmns$unit                     <- c("NS","WB")
  dmns$area                     <- c("A","B","C","D","F")
  fishery                       <- FLCatch(landings.n=FLQuant(NA,dimnames=dmns))
  fishery@landings.wt[,,"NS",,"A"]  <- fisheries[["NS"]]@landings.wt[,,"A"]
  fishery@landings.wt[,,"NS",,"B"]  <- fisheries[["NS"]]@landings.wt[,,"B"]
  fishery@landings.wt[,,"NS",,"C"]  <- fisheries[["NS"]]@landings.wt[,,"C"]
  fishery@landings.wt[,,"NS",,"D"]  <- fisheries[["NS"]]@landings.wt[,,"D"]
  fishery@landings.wt[,,"WB",,"A"]  <- fisheries[["WB"]]@landings.wt[,,"A"]
  fishery@landings.wt[,,"WB",,"C"]  <- fisheries[["WB"]]@landings.wt[,,"C"]
  fishery@landings.wt[,,"WB",,"D"]  <- fisheries[["WB"]]@landings.wt[,,"D"]
  fishery@landings.wt[,,"WB",,"F"]  <- fisheries[["WB"]]@landings.wt[,,"F"]
  units(fishery)[1:10]          <- as.list(c("tonnes","thousands","kg",NA,"tonnes","thousands","kg",NA,NA,"euro"))
  range(fishery)["plusgroup"]   <- range(fishery)["max"]
  
  #--------------------------------------------------------------
  #- Setup biol (2 populations)
  #--------------------------------------------------------------

  #- First NSAS
  NSASorig                      <- window(biols[[1]][,,,,,1:nits],start=histMinYr,end=futureMaxYr)
  biol                          <- NSASorig
  dmns                          <- dimnames(biol@n); dmns$unit <- "NS"
  biol                          <- setDimnames(biol,dmns)
  bunits                        <- c("NS","WB")
  biol                          <- expand(biol,unit=bunits)

  #- Then add the south
  WBSSorig                      <- window(biols[["WB"]][,,,,,1:nits],start=histMinYr,end=futureMaxYr)
  biols[["WB"]]                 <- WBSSorig
  dmns                          <- dimnames(biols[["WB"]]@m)
  dmns$unit                     <- "WB"
  biols[["WB"]]                 <- setDimnames(biols[["WB"]],dmns)
  for(iSlot in slotNames(biol)[1:5])
    slot(biol,iSlot)[,,"WB"]    <- slot(window(biols[[2]],start=histMinYr,end=futureMaxYr),iSlot)
  units(biol)[1:5]              <- as.list(c("thousands",NA,"kg",NA,NA))
  biol                          <- window(biol,end=futureMaxYr)
  dmns                          <- dimnames(biol@n); dmns$area <- "unique"
  biol                          <- setDimnames(biol,dmns)

  #--------------------------------------------------------------
  #- Stocks setup & catch residuals & devN & devF
  #--------------------------------------------------------------
  stocks[["NS"]]                <- window(stocks[["NS"]][,,,,,1:nits],start=histMinYr,end=futureMaxYr)
  stocks[["WB"]]                <- window(stocks[["WB"]][,,,,,1:nits],start=histMinYr,end=futureMaxYr)
  dmns                          <- dimnames(stocks[[1]]@stock.n)
  dmns$unit                     <- c("NS","WB")
  stock                         <- FLStock(stock.n=FLQuant(NA,dimnames=dmns))
  stock[,,"NS"]                 <- stocks[["NS"]]
  stock[,,"WB"]                 <- stocks[["WB"]]

  catchResids[["NS"]]           <- window(catchResids[["NS"]][,,,,,1:nits],start=histMinYr,end=futureMaxYr)
  catchResids[["WB"]]           <- window(catchResids[["WB"]][,,,,,1:nits],start=histMinYr,end=futureMaxYr)
  dmns                          <- dimnames(catchResids[[1]])
  dmns$unit                     <- c("NS","WB")
  catchResid                    <- FLQuant(NA,dimnames=dmns)
  catchResid[,,"NS"]            <- catchResids[["NS"]]
  catchResid[,,"WB"]            <- catchResids[["WB"]]

  devNs[["NS"]]                 <- window(devNs[["NS"]][,,,,,1:nits],start=histMinYr,end=futureMaxYr)
  devNs[["WB"]]                 <- window(devNs[["WB"]][,,,,,1:nits],start=histMinYr,end=futureMaxYr)
  dmns                          <- dimnames(devNs[[1]])
  dmns$unit                     <- c("NS","WB")
  devN                          <- FLQuant(NA,dimnames=dmns)
  devN[,,"NS"]                  <- devNs[["NS"]]
  devN[,,"WB"]                  <- devNs[["WB"]]
  
  devFs[["NS"]]                 <- window(devFs[["NS"]][,,,,,1:nits],start=histMinYr,end=futureMaxYr)
  devFs[["WB"]]                 <- window(devFs[["WB"]][,,,,,1:nits],start=histMinYr,end=futureMaxYr)
  dmns                          <- dimnames(devFs[[1]])
  dmns$unit                     <- c("NS","WB")
  devF                          <- FLQuant(NA,dimnames=dmns)
  devF[,,"NS"]                  <- devFs[["NS"]]
  devF[,,"WB"]                  <- devFs[["WB"]]

  #--------------------------------------------------------------
  #- Combine catches into one object
  #--------------------------------------------------------------

  #- Sync catches between fishery and original object
  fishery@landings.n[,ac(histPeriod),"NS",,"A"]     <- fisheries[["NS"]]@landings.n[,ac(histPeriod),"A"]
  fishery@landings.n[,ac(histPeriod),"NS",,"B"]     <- fisheries[["NS"]]@landings.n[,ac(histPeriod),"B"]
  fishery@landings.n[,ac(histPeriod),"NS",,"C"]     <- fisheries[["NS"]]@landings.n[,ac(histPeriod),"C"]
  fishery@landings.n[,ac(histPeriod),"NS",,"D"]     <- fisheries[["NS"]]@landings.n[,ac(histPeriod),"D"]
  fishery@landings.n[,ac(histPeriod),"WB",,"A"]     <- fisheries[["WB"]]@landings.n[,ac(histPeriod),"A"]
  fishery@landings.n[,ac(histPeriod),"WB",,"C"]     <- fisheries[["WB"]]@landings.n[,ac(histPeriod),"C"]
  fishery@landings.n[,ac(histPeriod),"WB",,"D"]     <- fisheries[["WB"]]@landings.n[,ac(histPeriod),"D"]
  fishery@landings.n[,ac(histPeriod),"WB",,"F"]     <- fisheries[["WB"]]@landings.n[,ac(histPeriod),"F"]

  fishery@landings.sel[]                            <- 0
  for(iUnit in c("NS","WB")){
    for(iFleet in c("A","B","C","D","F")){
      if(iFleet %in% dimnames(propNs[[iUnit]])$unit){
        fishery@landings.sel[,ac(histPeriod),iUnit,,iFleet] <- sweep(harvest(stocks[[iUnit]])[,ac(histMinYr:histMaxYr)],c(1,3:5),propNs[[iUnit]][,,iFleet],"*")
      }
    }
  }
  fishery@landings.sel[,ac(projPeriod),"NS",,"A"]   <- fisheries[["NS"]]@landings.sel[,ac(projPeriod),"A"]
  fishery@landings.sel[,ac(projPeriod),"NS",,"B"]   <- fisheries[["NS"]]@landings.sel[,ac(projPeriod),"B"]
  fishery@landings.sel[,ac(projPeriod),"NS",,"C"]   <- fisheries[["NS"]]@landings.sel[,ac(projPeriod),"C"]
  fishery@landings.sel[,ac(projPeriod),"NS",,"D"]   <- fisheries[["NS"]]@landings.sel[,ac(projPeriod),"D"]
  fishery@landings.sel[,ac(projPeriod),"WB",,"A"]   <- fisheries[["WB"]]@landings.sel[,ac(projPeriod),"A"]
  fishery@landings.sel[,ac(projPeriod),"WB",,"C"]   <- fisheries[["WB"]]@landings.sel[,ac(projPeriod),"C"]
  fishery@landings.sel[,ac(projPeriod),"WB",,"D"]   <- fisheries[["WB"]]@landings.sel[,ac(projPeriod),"D"]
  fishery@landings.sel[,ac(projPeriod),"WB",,"F"]   <- fisheries[["WB"]]@landings.sel[,ac(projPeriod),"F"]

  #- Set all units correct
  units(fishery)[1:10]  <- as.list(c("tonnes","thousands","kg",NA,"tonnes","thousands","kg",NA,NA,"euro"))
  units(biol)[c(1,3)]   <- as.list(c("thousands","kg"))

  #- Source recruitment scenario functions
  source(file.path(codePath,"functions.r"))

  #- Save starting conditions
  save.image(file.path(combPath,"startingConditions.RData"))
  load(file.path(combPath,"startingConditions.RData"))
  source(file.path(FLMETApath,"FLMETA.r"))
#------------------------------------------------------------------------------#
# 1) Define the HCR & historic TACs
#------------------------------------------------------------------------------#

  mpPoints                              <- list("NS"=FIAVBB$opt1,"WB"=NA)

  maxEff                                <- 1000
  TAC                                   <- FLQuant(NA,dimnames=list(age="all",year=histMinYr:(futureMaxYr+3),unit=c("NS","WB"),season="all",area=c("A","B","C","D","F"),iter=1:nits))
  TAC[,ac(2001:2014),"NS",,"A"]         <- c(265,265,400,460,535,455,341,201,171,164,200,405,478,470) *1000
  TAC[,ac(2001:2014),"NS",,"B"]         <- c(36,36,52,38,50,42,32,19,16,14,16,18,14,13)               *1000
  TAC[,ac(2001:2013),"NS",,"C"]         <- c(34,17,24.1,13.4,22.9,11.6,16.4,9.2,5.1,12.0,6.6,7.8,11.8)*1000
  TAC[,ac(2001:2013),"NS",,"D"]         <- c(12,9,8.4,10.8,9.0,3.4,3.4,3.7,1.5,1.8,1.8,4.4,1.6)       *1000
  TAC[,ac(2006:2014),"WB",,"C"]         <- c(54.5,19.5,6.7,10.5,11.2,14.2,24.1,29.2,27)               *1000
  TAC[,ac(2006:2014),"WB",,"D"]         <- NA
  TAC[,ac(2006:2014),"WB",,"F"]         <- c(47.5,49.5,45,27.2,22.7,15.8,20.9,25.8,19.8) *1000
  TACusage                              <- FLQuant(array(1,dim=c(1,length(histMinYr:futureMaxYr)+3,2,1,5,nits)),dimnames=dimnames(TAC[,ac(histMinYr:(futureMaxYr+3))]))
  HCRTAC                                <- TAC; HCRTAC[] <- NA; SSB <- HCRTAC[,,,,1]; HCRSSB <- SSB

#------------------------------------------------------------------------------#
# 2) Start running the MSE
#------------------------------------------------------------------------------#

  start.time          <- Sys.time()
  bunit               <- dimnames(biol@n)$unit
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
      #- Survivors
    survivors   <- n(biol)[,ac(iYr-1)] * exp(-z)
    n(biol)[ac((range(biol,"min")+1):range(biol,"max")),ac(iYr),,] <- survivors[-dim(survivors)[1],,,,,]@.Data

      #- Plusgroup
    if (!is.na(range(biol,"plusgroup"))){
      n(biol)[ac(range(biol,"max")),ac(iYr),] <- n(biol)[ac(range(biol,"max")),ac(iYr),] + survivors[ac(range(biol,"max"))]
    }

    cat("\n Finished biology \n")
    cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))

    #- Update fishery to year iYr-1
    landings.n(fishery)[,ac(iYr-1)]     <- sweep(sweep(landings.sel(fishery)[,ac(iYr-1),,,,],c(1:4,6),z,"/"),c(1:4,6),n(biol)[,ac(iYr-1)]*(1-exp(-z)),"*")

    #- Create stock object for assessment
    yrmin1      <- iYr -1
    TaY         <- yrmin1               #Terminal assessment year
    ImY         <- TaY+1                #Intermediate Year
    FcY         <- TaY+2                #Forecast year

    idxyrmin1   <- which(dimnames(biol@n)$year == yrmin1)
    tmp_biol    <- biol[,1:idxyrmin1]   #Same but faster as window(biol,histMinYr,yrmin1)
    tmp_fishery <- fishery[,1:idxyrmin1]#Same but faster as window(fishery,histMinYr,yrmin1)
    tmp_stocks  <- stock[,1:idxyrmin1] #Same but faster as window(stocks,histMinYr,yrmin1)

    #- Update stocks to year iYr -1
    tmp_stocks  <- updateStocks(tmp_stocks,tmp_fishery,yrmin1,tmp_biol,catchResid)

    #- Overwrite results from update to stock again (but not for 2011, as that result is already known)
    if(iYr > an(projPeriod[1]))
      stock     <- tmp2stocks(stock,tmp_stocks,TaY)

    cat("\n Finished update \n")
    cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))

    #-Do the assessment
    stock[,ac(histMinYr:TaY)]@stock.n   <- biol@n[,ac(histMinYr:TaY)]        * devN[,ac(histMinYr:TaY),,,ac(iYr),]
    stock[,ac(histMinYr:TaY)]@harvest   <- areaSums(landings.sel(fishery)[,ac(histMinYr:TaY)]) * devF[,ac(histMinYr:TaY),,,ac(iYr),]
    stock@stock[,ac(histMinYr:TaY)]     <- computeStock(stock[,ac(histMinYr:TaY)])
    survivors[ac(0),]                   <- biol@n[ac(0),ac(iYr)]            * devN[ac(0),ac(iYr),,,ac(iYr),]

    #Set plusgroup at 7 (which is true plusgroup - recruitment)
    survivors[-1,]    <- FLQuant(setPlusGroup(stock[,ac(TaY)]@stock.n * exp(-stock[,ac(TaY)]@harvest-stock[,ac(TaY)]@m),7)@.Data,
                                 dimnames=list(age=dimnames(stock@stock.n)$age[-1],year=ac(TaY),unit=dimnames(stock@stock.n)$unit,
                                               season=dimnames(stock@stock.n)$season,area=dimnames(stock@stock.n)$area,iter=1:nits))

    cat("\n Finished stock assessment \n")
    cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))

    #- Project 4-fleet setup

    projWBSS                  <- projectWBSS()
    projNSAS                  <- projectNSH(stock[,1:idxyrmin1,"NS"],survivors[,,"NS"],tmp_fishery[,,"NS",,c("A","B","C","D")],iYr,TAC[,,"NS",,c("A","B","C","D")],mpPoints[["NS"]]$scen,NULL,histMaxYr,mpPoints[["NS"]])

    TAC[,   ac(FcY),"NS",,c("A","B","C","D")]          <- projNSAS[["TAC"]]
    HCRTAC[,ac(FcY),"NS",,c("A","B","C","D")]          <- projNSAS[["HCRTAC"]]
    HCRSSB[,ac(FcY),"NS",,c("A","B","C","D")]          <- projNSAS[["SSB"]][["HCRSSB"]][,ac(FcY)]
    SSB[,   ac(FcY),"NS",,c("A","B","C","D")]          <- projNSAS[["SSB"]][["SSB"]][,ac(FcY)]

    TAC[,   ac(FcY),"WB",,c("A","C","D","F")]          <- projWBAS[["TAC"]]
    HCRTAC[,ac(FcY),"WB",,c("A","C","D","F")]          <- projWBAS[["HCRTAC"]]
    HCRSSB[,ac(FcY),"WB",,c("A","C","D","F")]          <- projWBAS[["SSB"]][["HCRSSB"]][,ac(FcY)]
    SSB[,   ac(FcY),"WB",,c("A","C","D","F")]          <- projWBAS[["SSB"]][["SSB"]][,ac(FcY)]
    
    #- 50% transfer of CD-Fleet TAC to A fleet
    TAC[,ac(FcY),"NS",,c("A")] <- TAC[,ac(FcY),"NS",,c("A")] + 0.5 * areaSums(TAC[,ac(FcY),"NS",,c("C","D")])

    cat("\n Finished forecast \n")
    cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))

    #-Calculate effort accordingly (assuming constant catchability)
    landings.sel(fishery[,ac(ImY),"NS",,c("A","B","C","D")]) <- sweep(landings.sel(fishery)[,ac(ImY)],c(3,6),pmin(maxEff,f31tF(TAC*TACusage,biol,ImY,fishery)),"*")
    landings.sel(fishery[,ac(ImY),"WB",,c("A","C","D","F")]) <- sweep(landings.sel(fishery)[,ac(ImY)],c(3,6),pmin(maxEff,f31tF(TAC*TACusage,biol,ImY,fishery)),"*")
    cat("\n Finished effort calc \n")
    cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
  }
  save.image(file=paste(outPath,runName,"_",settings$RecRegime,".RData",sep=""))


