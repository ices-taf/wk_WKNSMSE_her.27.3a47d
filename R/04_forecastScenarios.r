#-------------------------------------------------------------------------------
# WKHERMPII
#
# Author: Niels Hintzen
#         IMARES, The Netherland
#
# Performs an MSE of North Sea Herring under different TAC scenario's
#
# Date: 03-Oct-2011
#
# Build for R2.8.1
#-------------------------------------------------------------------------------
#projectNSH(window(stocks,start=histMinYr,end=TaY),survivors,tmp_fishery,iYr,TAC,scen,iFailRuns=NULL,histMaxYr,mpPoints)
projectNSH <- function(iStocks,iSurvivors,iFishery,iYr,iTAC,iScenario=scen,iFailRuns=NULL,iHistMaxYr,mpPoints){
  require(minpack.lm)

  lin     <- substr(R.Version()$os,1,3)== "lin"
  stk     <- iStocks
  stk.sur <- iSurvivors
  #===============================================================================
  # Setup control file
  #===============================================================================

  DtY         <- ac(range(stk)["maxyear"]) #Data year
  ImY         <- ac(an(DtY)+1) #Intermediate year
  FcY         <- ac(an(DtY)+2) #Forecast year
  CtY         <- ac(an(DtY)+3) #Continuation year
  CtY1        <- ac(an(DtY)+4)
  FuY         <- c(ImY,FcY,CtY,CtY1)#Future years
  pyears      <- an(DtY) - iHistMaxYr #Years away from original historic max year
  
  #- We substract pyears from the starting point, because with every year forward, we move the year ranges one forward too.
  RECS        <- FLQuants("ImY"=iSurvivors[1,],"FcY"=exp(apply(log(rec(stk)[,ac((range(stk)["maxyear"]-8-pyears):(range(stk)["maxyear"]))]),3:6,mean,na.rm=T)),
                                               "CtY"=exp(apply(log(rec(stk)[,ac((range(stk)["maxyear"]-8-pyears):(range(stk)["maxyear"]))]),3:6,mean,na.rm=T)))

  yrs1        <- list("m.spwn","harvest.spwn","stock.wt")
  yrs3        <- list("mat")
  yrs5        <- list("m")

  dsc         <- "North Sea Herring"
  nam         <- "NSH"
  dms         <- dimnames(stk@m)
  dms$year    <- c(rev(rev(dms$year)[1:3]),ImY,FcY,CtY,CtY1)
  dms$unit    <- dimnames(iFishery@landings.n)$area

  f01         <- ac(0:1)
  f26         <- ac(2:6)

  #===============================================================================
  # Setup stock file
  #===============================================================================

  stf         <- FLStock(name=nam,desc=dsc,m=FLQuant(NA,dimnames=dms))
  for(i in dms$unit) stf[,,i] <- window(stk,start=an(dms$year)[1],end=rev(an(dms$year))[1])
  units(stf)  <- units(stk)
  # Fill slots that are the same for all fleets
  for(i in c(unlist(yrs1),unlist(yrs3),unlist(yrs5))){
    if(i %in% unlist(yrs1)){ for(j in dms$unit){ slot(stf,i)[,FuY,j] <- slot(stk,i)[,DtY]}}
    if(i %in% unlist(yrs3)){ for(j in dms$unit){ slot(stf,i)[,FuY,j] <- yearMeans(slot(stk,i)[,ac((an(DtY)-2):an(DtY))])}}
    if(i %in% unlist(yrs5)){ for(j in dms$unit){ slot(stf,i)[,FuY,j] <- yearMeans(slot(stk,i)[,ac((an(DtY)-4):an(DtY))])}}
  }

  # Fill slots that are unique for the fleets
  for(iFuY in FuY){
    stf@harvest[,iFuY]       <- sweep(sweep(iFishery@landings.n[,DtY],c(1,6),areaSums(iFishery@landings.n[,DtY]),"/"),c(1,6),stk@harvest[,DtY],"*")
    stf@catch.wt[,iFuY]      <- iFishery@landings.wt[,ac(DtY)]
    stf@landings.wt[,iFuY]   <- iFishery@landings.wt[,ac(DtY)]
  }

  # Fill slots that have no meaning for NSAS
  stf@discards.n[]        <- 0
  stf@discards[]          <- 0
  stf@discards.wt[]       <- 0

  #===============================================================================
  # Intermediate year
  #===============================================================================

  for(i in dms$unit)        stf@stock.n[,ImY,i]     <- stk.sur
  stf@harvest[,ImY]         <- fleet.harvestFF(stk=stf,iYr=ImY,TACS=iTAC[,ImY])
  for(i in dms$unit){
    stf@catch.n[,ImY,i]     <- stf@stock.n[,ImY,i]*(1-exp(-unitSums(stf@harvest[,ImY])-stf@m[,ImY,i]))*(stf@harvest[,ImY,i]/(unitSums(stf@harvest[,ImY])+stf@m[,ImY,i]))
    stf@catch[,ImY,i]       <- computeCatch(stf[,ImY,i])
    stf@landings.n[,ImY,i]  <- stf@catch.n[,ImY,i]
    stf@landings[,ImY,i]    <- computeLandings(stf[,ImY,i])
  }

  #===============================================================================
  # Forecast year
  #===============================================================================

  for(i in dms$unit) stf@stock.n[1,FcY,i]                    <- RECS$FcY
  for(i in dms$unit) stf@stock.n[2:(dims(stf)$age-1),FcY,i]  <- (stf@stock.n[,ImY,1]*exp(-unitSums(stf@harvest[,ImY])-stf@m[,ImY,1]))[ac(range(stf)["min"]:(range(stf)["max"]-2)),]
  for(i in dms$unit) stf@stock.n[dims(stf)$age,FcY,i]        <- apply((stf@stock.n[,ImY,1]*exp(-unitSums(stf@harvest[,ImY])-stf@m[,ImY,1]))[ac((range(stf)["max"]-1):range(stf)["max"]),],2:6,sum,na.rm=T)

  ###--- Management options ---###

    #- First calculate the TAC according to the HCR (not the LTMP)
  res <- matrix(NA,nrow=dims(stf)$unit,ncol=dims(stf)$iter,dimnames=list(dimnames(stf@stock.n)$unit,dimnames(stf@stock.n)$iter))
  if(is.null(iFailRuns) == T){
    for(iTer in 1:dims(stf)$iter){
      if(lin)  res[,iTer]              <- nls.lm(par=rep(1,dims(stf)$unit),find.FABF,stk=stf[,FcY,,,,iTer],f01=f01,f26=f26,TACS=iTAC[,FcY,,,,iTer],
                                                 mpPoints=mpPoints,jac=NULL,lower=rep(0,dims(stf)$unit),upper=rep(1e5,dims(stf)$unit))$par
      if(!lin) res[,iTer]              <- nls.lm(par=rep(1,dims(stf)$unit),find.FABF,stk=stf[,FcY,,,,iTer],f01=f01,f26=f26,TACS=iTAC[,FcY,,,,iTer],
                                                 mpPoints=mpPoints,jac=NULL)$par
    }

  } else {
      for(iTer in (1:dims(stf)$iter)[-iFailRuns]){
        if(lin)  res[,iTer]            <- nls.lm(par=rep(1,dims(stf)$unit),find.FABF,stk=iter(stf[,FcY],iTer),f01=f01,f26=f26,TACS=iter(iTAC[,FcY],iTer),
                                                 mpPoints=mpPoints,jac=NULL,lower=rep(0,dims(stf)$unit),upper=rep(1e5,dims(stf)$unit))$par
        if(!lin) res[,iTer]            <- nls.lm(par=rep(1,dims(stf)$unit),find.FABF,stk=iter(stf[,FcY],iTer),f01=f01,f26=f26,TACS=iter(iTAC[,FcY],iTer),
                                                 mpPoints=mpPoints,jac=NULL)$par
      }
    }
    
  stf@harvest[,FcY]         <- sweep(stf@harvest[,FcY],c(3,6),res,"*")
  for(i in dms$unit){
    stf@catch.n[,FcY,i]     <- stf@stock.n[,FcY,i]*(1-exp(-unitSums(stf@harvest[,FcY])-stf@m[,FcY,i]))*(stf@harvest[,FcY,i]/(unitSums(stf@harvest[,FcY])+stf@m[,FcY,i]))
    stf@catch[,FcY,i]       <- computeCatch(stf[,FcY,i])
    stf@landings.n[,FcY,i]  <- stf@catch.n[,FcY,i]
    stf@landings[,FcY,i]    <- computeLandings(stf[,FcY,i])
  }
  HCRTAC                    <- round(harvestCatch(stf,FcY),0)
  SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
  HCRSSB                    <- SSB
  F26                       <- apply(unitSums(stf@harvest[f26,FcY]),2:6,mean)
  F01                       <- apply(unitSums(stf@harvest[f01,FcY]),2:6,mean)
  print(paste(round(c(apply(F26,1,mean)),5),round(c(apply(F01,1,mean)),5)))


  #-----------------------------------------------------------------------------
  #- Scenario no IAV
  #-----------------------------------------------------------------------------
  if(iScenario == "noIAV"){
    #- No change to TAC calculation
    iTAC[,FcY]              <- HCRTAC
  }

  #-----------------------------------------------------------------------------
  #- Scenario current Long Term Management Plan
  #-----------------------------------------------------------------------------
  if(iScenario == "LTMP"){
    iTAC[,FcY]                <- HCRTAC
    bidx                      <- which(HCRTAC[,,c("A")] > 1.15*iTAC[,ImY,,,c("A")]) #bigger than  1.15x ImY TAC
    sidx                      <- which(HCRTAC[,,c("A")] < 0.85*iTAC[,ImY,,,c("A")]) #smaller than 0.85x ImY TAC
    if(length(bidx)>0) iTAC[,FcY,,,c("A"),bidx] <- 1.15*iTAC[,ImY,,,c("A"),bidx]
    if(length(sidx)>0) iTAC[,FcY,,,c("A"),sidx] <- 0.85*iTAC[,ImY,,,c("A"),sidx]

    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY])
    SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
                                           exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
    idxssb                    <- which(SSB < ifelse(mpPoints$stabilityBreak == "Blim",mpPoints$Blim,mpPoints$Btrigger))
    if(length(idxssb)>0) iTAC[,FcY,,,c("A"),idxssb] <- HCRTAC[,,"A",,,idxssb]
  }

  #-----------------------------------------------------------------------------
  #- Scenario: if HCR TAC results in F outside of 15% of target F, remove TAC contraint, else use it
  #-----------------------------------------------------------------------------
  if(iScenario == "FIAVBB"){
    tmpTAC                    <- iTAC[,FcY,]
    tmpTAC[]                  <- HCRTAC
    iTAC[,FcY]                <- HCRTAC
    #- First calculate what the F would be if the TAC would be constraint by 15%
    bidx                      <- which(HCRTAC[,,c("A")] > 1.15*iTAC[,ImY,,,c("A")])
    sidx                      <- which(HCRTAC[,,c("A")] < 0.85*iTAC[,ImY,,,c("A")])
    if(length(bidx)>0) tmpTAC[,,,,c("A"),bidx]   <- 1.15*iTAC[,ImY,,,c("A"),bidx]
    if(length(sidx)>0) tmpTAC[,,,,c("A"),sidx]   <- 0.85*iTAC[,ImY,,,c("A"),sidx]
    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=tmpTAC)

    for(i in dms$unit){
      stf@catch.n[,FcY,i]     <- stf@stock.n[,FcY,i]*(1-exp(-unitSums(stf@harvest[,FcY])-stf@m[,FcY,i]))*(stf@harvest[,FcY,i]/(unitSums(stf@harvest[,FcY])+stf@m[,FcY,i]))
      stf@catch[,FcY,i]       <- computeCatch(stf[,FcY,i])
      stf@landings.n[,FcY,i]  <- stf@catch.n[,FcY,i]
      stf@landings[,FcY,i]    <- computeLandings(stf[,FcY,i])
    }
    #- Determine which iters are outside the +/- 15% from target F or inside this range
    outidx                    <- which(apply(unitSums(stf@harvest[f26,FcY]),2:6,mean,na.rm=T) <  0.85 * F26 |
                                       apply(unitSums(stf@harvest[f26,FcY]),2:6,mean,na.rm=T) >  1.15 * F26)
    inidx                     <- which(apply(unitSums(stf@harvest[f26,FcY]),2:6,mean,na.rm=T) >= 0.85 * F26 &
                                       apply(unitSums(stf@harvest[f26,FcY]),2:6,mean,na.rm=T) <= 1.15 * F26)
    #- If outside this range, do not apply TAC constraint and use 15% cap on F
    if(length(outidx)>0){
      sidx                    <- which((apply(unitSums(stf@harvest[f26,FcY,,,,outidx]),2:6,mean,na.rm=T)) <  0.9 * F26[,,,,,outidx])
      bidx                    <- which((apply(unitSums(stf@harvest[f26,FcY,,,,outidx]),2:6,mean,na.rm=T)) >  1.1 * F26[,,,,,outidx])
      
      if(length(sidx)>0) stf@harvest[,FcY,"A",,,outidx[sidx]]       <- sweep(stf@harvest[,FcY,"A",,,outidx[sidx]],2:6,
                                                                        ((F26[,,,,,outidx[sidx]] * 0.9) - apply(unitSums(stf@harvest[f26,FcY,c("B","C","D"),,,outidx[sidx]]),2:6,mean,na.rm=T)) /
                                                                          apply(unitSums(stf@harvest[f26,FcY,"A",,,outidx[sidx]]),2:6,mean,na.rm=T),"*")
      if(length(bidx)>0) stf@harvest[,FcY,"A",,,outidx[bidx]]       <- sweep(stf@harvest[,FcY,"A",,,outidx[bidx]],2:6,
                                                                        ((F26[,,,,,outidx[bidx]] * 1.1) - apply(unitSums(stf@harvest[f26,FcY,c("B","C","D"),,,outidx[bidx]]),2:6,mean,na.rm=T)) /
                                                                          apply(unitSums(stf@harvest[f26,FcY,"A",,,outidx[bidx]]),2:6,mean,na.rm=T),"*")

      iTAC[,FcY,,,c("A"),outidx]                <-  quantSums(stf@catch.wt[,FcY,"A",,,outidx] * stf@stock.n[,FcY,"A",,,outidx] *
                                                    (1-exp(-unitSums(stf@harvest[,FcY,,,,outidx])-stf@m[,FcY,"A",,,outidx])) *
                                                    (stf@harvest[,FcY,"A",,,outidx]/(unitSums(stf@harvest[,FcY,,,,outidx])+stf@m[,FcY,"A",,,outidx])))
    }
    if(length(inidx)>0){
      bidx                    <- which(HCRTAC[,,c("A"),,,inidx] > 1.15*iTAC[,ImY,,,c("A"),inidx])
      sidx                    <- which(HCRTAC[,,c("A"),,,inidx] < 0.85*iTAC[,ImY,,,c("A"),inidx])
      if(length(bidx)>0) iTAC[,FcY,,,c("A"),inidx[bidx]] <- 1.15*iTAC[,ImY,,,c("A"),inidx[bidx]]
      if(length(sidx)>0) iTAC[,FcY,,,c("A"),inidx[sidx]] <- 0.85*iTAC[,ImY,,,c("A"),inidx[sidx]]
    }
    
    #- third banking and borrowing
    if(iYr == "2014"){             #Bank
      iTAC[,FcY,,,c("A")]         <- 0.9 * iTAC[,FcY,,,"A"]
    }
    if(iYr == "2014"){             #Repay banked part                               Borrow
      iTAC[,FcY,,,c("A")]         <- (iTAC[,ImY,,,c("A")] / 0.9 - iTAC[,ImY,,,c("A")]) + (1.1 * iTAC[,FcY,,,c("A")])
    }
    if(an(iYr) > 2015){            #Repay borrowed part                             Borrow
      iTAC[,FcY,,,c("A")]         <- (iTAC[,ImY,,,c("A")] / 1.1 - iTAC[,ImY,,,c("A")]) + (1.1 * iTAC[,FcY,,,c("A")])
    }

    #- If SSB below Blim, then no cap at all
    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY])
    SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
                                           exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
    idxssb                    <- which(SSB < ifelse(mpPoints$stabilityBreak == "Blim",mpPoints$Blim,mpPoints$Btrigger))
    if(length(idxssb)>0) iTAC[,FcY,,,c("A"),idxssb] <- HCRTAC[,,"A",,,idxssb]
  }
  
  #-----------------------------------------------------------------------------
  #- Scenario take mean of HCR TAC and ImY TAC
  #-----------------------------------------------------------------------------
  if(iScenario == "meanTAC"){
    iTAC[,FcY,c("B","C","D")] <- HCRTAC[,,c("B","C","D")]
    iTAC[,FcY,"A"]            <- (HCRTAC[,,"A"] + iTAC[,ImY,"A"]) / 2
    
    #- If SSB below Blim, then no cap at all
    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY])
    SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
                                           exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
    idxssb                    <- which(SSB < ifelse(mpPoints$stabilityBreak == "Blim",mpPoints$Blim,mpPoints$Btrigger))
    if(length(idxssb)>0) iTAC[,FcY,c("A"),,,idxssb] <- HCRTAC[,,"A",,,idxssb]
  }
  
  #-----------------------------------------------------------------------------
  #- Scenario bank and borrow: 1st year, bank, other years repay + borrow
  #-----------------------------------------------------------------------------
  if(iScenario == "BB"){

    #- First HCR with 15% IAV
    iTAC[,FcY]                <- HCRTAC
    bidx                      <- which(HCRTAC[,,c("A")] > 1.15*iTAC[,ImY,c("A")]) #bigger than  1.15x ImY TAC
    sidx                      <- which(HCRTAC[,,c("A")] < 0.85*iTAC[,ImY,c("A")]) #smaller than 0.85x ImY TAC
    if(length(bidx)>0) iTAC[,FcY,c("A"),,,bidx] <- 1.15*iTAC[,ImY,c("A"),,,bidx]
    if(length(sidx)>0) iTAC[,FcY,c("A"),,,sidx] <- 0.85*iTAC[,ImY,c("A"),,,sidx]

    #- Second banking and borrowing
    if(iYr == "2012"){             #Bank
      iTAC[,FcY,c("A")]         <- 0.9 * iTAC[,FcY,"A"]
    }
    if(iYr == "2013"){             #Repay banked part                               Borrow
      iTAC[,FcY,c("A")]         <- (iTAC[,ImY,c("A")] / 0.9 - iTAC[,ImY,c("A")]) + (1.1 * iTAC[,FcY,c("A")])
    }
    if(an(iYr) > 2013){            #Repay borrowed part                             Borrow
      iTAC[,FcY,c("A")]         <- (iTAC[,ImY,c("A")] / 1.1 - iTAC[,ImY,c("A")]) + (1.1 * iTAC[,FcY,c("A")])
    }

    #- Make sure no SSB is below Blim
    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY])
    SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
                                           exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
    idxssb                    <- which(SSB < ifelse(mpPoints$stabilityBreak == "Blim",mpPoints$Blim,mpPoints$Btrigger))
    if(length(idxssb)>0) iTAC[,FcY,c("A"),,,idxssb] <- HCRTAC[,,"A",,,idxssb]
  }
  
    #-----------------------------------------------------------------------------
  #- Scenario bank and borrow: 1st year, bank, other years repay + borrow
  #-----------------------------------------------------------------------------
  if(iScenario == "Bank"){

    #- First HCR with 15% IAV
    iTAC[,FcY]                <- HCRTAC
    bidx                      <- which(HCRTAC[,,c("A")] > 1.15*iTAC[,ImY,c("A")]) #bigger than  1.15x ImY TAC
    sidx                      <- which(HCRTAC[,,c("A")] < 0.85*iTAC[,ImY,c("A")]) #smaller than 0.85x ImY TAC
    if(length(bidx)>0) iTAC[,FcY,c("A"),,,bidx] <- 1.15*iTAC[,ImY,c("A"),,,bidx]
    if(length(sidx)>0) iTAC[,FcY,c("A"),,,sidx] <- 0.85*iTAC[,ImY,c("A"),,,sidx]

    #- Second banking
    if(iYr == "2012"){             #Bank
      iTAC[,FcY,c("A")]         <- 0.9 * iTAC[,FcY,"A"]
    }
    if(an(iYr) >= 2013){             #Repay banked part                               Bank
      iTAC[,FcY,c("A")]         <- (iTAC[,ImY,c("A")] / 0.9 - iTAC[,ImY,c("A")]) + (0.9 * iTAC[,FcY,c("A")])
    }

    #- Make sure no SSB is below Blim
    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY])
    SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
                                           exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
    idxssb                    <- which(SSB < ifelse(mpPoints$stabilityBreak == "Blim",mpPoints$Blim,mpPoints$Btrigger))
    if(length(idxssb)>0) iTAC[,FcY,c("A"),,,idxssb] <- HCRTAC[,,"A",,,idxssb]
  }

  
  #-----------------------------------------------------------------------------
  #- Scenario bank and borrow: 1st year, bank, other years repay + borrow
  #-----------------------------------------------------------------------------
  if(iScenario == "F3years"){
    iTAC[,FcY]                <- HCRTAC
    #- Project 2 years forward excluding FcY (which in total = 3 years ahead)
    for(i in ac(c(CtY,CtY1))){
      for(j in dimnames(stf@stock.n)$unit){
        stf@stock.n[1,i,j]    <- RECS$CtY
        stf@harvest[,i,j]     <- stf@harvest[,FcY,j]
        survivors             <- stock.n(stf)[,ac(an(i)-1),j] * exp(-unitSums(stf@harvest[,ac(an(i)-1)])-stf@m[,ac(an(i)-1),j])
        stf@stock.n[-1,i,j]   <- survivors[-dim(survivors)[1],,,,,]@.Data
        if (!is.na(range(stf,"plusgroup"))){
          stock.n(stf)[ac(range(stf,"max")),i,j] <- stock.n(stf)[ac(range(stf,"max")),i,j] + survivors[dim(survivors)[1],,,,,]@.Data
        }
      }
    }
    
    TACs                  <- iTAC[,FuY]; TACs[] <- NA
    for(i in FuY)
      TACs[,i]            <- round(harvestCatch(stf,i),0)
      
    #- Three year average
    iTAC[,FcY,c("A")]     <- apply(TACs[,ac(an(FcY):an(CtY1)),"A"],3:6,mean,na.rm=T)
    
    #- Make sure no SSB is below Blim
    stf@harvest[,FcY]         <- fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY])
    SSB                       <- quantSums(stf@stock.n[,FcY,1]*stf@stock.wt[,FcY,1]*stf@mat[,FcY,1]*
                                           exp(-unitSums(stf@harvest[,FcY])*stf@harvest.spwn[,FcY,1]-stf@m[,FcY,1]*stf@m.spwn[,FcY,1]))
    idxssb                    <- which(SSB < ifelse(mpPoints$stabilityBreak == "Blim",mpPoints$Blim,mpPoints$Btrigger))
    if(length(idxssb)>0) iTAC[,FcY,c("A"),,,idxssb] <- HCRTAC[,,"A",,,idxssb]
  }

return(list(TAC=iTAC[,FcY],fSTF=fleet.harvestFF(stk=stf,iYr=FcY,TACS=iTAC[,FcY]),HCRTAC=HCRTAC,SSB=list(HCRSSB=HCRSSB,SSB=SSB)))}




