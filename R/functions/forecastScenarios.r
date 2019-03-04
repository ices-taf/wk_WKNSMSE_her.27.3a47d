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
#projectNSH(stkAssessment,fishery,iYr,TAC,histMaxYr,referencePoints,managementRule,catchVar)
projectNSH <- function(iStocks,iFishery,iYr,iTAC,iHistMaxYr,mpPoints,managementRule){
  require(minpack.lm)

  lin     <- substr(R.Version()$os,1,3)== "lin"
  stk     <- iStocks
  #===============================================================================
  # Setup control file
  #===============================================================================

  DtY         <- ac(range(stk)["maxyear"]-1) #Data year
  ImY         <- ac(an(DtY)+1) #Intermediate year
  FcY         <- ac(an(DtY)+2) #Forecast year
  CtY         <- ac(an(DtY)+3) #Continuation year
  CtY1        <- ac(an(DtY)+4)
  FuY         <- c(ImY,FcY,CtY,CtY1)#Future years
  pyears      <- an(DtY) - iHistMaxYr #Years away from original historic max year
  
  #- We substract pyears from the starting point, because with every year forward, we move the year ranges one forward too.
  RECS        <- FLQuants("ImY"=stk@stock.n[1,ac(ImY)],"FcY"=exp(apply(log(rec(stk)[,ac((range(stk)["maxyear"]-8-pyears):(range(stk)["maxyear"]))]),3:6,mean,na.rm=T)),
                                                       "CtY"=exp(apply(log(rec(stk)[,ac((range(stk)["maxyear"]-8-pyears):(range(stk)["maxyear"]))]),3:6,mean,na.rm=T)))

  yrs1        <- list("m.spwn","harvest.spwn","stock.wt")
  yrs3        <- list("mat")
  yrs5        <- list("m")

  dsc         <- "North Sea Herring"
  nam         <- "NSH"
  dms         <- dimnames(stk@m)
  dms$year    <- c(rev(rev(dms$year)[1:3]),FcY,CtY,CtY1)
  dms$unit    <- c("A","B","C","D")

  f01         <- ac(0:1)
  f26         <- ac(2:6)

  #===============================================================================
  # Setup stock file
  #===============================================================================

  stf         <- FLStock(FLQuant(NA,dimnames=dms),name=nam,desc=dsc)
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

  iTAC[,ImY]                <- iTAC[,ImY]
  iCATCH                    <- iTAC
  iCATCH[,ImY,,,"A"]        <- iTAC[,ImY,,,"A"] + TAC_var[ImY,,'Ctransfer'] * iTAC[,ImY,,,"C"]
  iCATCH[,ImY,,,"B"]        <- iTAC[,ImY,,,"B",drop=T] * TAC_var[ImY,,'Buptake',drop=T]
  iCATCH[,ImY,,,"C"]        <- (1-TAC_var[ImY,,'Ctransfer']) * iTAC[,ImY,,,"C"]
  iCATCH[,ImY,,,"D"]        <- iTAC[,ImY,,,"D",drop=T] * TAC_var[ImY,,'Duptake',drop=T] * TAC_var[ImY,,'Dsplit',drop=T]
  stf@harvest[,ImY]         <- fleet.harvestFF(stk=stf,iYr=ImY,TACS=iCATCH[,ImY])

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
  stf@harvest[,FcY]                                          <- stf@harvest[,ImY]

  ###--- Management options ---###

  #- HCR A
  res <- matrix(NA,nrow=2,ncol=dims(stf)$iter,dimnames=list(c("A","B"),dimnames(stf@stock.n)$iter))

  if(managementRule$HCR == "A"){
    for(iTer in 1:dims(stf)$iter){
      res[,iTer]              <- nls.lm(par=rep(1,2),find.FAB_HCRA,
                                        stk=stf[,FcY,c("A","B"),,,iTer],
                                        f01=f01,
                                        f26=f26,
                                        TACS=iTAC[,FcY,c("A","B"),,,iTer],
                                        mpPoints=mpPoints,
                                        jac=NULL,lower=rep(1e-8,2),upper=rep(1e5,2),control=nls.lm.control(maxiter=1000))$par
    }
  }
  #- HCR B
  if(managementRule$HCR == "B"){
    for(iTer in 1:dims(stf)$iter){
      res[,iTer]              <- nls.lm(par=rep(1,2),find.FAB_HCRB,stk=stf[,FcY,c("A","B"),,,iTer],f01=f01,f26=f26,TACS=iTAC[,FcY,c("A","B"),,,iTer],
                                        mpPoints=mpPoints,jac=NULL,lower=rep(1e-8,2),upper=rep(1e5,2),control=nls.lm.control(maxiter=1000))$par
    }
  }


  #- Calculate TACs
  stf@harvest[,FcY,c("A","B")]         <- sweep(stf@harvest[,FcY,c("A","B")],c(3,6),res,"*")
  for(i in dms$unit){
    stf@catch.n[,FcY,i]     <- stf@stock.n[,FcY,i]*(1-exp(-unitSums(stf@harvest[,FcY])-stf@m[,FcY,i]))*(stf@harvest[,FcY,i]/(unitSums(stf@harvest[,FcY])+stf@m[,FcY,i]))
    stf@catch[,FcY,i]       <- computeCatch(stf[,FcY,i])
    stf@landings.n[,FcY,i]  <- stf@catch.n[,FcY,i]
    stf@landings[,FcY,i]    <- computeLandings(stf[,FcY,i])
  }
  iTAC[,FcY,,,c("A","B")]   <- stf@catch[,FcY,c("A","B")]
  #- Get SSB as a results from the TACS
  totF                      <- unitSums(stf@harvest[,FcY,c("A","B")])
  SSBHCR                    <- quantSums(stf@stock.n[,FcY,1] * stf@stock.wt[,FcY,1] * exp(-totF*stf@harvest.spwn[,FcY,1] - stf@m[,FcY,1]) * stf@mat[,FcY,1])
  ssb.CtY                   <- NA
  print(SSBHCR)
  print(iTAC[,FcY,,,c("A","B")])
  print(quantMeans(unitSums(stf@harvest[f26,FcY,c("A","B")])))
  print(quantMeans(unitSums(stf@harvest[f01,FcY,c("A","B")])))

  #-----------------------------------------------------------------------------
  #- TAC IAV
  #-----------------------------------------------------------------------------
  if(is.null(managementRule$TACIAV) == F){
    mrF                             <- managementRule$TACIAV
    if(managementRule$TACIAV[1] == "A"){
      idx                             <- which(SSBHCR > mpPoints$Btrigger)
      if(length(idx)>0){
        for(imrF in mrF){
          bidx                          <- which(iTAC[,FcY,,,imrF,idx] > (1.25 * iTAC[,ImY,,,imrF,idx]))
          iTAC[,FcY,,,imrF,idx[bidx]]   <- 1.25 * iTAC[,ImY,,,imrF,idx[bidx]]
          sidx                          <- which(iTAC[,FcY,,,imrF,idx] < (0.8 * iTAC[,ImY,,,imrF,idx]))
          iTAC[,FcY,,,imrF,idx[sidx]]   <- 0.8 * iTAC[,ImY,,,imrF,idx[sidx]]
        }
      }
    }
    if(managementRule$TACIAV[1] == "E"){
      mrF                             <- c("A","B")
      idx                             <- 1:dim(iTAC)[6]
      for(imrF in mrF){
        bidx                          <- which(iTAC[,FcY,,,imrF,idx] > (1.25 * iTAC[,ImY,,,imrF,idx]))
        iTAC[,FcY,,,imrF,idx[bidx]]   <- 1.25 * iTAC[,ImY,,,imrF,idx[bidx]]
        sidx                          <- which(iTAC[,FcY,,,imrF,idx] < (0.8 * iTAC[,ImY,,,imrF,idx]))
        iTAC[,FcY,,,imrF,idx[sidx]]   <- 0.8 * iTAC[,ImY,,,imrF,idx[sidx]]
      }
    }
  }

  #-----------------------------------------------------------------------------
  #- Scenario bank and borrow: 1st year, bank, other years repay + borrow
  #-----------------------------------------------------------------------------
  if(is.null(managementRule$BB) == F){
    if(managementRule$BB[1] == "A"){
      bidx                        <- which(SSBHCR > mpPoints$Btrigger)
      mrF                         <- managementRule$BB
      if(length(bidx)>0){
        for(imrF in mrF){
          if(iYr == "2018"){             #Bank
            iTAC[,FcY,,,imrF,bidx]  <- 0.9 * iTAC[,FcY,,,imrF,bidx]
          }
          if(iYr == "2019"){             #Repay banked part                               Borrow
            iTAC[,FcY,,,imrF,bidx]  <- (iTAC[,ImY,,,imrF,bidx] / 0.9 - iTAC[,ImY,,,imrF,bidx]) + (1.1 * iTAC[,FcY,,,imrF,bidx])
          }
          if(an(iYr) > 2019){            #Repay borrowed part                             Borrow
            iTAC[,FcY,,,imrF,bidx]  <- (iTAC[,ImY,,,imrF,bidx] / 1.1 - iTAC[,ImY,,,imrF,bidx]) + (1.1 * iTAC[,FcY,,,imrF,bidx])
          }
        }
      }
    }
    if(managementRule$BB[1] == "E"){
    
      #Update to continuation year
      stf@harvest[,CtY,c("A","B")]                                <- fleet.harvestFF(stk=stf[,,c("A","B")],iYr=CtY,TACS=iTAC[,FcY,,,c("A","B")])
      for(i in dms$unit) stf@stock.n[1,CtY,i]                     <- RECS$CtY
      for(i in dms$unit) stf@stock.n[2:(dims(stf)$age-1),CtY,i]   <- (stf@stock.n[,FcY,1]*exp(-unitSums(stf@harvest[,FcY,c("A","B")])-stf@m[,FcY,1]))[ac(range(stf)["min"]:(range(stf)["max"]-2)),]
      for(i in dms$unit) stf@stock.n[dims(stf)$age,CtY,i]         <- apply((stf@stock.n[,FcY,1]*exp(-unitSums(stf@harvest[,FcY,c("A","B")])-stf@m[,FcY,1]))[ac((range(stf)["max"]-1):range(stf)["max"]),],2:6,sum,na.rm=T)
      ssb.CtY                                                     <- quantSums(stf@stock.n[,CtY,1] * stf@stock.wt[,CtY,1]*stf@mat[,CtY,1]*exp(-unitSums(stf@harvest[,FcY,c("A","B")])*stf@harvest.spwn[,CtY,1]-stf@m[,CtY,1]*stf@m.spwn[,CtY,1])) #assume same harvest as in FcY
      totF                                                        <- unitSums(stf@harvest[,FcY,c("A","B")])

      mrF                         <- c("A","B")
      idx                         <- which(SSBHCR > mpPoints$Bpa | quantMeans(totF[ac(2:6),]) < mpPoints$Fpa)
      bidx                        <- idx[which(SSBHCR[,,,,,idx] > mpPoints$Bpa & ssb.CtY[,,,,,idx] > mpPoints$Bpa)]
      if(length(bidx)>0){
        for(imrF in mrF){
          if(iYr == "2018"){             #Bank
            iTAC[,FcY,,,imrF,bidx]  <- 0.9 * iTAC[,FcY,,,imrF,bidx]
          }
          if(iYr == "2019"){             #Repay banked part                               Borrow
            iTAC[,FcY,,,imrF,bidx]  <- (iTAC[,ImY,,,imrF,bidx] / 0.9 - iTAC[,ImY,,,imrF,bidx]) + (1.1 * iTAC[,FcY,,,imrF,bidx])
          }
          if(an(iYr) > 2019){            #Repay borrowed part                             Borrow
            iTAC[,FcY,,,imrF,bidx]  <- (iTAC[,ImY,,,imrF,bidx] / 1.1 - iTAC[,ImY,,,imrF,bidx]) + (1.1 * iTAC[,FcY,,,imrF,bidx])
          }
        }
      }
    }
  }
  

return(list(TAC=iTAC[,FcY],SSB=list("FcY"=SSBHCR,"CtY"=ssb.CtY),Fbar=totF))}




