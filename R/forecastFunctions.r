#-------------------------------------------------------------------------------
#- Calculate the harvest by fleet given a TAC for a stock
#-------------------------------------------------------------------------------
fleet.harvestFF<- function(stk,iYr,TACS){
                    lin     <- substr(R.Version()$os,1,3)== "lin"
                    nUnits                          <- dims(stk)$unit
                    nIter                           <- dims(TACS)$iter
                    res                             <- matrix(NA,ncol=nIter,nrow=nUnits,dimnames=list(units=dimnames(stk@stock.n)$unit,iter=1:nIter))
                    for(iTer in 1:nIter){
                      Ns      = c(stk@stock.n[,iYr,1,,,iTer]@.Data)
                      Fs      = stk@harvest[,iYr,,,,iTer,drop=T]@.Data
                      Cwts    = stk@catch.wt[,iYr,,,,iTer,drop=T]@.Data
                      Ms      = c(stk@m[,iYr,1,,,iTer]@.Data)
                      if(lin) res[,iTer] <- nls.lm(par=rep(1,nUnits),rescaleFFF,Ns=Ns,Fs=Fs,Cwts=Cwts,Ms=Ms,TACS=c(TACS[,,,,,iTer]),lower=rep(0,nUnits),upper=rep(1e5,nUnits),jac=NULL,nls.lm.control(ftol = (.Machine$double.eps)))$par
                      if(!lin)res[,iTer] <- nls.lm(par=rep(1,nUnits),rescaleFFF,Ns=Ns,Fs=Fs,Cwts=Cwts,Ms=Ms,TACS=c(TACS[,,,,,iTer]),lower=rep(0,nUnits),upper=rep(1e5,nUnits),jac=NULL,nls.lm.control(ftol = (.Machine$double.eps)))$par
                    }
                    stk@harvest[,iYr]               <- sweep(stk@harvest[,iYr],c(3,6),res,"*")
                 return(stk@harvest[,iYr])}

#-------------------------------------------------------------------------------
#- Function to scale the F pattern for stock
#-------------------------------------------------------------------------------
rescaleFFF     <- function(mult,Ns.=Ns,Fs.=Fs,Cwts.=Cwts,Ms.=Ms,TACS.=TACS){
                    Fs.    <- t(t(Fs.) * mult)
                    stkZ  <- rowSums(Fs.) + Ms.
                    res   <- sqrt(c((TACS. - c(colSums(sweep(sweep(Fs.,1,stkZ,"/")* Cwts.,1,Ns. * (1-exp(-stkZ)),"*")))))^2)
                 return(res)}
                 
#-------------------------------------------------------------------------------
#- Management plan: calculate A and B TAC
#-------------------------------------------------------------------------------

                 
find.FAB_HCRA  <- function(mult,stk.=stk,f01.=f01,f26.=f26,TACS.=TACS,mpPoints.=mpPoints){
                    Fs      <- sweep(stk.@harvest@.Data,3,mult,"*")
                    Ns      <- stk.@stock.n@.Data[,,1,,,]
                    Ms      <- stk.@m@.Data[,,1,,,]
                    Hspwns  <- stk.@harvest.spwn@.Data[,,1,,,]
                    Mspwns  <- stk.@m.spwn@.Data[,,1,,,]
                    Mats    <- stk.@mat@.Data[,,1,,,]
                    Swghts  <- stk.@stock.wt@.Data[,,1,,,]
                    Cwghts  <- stk.@catch.wt@.Data

                    bigF              <- apply(Fs,1,sum)
                    ssb               <- sum(Ns * Swghts * exp(-bigF*Hspwns - Ms*Mspwns) * Mats)
                    if(ssb < mpPoints.$Btrigger){
                      resA <- mpPoints.$Ftarget*ssb/mpPoints.$Btrigger
                      resB <- mpPoints.$F01*ssb/mpPoints.$Btrigger
                    }
                    if(ssb > mpPoints.$Btrigger){
                      resA <- mpPoints.$Ftarget
                      resB <- mpPoints.$F01
                    }
                    fbarB     <- mean(bigF[f01.])
                    fbarA     <- mean(bigF[f26.])
                    ret       <- c(sqrt(c((fbarA-resA)^2,(fbarB-resB)^2)^2))
                 return(ret)}

find.FAB_HCRB  <- function(mult,stk.=stk,f01.=f01,f26.=f26,TACS.=TACS,mpPoints.=mpPoints){
                    Fs      <- sweep(stk.@harvest@.Data,3,mult,"*")
                    Ns      <- stk.@stock.n@.Data[,,1,,,]
                    Ms      <- stk.@m@.Data[,,1,,,]
                    Hspwns  <- stk.@harvest.spwn@.Data[,,1,,,]
                    Mspwns  <- stk.@m.spwn@.Data[,,1,,,]
                    Mats    <- stk.@mat@.Data[,,1,,,]
                    Swghts  <- stk.@stock.wt@.Data[,,1,,,]
                    Cwghts  <- stk.@catch.wt@.Data

                    bigF              <- apply(Fs,1,sum)
                    ssb               <- sum(Ns * Swghts * exp(-bigF*Hspwns - Ms*Mspwns) * Mats)
                    if(ssb < mpPoints.$Btrigger & ssb > mpPoints.$Blim){
                      resA <- mpPoints.$Ftarget*ssb/mpPoints.$Btrigger
                      resB <- mpPoints.$F01*ssb/mpPoints.$Btrigger
                    }
                    if(ssb <= mpPoints.$Blim){
                      resA <- 0.1
                      resB <- 0.04
                    }
                    if(ssb >= mpPoints.$Btrigger){
                      resA <- mpPoints.$Ftarget
                      resB <- mpPoints.$F01
                    }
                    fbarB     <- mean(bigF[f01.])
                    fbarA     <- mean(bigF[f26.])
                    ret       <- c(sqrt(c((fbarA-resA)^2,(fbarB-resB)^2)^2))
                 return(ret)}

#-------------------------------------------------------------------------------
#- Convert TACs from forecast to realised F by the fishing fleet
#-------------------------------------------------------------------------------

TAC2sel <- function(mult,iYr,iBiol,iFishery,iTAC,catchVar,TAC_var,iTer){
  Ns  <- iBiol@stock.n[,iYr,,,,iTer,drop=T];
  Fs  <- sweep(iFishery@landings.sel[,iYr,,,,iTer],3:6,mult,"*")[,drop=T]
  Ms  <- iBiol@m[,iYr,,,,iTer,drop=T]
  Wts <- iFishery@landings.wt[,iYr,,,,iTer,drop=T]

  Cs  <- colSums(sweep(sweep(Fs,1,rowSums(Fs)+Ms,"/") * Wts,1,Ns * (1-exp(-rowSums(Fs)-Ms)),"*"))

  Atarget <- iTAC[,FcY,,,"A",iTer,drop=T]
  Btarget <- iTAC[,FcY,,,"B",iTer,drop=T]
  Ctarget <- sum(rowSums(Fs)*catchVar[1,FcY,,"FCprop",,iTer,drop=T])
  Dtarget <- iTAC[,FcY,,,"D",iTer,drop=T] * TAC_var[FcY,iTer,"Dsplit"] * TAC_var[FcY,iTer,"Duptake"]

  ret     <- sqrt((c(Atarget,Btarget,Ctarget,Dtarget) - c(Cs[1],Cs[2],sum(Fs[,3]),Cs[4]))^2)
return(ret)}