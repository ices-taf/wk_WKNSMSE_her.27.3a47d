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
#
# Note: The Cfleet is defined as a proportion of F
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
TAC[,ac((max(TAC_C[,1])+1):(futureMaxYr+3)),,,"C"] <- rnorm(length((max(TAC_C[,1])+1):(futureMaxYr+3)),
                                                            mean(TAC_C[,2]),
                                                            sd(TAC_C[,2])) # set random TAC for the C fleet in WB


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
# 4) prepare stf object for intermediate year prior to entering MSE
#
#------------------------------------------------------------------------------#

# standard deviation for the recruitment for the first projected year. This is 
# because the recruitment is infered from the SAM object and residulas are not 
# computed for the first year as this was used to create random samples
# In the subsequent years, the standard deviation can be infered from the SAM object
sd_rec_init <- sd(log(subset(rec(NSH.sim),year==histMaxYr+1)$value))^2

# create stf object and calculate ImY
stf <- stf_ImY( NSH.sim,
                stocks,
                fisheryFuture,
                TAC,
                TAC_var,
                FCPropIts,
                c('2018','2019','2020'))

#------------------------------------------------------------------------------#
# 5) Start running the MSE
#------------------------------------------------------------------------------#

# loop on iterations. Now just testing 1 iteration
idxIter <- 1
iYr     <-2018

start.time          <- Sys.time()
bunit               <- dimnames(biol@n)$unit
for (iYr in an(projPeriod)){
  cat(iYr,"\n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))
  
  #-------------------------------------------------------------------------------
  # Biology
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # Assessment
  #-------------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  # Fishery, i.e. Short term forecast
  #-------------------------------------------------------------------------------
  
  # define years for fisheries
  DtY <- ac(iYr-1)  # terminal year in the assessment
  ImY <- ac(iYr)    # intermediate year in the short term forecast, ie, current year
  FcY <- ac(iYr+1)  # year for which the advice is given, ie, forecast two years ahead the last year in the assessment
  CtY <- ac(an(DtY)+3)             #Continuation year
  FuY <- c(ImY,FcY,CtY) # combination of future years
  
  ############# Recruitment #############
  # Recruitment for intermediate year, forecast year and continuation year
  # intermediate year: taken from SAM
  # forecast and continuation years: 10 years weighted average
  if(iYr == histMaxYr+1){
    recWeights <- sd_rec_init
  }else{
    recWeights  <- subset(NSH.sam@params, name=="logR")$std.dev^2
  }
  
  RECS        <- FLQuants("ImY" =FLQuant(subset(rec(NSH.sim[[idxIter]]),
                                                year==ImY)$value,
                                         dimnames=list(age="0",
                                                       year=ImY,
                                                       unit="unique",
                                                       season="all",
                                                       area="unique",
                                                       iter=1:nits)),
                          "FcY" =exp(apply(log(rec(stocks[[idxIter]])[,ac((an(DtY)-10):DtY)]),3:6,weighted.mean,w=1/rev(rev(recWeights)[2:12]),na.rm=T)),
                          "CtY" =exp(apply(log(rec(stk)[,ac((an(DtY)-10):DtY)]),3:6,weighted.mean,w=1/rev(rev(recWeights)[2:12]),na.rm=T)))
  
  # computation recruitment for stf
  
  # update fishery object for forecast year
  for(i in dms$unit) stf@stock.n[1,FcY,i]                     <- RECS$FcY
  for(i in dms$unit) stf@stock.n[2:(dims(stf)$age-1),FcY,i]   <- (stf@stock.n[,ImY,1]*exp(-unitSums(stf@harvest[,ImY])-stf@m[,ImY,1]))[ac(range(stf)["min"]:(range(stf)["max"]-2)),]
  for(i in dms$unit) stf@stock.n[dims(stf)$age,FcY,i]         <- apply((stf@stock.n[,ImY,1]*exp(-unitSums(stf@harvest[,ImY])-stf@m[,ImY,1]))[ac((range(stf)["max"]-1):range(stf)["max"]),],2:6,sum,na.rm=T)
  
  # apply HCR to obtain F2-6 and F0-1
  
  
  #----------------------------------------
  # Assessment
  #----------------------------------------
  
  #----------------------------------------
  # Operating model
  #----------------------------------------
  
  
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


