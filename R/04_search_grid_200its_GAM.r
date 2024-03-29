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

library(mgcv)
library(lattice)

#path <- 'D:/git/wk_WKNSMSE_her.27.3a47d/R'
path <- 'C:/git/wk_WKNSMSE_her.27.3a47d/R'
assessment_name   <- "NSAS_WKNSMSE2018"
try(setwd(path),silent=TRUE)

# paths to different subfolders
dataPath      <- file.path(".","data/")
outPath       <- file.path(".","results/")
scriptPath    <- file.path(".","side_scripts/")
functionPath  <- file.path(".","functions/")

outputName <- 'grid_search_interp'

PDF <- FALSE
PNG <- ifelse(PDF,F,T)
if(PDF) pdf(file.path(outPath,'plots',paste0(outputName,".pdf")))
if(PNG) png(file.path(outPath,'plots',paste0(outputName,"_%02d.png")),
            units = "px",
            height=800,
            width=672,
            bg = "white")

#-------------------------------------------------------------------------------
# 2) Reading data, fitting gam model to grid and plot fitted grid
#-------------------------------------------------------------------------------

# interpolated sequences
FtarInter   <- seq(0.15,0.30,0.01)
BtrigInter  <- seq(0.8e6,2.2e6,0.1e6)


list_HCR <- c('A',
              'B',
              'A_IAV_AB_BB_AB',
              'A_IAV_A_BB_A',
              'B_IAV_E_BB_E')

optPoints <- data.frame(matrix(ncol = 3, nrow = length(list_HCR)))
colnames(optPoints) <- c('case','Ftar','Btrig')

# loop through the grids

for(idxHCR in 1:length(list_HCR)){
  HCRString <- list_HCR[idxHCR]
  
  LTR     <- read.csv(file.path(outPath,paste0('grid_HCR_',HCRString),paste0('grid_search_',HCRString,'_LTR3_vec','.csv')))
  LTY     <- read.csv(file.path(outPath,paste0('grid_HCR_',HCRString),paste0('grid_search_',HCRString,'_LTY_vec','.csv')))
  IAV     <- read.csv(file.path(outPath,paste0('grid_HCR_',HCRString),paste0('grid_search_',HCRString,'_IAV_vec','.csv')))
  
  gamData           <- expand.grid(Ftarget = FtarInter,
                                   Btrigger=BtrigInter)
  
  # fit LTR
  gamData$LTR_gam       <- predict.gam(gam(Risk3 ~  s(Ftarget,Btrigger),data=LTR),newdata=gamData)
  
  # fit LTY
  gamData$LTY_gam       <- predict.gam(gam(yield ~  s(Ftarget,Btrigger),data=LTY),newdata=gamData)
  
  # fit IAV
  gamData$IAV_gam       <- predict.gam(gam(IAV ~  s(Ftarget,Btrigger),data=IAV),newdata=gamData)
  
  idxCoarse <- (which(!is.na(match(gamData$Ftarget,LTR$Ftarget)) & !is.na(match(gamData$Btrigger,LTR$Btrigger))))
  
  gamDataCoarse <- gamData[idxCoarse,]
  
  for(idxCoarse in 1:dim(LTR)[1]){
    idxInit <- which(gamDataCoarse$Ftarget[idxCoarse] == LTR$Ftarget & gamDataCoarse$Btrigger[idxCoarse] == LTR$Btrigger)
    gamDataCoarse$LTR_sim[idxCoarse]    <- LTR$Risk3[idxInit]
    gamDataCoarse$LTY_sim[idxCoarse]    <- LTY$yield[idxInit]
    gamDataCoarse$IAV_sim[idxCoarse]    <- IAV$IAV[idxInit]
    
    # calculate differences in % between the simulations and the GAM fit
    gamDataCoarse$LTR_gamFit[idxCoarse] <- abs(gamDataCoarse$LTR_sim[idxCoarse] - gamDataCoarse$LTR_gam[idxCoarse])
  }
  
  filtRes <- gamData[which(gamData$LTR_gam <= 0.05 & gamData$LTR_gam > 0),] # values above 0.05
  
  optIdx <- which(filtRes$LTY_gam == max(filtRes$LTY_gam))
  
  optPoints$Ftar[idxHCR]  <- filtRes$Ftarget[optIdx]
  optPoints$Btrig[idxHCR] <- filtRes$Btrigger[optIdx]
  optPoints$case[idxHCR]  <- HCRString
  
  # plot LTR
  gamLTRPlot        <- gamData[which(gamData$LTR_gam <= 0.05),] # values above 0.05
  
  print(levelplot(LTR_gam ~  Btrigger*Ftarget, 
                  data=filtRes,
                  xlim = c(min(BtrigInter),max(BtrigInter)),
                  ylim = c(min(FtarInter),max(FtarInter)),
                  main=paste0(HCRString," - LTR")))
  
  # plot LTY
  gamLTYPlot        <- gamData[which(gamData$LTR_gam <= 0.05),] # values above 0.05
  print(levelplot(LTY_gam ~  Btrigger*Ftarget, 
                  data=filtRes,
                  xlim = c(min(BtrigInter),max(BtrigInter)),
                  ylim = c(min(FtarInter),max(FtarInter)),
                  main=paste0(HCRString," - LTY")))
  
  # plot IAV
  #print(levelplot(IAV_gam ~  Ftarget * Btrigger, 
  #                data=gamData,
  #                main=paste0(HCRString," - IAV")))
  
  # plot LTR difference with simulations
  print(levelplot(LTR_gamFit ~  Btrigger*Ftarget, 
                  data=gamDataCoarse,
                  xlim = c(min(BtrigInter),max(BtrigInter)),
                  ylim = c(min(FtarInter),max(FtarInter)),
                  main=paste0(HCRString," - diff")))
}

dev.off()
