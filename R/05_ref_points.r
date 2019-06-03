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
path <- 'E:/git/wk_WKNSMSE_her.27.3a47d/R'
assessment_name   <- "NSAS_WKNSMSE2018"
try(setwd(path),silent=TRUE)

# paths to different subfolders
dataPath      <- file.path(".","data/")
outPath       <- file.path(".","results/")
scriptPath    <- file.path(".","side_scripts/")
functionPath  <- file.path(".","functions/")


list_HCR <- c('A',
              'B',
              'A_IAV_AB_BB_AB',
              'A_IAV_A_BB_A')

optPoints <- data.frame(matrix(ncol = 3, nrow = length(list_HCR)))
colnames(optPoints) <- c('case','Ftar','Btrig')

for(idxHCR in 1:length(list_HCR)){
  HCRString <- list_HCR[idxHCR]
  
  LTR     <- read.csv(file.path(outPath,paste0('grid_HCR_',HCRString),paste0('grid_search_1000its_',HCRString,'_LTR3_vec_1000its','.csv')))
  LTY     <- read.csv(file.path(outPath,paste0('grid_HCR_',HCRString),paste0('grid_search_1000its_',HCRString,'_LTY_vec_1000its','.csv')))
  
  idxFiltRisk <- which(LTR$Risk3 <= 0.05)
  
  LTR <- LTR[idxFiltRisk,]
  LTY <- LTY[idxFiltRisk,]
  
  idxOpt <- which(LTY$yield == max(LTY$yield))
  
  optPoints$Ftar[idxHCR]  <- LTY$Ftarget[idxOpt]
  optPoints$Btrig[idxHCR] <- LTY$Btrigger[idxOpt]
  optPoints$case[idxHCR]  <- HCRString
}