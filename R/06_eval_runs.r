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

library(minpack.lm)  # install.packages("minpack.lm")
library(stats)
library(FLSAM)

# define path to directory
#path          <- "D:/Work/Herring MSE/NSAS/"
#path              <- "D:/git/wk_WKNSMSE_her.27.3a47d/R/"
#path              <- "F:/WKNSMSE/wk_WKNSMSE_her.27.3a47d/R"
#path <- 'E:/git/wk_WKNSMSE_her.27.3a47d/R'
path <- 'E:/git/wk_WKNSMSE_her.27.3a47d/R'
#path <- 'D:/git/wk_WKNSMSE_her.27.3a47d/R'
assessment_name   <- "NSAS_WKNSMSE2018"
try(setwd(path),silent=TRUE)

# paths to different subfolders
dataPath      <- file.path(".","data/")
outPath       <- file.path(".","results/")
scriptPath    <- file.path(".","side_scripts/")
functionPath  <- file.path(".","functions/")

outputName <- 'eval_runs'

referencePoints <- list(Fmsy = 0.26,
                        Fsq  = NA,
                        Flim = 0.34,
                        Fpa  = 0.30,
                        Blim = 800000,
                        Bpa  = 900000,
                        MSYBtrigger = 1400000)
nits <- 1000

# initialaze tables
eval.options  <- c( 'Fopt',
                    'F0.9',
                    'F1.1',
                    'Fmsy')

eval.case     <- c('A',
                   'B',
                   'A_IAV_A_BB_A',
                   'A_IAV_AB_BB_AB')

eval.metrics  <- c('yield',
                   'SSB',
                   'risk',
                   'IAV')

eval.period   <- list(c('ST',ac(2019:2021)),
                      c('MT',ac(2022:2026)),
                      c('LT',ac(2027:2036)))

eval.table.colNames <- c('HCR','MS','F_option','Ftar','Btrig')
for(idx1 in 1:length(eval.metrics)){
  for(idx2 in 1:length(eval.period)){
    myColName <- paste0(eval.metrics[idx1],'_',eval.period[[idx2]][1])
    eval.table.colNames <- append(eval.table.colNames,myColName)
  }
}

eval.table <- data.frame(matrix(ncol = length(eval.table.colNames),
                                nrow = length(eval.options)*length(eval.case)))

colnames(eval.table) <- eval.table.colNames

idxMain <- 1
for(idx1 in 1:length(eval.case)){
  for(idx2 in 1:length(eval.options)){
    if(eval.case[idx1] =='A') HCRCode <- 'A'
    if(eval.case[idx1] == 'B') HCRCode <- 'B'
    if(eval.case[idx1] =='A_IAV_A_BB_A') HCRCode <- 'A+C'
    if(eval.case[idx1] =='A_IAV_AB_BB_AB') HCRCode <- 'A+D'
    if(eval.case[idx1] =='B_IAV_E_BB_E') HCRCode <- 'B+E'
    
    eval.table[idxMain,'HCR']       <- eval.case[idx1]
    eval.table[idxMain,'F_option']  <- eval.options[idx2]
    eval.table[idxMain,'MS']        <- HCRCode
    idxMain <- idxMain + 1
  }
}

# read reference points for each case
eval.ref <- read.table(file.path(outPath,'ref_points.csv'),sep=',')

#-------------------------------------------------------------------------------
# 2) fill in evaluation run table
#-------------------------------------------------------------------------------

# define string array for f01 and f26
f01         <- ac(0:1)
f26         <- ac(2:6)

# loop on HCRs
for(idxHCR in 1:length(eval.case)){
  
  HCRLoad <- eval.case[idxHCR]
  Btrig   <- as.numeric(as.vector(eval.ref[eval.ref[,1] == HCRLoad,3]))
  Ftar    <- as.numeric(as.vector(eval.ref[eval.ref[,1] == HCRLoad,2]))
  
  # loop on F cases
  for(idxFOption in 1:length(eval.options)){
    if(eval.options[idxFOption] == 'F0.9') Foption <- Ftar*0.9
    if(eval.options[idxFOption] == 'F1.1') Foption <- Ftar*1.1
    if(eval.options[idxFOption] == 'Fmsy') Foption <- referencePoints$Fmsy
    if(eval.options[idxFOption] == 'Fopt') Foption <- Ftar
    
    idxTable <- which(eval.table['HCR'] == HCRLoad & eval.table['F_option'] == eval.options[idxFOption])
    
    fileName <- paste0('NSAS_Ftar_',Foption,'_Btrig_',Btrig,'.*')
    
    fileList <- list.files(file.path(outPath,paste0('grid_HCR_',HCRLoad),'eval_run'),pattern = fileName)
    
    print(fileList)
    load(file.path(outPath,paste0('grid_HCR_',HCRLoad),'eval_run',fileList))
    
    biol@catch    <- computeCatch(biol)
    biol@stock    <- computeStock(biol)
    biol@landings <- computeLandings(biol)
    
    eval.table[idxTable,'Ftar']   <- Foption
    eval.table[idxTable,'Btrig']  <- Btrig
    
    # loop on metric period
    for(idxMetric in 1:length(eval.period)){
      metricsPeriod <- eval.period[[idxMetric]][2:length(eval.period[[idxMetric]])]
      
      SSB <- ssb(biol[,metricsPeriod])
      SSB <- drop(SSB)
      
      yieldFieldName  <- paste0('yield_',eval.period[[idxMetric]][1])
      SSBFieldName    <- paste0('SSB_',eval.period[[idxMetric]][1])
      riskFieldName   <- paste0('risk_',eval.period[[idxMetric]][1])
      IAVFieldName    <- paste0('IAV_',eval.period[[idxMetric]][1])
      
      # yield
      catchQuant  <- apply(drop(biol[,metricsPeriod]@catch), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
      
      eval.table[idxTable,yieldFieldName] <- mean(catchQuant['50%',])
      
      # SSB
      SSBQuant  <- apply(SSB, 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
      
      eval.table[idxTable,SSBFieldName] <- mean(SSBQuant['50%',])
      
      # risk of SSB < Blim
      SSB_riskMat <- array(FALSE,dim=dim(SSB))
      SSB_bool    <- array(FALSE,dim=c(1,nits))
      
      
      for(idxIter in 1:nits){
        # store value per year
        SSB_riskMat[which(SSB[,idxIter] < referencePoints$Blim),idxIter] <- TRUE
        
        # TRUE/FALSE for each iteration
        if(length(which(SSB[,idxIter] < referencePoints$Blim)!=0))
          SSB_bool[idxIter] <- TRUE
      }
      
      SSB_prob <- array(NA,dim=c(1,length(metricsPeriod)))
      
      for(idxProb in 1:length(metricsPeriod)){
        SSB_prob[idxProb] <- length(which(SSB_riskMat[idxProb,] == TRUE))/nits
      }
      
      eval.table[idxTable,riskFieldName] <- max(SSB_prob) # risk 3 here
      
      # IAV
      IAVMat  <- apply(drop(biol[,metricsPeriod]@catch), 2, diff, na.rm=TRUE) # get the differences between years
      IAVMat  <- abs(IAVMat/drop(biol[,metricsPeriod[1:(length(metricsPeriod)-1)]]@catch))# difference relative to previous year
      
      IAVQuant    <- apply(IAVMat, 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
      
      eval.table[idxTable,IAVFieldName] <- mean(IAVQuant['50%',])
      
    }
    
  }
}

