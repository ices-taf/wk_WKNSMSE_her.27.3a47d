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
library(stats)

# define path to directory
#path          <- "D:/Work/Herring MSE/NSAS/"
#path              <- "D:/git/wk_WKNSMSE_her.27.3a47d/R/"
#path              <- "F:/WKNSMSE/wk_WKNSMSE_her.27.3a47d/R"
path <- 'E:/git/wk_WKNSMSE_her.27.3a47d/R'
#path <- 'D:/git/wk_WKNSMSE_her.27.3a47d/R'
assessment_name   <- "NSAS_WKNSMSE2018"
try(setwd(path),silent=TRUE)

# paths to different subfolders
dataPath      <- file.path(".","data/")
outPath       <- file.path(".","results/")
scriptPath    <- file.path(".","side_scripts/")
functionPath  <- file.path(".","functions/")

#-------------------------------------------------------------------------------
# 2) plotting grid search for HCR A+C
#-------------------------------------------------------------------------------

HCR <- 'A_IAV_A'
#HCR <- 'B'

fileList <- list.files(file.path(outPath,paste0('grid_HCR_',HCR)))

load(file.path(outPath,paste0('grid_HCR_',HCR),fileList[1])) 

Ftar      <- array(NA, dim=c(1,length(fileList)))
Btrig     <- array(NA, dim=c(1,length(fileList)))
LTY       <- array(NA, dim=c(1,length(fileList)))
LTR       <- array(NA, dim=c(1,length(fileList)))
IAV       <- array(NA, dim=c(1,length(fileList)))
idxFile   <- 1

metricsPeriod <- projPeriod[(length(projPeriod)-10):(length(projPeriod)-1)]

for(fileName in fileList){
  print(fileName)
  load(file.path(outPath,paste0('grid_HCR_',HCR),fileName)) 
  
  
  fileNameSplit <- strsplit(fileName,'_')
  
  Ftar[idxFile]   <- as.numeric(fileNameSplit[[1]][3])
  Btrig[idxFile]  <- as.numeric(fileNameSplit[[1]][5])
  
  biol@catch    <- computeCatch(biol)
  biol@stock    <- computeStock(biol)
  biol@landings <- computeLandings(biol)
  
  # risk of SSB < Blim
  SSB <- ssb(biol[,metricsPeriod])
  SSB <- drop(SSB)
  
  SSB_riskMat <- array(FALSE,dim=dim(SSB))
  SSB_bool    <- array(FALSE,dim=c(1,nits))
  
  for(idxIter in 1:nits){
    # store value per year
    SSB_riskMat[which(SSB[,idxIter] < referencePoints$Blim),idxIter] <- TRUE
    # TRUE/FALSE for each iteration
    if(length(which(SSB[,idxIter] < referencePoints$Blim)!=0))
      SSB_bool[idxIter] <- TRUE
  }

  # IAV
  IAVMat  <- apply(drop(biol[,metricsPeriod]@catch), 2, diff, na.rm=TRUE) # get the differences between years
  IAVMat  <- IAVMat/drop(biol[,metricsPeriod[1:(length(metricsPeriod)-1)]]@catch)# difference relative to previous year

  IAVQuant    <- apply(IAVMat, 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
  
  # LTY
  catchQuant  <- apply(drop(biol[,metricsPeriod]@catch), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
  
  # store value for each file
  LTR[idxFile]      <- length(which(SSB_bool))/nits
  LTY[idxFile]      <- mean(catchQuant['50%',])
  IAV[idxFile]      <- mean(IAVQuant['50%',])
  
  
  idxFile <- idxFile + 1
}

FtarUnique    <- unique(t(Ftar))
FtarUnique    <- sort(FtarUnique)
BtrigUnique   <- unique(t(Btrig))
BtrigUnique   <- sort(BtrigUnique)

LTYMat  <- array(NA, dim=c(length(FtarUnique),length(BtrigUnique)),dimnames = list(FtarUnique,BtrigUnique)) # long term yield
IAVMat  <- array(NA, dim=c(length(FtarUnique),length(BtrigUnique)),dimnames = list(FtarUnique,BtrigUnique)) # IAV
LTRMat  <- array(NA, dim=c(length(FtarUnique),length(BtrigUnique)),dimnames = list(FtarUnique,BtrigUnique)) # long term risk

for(idxFtar in 1:length(FtarUnique)){
  for(idxBtrig in 1:length(BtrigUnique)){
    idxMatch <- which((Ftar %in% FtarUnique[idxFtar]) & (Btrig %in% BtrigUnique[idxBtrig]))
    
    if(length(idxMatch)!=0){
      LTYMat[idxFtar,idxBtrig] <- LTY[idxMatch]
      IAVMat[idxFtar,idxBtrig] <- IAV[idxMatch]
      LTRMat[idxFtar,idxBtrig] <- LTR[idxMatch]
    }
  }
}

write.table(LTYMat,file.path(outPath,paste0('grid_HCR_',HCR,'_LTY.csv')),sep = ",",col.names=NA)
write.table(IAVMat,file.path(outPath,paste0('grid_HCR_',HCR,'_IAV.csv')),sep = ",",col.names=NA)
write.table(LTRMat,file.path(outPath,paste0('grid_HCR_',HCR,'_LTR.csv')),sep = ",",col.names=NA)

################### Plotting
require(reshape2)
require(ggplot2)
library(gridExtra)
library(RColorBrewer)


# create data frame for plotting
plotMat=cbind(melt(LTYMat),melt(abs(IAVMat))$value,melt(abs(LTRMat))$value)
names(plotMat)=c("Ftarget","Btrigger","LTY","IAV","LTR")

myPalette <- colorRampPalette(brewer.pal(11, "RdYlGn"))

# long term yield
p1 <- ggplot(plotMat)
p1 <- p1 + geom_tile(aes(x=plotMat$Btrigger, y=plotMat$Ftarget, fill=plotMat$LTY)) 
p1 <- p1 + scale_fill_gradientn(name='Long term yield',colours = myPalette(4),na.value="white")
p1 <- p1 + xlab('Btrigger') + ylab('Ftarget')
p1 <- p1 + scale_x_continuous(breaks=BtrigUnique) + scale_y_continuous(breaks=FtarUnique)

# IAV
p2 <- ggplot(plotMat)
p2 <- p2 + geom_tile(aes(x=plotMat$Btrigger, y=plotMat$Ftarget, fill=plotMat$IAV)) 
p2 <- p2 + scale_fill_gradientn(name='IAV',colours = rev(myPalette(4)),na.value="white")
p2 <- p2 + xlab('Btrigger') + ylab('Ftarget')
p2 <- p2 + scale_x_continuous(breaks=BtrigUnique) + scale_y_continuous(breaks=FtarUnique)

# Long term risk
p3 <- ggplot(plotMat)
p3 <- p3 + geom_tile(aes(x=plotMat$Btrigger, y=plotMat$Ftarget, fill=plotMat$LTR)) 
p3 <- p3 + scale_fill_gradientn(name='Long term risk',colours = rev(myPalette(4)),
                                limits=c(0,0.55),na.value="white")
p3 <- p3 + xlab('Btrigger') + ylab('Ftarget')
p3 <- p3 + scale_x_continuous(breaks=BtrigUnique) + scale_y_continuous(breaks=FtarUnique)

# plot matrices
grid.arrange(p1, p2, p3)



## dump

for(idxIter in 1:nits){
  print(idxIter)
  plot(iter(biol,idxIter))
}