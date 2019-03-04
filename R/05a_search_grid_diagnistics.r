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
require(reshape2)
require(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(grid)

# define path to directory
#path          <- "D:/Work/Herring MSE/NSAS/"
#path              <- "D:/git/wk_WKNSMSE_her.27.3a47d/R/"
#path              <- "F:/WKNSMSE/wk_WKNSMSE_her.27.3a47d/R"
#path <- 'E:/git/wk_WKNSMSE_her.27.3a47d/R'
path <- 'D:/git/wk_WKNSMSE_her.27.3a47d/R'
assessment_name   <- "NSAS_WKNSMSE2018"
try(setwd(path),silent=TRUE)

# paths to different subfolders
dataPath      <- file.path(".","data/")
outPath       <- file.path(".","results/")
scriptPath    <- file.path(".","side_scripts/")
functionPath  <- file.path(".","functions/")

outputName <- 'grid_search'

PDF <- FALSE
PNG <- ifelse(PDF,F,T)
if(PDF) pdf(file.path(outPath,'plots',paste0(outputName,".pdf")))
if(PNG) png(file.path(outPath,'plots',paste0(outputName,"_%02d.png")),
            units = "px",
            height=800,
            width=672,
            bg = "white")

#-------------------------------------------------------------------------------
# 2) plotting grid search for different cases
#-------------------------------------------------------------------------------

#list_HCR <- c('A','B')
list_HCR <- c('A',
              'B',
              'B_IAV_AB_BB_AB',
              'A_IAV_AB_BB_AB',
              'A_IAV_A_BB_A')

#'B',
#'B_IAV_AB_BB_AB',
#'A_IAV_AB_BB_AB',
#'A_IAV_A_BB_A'

for(idxHCR in 1:length(list_HCR)){

# choose HCR combination to plot
#HCRPlot <- 'A_IAV_A_BB_A'
#HCRPlot <- 'A_IAV_AB_BB_AB'
#HCRPlot <- 'A'
#HCRPlot <- 'B_IAV_A_BB_A'
#HCRPlot <- 'B_IAV_AB_BB_AB'
#HCRPlot <- 'B'

HCRPlot <- list_HCR[idxHCR]

fileList <- list.files(file.path(outPath,paste0('grid_HCR_',HCRPlot)),pattern = "\\.RData$")

load(file.path(outPath,paste0('grid_HCR_',HCRPlot),fileList[1])) 

Ftar      <- array(NA, dim=c(1,length(fileList)))
Btrig     <- array(NA, dim=c(1,length(fileList)))
LTY       <- array(NA, dim=c(1,length(fileList)))
LTR       <- array(NA, dim=c(1,length(fileList)))
LTIAV     <- array(NA, dim=c(1,length(fileList)))
LTFbar    <- array(NA, dim=c(1,length(fileList)))
idxFile   <- 1

metricsPeriod <- projPeriod[(length(projPeriod)-10):(length(projPeriod)-1)]

for(fileName in fileList){
  print(fileName)
  load(file.path(outPath,paste0('grid_HCR_',HCRPlot),fileName)) 
  
  
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
  
  SSB_prob <- array(NA,dim=c(1,length(metricsPeriod)))
  
  for(idxProb in 1:length(metricsPeriod)){
    SSB_prob[idxProb] <- length(which(SSB_riskMat[idxProb,] == TRUE))/nits
  }
  
  # Fbar
  fbarMat <- drop(fbar(biol[,metricsPeriod]))
  
  fbarQuant <- apply(fbarMat, 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

  # IAV
  IAVMat  <- apply(drop(biol[,metricsPeriod]@catch), 2, diff, na.rm=TRUE) # get the differences between years
  IAVMat  <- abs(IAVMat/drop(biol[,metricsPeriod[1:(length(metricsPeriod)-1)]]@catch))# difference relative to previous year

  IAVQuant    <- apply(IAVMat, 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
  
  # LTY
  catchQuant  <- apply(drop(biol[,metricsPeriod]@catch), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
  
  # store value for each file
  LTR[idxFile]      <- max(SSB_prob) #length(which(SSB_bool))/nits
  LTY[idxFile]      <- mean(catchQuant['50%',])
  LTIAV[idxFile]    <- mean(IAVQuant['50%',])
  LTFbar[idxFile]   <- mean(fbarQuant['50%',])
  
  
  idxFile <- idxFile + 1
}

FtarUnique    <- unique(t(Ftar))
FtarUnique    <- sort(FtarUnique)
BtrigUnique   <- unique(t(Btrig))
BtrigUnique   <- sort(BtrigUnique)

LTYMat      <- array(NA, dim=c(length(FtarUnique),length(BtrigUnique)),dimnames = list(FtarUnique,BtrigUnique)) # long term yield
IAVMat      <- array(NA, dim=c(length(FtarUnique),length(BtrigUnique)),dimnames = list(FtarUnique,BtrigUnique)) # IAV
LTRMat      <- array(NA, dim=c(length(FtarUnique),length(BtrigUnique)),dimnames = list(FtarUnique,BtrigUnique)) # long term risk
LTFbarMat   <- array(NA, dim=c(length(FtarUnique),length(BtrigUnique)),dimnames = list(FtarUnique,BtrigUnique)) # long term risk

for(idxFtar in 1:length(FtarUnique)){
  for(idxBtrig in 1:length(BtrigUnique)){
    idxMatch <- which((Ftar %in% FtarUnique[idxFtar]) & (Btrig %in% BtrigUnique[idxBtrig]))
    
    if(length(idxMatch)!=0){
      LTYMat[idxFtar,idxBtrig]    <- LTY[idxMatch]
      IAVMat[idxFtar,idxBtrig]    <- abs(LTIAV[idxMatch])
      LTRMat[idxFtar,idxBtrig]    <- LTR[idxMatch]
      LTFbarMat[idxFtar,idxBtrig] <- LTFbar[idxMatch]
    }
  }
}

################## LTR ##################

# Turn data into data.frame for writing
LTR_vec             <- as.data.frame(as.table(LTRMat))
colnames(LTR_vec)   <- c("Ftarget","Btrigger","Risk")
LTR_vec$Ftarget     <- as.numeric(as.character(LTR_vec$Ftarget))
LTR_vec$Btrigger    <- as.numeric(as.character(LTR_vec$Btrigger))

write.table(LTR_vec,
            file.path(outPath,paste0('grid_HCR_',HCRPlot),paste0(outputName,'_',HCRPlot,'_LTR_vec','.csv')),
            sep = ",",
            row.names = FALSE)

################## LTY ##################

#- LTY - Turn data into data.frame for writing
LTY_vec             <- as.data.frame(as.table(LTYMat))
colnames(LTY_vec)   <- c("Ftarget","Btrigger","yield")
LTY_vec$Ftarget     <- as.numeric(as.character(LTY_vec$Ftarget))
LTY_vec$Btrigger    <- as.numeric(as.character(LTY_vec$Btrigger))

write.table(LTY_vec,
            file.path(outPath,paste0('grid_HCR_',HCRPlot),paste0(outputName,'_',HCRPlot,'_LTY_vec','.csv')),
            sep = ",",
            row.names = FALSE)

################## IAV ##################

#- IAV - Turn data into data.frame for writing
IAV_vec            <- as.data.frame(as.table(IAVMat))
colnames(IAV_vec)  <- c("Ftarget","Btrigger","IAV")
IAV_vec$Ftarget    <- as.numeric(as.character(IAV_vec$Ftarget))
IAV_vec$Btrigger   <- as.numeric(as.character(IAV_vec$Btrigger))

write.table(IAV_vec,
            file.path(outPath,paste0('grid_HCR_',HCRPlot),paste0(outputName,'_',HCRPlot,'_IAV_vec','.csv')),
            sep = ",",
            row.names = FALSE)

################## Fbar ##################

#- Fbar - Turn data into data.frame for writing
Fbar_vec            <- as.data.frame(as.table(LTFbarMat))
colnames(IAV_vec)   <- c("Ftarget","Btrigger","Fbar")
IAV_vec$Ftarget     <- as.numeric(as.character(IAV_vec$Ftarget))
IAV_vec$Btrigger    <- as.numeric(as.character(IAV_vec$Btrigger))

write.table(IAV_vec,
            file.path(outPath,paste0('grid_HCR_',HCRPlot),paste0(outputName,'_',HCRPlot,'_Fbar_vec','.csv')),
            sep = ",",
            row.names = FALSE)


################### write tables

if(idxHCR == 1){
  # write title
  write.table(paste0('grid_HCR_',HCRPlot),
              file.path(outPath,paste0(outputName,'.csv')),
              sep = ",",
              col.names=FALSE,
              row.names = FALSE,
              quote = FALSE)
}else{
  # write title
  write.table(paste0('grid_HCR_',HCRPlot),
              file.path(outPath,paste0(outputName,'.csv')),
              sep = ",",
              col.names=FALSE,
              row.names = FALSE,
              quote = FALSE,
              append=TRUE)
}

# write LTY results
BtrigString <- c('LTY',BtrigUnique)
write.table(rbind(BtrigString,cbind(FtarUnique,LTYMat)),
            file.path(outPath,paste0(outputName,'.csv')),
            sep = ",",
            col.names=FALSE,
            row.names = FALSE,
            quote = FALSE,
            append=TRUE)

# write IAV results
BtrigString <- c('IAV',BtrigUnique)
write.table(rbind(BtrigString,cbind(FtarUnique,IAVMat)),
            file.path(outPath,paste0(outputName,'.csv')),
            sep = ",",
            col.names=FALSE,
            row.names = FALSE,
            quote = FALSE,
            append=TRUE)

# write LTR results
BtrigString <- c('LTR',BtrigUnique)
write.table(rbind(BtrigString,cbind(FtarUnique,LTRMat)),
            file.path(outPath,paste0(outputName,'.csv')),
            sep = ",",
            col.names=FALSE,
            row.names = FALSE,
            quote = FALSE,
            append=TRUE)

# write fbar results
BtrigString <- c('Fbar',BtrigUnique)
write.table(rbind(BtrigString,cbind(FtarUnique,LTFbarMat)),
            file.path(outPath,paste0(outputName,'.csv')),
            sep = ",",
            col.names=FALSE,
            row.names = FALSE,
            quote = FALSE,
            append=TRUE)

#write.table(LTYMat,file.path(outPath,paste0('grid_HCR_',HCRPlot,'_LTY.csv')),sep = ",",col.names=NA,append=TRUE)
#write.table(IAVMat,file.path(outPath,paste0('grid_HCR_',HCRPlot,'_IAV.csv')),sep = ",",col.names=NA,append=TRUE)
#write.table(LTRMat,file.path(outPath,paste0('grid_HCR_',HCRPlot,'_LTR.csv')),sep = ",",col.names=NA,append=TRUE)

################### Plotting

# create data frame for plotting
plotMat <- as.data.frame(cbind(melt(LTYMat),melt(abs(IAVMat))$value,melt(LTRMat)$value,melt(LTFbarMat)$value))
colnames(plotMat)   <- c("Ftarget","Btrigger","LTY","IAV","LTR",'Fbar')
plotRowNames        <- ac(round(plotMat$LTY))
plotRowNames[which(is.na(plotRowNames))] <- paste0(plotRowNames[which(is.na(plotRowNames))],ac(1:length(which(is.na(plotRowNames)))))
rownames(plotMat)   <- plotRowNames

myPalette <- colorRampPalette(brewer.pal(11, "RdYlGn"))

# long term yield
plotLabelsLTY    <- ac(round(plotMat$LTY))
treshold_low  <- 3.2*1e5
treshold_high <- 3.8*1e5
plotMat$LTY[plotMat$LTY < treshold_low]   <- treshold_low
plotMat$LTY[plotMat$LTY > treshold_high]  <- treshold_high

p1 <- ggplot(plotMat,aes(x=plotMat$Btrigger, y=plotMat$Ftarget, fill=plotMat$LTY,label=plotLabelsLTY))
p1 <- p1 + geom_tile(aes(x=plotMat$Btrigger, y=plotMat$Ftarget, fill=plotMat$LTY)) 
p1 <- p1 + scale_fill_gradientn(name='Long term yield',colours = myPalette(4),
                                limits=c(treshold_low,treshold_high),na.value="white")
p1 <- p1 + xlab('Btrigger') + ylab('Ftarget')
p1 <- p1 + scale_x_continuous(breaks=BtrigUnique) + scale_y_continuous(breaks=FtarUnique)
p1 <- p1 + geom_text()
p1 <- p1 + theme( panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  panel.border = element_blank())

# IAV
plotLabelsIAV    <- ac(round(plotMat$IAV*1e4)/1e4)
treshold_low  <- 0.10
treshold_high <- 0.25
plotMat$IAV[plotMat$IAV < treshold_low]   <- treshold_low
plotMat$IAV[plotMat$IAV > treshold_high]  <- treshold_high

p2 <- ggplot(plotMat,aes(x=plotMat$Btrigger, y=plotMat$Ftarget, fill=plotMat$IAV,label=plotLabelsIAV))
p2 <- p2 + geom_tile(aes(x=plotMat$Btrigger, y=plotMat$Ftarget, fill=plotMat$IAV)) 
p2 <- p2 + scale_fill_gradientn(name='IAV',colours = rev(myPalette(4)),
                                limits=c(treshold_low,treshold_high),na.value="white")
p2 <- p2 + xlab('Btrigger') + ylab('Ftarget')
p2 <- p2 + scale_x_continuous(breaks=BtrigUnique) + scale_y_continuous(breaks=FtarUnique)
p2 <- p2 + geom_text()
p2 <- p2 + theme( panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  panel.border = element_blank())

# Long term risk
plotLabelsLTR    <- ac(round(plotMat$LTR*1e4)/1e4)
treshold_low  <- 0
treshold_high <- 0.25
plotMat$LTR[plotMat$LTR < treshold_low]   <- treshold_low
plotMat$LTR[plotMat$LTR > treshold_high]  <- treshold_high

p3 <- ggplot(plotMat,aes(x=plotMat$Btrigger, y=plotMat$Ftarget, fill=plotMat$LTR,label=plotLabelsLTR))
p3 <- p3 + geom_tile(aes(x=plotMat$Btrigger, y=plotMat$Ftarget, fill=plotMat$LTR)) 
p3 <- p3 + scale_fill_gradientn(name='Long term risk',colours = rev(myPalette(4)),
                                limits=c(treshold_low,treshold_high),na.value="white")
p3 <- p3 + xlab('Btrigger') + ylab('Ftarget')
p3 <- p3 + scale_x_continuous(breaks=BtrigUnique) + scale_y_continuous(breaks=FtarUnique)
p3 <- p3 + geom_text()
p3 <- p3 + theme( panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  panel.border = element_blank())

# Long term risk
plotLabelsLTR    <- ac(round(plotMat$LTR*1e2)/1e2)
treshold_low  <- 0
treshold_high <- 0.25
plotMat$LTR[plotMat$LTR < treshold_low]   <- treshold_low
plotMat$LTR[plotMat$LTR > treshold_high]  <- treshold_high

p3 <- ggplot(plotMat,aes(x=plotMat$Btrigger, y=plotMat$Ftarget, fill=plotMat$LTR,label=plotLabelsLTR))
p3 <- p3 + geom_tile(aes(x=plotMat$Btrigger, y=plotMat$Ftarget, fill=plotMat$LTR)) 
p3 <- p3 + scale_fill_gradientn(name='Long term risk',colours = rev(myPalette(4)),
                                limits=c(treshold_low,treshold_high),na.value="white")
p3 <- p3 + xlab('Btrigger') + ylab('Ftarget')
p3 <- p3 + scale_x_continuous(breaks=BtrigUnique) + scale_y_continuous(breaks=FtarUnique)
p3 <- p3 + geom_text()
p3 <- p3 + theme( panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  panel.border = element_blank())

# Fbar
plotLabelsLTFbar    <- ac(round(plotMat$Fbar*1e4)/1e4)
treshold_low  <- 0.15
treshold_high <- 0.35
plotMat$Fbar[plotMat$Fbar < treshold_low]   <- treshold_low
plotMat$Fbar[plotMat$Fbar > treshold_high]  <- treshold_high

p4 <- ggplot(plotMat,aes(x=plotMat$Btrigger, y=plotMat$Ftarget, fill=plotMat$Fbar,label=plotLabelsLTFbar))
p4 <- p4 + geom_tile(aes(x=plotMat$Btrigger, y=plotMat$Ftarget, fill=plotMat$Fbar)) 
p4 <- p4 + scale_fill_gradientn(name='Fbar',colours = rev(myPalette(4)),
                                limits=c(treshold_low,treshold_high),na.value="white")
p4 <- p4 + xlab('Btrigger') + ylab('Ftarget')
p4 <- p4 + scale_x_continuous(breaks=BtrigUnique) + scale_y_continuous(breaks=FtarUnique)
p4 <- p4 + geom_text()
p4 <- p4 + theme( panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  panel.border = element_blank())

# plot matrices
HCRstr  <- strsplit(HCRPlot,'IAV')
IAVstr  <- strsplit(HCRstr[[1]][2],'BB')
BBstr   <- strsplit(IAVstr[[1]][2],'_')[[1]][2]
IAVstr  <- strsplit(IAVstr[[1]][1],'_')[[1]][2]
HCRstr  <- strsplit(HCRstr[[1]][1],'_')[[1]][1]

if(HCRPlot =='A') HCRCode <- 'A'
if(HCRPlot == 'B') HCRCode <- 'B'
if(HCRPlot =='A_IAV_A_BB_A') HCRCode <- 'A+C'
if(HCRPlot =='B_IAV_A_BB_A') HCRCode <- 'B+C'
if(HCRPlot =='A_IAV_AB_BB_AB') HCRCode <- 'A+D'
if(HCRPlot =='B_IAV_AB_BB_AB') HCRCode <- 'B+D'

titlePlot <- paste0(HCRCode,' - IAV fleet ',IAVstr,' - BB fleet ',BBstr)

p <- grid.arrange(p1, p2, p3, p4,top = textGrob(titlePlot,gp=gpar(fontsize=20,
                                                               font=2)),ncol=1)

print(p)
}

dev.off()

