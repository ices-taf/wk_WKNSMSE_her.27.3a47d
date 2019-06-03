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

nits <- 1000

dFtar   <- 0.01
dBtrig  <- 0.1e6
#FtarSeq   <- seq(0.14,0.28,dFtar)
#BtrigSeq  <- seq(1e6,2.4e6,dBtrig)

FtarSeq   <- seq(0.19,0.26,dFtar)
BtrigSeq  <- seq(1.1e6,1.7e6,dBtrig)

outputName <- paste0('grid_search_',nits,'its')

plotIndRIsk     <- FALSE
writeIndTables  <- TRUE

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

# define string array for f01 and f26
f01         <- ac(0:1)
f26         <- ac(2:6)

#list_HCR <- c('A','B')
list_HCR <- c('A')

#list_HCR <- c('B_IAV_E_BB_E')

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
  
  fileList <- list.files(file.path(outPath,paste0('grid_HCR_',HCRPlot),paste0(nits,'its')),pattern = "\\.RData$")
  
  load(file.path(outPath,paste0('grid_HCR_',HCRPlot),paste0(nits,'its'),fileList[1])) 
  
  Ftar                <- array(NA, dim=c(1,length(fileList)))
  Btrig               <- array(NA, dim=c(1,length(fileList)))
  LTY                 <- array(NA, dim=c(1,length(fileList)))
  LTR1                <- array(NA, dim=c(1,length(fileList)))
  LTR3                <- array(NA, dim=c(1,length(fileList)))
  LTIAV               <- array(NA, dim=c(1,length(fileList)))
  LTF01               <- array(NA, dim=c(1,length(fileList)))
  LTF26               <- array(NA, dim=c(1,length(fileList)))
  LT_RiskTrend        <- array(NA, dim=c(1,length(fileList)))
  LT_RiskTrendSig     <- array(NA, dim=c(1,length(fileList)))
  LT_SSBTailTrend     <- array(NA, dim=c(1,length(fileList)))
  LT_SSBTailTrendSig  <- array(NA, dim=c(1,length(fileList)))
  idxFile   <- 1
  
  metricsPeriod <- projPeriod[(length(projPeriod)-10):(length(projPeriod)-1)]
  
  for(fileName in fileList){
    print(fileName)
    load(file.path(outPath,paste0('grid_HCR_',HCRPlot),paste0(nits,'its'),fileName)) 
    
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
    
    SSBTail     <- apply(SSB, 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
    lmFrame <- as.data.frame(t(rbind(an(metricsPeriod),SSB_prob,SSBTail[1,])))
    colnames(lmFrame) <- c('years','risk','SSB5per')
    
    # trend in annual risk
    riskLR <- lm(risk~years,data=lmFrame)
    
    # trend in SSB 5% tail
    SSBTailLR <- lm(SSB5per~years,data=lmFrame)
    
    if(plotIndRIsk == TRUE){
      
      sink(file.path(outPath,paste0('lm_',HCRPlot,'_ftar_',ftarget,'_btrig_',btrigger,'.txt')))
      print('annual risk')
      print(summary(riskLR))
      print('5% tail')
      print(summary(SSBTailLR))
      print('slope as % of average')
      print(SSBTailLR$coefficients[2]/mean(lmFrame$SSB5per))
      print('Ftar/Btrig')
      print(c(ftarget,btrigger))
      sink()
      
      
      # plot
      par(mfrow=c(2,1))
      plot(metricsPeriod,SSBTail[2,],type='l',
           ylim=c(0,max(SSBTail[3,])),
           ylab='SSB',
           xlab='years',main=paste0(HCRPlot,'_ftar_',ftarget,'_btrig_',btrigger))
      lines(metricsPeriod,SSBTail[1,],type='l',col='blue',lty=2)
      lines(metricsPeriod,SSBTail[3,],type='l',col='blue',lty=2)
      
      plot(metricsPeriod,SSB_prob,type='l',
           ylim=c(0,0.08),
           ylab='Annual risk',
           xlab='years',
           xlim=c(an(min(metricsPeriod)),an(max(metricsPeriod))))
      lines(c(2010,2040),c(0.05,0.05),lty=2,col='green')
    }
    
    
    # Fbar
    FHCR26      <- FHCR[f26,ac(2019:2038)]
    FHCR01      <- FHCR[f01,ac(2019:2038)]
    
    f01Mat <- apply(drop(FHCR01[,metricsPeriod]),c(2,3),'mean')
    f26Mat <- apply(drop(FHCR26[,metricsPeriod]),c(2,3),'mean')
    
    f01Quant <- apply(f01Mat, 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
    f26Quant <- apply(f26Mat, 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
    
    
    #
    TACBMat <- drop(TAC[,metricsPeriod,,,"B"])
    
    TACBQuant <- apply(TACBMat, 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
    
    mean(TACBQuant['50%',])
    
    mean(TAC[,,,,c("B")], na.rm=1)
    
    # IAV
    IAVMat  <- apply(drop(biol[,metricsPeriod]@catch), 2, diff, na.rm=TRUE) # get the differences between years
    IAVMat  <- abs(IAVMat/drop(biol[,metricsPeriod[1:(length(metricsPeriod)-1)]]@catch))# difference relative to previous year
    
    IAVQuant    <- apply(IAVMat, 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
    
    # LTY
    catchQuant  <- apply(drop(biol[,metricsPeriod]@catch), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
    
    # store value for each file
    LTR1[idxFile]                   <- mean(SSB_prob) #length(which(SSB_bool))/nits
    LTR3[idxFile]                   <- max(SSB_prob) #length(which(SSB_bool))/nits
    LTY[idxFile]                    <- mean(catchQuant['50%',])
    LTIAV[idxFile]                  <- mean(IAVQuant['50%',])
    LTF01[idxFile]                  <- mean(f01Quant['50%',])
    LTF26[idxFile]                  <- mean(f26Quant['50%',])
    LT_RiskTrend[idxFile]           <- summary(riskLR)$coefficients[2,1]
    LT_RiskTrendSig[idxFile]        <- ifelse(summary(riskLR)$coefficients[2,4] > 0.05,0,1) # 0 is not significant. 1 is significant
    LT_SSBTailTrend[idxFile]        <- summary(SSBTailLR)$coefficients[2,1]
    LT_SSBTailTrendSig[idxFile]     <- ifelse(summary(SSBTailLR)$coefficients[2,4] > 0.05,0,1) # 0 is not significant. 1 is significant
    
    
    idxFile <- idxFile + 1
  }
  
  FtarUnique    <- unique(t(Ftar))
  FtarUnique    <- sort(FtarUnique)
  BtrigUnique   <- unique(t(Btrig))
  BtrigUnique   <- sort(BtrigUnique)
  
  LTYMat      <- array(NA, dim=c(length(FtarSeq),length(BtrigSeq)),dimnames = list(FtarSeq,BtrigSeq)) # long term yield
  IAVMat      <- array(NA, dim=c(length(FtarSeq),length(BtrigSeq)),dimnames = list(FtarSeq,BtrigSeq)) # IAV
  LTR1Mat     <- array(NA, dim=c(length(FtarSeq),length(BtrigSeq)),dimnames = list(FtarSeq,BtrigSeq)) # long term risk
  LTR3Mat     <- array(NA, dim=c(length(FtarSeq),length(BtrigSeq)),dimnames = list(FtarSeq,BtrigSeq)) # long term risk
  LTF01Mat    <- array(NA, dim=c(length(FtarSeq),length(BtrigSeq)),dimnames = list(FtarSeq,BtrigSeq)) # long realised F01
  LTF26Mat        <- array(NA, dim=c(length(FtarSeq),length(BtrigSeq)),dimnames = list(FtarSeq,BtrigSeq)) # long realised F26
  LT_RiskTrendMat    <- array(NA, dim=c(length(FtarSeq),length(BtrigSeq)),dimnames = list(FtarSeq,BtrigSeq)) # long realised F26
  LT_RiskTrendSigMat    <- array(NA, dim=c(length(FtarSeq),length(BtrigSeq)),dimnames = list(FtarSeq,BtrigSeq)) # long realised F26
  LT_SSBTailTrendMat    <- array(NA, dim=c(length(FtarSeq),length(BtrigSeq)),dimnames = list(FtarSeq,BtrigSeq)) # long realised F26
  LT_SSBTailTrendSigMat    <- array(NA, dim=c(length(FtarSeq),length(BtrigSeq)),dimnames = list(FtarSeq,BtrigSeq)) # long realised F26
  
  for(idxFtar in 1:length(FtarUnique)){
    for(idxBtrig in 1:length(BtrigUnique)){
      idxMatch <- which((Ftar %in% FtarUnique[idxFtar]) & (Btrig %in% BtrigUnique[idxBtrig]))
      
      if(length(idxMatch)!=0){
        LTYMat[ac(FtarUnique[idxFtar]),
               ac(BtrigUnique[idxBtrig])]    <- LTY[idxMatch]
        IAVMat[ac(FtarUnique[idxFtar]),
               ac(BtrigUnique[idxBtrig])]    <- abs(LTIAV[idxMatch])
        LTR1Mat[ac(FtarUnique[idxFtar]),
                ac(BtrigUnique[idxBtrig])]   <- LTR1[idxMatch]
        LTR3Mat[ac(FtarUnique[idxFtar]),
                ac(BtrigUnique[idxBtrig])]   <- LTR3[idxMatch]
        LTF01Mat[ac(FtarUnique[idxFtar]),
                 ac(BtrigUnique[idxBtrig])]  <- LTF01[idxMatch]
        LTF26Mat[ac(FtarUnique[idxFtar]),
                 ac(BtrigUnique[idxBtrig])]  <- LTF26[idxMatch]
        # trend matrices
        LT_RiskTrendMat[ac(FtarUnique[idxFtar]),
                        ac(BtrigUnique[idxBtrig])]        <- LT_RiskTrend[idxMatch]
        LT_RiskTrendSigMat[ac(FtarUnique[idxFtar]),
                           ac(BtrigUnique[idxBtrig])]     <- LT_RiskTrendSig[idxMatch]
        LT_SSBTailTrendMat[ac(FtarUnique[idxFtar]),
                           ac(BtrigUnique[idxBtrig])]     <- LT_SSBTailTrend[idxMatch]
        LT_SSBTailTrendSigMat[ac(FtarUnique[idxFtar]),
                              ac(BtrigUnique[idxBtrig])]  <- LT_SSBTailTrendSig[idxMatch]
      }
    }
  }
  
  ################### Plotting
  
  # create data frame for plotting
  plotMat <- as.data.frame(cbind(melt(LTYMat),melt(abs(IAVMat))$value,melt(LTR1Mat)$value,melt(LTR3Mat)$value,melt(LTF01Mat)$value,melt(LTF26Mat)$value))
  colnames(plotMat)   <- c("Ftarget","Btrigger","LTY","IAV","LTR1","LTR3",'LTRF01','LTRF26')
  plotRowNames        <- ac(round(plotMat$LTY))
  plotRowNames[which(is.na(plotRowNames))] <- paste0(plotRowNames[which(is.na(plotRowNames))],ac(1:length(which(is.na(plotRowNames)))))
  rownames(plotMat)   <- plotRowNames
  
  #myPalette <- colorRampPalette(brewer.pal(11, "RdYlGn"))
  
  myPalette1 <- colorRampPalette(c("red", "green", "blue"),
                                  space = "rgb")
  
  myPalette2 <- colorRampPalette(c("blue", "green",'yellow'),
                                 space = "rgb")
  
  # long term yield
  plotLabelsLTY    <- ac(round(plotMat$LTY/1e03))
  treshold_low  <- 3.3*1e5#/1e03
  treshold_high <- 3.6*1e5#/1e03
  plotMat$LTY[plotMat$LTY < treshold_low]   <- treshold_low
  plotMat$LTY[plotMat$LTY > treshold_high]  <- treshold_high
  
  p1 <- ggplot(plotMat,aes(x=plotMat$Btrigger, y=plotMat$Ftarget, fill=plotMat$LTY,label=plotLabelsLTY))
  p1 <- p1 + geom_tile(aes(x=plotMat$Btrigger, y=plotMat$Ftarget, fill=plotMat$LTY)) 
  p1 <- p1 + scale_fill_gradientn(name='Long term yield',colours = myPalette2(4),
                                  limits=c(treshold_low,treshold_high),na.value="white")
  p1 <- p1 + xlab('Btrigger') + ylab('Ftarget')
  p1 <- p1 + scale_x_continuous(breaks=BtrigSeq) + scale_y_continuous(breaks=FtarSeq)
  p1 <- p1 + geom_text()
  p1 <- p1 + theme( panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(), 
                    panel.border = element_blank())
  p1 <- p1 + geom_vline(xintercept=BtrigSeq-dBtrig/2)
  p1 <- p1 + geom_hline(yintercept=FtarSeq-dFtar/2)
  
  # Long term risk 3
  plotLabelsLTR3    <- ac(round(plotMat$LTR3*1e4)/1e4)
  treshold_low  <- 0.04
  treshold_high <- 0.06
  plotMat$LTR3[plotMat$LTR3 < treshold_low]   <- treshold_low
  plotMat$LTR3[plotMat$LTR3 > treshold_high]  <- treshold_high
  
  p4 <- ggplot(plotMat,aes(x=plotMat$Btrigger, y=plotMat$Ftarget, fill=plotMat$LTR3,label=plotLabelsLTR3))
  p4 <- p4 + geom_tile(aes(x=plotMat$Btrigger, y=plotMat$Ftarget, fill=plotMat$LTR3)) 
  p4 <- p4 + scale_fill_gradientn(name='Risk 3',colours = rev(myPalette1(6)),
                                  limits=c(treshold_low,treshold_high),na.value="white")
  p4 <- p4 + xlab('Btrigger') + ylab('Ftarget')
  p4 <- p4 + scale_x_continuous(breaks=BtrigSeq) + scale_y_continuous(breaks=FtarSeq)
  p4 <- p4 + geom_text()
  p4 <- p4 + theme( panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(), 
                    panel.border = element_blank())
  p4 <- p4 + geom_vline(xintercept=BtrigSeq-dBtrig/2)
  p4 <- p4 + geom_hline(yintercept=FtarSeq-dFtar/2)
  
  # plot matrices
  HCRstr  <- strsplit(HCRPlot,'IAV')
  IAVstr  <- strsplit(HCRstr[[1]][2],'BB')
  BBstr   <- strsplit(IAVstr[[1]][2],'_')[[1]][2]
  IAVstr  <- strsplit(IAVstr[[1]][1],'_')[[1]][2]
  HCRstr  <- strsplit(HCRstr[[1]][1],'_')[[1]][1]
  
  if(HCRPlot =='A') HCRCode <- 'A'
  if(HCRPlot == 'B') HCRCode <- 'B'
  if(HCRPlot =='A_IAV_A_BB_A') HCRCode <- 'A+C'
  if(HCRPlot =='A_IAV_AB_BB_AB') HCRCode <- 'A+D'
  if(HCRPlot =='B_IAV_E_BB_E') HCRCode <- 'B+E'
  if(HCRPlot =='B_IAV__BB_E') HCRCode <- 'B+E_alt'
  
  p <- grid.arrange(p1, p4)
  
  print(p)
}

dev.off()

