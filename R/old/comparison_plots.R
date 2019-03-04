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
require(reshape2)
require(ggplot2)
library(gridExtra)
library(RColorBrewer)

# define path to directory
#path          <- "D:/Work/Herring MSE/NSAS/"
path              <- "D:/git/wk_WKNSMSE_her.27.3a47d/R/"
#path              <- "F:/WKNSMSE/wk_WKNSMSE_her.27.3a47d/R"
#path <- 'E:/wk_WKNSMSE_her.27.3a47d/R'
#path <- 'D:/git/wk_WKNSMSE_her.27.3a47d/R'
assessment_name   <- "NSAS_WKNSMSE2018"
try(setwd(path),silent=TRUE)

# paths to different subfolders
dataPath      <- file.path(".","data/")
outPath       <- file.path(".","results/")
scriptPath    <- file.path(".","side_scripts/")
functionPath  <- file.path(".","functions/")

outputName <- 'comparison_alt_model'

#PDF <- FALSE
#PNG <- ifelse(PDF,F,T)
#if(PDF) pdf(file.path(outPath,'plots',paste0(outputName,".pdf")))
#if(PNG) png(file.path(outPath,'plots',paste0(outputName,"_%02d.png")),
#            units = "px",
#            height=800,
#            width=672,
#            bg = "white")

##############################
caseName <- 'HCRA - 0.24/1.4e06'

# load object
load(file.path(outPath,'comp','NSAS_Ftar_0.24_Btrig_1400000_HCR_A_TACIAV__BB__200iters.RData.RData'))

metricsPeriod <- projPeriod[(length(projPeriod)-10):(length(projPeriod)-1)]

SSB_class1   <- ssb(biol)
fbar_class1  <- fbar(biol)
catch_class1 <- computeCatch(biol)

SSB_class <- SSB_class1[,metricsPeriod]
SSB_class <- drop(SSB_class)

SSB_riskMat <- array(FALSE,dim=dim(SSB_class))
SSB_bool    <- array(FALSE,dim=c(1,nits))

for(idxIter in 1:nits){
  # store value per year
  SSB_riskMat[which(SSB_class[,idxIter] < referencePoints$Blim),idxIter] <- TRUE
  
  # TRUE/FALSE for each iteration
  if(length(which(SSB_class[,idxIter] < referencePoints$Blim)!=0))
    SSB_bool[idxIter] <- TRUE
}

SSB_prob <- array(NA,dim=c(1,length(metricsPeriod)))

for(idxProb in 1:length(metricsPeriod)){
  SSB_prob[idxProb] <- length(which(SSB_riskMat[idxProb,] == TRUE))/nits
}

print(paste0('risk=',max(SSB_prob)))


load(file.path(outPath,'comp','NSAS_Ftar_0.24_Btrig_1400000_HCR_A_TACIAV__BB__200iters_altC.RData.RData'))

metricsPeriod <- projPeriod[(length(projPeriod)-10):(length(projPeriod)-1)]

SSB_alt1     <- ssb(biol)
fbar_alt1    <- fbar(biol)
catch_alt1   <- computeCatch(biol)

SSB_alt <- SSB_alt1[,metricsPeriod]
SSB_alt <- drop(SSB_alt)

SSB_riskMat <- array(FALSE,dim=dim(SSB_alt))
SSB_bool    <- array(FALSE,dim=c(1,nits))

for(idxIter in 1:nits){
  # store value per year
  SSB_riskMat[which(SSB_alt[,idxIter] < referencePoints$Blim),idxIter] <- TRUE
  
  # TRUE/FALSE for each iteration
  if(length(which(SSB_alt[,idxIter] < referencePoints$Blim)!=0))
    SSB_bool[idxIter] <- TRUE
}

SSB_prob <- array(NA,dim=c(1,length(metricsPeriod)))

for(idxProb in 1:length(metricsPeriod)){
  SSB_prob[idxProb] <- length(which(SSB_riskMat[idxProb,] == TRUE))/nits
}

print(paste0('risk=',max(SSB_prob)))

### plotting
years <- 2010:max(an(fullPeriod[1:length(fullPeriod)-1]))
nIndPlot  <- 3
plotSel   <- round(runif(nIndPlot,min=1,max=nits))

# SSB
plotQuant   <- SSB_class1
plotQuant2  <- SSB_alt1
Cwt <- apply(drop(plotQuant[,ac(years)]), 1, quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE)
Cwt2 <- apply(drop(plotQuant2[,ac(years)]), 1, quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE)

plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='SSB')
lines(years, Cwt[2,], type="l",lwd=2)
polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)

lines(years, Cwt2[2,], type="l",lwd=2)
polygon(c(years,rev(years)),c(Cwt2[1,],rev(Cwt2[3,])),col=rgb(0,0,1,0.5),lty=0)
legend("topright",
       c("default","alternative"),
       col=c(rgb(1,0,0,1),rgb(0,0,1,1)),lty=1,lwd=3)

#mtext(caseName, line = -2, cex=1.5, font=2,outer = TRUE)

lines(c(2017,2017),c(0,20e08), type="l",lwd=5)

# fbar
plotQuant   <- fbar_class1
plotQuant2  <- fbar_alt1
Cwt <- apply(drop(plotQuant[,ac(years)]), 1, quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE)
Cwt2 <- apply(drop(plotQuant2[,ac(years)]), 1, quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE)

plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='fbar')
lines(years, Cwt[2,], type="l",lwd=2)
polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)

lines(years, Cwt2[2,], type="l",lwd=2)
polygon(c(years,rev(years)),c(Cwt2[1,],rev(Cwt2[3,])),col=rgb(0,0,1,0.5),lty=0)
legend("topright",
       c("default","alternative"),
       col=c(rgb(1,0,0,1),rgb(0,0,1,1)),lty=1,lwd=3)

#mtext(caseName, line = -2, cex=1.5, font=2,outer = TRUE)


lines(c(2017,2017),c(0,20e08), type="l",lwd=5)




