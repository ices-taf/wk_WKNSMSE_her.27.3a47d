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
library(grid)
library(RColorBrewer)
library(tidyr)


# define path to directory
#path          <- "D:/Work/Herring MSE/NSAS/"
#path              <- "D:/git/wk_WKNSMSE_her.27.3a47d/R/"
#path              <- "F:/WKNSMSE/wk_WKNSMSE_her.27.3a47d/R"
#path <- 'E:/wk_WKNSMSE_her.27.3a47d/R'
#path <- 'D:/git/wk_WKNSMSE_her.27.3a47d/R'
path <- 'E:/git/wk_WKNSMSE_her.27.3a47d/R'
assessment_name   <- "NSAS_WKNSMSE2018"
try(setwd(path),silent=TRUE)

# paths to different subfolders
dataPath      <- file.path(".","data/")
outPath       <- file.path(".","results/")
scriptPath    <- file.path(".","side_scripts/")
functionPath  <- file.path(".","functions/")

outputName <- 'output_diagnostics_BD0'

PDF <- FALSE
PNG <- ifelse(PDF,F,T)
if(PDF) pdf(file.path(outPath,'plots',paste0(outputName,".pdf")))
if(PNG) png(file.path(outPath,'plots',paste0(outputName,"_%02d.png")),
            units = "px",
            height=800,
            width=672,
            bg = "white")

#-------------------------------------------------------------------------------
# 2) Load objects
#-------------------------------------------------------------------------------

nits <- 1000
# run to load
HCR   <- 'A'
IAV   <- NULL
BB    <- NULL
if(length(IAV) != 0){
  runFolder <- paste0('grid_','HCR_',HCR,'_IAV_',IAV,'_BB_',BB)
}else{
  runFolder <- paste0('grid_','HCR_',HCR)
}

Ftar  <- 0.22
Btrig <- 1.4e06

runName         <- paste0("NSAS_Ftar_",Ftar,
                          "_Btrig_",Btrig,
                          "_HCR_",HCR,
                          "_TACIAV_",IAV,
                          "_BB_",BB,
                          "_",nits,"iters.RData")



# load MSE parameters
load(file.path(outPath,paste0(assessment_name,'_parameters_MSE_',ac(nits),'.RData')))

load(file.path(outPath,paste0(assessment_name,'_sf_noLAI.Rdata')))
load(file.path(outPath,paste0(assessment_name,'_mf_noLAI.Rdata')))

#-------------------------------------------------------------------------------
# 2) plotting
#-------------------------------------------------------------------------------

load('E:/git/wk_WKNSMSE_her.27.3a47d/R/results/grid_HCR_A_IAV_A_BB_A/eval_run/NSAS_Ftar_0.22_Btrig_1400000_HCR_A_TACIAV_A_BB_A_1000iters.RData')

strFleet    <- c('A','B','C','D')
nFleets     <- length(strFleet)
nAges       <- dim(biol)[1]
surveyNames <- names(surveys)

recAll1    <- rec(biol)
ssbAll1    <- ssb(biol)
biol@catch <- computeCatch(biol)
catchAll1  <- biol@catch
fbarAll1   <- fbar(biol)

load('E:/git/wk_WKNSMSE_her.27.3a47d/R/results/grid_HCR_A_IAV_A_BB_A/NSAS_Ftar_0.22_Btrig_1400000_HCR_A_TACIAV_A_BB_A_1000_noBD_iters.RData')

recAll2    <- rec(biol)
ssbAll2    <- ssb(biol)
biol@catch <- computeCatch(biol)
catchAll2  <- biol@catch
fbarAll2   <- fbar(biol)

################################################################################
# stock time series
################################################################################

#years <- an(fullPeriod[1:length(fullPeriod)-1])
years <- 1980:max(an(fullPeriod[1:length(fullPeriod)-1]))
nIndPlot  <- 3
plotSel   <- round(runif(nIndPlot,min=1,max=nits))

par(mfrow=c(4,1))
# plot recruitment
plotQuant1 <- recAll1
plotQuant2 <- recAll2
Cwt1 <- apply(drop(plotQuant1[,ac(years)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
Cwt2 <- apply(drop(plotQuant2[,ac(years)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

plot(c(years,rev(years)),c(Cwt1[1,],rev(Cwt1[3,])),xlab='year',pch = ".",ylab='Recruitment (thousands)',cex.lab=1.5)
lines(years, Cwt1[2,], type="l",lwd=2,col=rgb(1,0,0),cex.lab=1.5)
polygon(c(years,rev(years)),c(Cwt1[1,],rev(Cwt1[3,])),col=rgb(1,0,0,0.5),lty=0,cex.lab=1.5)
lines(years, Cwt2[2,], type="l",lwd=2,col=rgb(0,1,0),cex=3)
polygon(c(years,rev(years)),c(Cwt2[1,],rev(Cwt2[3,])),col=rgb(0,1,0,0.5),lty=0,cex.lab=1.5)

for(idxPlot in 1:nIndPlot){
  lines(an(projPeriod), drop(iter(recAll2[1,projPeriod],plotSel[idxPlot])), type="l",lwd=0.5,cex.lab=1.5)
}

lines(c(2019,2019),c(0,20e08), type="l",lwd=5)


# plot ssb
plotQuant1 <- ssbAll1
plotQuant2 <- ssbAll2
Cwt1 <- apply(drop(plotQuant1[,ac(years)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
Cwt2 <- apply(drop(plotQuant2[,ac(years)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

plot(c(years,rev(years)),c(Cwt1[1,],rev(Cwt1[3,])),xlab='year',pch = ".",ylab='SSB (tonnes)',cex.lab=1.5)
lines(years, Cwt1[2,], type="l",lwd=2,col=rgb(1,0,0),cex.lab=1.5)
polygon(c(years,rev(years)),c(Cwt1[1,],rev(Cwt1[3,])),col=rgb(1,0,0,0.5),lty=0,cex.lab=1.5)
lines(years, Cwt2[2,], type="l",lwd=2,col=rgb(0,1,0),cex.lab=1.5)
polygon(c(years,rev(years)),c(Cwt2[1,],rev(Cwt2[3,])),col=rgb(0,1,0,0.5),lty=0,cex.lab=1.5)

for(idxPlot in 1:nIndPlot){
  lines(an(projPeriod), drop(iter(ssbAll2[1,projPeriod],plotSel[idxPlot])), type="l",lwd=0.5,cex.lab=1.5)
}

lines(c(2019,2019),c(0,20e08), type="l",lwd=5)

# plot catch
plotQuant1 <- catchAll1
plotQuant2 <- catchAll2
Cwt1 <- apply(drop(plotQuant1[,ac(years)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
Cwt2 <- apply(drop(plotQuant2[,ac(years)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

plot(c(years,rev(years)),c(Cwt1[1,],rev(Cwt1[3,])),xlab='year',pch = ".",ylab='Catch (tonnes)',cex.lab=1.5)
lines(years, Cwt1[2,], type="l",lwd=2,col=rgb(1,0,0),cex.lab=1.5)
polygon(c(years,rev(years)),c(Cwt1[1,],rev(Cwt1[3,])),col=rgb(1,0,0,0.5),lty=0,cex.lab=1.5)
lines(years, Cwt2[2,], type="l",lwd=2,col=rgb(0,1,0),cex.lab=1.5)
polygon(c(years,rev(years)),c(Cwt2[1,],rev(Cwt2[3,])),col=rgb(0,1,0,0.5),lty=0,cex.lab=1.5)

for(idxPlot in 1:nIndPlot){
  lines(an(projPeriod), drop(iter(catchAll2[1,projPeriod],plotSel[idxPlot])), type="l",lwd=0.5,cex.lab=1.5)
}

lines(c(2019,2019),c(0,20e08), type="l",lwd=5)

# plot fbar
plotQuant1 <- fbarAll1
plotQuant2 <- fbarAll2
Cwt1 <- apply(drop(plotQuant1[,ac(years)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
Cwt2 <- apply(drop(plotQuant2[,ac(years)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

plot(c(years,rev(years)),c(Cwt1[1,],rev(Cwt1[3,])),xlab='year',pch = ".",ylab='F2-6',cex.lab=1.5)
lines(years, Cwt1[2,], type="l",lwd=2,col=rgb(1,0,0),cex.lab=1.5)
polygon(c(years,rev(years)),c(Cwt1[1,],rev(Cwt1[3,])),col=rgb(1,0,0,0.5),lty=0,cex.lab=1.5)
lines(years, Cwt2[2,], type="l",lwd=2,col=rgb(0,1,0),cex.lab=1.5)
polygon(c(years,rev(years)),c(Cwt2[1,],rev(Cwt2[3,])),col=rgb(0,1,0,0.5),lty=0,cex.lab=1.5)

for(idxPlot in 1:nIndPlot){
  lines(an(projPeriod), drop(iter(fbarAll2[1,projPeriod],plotSel[idxPlot])), type="l",lwd=0.5,cex.lab=1.5)
}

lines(c(2019,2019),c(0,20e08), type="l",lwd=5)

#mtext(runName, line = -2, cex=1.5, font=2,outer = TRUE)

################################################################################
# stock time series BD=0
################################################################################

#years <- an(fullPeriod[1:length(fullPeriod)-1])
years <- 1980:max(an(fullPeriod[1:length(fullPeriod)-1]))
nIndPlot  <- 3
plotSel   <- c(103,594,804)#round(runif(nIndPlot,min=1,max=nits))

par(mfrow=c(4,1))
# plot recruitment
plotQuant <- recAll2
Cwt <- apply(drop(plotQuant[,ac(years)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='Recruitment (thousands)',cex.lab=1.5)
lines(years, Cwt[2,], type="l",lwd=2,cex.lab=1.5)
polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0,cex.lab=1.5)

for(idxPlot in 1:nIndPlot){
  lines(an(projPeriod), drop(iter(recAll2[1,projPeriod],plotSel[idxPlot])), type="l",lwd=0.5,cex.lab=1.5)
}

lines(c(2019,2019),c(0,20e08), type="l",lwd=5)


# plot ssb
plotQuant <- ssbAll2
Cwt <- apply(drop(plotQuant[,ac(years)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='SSB (tonnes)',cex.lab=1.5)
lines(years, Cwt[2,], type="l",lwd=2,cex.lab=1.5)
polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0,cex.lab=1.5)

for(idxPlot in 1:nIndPlot){
  lines(an(projPeriod), drop(iter(ssbAll2[1,projPeriod],plotSel[idxPlot])), type="l",lwd=0.5,cex.lab=1.5)
}

lines(c(2019,2019),c(0,20e08), type="l",lwd=5)

# plot catch
plotQuant <- catchAll2
plotQuant[is.na(plotQuant)] <- 0
Cwt <- apply(drop(plotQuant[,ac(years)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='Catch (tonnes)',cex.lab=1.5)
lines(years, Cwt[2,], type="l",lwd=2,cex.lab=1.5)
polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0,cex.lab=1.5)

for(idxPlot in 1:nIndPlot){
  lines(an(projPeriod), drop(iter(catchAll2[1,projPeriod],plotSel[idxPlot])), type="l",lwd=0.5,cex.lab=1.5)
}

lines(c(2019,2019),c(0,20e08), type="l",lwd=5)

# plot fbar
plotQuant <- fbarAll2
Cwt <- apply(drop(plotQuant[,ac(years)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='F2-6',cex.lab=1.5)
lines(years, Cwt[2,], type="l",lwd=2,cex.lab=1.5)
polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0,cex.lab=1.5)

for(idxPlot in 1:nIndPlot){
  lines(an(projPeriod), drop(iter(fbarAll2[1,projPeriod],plotSel[idxPlot])), type="l",lwd=0.5,cex.lab=1.5)
}

lines(c(2019,2019),c(0,20e08), type="l",lwd=5)

#mtext(runName, line = -2, cex=1.5, font=2,outer = TRUE)



dev.off()
