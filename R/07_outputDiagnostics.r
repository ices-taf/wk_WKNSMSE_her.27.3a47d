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
path <- 'D:/git/wk_WKNSMSE_her.27.3a47d/R'
assessment_name   <- "NSAS_WKNSMSE2018"
try(setwd(path),silent=TRUE)

# paths to different subfolders
dataPath      <- file.path(".","data/")
outPath       <- file.path(".","results/")
scriptPath    <- file.path(".","side_scripts/")
functionPath  <- file.path(".","functions/")

outputName <- 'output_diagnostics'

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

nits <- 200
# run to load
HCR   <- 'B'
IAV   <- 'AB'
BB    <- 'AB'
runFolder <- paste0('grid_','HCR_',HCR,'_IAV_',IAV,'_BB_',BB)
Ftar  <- 0.26
Btrig <- 1.4e06

runName         <- paste0("NSAS_Ftar_",Ftar,
                          "_Btrig_",Btrig,
                          "_HCR_",HCR,
                          "_TACIAV_",IAV,
                          "_BB_",BB,
                          "_",nits,"iters.RData")

load(file.path(outPath,runFolder,paste0(runName,'.RData')))




# load MSE parameters
load(file.path(outPath,paste0(assessment_name,'_parameters_MSE_',ac(nits),'.RData')))

strFleet    <- c('A','B','C','D')
nFleets     <- length(strFleet)
nAges       <- dim(biol)[1]
surveyNames <- names(surveys)

load(file.path(outPath,paste0(assessment_name,'_sf_noLAI.Rdata')))
load(file.path(outPath,paste0(assessment_name,'_mf_noLAI.Rdata')))

#-------------------------------------------------------------------------------
# 2) plotting
#-------------------------------------------------------------------------------

################################################################################
# stock time series
################################################################################

recAll    <- rec(biol)
ssbAll    <- ssb(biol)
biol@catch <- computeCatch(biol)
catchAll  <- biol@catch
fbarAll   <- fbar(biol)

#years <- an(fullPeriod[1:length(fullPeriod)-1])
years <- 2000:max(an(fullPeriod[1:length(fullPeriod)-1]))
nIndPlot  <- 3
plotSel   <- round(runif(nIndPlot,min=1,max=nits))

par(mfrow=c(2,2))
# plot recruitment
plotQuant <- recAll
Cwt <- apply(drop(plotQuant[,ac(years)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='recruitment')
lines(years, Cwt[2,], type="l",lwd=2)
polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)

for(idxPlot in 1:nIndPlot){
  lines(an(projPeriod), drop(iter(recAll[1,projPeriod],plotSel[idxPlot])), type="l",lwd=0.5)
}

lines(c(2018,2018),c(0,20e08), type="l",lwd=5)


# plot ssb
plotQuant <- ssbAll
Cwt <- apply(drop(plotQuant[,ac(years)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='SSB')
lines(years, Cwt[2,], type="l",lwd=2)
polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)

for(idxPlot in 1:nIndPlot){
  lines(an(projPeriod), drop(iter(ssbAll[1,projPeriod],plotSel[idxPlot])), type="l",lwd=0.5)
}

lines(c(2018,2018),c(0,20e08), type="l",lwd=5)

# plot catch
plotQuant <- catchAll
plotQuant[is.na(plotQuant)] <- 0
Cwt <- apply(drop(plotQuant[,ac(years)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='catch')
lines(years, Cwt[2,], type="l",lwd=2)
polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)

for(idxPlot in 1:nIndPlot){
  lines(an(projPeriod), drop(iter(catchAll[1,projPeriod],plotSel[idxPlot])), type="l",lwd=0.5)
}

lines(c(2018,2018),c(0,20e08), type="l",lwd=5)

# plot fbar
plotQuant <- fbarAll
Cwt <- apply(drop(plotQuant[,ac(years)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='F2-6')
lines(years, Cwt[2,], type="l",lwd=2)
polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)

for(idxPlot in 1:nIndPlot){
  lines(an(projPeriod), drop(iter(fbarAll[1,projPeriod],plotSel[idxPlot])), type="l",lwd=0.5)
}

lines(c(2018,2018),c(0,20e08), type="l",lwd=5)

mtext(runName, line = -2, cex=1.5, font=2,outer = TRUE)

################################################################################
# assessment biol comparison
################################################################################

par(mfrow=c(2,1))
# SSB
years <- an(projPeriod[1:(length(projPeriod)-1)])
plotQuant <- ssb(biol[,projPeriod])/ssb(stkAssessment[,projPeriod])
Cwt <- apply(drop(plotQuant[,ac(years)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='SSB')
lines(years, Cwt[2,], type="l",lwd=2)
polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)
# fbar
years <- an(projPeriod[1:(length(projPeriod)-1)])
plotQuant <- fbar(biol[,projPeriod])/fbar(stkAssessment[,projPeriod])
Cwt <- apply(drop(plotQuant[,ac(years)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='fbar')
lines(years, Cwt[2,], type="l",lwd=2)
polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)

################################################################################
# fishing selectivities
################################################################################

# A fleet
fishSel       <- as.data.frame(fishery@landings.sel[,projPeriod])
#fishSel       <- as.data.frame(sweep(fishery@landings.sel[,projPeriod],c(2:6),quantMeans(fishery@landings.sel[,projPeriod]),"/"))
fishSel$run   <- 'proj'

fishSelHist     <-as.data.frame(NSH3f.sam@harvest[,ac(2005:2017)])
#fishSelHist     <-as.data.frame(sweep(NSH3f.sam@harvest,c(2:6),quantMeans(NSH3f.sam@harvest),"/"))
fishSelHist$run <- 'hist'
fishSel$age     <- ac(fishSel$age)

fishSel <- rbind(fishSel,fishSelHist)

# fleet A
p1 <-  ggplot(fishSel[fishSel$area == 'A',],aes(x=age, y=data,fill=run)) + 
              geom_boxplot(outlier.alpha = 0.1) + ylab('F') + ylim(0,0.7) + ggtitle('F A fleet') + 
              theme(plot.title = element_text(size = 15, face = "bold"))

# fleet B
p2 <-  ggplot(fishSel[fishSel$area == 'B' | fishSel$area == 'BD',],aes(x=age, y=data,fill=run)) + 
              geom_boxplot(outlier.alpha = 0.1) + ylab('F') + ylim(0,3*1e-2) + ggtitle('F B fleet') + 
              theme(plot.title = element_text(size = 15, face = "bold"))

# fleet C
p3 <-  ggplot(fishSel[fishSel$area == 'C',],aes(x=age, y=data,fill=run)) + 
              geom_boxplot(outlier.alpha = 0.1) + ylab('F') + ylim(0,3*1e-2) + ggtitle('F C fleet') + 
              theme(plot.title = element_text(size = 15, face = "bold"))

# fleet D
p4 <-  ggplot(fishSel[fishSel$area == 'D' | fishSel$area == 'BD',],aes(x=age, y=data,fill=run)) + 
              geom_boxplot(outlier.alpha = 0.1) + ylab('F') + ylim(0,3*1e-2) + ggtitle('F D fleet') + 
              theme(plot.title = element_text(size = 15, face = "bold"))

p <- grid.arrange(p1, p2, p3, p4)


dev.off()
