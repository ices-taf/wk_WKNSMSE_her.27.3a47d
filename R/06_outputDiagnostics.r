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
path <- 'C:/git/wk_WKNSMSE_her.27.3a47d/R'
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
HCR   <- 'A'
IAV   <- NULL
BB    <- NULL
if(length(IAV) != 0){
  runFolder <- paste0('grid_','HCR_',HCR,'_IAV_',IAV,'_BB_',BB)
}else{
  runFolder <- paste0('grid_','HCR_',HCR)
}

Ftar  <- 0.26
Btrig <- 1.4e06

runName         <- paste0("NSAS_Ftar_",Ftar,
                          "_Btrig_",Btrig,
                          "_HCR_",HCR,
                          "_TACIAV_",IAV,
                          "_BB_",BB,
                          "_",nits,"iters.RData")

load(file.path(outPath,runFolder,paste0(runName)))



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

recAll    <- rec(biol)
ssbAll    <- ssb(biol)
biol@catch <- computeCatch(biol)
catchAll  <- biol@catch
fbarAll   <- fbar(biol)

################################################################################
# stock time series
################################################################################

#years <- an(fullPeriod[1:length(fullPeriod)-1])
years <- 1980:max(an(fullPeriod[1:length(fullPeriod)-1]))
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

lines(c(2019,2019),c(0,20e08), type="l",lwd=5)


# plot ssb
plotQuant <- ssbAll
Cwt <- apply(drop(plotQuant[,ac(years)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='SSB')
lines(years, Cwt[2,], type="l",lwd=2)
polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)

for(idxPlot in 1:nIndPlot){
  lines(an(projPeriod), drop(iter(ssbAll[1,projPeriod],plotSel[idxPlot])), type="l",lwd=0.5)
}

lines(c(2019,2019),c(0,20e08), type="l",lwd=5)

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

lines(c(2019,2019),c(0,20e08), type="l",lwd=5)

# plot fbar
plotQuant <- fbarAll
Cwt <- apply(drop(plotQuant[,ac(years)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='F2-6')
lines(years, Cwt[2,], type="l",lwd=2)
polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)

for(idxPlot in 1:nIndPlot){
  lines(an(projPeriod), drop(iter(fbarAll[1,projPeriod],plotSel[idxPlot])), type="l",lwd=0.5)
}

lines(c(2019,2019),c(0,20e08), type="l",lwd=5)

mtext(runName, line = -2, cex=1.5, font=2,outer = TRUE)

################################################################################
# recruitment relationships
################################################################################

par(mfrow=c(1,2))

SSBVec <- seq(0,2e06,100000)

# Ricker
CRI <- array(NA,dim=c(length(SSBVec),length(itersRI)))

for(idx in 1:length(itersRI)){
  CRI[,idx] <- SSBVec*as.vector(paramRec['a',itersRI[idx]])*exp(-as.vector(paramRec['b',itersRI[idx]])*SSBVec)
}

plot(SSBVec,CRI[,1],type='l',ylab='Recruitment',xlab='SSB',ylim=c(0,5e07),main='Ricker')

for(idx in 2:length(itersRI)){
  lines(SSBVec,CRI[,idx])
}
panel.grid()

# Segmented regression
CSSR <- array(NA,dim=c(length(SSBVec),length(itersSR)))

for(idx in 1:length(itersSR)){
  idxSSR1 <- which(SSBVec <= as.vector(paramRec['b',itersSR[idx]]))
  
  CSSR[idxSSR1,idx]    <- as.vector(paramRec['a',itersSR[idx]])*SSBVec[idxSSR1] # SSB < b (slope)
  if(max(idxSSR1) < dim(CSSR)[1]){
    CSSR[(max(idxSSR1)+1):dim(CSSR)[1],idx]   <- replicate(dim(CSSR)[1]-length(idxSSR1),as.vector(paramRec['a',itersSR[idx]])*as.vector(paramRec['b',itersSR[idx]])) # SSB > b (plateau)
  }
}

plot(SSBVec,CSSR[,1],type='l',ylab='Recruitment',xlab='SSB',ylim=c(0,5e07),main='Segmented regression')

for(idx in 2:length(itersSR)){
  lines(SSBVec,CSSR[,idx])
}
panel.grid()

################################################################################
# recruitment scatter plot
################################################################################

par(mfrow=c(1,1))

ssbProj <- ssbAll[,projPeriod[1:(length(projPeriod)-2)]]
ssbHist <- ssbAll[,histPeriod[1:(length(histPeriod)-1)]]

reProj  <- as.vector(biol@stock.n[1,projPeriod[2:(length(projPeriod)-1)]])
recHist <- as.vector(biol@stock.n[1,histPeriod[2:length(histPeriod)]])

plot(ssbProj,reProj,xlab='SSB',ylab='Recruitment',pch =16,col='red',cex = .1)
points(ssbHist,recHist,pch =16,cex = .1)

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
