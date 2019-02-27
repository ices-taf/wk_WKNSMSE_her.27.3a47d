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

outputName <- 'input_diagnostics'

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
# load object
load(file.path(outPath,paste0(assessment_name,'_init_MSE_',ac(nits),'.RData')))
stkAssessement.ctrl <- NSH.ctrl

# load MSE parameters
load(file.path(outPath,paste0(assessment_name,'_parameters_MSE_',ac(nits),'.RData')))

strFleet    <- c('A','B','C','D')
nFleets     <- length(strFleet)
nAges       <- dim(biol)[1]
surveyNames <- names(surveys)

load(file.path(outPath,paste0(assessment_name,'_sf_noLAI.Rdata')))

# 100 iteration object not here yet
load(file.path(outPath,'SplitUptakes200.RData'))

#-------------------------------------------------------------------------------
# 2) Plotting
#-------------------------------------------------------------------------------

################################################################################
# catch weight at age
################################################################################
plotQuant <- biol@catch.wt
years <- an(colnames(plotQuant))
par(mfrow=c(3,3))

for(idxAge in 1:dim(plotQuant)[1]){
  Cwt <- apply(drop(plotQuant[idxAge,]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

  plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab=paste0('age ',idxAge-1))#,',ylim=c(min(a$lbnd),max(min(a$ubnd))),main='',)
  lines(years, Cwt[2,], type="l",lwd=2)
  polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)
}
mtext("catch weight at age", line = -2, cex=1.5, font=2,outer = TRUE)

################################################################################
# catch weight at age multi-fleet
################################################################################
#plot(iter(fishery@landings.wt[1,,2],1))

################################################################################
# stock weight at age
################################################################################
plotQuant <- biol@stock.wt
years <- an(colnames(plotQuant))
par(mfrow=c(3,3))
for(idxAge in 1:dim(plotQuant)[1]){
  Cwt <- apply(drop(plotQuant[idxAge,]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
  
  plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab=paste0('age ',idxAge-1))#,',ylim=c(min(a$lbnd),max(min(a$ubnd))),main='',)
  lines(years, Cwt[2,], type="l",lwd=2)
  polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)
}
mtext("stock weight at age", line = -2, cex=1.5, font=2, outer = TRUE)

################################################################################
# maturity at age
################################################################################
plotQuant <- biol@mat
years <- an(colnames(plotQuant))
par(mfrow=c(3,1))
for(idxAge in 3:5){
  Cwt <- apply(drop(plotQuant[idxAge,]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
  
  plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab=paste0('age ',idxAge-1))#,',ylim=c(min(a$lbnd),max(min(a$ubnd))),main='',)
  lines(years, Cwt[2,], type="l",lwd=2)
  polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)
}
mtext("maturity at age", line = -2, cex=1.5, font=2, outer = TRUE)

################################################################################
# M
################################################################################
plotQuant <- biol@m
years <- an(colnames(plotQuant))
par(mfrow=c(3,3))
for(idxAge in 1:dim(plotQuant)[1]){
  Cwt <- apply(drop(plotQuant[idxAge,]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
  
  plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab=paste0('age ',idxAge-1))#,',ylim=c(min(a$lbnd),max(min(a$ubnd))),main='',)
  lines(years, Cwt[2,], type="l",lwd=2)
  polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)
}
mtext("Natural mortality", line = -2, cex=1.5, font=2, outer = TRUE)

################################################################################
# random samples survey indices
################################################################################

# survey loop
for(surveyName in surveys@names){
  nAges <- length(surveys[[surveyName]]@range[1]:surveys[[surveyName]]@range[2])
  par(mfrow=n2mfrow(nAges))
  # age loop
  for(ageIdx in surveys[[surveyName]]@range[1]:surveys[[surveyName]]@range[2]){
    a     <- NSH.tun[[surveyName]]@index[,drop=FALSE]
    years <- as.numeric(colnames(NSH.tun[[surveyName]]@index))
    plot(years,a[ac(ageIdx),],type='l',ylab=c('age',ac(ageIdx)))
    
    for(idxIter in 1:100){
      b     <-surveys[[surveyName]]@index[ac(ageIdx),,,,,idxIter]
      b     <- b[,match(years,colnames(b))]
      lines(years,b,col='green')
    }
    lines(years,a[ac(ageIdx),],type='l',ylab=c('age',ac(ageIdx)),lwd=3)
  }
  mtext(paste0('Random samples - ',surveyName), line = -2, cex=1.5, font=2, outer = TRUE)
}

################################################################################
# random samples catches
################################################################################

nAges <- dim(biol@catch.n)[1]
par(mfrow=n2mfrow(nAges))
for(ageIdx in 1:dim(biol@catch.n)[1]){
  a     <- NSH@catch.n[,drop=FALSE]
  years <- histMinYr:histMaxYr # vector the years
  plot(years,a[ac(ageIdx-1),],type='l',ylab=c('age',ac(ageIdx-1)))
  for(idxIter in 1:100){
    b     <- iter(biol@catch.n[ageIdx],idxIter)
    b     <- b[,match(years,colnames(b))]
    lines(years,b,col='green')
  }
  lines(years,a[ac(ageIdx-1),],type='l',ylab=c('age',ac(ageIdx)),lwd=3)
}
mtext(paste0('Random samples - catches'), line = -2, cex=1.5, font=2,outer = TRUE)

################################################################################
# Recruitment 
################################################################################

# residuals
par(mfrow=c(1,1))
recPeriod <- ac(2002:2016)

histResiRec <- residuals(biol.sr)[,ac(an(recPeriod)[2]:(max(an(recPeriod))))]
projResiRec <- sr.res


plotQuant   <- projResiRec
years       <- an(projPeriod)
yearsHist   <- an(colnames(histResiRec))
projData    <- apply(drop(plotQuant[,ac(years)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
histData    <- apply(drop(histResiRec[,ac(yearsHist)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

dataAll     <- cbind(histData,projData)
yearsAll    <- c(yearsHist,years)

plot(c(yearsAll,rev(yearsAll)),c(dataAll[1,],rev(dataAll[3,])),xlab='year',pch = ".",ylab='Recruitment residuals')
lines(yearsAll, dataAll[2,], type="l",lwd=2)
polygon(c(yearsAll,rev(yearsAll)),c(dataAll[1,],rev(dataAll[3,])),col=rgb(1,0,0,0.5),lty=0)

# ACF residuals
par(mfrow=c(1,1))
idxIter <- round(runif(1,min=1,max=nits))
recPeriod <- ac(2002:2016)
a<-acf(iter(residuals(biol.sr)[,ac(an(recPeriod)[2]:(max(an(recPeriod)))-1)],idxIter))
plot(a, main = "recruitment autocorrelation")

################################################################################
# Proportion of F
################################################################################
plotQuant <- catchVar[1,,,'FCprop']
years     <- an(c(fecYears,projPeriod))
Cwt       <- apply(drop(plotQuant[,ac(years)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='Proportion of total F')
lines(years, Cwt[2,], type="l",lwd=2)
polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)

mtext(paste0('Proportion of Fc'), line = -2, cex=1.5, font=2,outer = TRUE)

################################################################################
# TAC variables
################################################################################
par(mfrow=c(3,1))
# C transfer
plotQuant <- TAC_var[,,'Ctransfer']
years     <- an(projPeriod)
Cwt       <- apply(drop(plotQuant[ac(years),]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='transfer TAC C fleet')
lines(years, Cwt[2,], type="l",lwd=2)
polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)

# D split
plotQuant   <- TAC_var[,,'Dsplit']
years       <- an(projPeriod)
Cwt         <- apply(drop(plotQuant[ac(years),]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
histDsplit  <-t(replicate(3,DSplitHist[,2]))

Cwt         <- cbind(histDsplit,Cwt)
years       <- min(DSplitHist[,1]):max(years)

plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='NSAS/WBSS split D fleet')
lines(years, Cwt[2,], type="l",lwd=2)
polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)

# uptake B fleet
plotQuant   <- TAC_var[,,'Buptake']
years       <- an(projPeriod)
Cwt         <- apply(drop(plotQuant[ac(years),]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

uptakeTab   <- read.table(file.path(dataPath,'over_underfishing2017.csv'),sep = ",")
uptakeBHist <- an(as.vector(uptakeTab[2:16,3]))

plot(uptakeTab[2:16,1], uptakeBHist, type="l",lwd=2,xlim=c(2003,2040))
lines(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='uptake B fleet')
lines(years, Cwt[2,], type="l",lwd=2)

polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)

mtext(paste0('TAC variables'), line = -2, cex=1.5, font=2,outer = TRUE)

################################################################################
# Process error
################################################################################

surv <- biol@stock.n[,histPeriod]*exp(-biol@harvest[,histPeriod]-biol@m[,histPeriod]) # effectively, this is age 1 to 8 in year + 1
surv[dim(surv)[1]-1] <- surv[dim(surv)[1]-1] + surv[dim(surv)[1]]
dimnames(surv)$age <- ac(1:9)

# process error
procError <-  surv[ac(1:8),histPeriod[1:(length(histPeriod)-1)]]/ # survivors age 1 to 8 (0 to 7 in surv object) year 1948 to 2017
              biol@stock.n[ac(1:8),histPeriod[2:length(histPeriod)]] # numbers at age, age 1 to 8

par(mfrow=n2mfrow(dim(procError)[1]))

for(idxAge in 1:dim(procError)[1]){
  yearsProc <- 1970:2008
  Cwt       <- apply(drop(procError[idxAge,ac(yearsProc)]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
  Cwt2      <- apply(drop(varProccError[idxAge,]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
  
  Cwt <- cbind(Cwt,Cwt2)
  
  yearsProc <- an(colnames(Cwt))
  plot(c(yearsProc,rev(yearsProc)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab=paste0('age ',idxAge))
  lines(yearsProc, Cwt[2,], type="l",lwd=2)
  polygon(c(yearsProc,rev(yearsProc)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)
}

################################################################################
# selectivities
################################################################################

idxIter <- round(runif(1,min=1,max=nits))
nYearsWin <- 5
nWin <- length(projPeriod)/nYearsWin
sel <- drop(iter(fishery@landings.sel,idxIter))

# A fleet
par(mfrow=n2mfrow(nWin))
idxFleet <- 1

yearCurrent <- an(projPeriod[1])

for(idxWins in 1:nWin){
  plot(sel[,ac(yearCurrent),idxFleet],type="l",lwd=2,ylim=c(0,1),xlab='year',ylab='Normalized F')
  lines(sel[,ac(yearCurrent+1),idxFleet],type="l",lwd=2,col='red')
  lines(sel[,ac(yearCurrent+2),idxFleet],type="l",lwd=2,col='blue')
  lines(sel[,ac(yearCurrent+3),idxFleet],type="l",lwd=2,col='green')
  lines(sel[,ac(yearCurrent+4),idxFleet],type="l",lwd=2,col='yellow')

  legend('topleft',legend = c(ac(yearCurrent),
                               ac(yearCurrent+1),
                               ac(yearCurrent+2),
                               ac(yearCurrent+3),
                               ac(yearCurrent+4)),
         lty = 1,
         col=c('black','red','blue','green','yellow'),
         cex = 1.5)

  yearCurrent <- yearCurrent + 5
}

mtext(paste0('Selectivity A fleet'), line = -2, cex=1.5, font=2,outer = TRUE)

# BD fleet
par(mfrow=n2mfrow(nWin))
idxFleet <- 2

yearCurrent <- an(projPeriod[1])

for(idxWins in 1:nWin){
  plot(sel[,ac(yearCurrent),idxFleet],type="l",lwd=2,ylim=c(0,1),xlab='year',ylab='Normalized F')
  lines(sel[,ac(yearCurrent+1),idxFleet],type="l",lwd=2,col='red')
  lines(sel[,ac(yearCurrent+2),idxFleet],type="l",lwd=2,col='blue')
  lines(sel[,ac(yearCurrent+3),idxFleet],type="l",lwd=2,col='green')
  lines(sel[,ac(yearCurrent+4),idxFleet],type="l",lwd=2,col='yellow')
  
  legend('topright',legend = c(ac(yearCurrent),
                              ac(yearCurrent+1),
                              ac(yearCurrent+2),
                              ac(yearCurrent+3),
                              ac(yearCurrent+4)),
         lty = 1,
         col=c('black','red','blue','green','yellow'),
         cex = 1.5)
  
  yearCurrent <- yearCurrent + 5
}
mtext(paste0('Selectivity BD fleet'), line = -2, cex=1.5, font=2,outer = TRUE)

# C fleet
par(mfrow=n2mfrow(nWin))
idxFleet <- 3

yearCurrent <- an(projPeriod[1])

for(idxWins in 1:nWin){
  plot(sel[,ac(yearCurrent),idxFleet],type="l",lwd=2,ylim=c(0,1),xlab='year',ylab='Normalized F')
  lines(sel[,ac(yearCurrent+1),idxFleet],type="l",lwd=2,col='red')
  lines(sel[,ac(yearCurrent+2),idxFleet],type="l",lwd=2,col='blue')
  lines(sel[,ac(yearCurrent+3),idxFleet],type="l",lwd=2,col='green')
  lines(sel[,ac(yearCurrent+4),idxFleet],type="l",lwd=2,col='yellow')
  
  legend('topright',legend = c(ac(yearCurrent),
                              ac(yearCurrent+1),
                              ac(yearCurrent+2),
                              ac(yearCurrent+3),
                              ac(yearCurrent+4)),
         lty = 1,
         col=c('black','red','blue','green','yellow'),
         cex = 1.5)
  
  yearCurrent <- yearCurrent + 5
}

mtext(paste0('Selectivity C fleet'), line = -2, cex=1.5, font=2,outer = TRUE)

################################################################################
# selectivities over the years
################################################################################

fishSel       <- as.data.frame(fishery@landings.sel[,projPeriod])
#fishSel       <- as.data.frame(sweep(fishery@landings.sel[,projPeriod],c(2:6),quantMeans(fishery@landings.sel[,projPeriod]),"/"))
fishSel$run   <- 'proj'
fishSel$year  <- ac(fishSel$year)

fishSel <- rbind(fishSel,fishSelHist)

fleetList <- c('A','B','C','D')

for(mFleet in fleetList){

# fleet A
ageSel <- 0
p1 <-  ggplot(fishSel[fishSel$area == mFleet & fishSel$age == ageSel,],aes(x=year, y=data)) + 
              geom_boxplot(outlier.alpha = 0.1) + ylab(paste0('F',' - age ',ageSel)) + 
              ggtitle(paste0('F ',mFleet,' fleet')) + 
              theme(plot.title = element_text(size = 15, face = "bold"))

ageSel <- 1
p2 <-  ggplot(fishSel[fishSel$area == mFleet & fishSel$age == ageSel,],aes(x=year, y=data)) + 
  geom_boxplot(outlier.alpha = 0.1) + ylab(paste0('F',' - age ',ageSel)) + 
  theme(plot.title = element_text(size = 15, face = "bold"))

ageSel <- 2
p3 <-  ggplot(fishSel[fishSel$area == mFleet & fishSel$age == ageSel,],aes(x=year, y=data)) + 
  geom_boxplot(outlier.alpha = 0.1) + ylab(paste0('F',' - age ',ageSel)) + 
  theme(plot.title = element_text(size = 15, face = "bold"))

ageSel <- 3
p4 <-  ggplot(fishSel[fishSel$area == 'A' & fishSel$age == ageSel,],aes(x=year, y=data)) + 
  geom_boxplot(outlier.alpha = 0.1) + ylab(paste0('F',' - age ',ageSel)) + 
  theme(plot.title = element_text(size = 15, face = "bold"))

ageSel <- 4
p5 <-  ggplot(fishSel[fishSel$area == mFleet & fishSel$age == ageSel,],aes(x=year, y=data)) + 
  geom_boxplot(outlier.alpha = 0.1) + ylab(paste0('F',' - age ',ageSel)) + 
  theme(plot.title = element_text(size = 15, face = "bold"))

ageSel <- 5
p6 <-  ggplot(fishSel[fishSel$area == mFleet & fishSel$age == ageSel,],aes(x=year, y=data)) + 
  geom_boxplot(outlier.alpha = 0.1) + ylab(paste0('F',' - age ',ageSel)) + 
  theme(plot.title = element_text(size = 15, face = "bold"))

ageSel <- 6
p7 <-  ggplot(fishSel[fishSel$area == mFleet & fishSel$age == ageSel,],aes(x=year, y=data)) + 
  geom_boxplot(outlier.alpha = 0.1) + ylab(paste0('F',' - age ',ageSel)) + 
  theme(plot.title = element_text(size = 15, face = "bold"))

ageSel <- 7
p8 <-  ggplot(fishSel[fishSel$area == mFleet & fishSel$age == ageSel,],aes(x=year, y=data)) + 
  geom_boxplot(outlier.alpha = 0.1) + ylab(paste0('F',' - age ',ageSel)) + 
  theme(plot.title = element_text(size = 15, face = "bold"))

ageSel <- 8
p9 <-  ggplot(fishSel[fishSel$area == mFleet & fishSel$age == ageSel,],aes(x=year, y=data)) + 
  geom_boxplot(outlier.alpha = 0.1) + ylab(paste0('F',' - age ',ageSel)) + 
  theme(plot.title = element_text(size = 15, face = "bold"))


p <- grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9)
print(p)
}


dev.off()

