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
path <- 'E:/wk_WKNSMSE_her.27.3a47d/R'
#path <- 'D:/git/wk_WKNSMSE_her.27.3a47d/R'
assessment_name   <- "NSAS_WKNSMSE2018"
try(setwd(path),silent=TRUE)

# paths to different subfolders
dataPath      <- file.path(".","data/")
outPath       <- file.path(".","results/")
scriptPath    <- file.path(".","side_scripts/")
functionPath  <- file.path(".","functions/")

#-------------------------------------------------------------------------------
# 2) Load objects
#-------------------------------------------------------------------------------

nits <- 1000
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

  plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='catch at age')#,',ylim=c(min(a$lbnd),max(min(a$ubnd))),main='',)
  lines(years, Cwt[2,], type="l",lwd=2)
  polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)
}

################################################################################
# catch weight at age multi-fleet
################################################################################
plot(iter(fishery@landings.wt[1,,2],1))

################################################################################
# stock weight at age
################################################################################
plotQuant <- biol@stock.wt
years <- an(colnames(plotQuant))
par(mfrow=c(3,3))
for(idxAge in 1:dim(plotQuant)[1]){
  Cwt <- apply(drop(plotQuant[idxAge,]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
  
  plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='catch at age')#,',ylim=c(min(a$lbnd),max(min(a$ubnd))),main='',)
  lines(years, Cwt[2,], type="l",lwd=2)
  polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)
}

################################################################################
# maturity at age
################################################################################
plotQuant <- biol@mat
years <- an(colnames(plotQuant))
par(mfrow=c(3,3))
for(idxAge in 1:dim(plotQuant)[1]){
  Cwt <- apply(drop(plotQuant[idxAge,]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
  
  plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='catch at age')#,',ylim=c(min(a$lbnd),max(min(a$ubnd))),main='',)
  lines(years, Cwt[2,], type="l",lwd=2)
  polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)
}

################################################################################
# M
################################################################################
plotQuant <- biol@m
years <- an(colnames(plotQuant))
par(mfrow=c(3,3))
for(idxAge in 1:dim(plotQuant)[1]){
  Cwt <- apply(drop(plotQuant[idxAge,]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
  
  plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='catch at age')#,',ylim=c(min(a$lbnd),max(min(a$ubnd))),main='',)
  lines(years, Cwt[2,], type="l",lwd=2)
  polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)
}

