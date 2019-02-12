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
#path <- 'E:/wk_WKNSMSE_her.27.3a47d/R'
path <- 'D:/git/wk_WKNSMSE_her.27.3a47d/R'
assessment_name   <- "NSAS_WKNSMSE2018"
try(setwd(path),silent=TRUE)

# paths to different subfolders
dataPath      <- file.path(".","data/")
outPath       <- file.path(".","results/")
scriptPath    <- file.path(".","side_scripts/")
functionPath  <- file.path(".","functions/")

source(file.path(functionPath,"stf_ImY.R"))
source(file.path(functionPath,"optF_ImY.R"))
source(file.path(functionPath,"stf_FcY.R"))
source(file.path(functionPath,"optF_FcY.R"))
source(file.path(functionPath,"MSE_assessment.R"))

#-------------------------------------------------------------------------------
# 2) Load objects
#-------------------------------------------------------------------------------

# load object
load(file.path(outPath,paste0(assessment_name,'_init_MSE_full.RData')))
stkAssessement.ctrl <- NSH.ctrl
biol@m.spwn[,ac(2018:2040)] <- 0.67

# load MSE parameters
load(file.path(outPath,paste0(assessment_name,'_parameters_MSE_full.RData')))

strFleet    <- c('A','B','C','D')
nFleets     <- length(strFleet)
nAges       <- dim(biol)[1]
surveyNames <- names(surveys)

load(file.path(outPath,paste0(assessment_name,'_sf_noLAI.Rdata')))

#-------------------------------------------------------------------------------
# 2) Plotting
#-------------------------------------------------------------------------------

# catch weight at age
years <- an(colnames(biol@catch.wt))
par(mfrow=c(3,3))
for(idxAge in 1:dim(biol@catch.wt)[1]){
  Cwt <- apply(drop(biol@catch.wt[idxAge,]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)

  plot(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),xlab='year',pch = ".",ylab='catch at age')#,',ylim=c(min(a$lbnd),max(min(a$ubnd))),main='',)
  lines(years, Cwt[2,], type="l",lwd=2)
  polygon(c(years,rev(years)),c(Cwt[1,],rev(Cwt[3,])),col=rgb(1,0,0,0.5),lty=0)
}



# stock weight at age
years <- an(colnames(biol@catch.wt))
par(mfrow=c(3,3))
for(idxAge in 1:dim(biol@catch.wt)[1]){
  Swt <- apply(drop(biol@stock.wt[idxAge,]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
  
  plot(c(years,rev(years)),c(Swt[1,],rev(Swt[3,])),xlab='year',pch = ".",ylab='catch at age')#,',ylim=c(min(a$lbnd),max(min(a$ubnd))),main='',)
  lines(years, Swt[2,], type="l",lwd=2)
  polygon(c(years,rev(years)),c(Swt[1,],rev(Swt[3,])),col=rgb(1,0,0,0.5),lty=0)
}


# stock weight at age
quant <- biol@stock.wt
years <- an(colnames(quant))
par(mfrow=c(3,3))
for(idxAge in 1:dim(biol@stock.wt)[1]){
  Swt <- apply(drop(biol@stock.wt[idxAge,]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
  
  plot(c(years,rev(years)),c(Swt[1,],rev(Swt[3,])),xlab='year',pch = ".",ylab='catch at age')#,',ylim=c(min(a$lbnd),max(min(a$ubnd))),main='',)
  lines(years, Swt[2,], type="l",lwd=2)
  polygon(c(years,rev(years)),c(Swt[1,],rev(Swt[3,])),col=rgb(1,0,0,0.5),lty=0)
}

# maturity at age
years <- an(colnames(biol@mat))
par(mfrow=c(3,3))
for(idxAge in 1:dim(biol@catch.wt)[1]){
  Swt <- apply(drop(biol@stock.wt[idxAge,]), 1, quantile, probs=c(0.05, 0.5, 0.95), na.rm=TRUE)
  
  plot(c(years,rev(years)),c(Swt[1,],rev(Swt[3,])),xlab='year',pch = ".",ylab='catch at age')#,',ylim=c(min(a$lbnd),max(min(a$ubnd))),main='',)
  lines(years, Swt[2,], type="l",lwd=2)
  polygon(c(years,rev(years)),c(Swt[1,],rev(Swt[3,])),col=rgb(1,0,0,0.5),lty=0)
}
