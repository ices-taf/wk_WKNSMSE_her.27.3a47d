#-------------------------------------------------------------------------------
# WKNSMSE
#
# Author: Benoit Berges
#         WMR, The Netherland
# email: benoit.berges@wur.nl
#
# script to test the effect of dropping LAI indices
#
# Date: 2018/12/11
#
# Build for R3.5.1, 64bits
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Running assessment
#-------------------------------------------------------------------------------

rm(list=ls())

library(FLSAM)
library(FLEDA)
#library(R.utils)

#path          <- "D:/Work/Herring MSE/NSAS/"
path              <- "D:/git/wk_WKNSMSE_her.27.3a47d/R/"
#path              <- "F:/WKNSMSE/wk_WKNSMSE_her.27.3a47d/R"
#path <- '/home/berge057/ICES/wk_WKNSMSE_her.27.3a47d/R/'
path <- 'E:/wk_WKNSMSE_her.27.3a47d/R'
assessment_name_noLAI   <- "NSAS_WKNSMSE2018"
try(setwd(path),silent=TRUE)

# paths to different subfolders
dataPath      <- file.path(".","data/")
outPath       <- file.path(".","results/")
scriptPath    <- file.path(".","side_scripts/")

# running assessment with no LAI
source(file.path(scriptPath,"setupAssessmentObjects_sf_noLAI.R"))
source(file.path(scriptPath,"setupControlObject_sf_noLAI.R"))

NSHnoLAI.ctrl   <- NSH.ctrl
NSHnoLAI.tun    <- NSH.tun
NSHnoLAI        <- NSH

NSHnoLAI.sam            <- FLSAM(NSHnoLAI,
                                 NSHnoLAI.tun,
                                 NSHnoLAI.ctrl)
NSHnoLAI@stock.n           <- NSHnoLAI.sam@stock.n[,ac(range(NSH)["minyear"]:range(NSH)["maxyear"])]
NSHnoLAI@harvest           <- NSHnoLAI.sam@harvest[,ac(range(NSH)["minyear"]:range(NSH)["maxyear"])]

NSHnoLAI.tun@residuals  <- F
n.retro.years           <- 7
NSHnoLAI.retro          <- retro( NSHnoLAI, 
                                  NSHnoLAI.tun,
                                  NSHnoLAI.ctrl,
                                  retro=n.retro.years)

#-------------------------------------------------------------------
# saving 
#-------------------------------------------------------------------

save(NSHnoLAI,
     NSHnoLAI.tun,
     NSHnoLAI.ctrl,
     NSHnoLAI.sam,
     NSHnoLAI.retro,
     file=file.path(outPath,paste0(assessment_name_noLAI,'_sf_retro.Rdata')))

# modify object names so it fits those in the other scripts
NSH       <- NSHnoLAI
NSH.tun   <- NSHnoLAI.tun
NSH.ctrl  <- NSHnoLAI.ctrl
NSH.sam   <- NSHnoLAI.sam
save(NSH,
     NSHnoLAI.tun,
     NSHnoLAI.ctrl,
     NSHnoLAI.sam,
     file=file.path(outPath,paste0(assessment_name_noLAI,'_sf.Rdata')))

#-------------------------------------------------------------------
# loading assessments with and without LAI index
#-------------------------------------------------------------------

load(file.path(outPath,paste0('NSAS_WKNSMSE2018_sf_retro.Rdata')))

load(file.path(outPath,paste0('NSAS_WKNSMSE2018_sf_noLAI_retro.Rdata')))

#-------------------------------------------------------------------
# Compute mohn's rhos
#-------------------------------------------------------------------

# SSB
SSB_mrNoLAI <- mohns.rho(NSHnoLAI.retro,span=n.retro.years,ref.year=2017,type="ssb")
mean(SSB_mrNoLAI$rho)
SSB_mrLAI   <- mohns.rho(NSH.retro,span=n.retro.years,ref.year=2017,type="ssb")
mean(SSB_mrLAI$rho)
# fbar
fbar_mrNoLAI <- mohns.rho(NSHnoLAI.retro,span=n.retro.years,ref.year=2017,type="fbar")
mean(fbar_mrNoLAI$rho)
fbar_mrLAI   <- mohns.rho(NSH.retro,span=n.retro.years,ref.year=2017,type="fbar")
mean(fbar_mrLAI$rho)
# recruitment
rec_mrNoLAI <- mohns.rho(NSHnoLAI.retro,span=n.retro.years,ref.year=2017,type="rec")
mean(rec_mrNoLAI$rho)
rec_mrLAI   <- mohns.rho(NSH.retro,span=n.retro.years,ref.year=2017,type="rec")
mean(rec_mrLAI$rho)

#-------------------------------------------------------------------
# plotting
#-------------------------------------------------------------------

ssbNoLAI      <- ssb(NSHnoLAI.sam)
recruitNoLAI  <- rec(NSHnoLAI.sam)
fbarNoLAI     <- fbar(NSHnoLAI.sam)

ssbLAI      <- ssb(NSH.sam)
recruitLAI  <- rec(NSH.sam)
fbarLAI     <- fbar(NSH.sam)

yrs   <- ssbLAI$year

par(mfrow=c(3,1))
# SSB
# with LAI
plot(yrs,
     ssbLAI$value,
     xlab='year',ylab='SSB',main='SSB',type='l',col=rgb(0,0,1,1),lwd=3)
polygon(c(yrs,rev(yrs)),
        c(ssbLAI$lbnd,rev(ssbLAI$ubnd)),
        col=rgb(0,0,1,0.1),lty=0)
# without LAI
lines(yrs,
      ssbNoLAI$value,
      col=rgb(1,0,0,1),lwd=3)
polygon(c(yrs,rev(yrs)),
        c(ssbNoLAI$lbnd,rev(ssbNoLAI$ubnd)),
        col=rgb(1,0,0,0.1),lty=0)

legend("topright",
       c("with LAI","without LAI"),
       col=c(rgb(0,0,1,1),rgb(1,0,0,1)),lty=1,lwd=3)

# Fbar
# with LAI
plot(yrs,
     fbarLAI$value,
     xlab='year',ylab='Fbar',main='Fbar',type='l',col=rgb(0,0,1,1),lwd=3)
polygon(c(yrs,rev(yrs)),
        c(fbarLAI$lbnd,rev(fbarLAI$ubnd)),
        col=rgb(0,0,1,0.1),lty=0)
# without LAI
lines(yrs,
      fbarNoLAI$value,
      col=rgb(1,0,0,1),lwd=3)
polygon(c(yrs,rev(yrs)),
        c(fbarNoLAI$lbnd,rev(fbarNoLAI$ubnd)),
        col=rgb(1,0,0,0.1),lty=0)

# Recruitment
# with LAI
plot(yrs,
     recruitLAI$value,
     xlab='year',ylab='rec',main='rec',type='l',col=rgb(0,0,1,1),lwd=3)
polygon(c(yrs,rev(yrs)),
        c(recruitLAI$lbnd,rev(recruitLAI$ubnd)),
        col=rgb(0,0,1,0.1),lty=0)
# without LAI
lines(yrs,
      recruitNoLAI$value,
      col=rgb(1,0,0,1),lwd=3)
polygon(c(yrs,rev(yrs)),
        c(recruitNoLAI$lbnd,rev(recruitNoLAI$ubnd)),
        col=rgb(1,0,0,0.1),lty=0)

# retrospective comparison
yrs_retro <- SSB_mrNoLAI$year

par(mfrow=c(3,1))
# SSB
plot(yrs_retro,
     SSB_mrLAI$rho,
     xlab='year',ylab='mohn rho SSB',main='mohn rho SSB',type='l',col=rgb(0,0,1,1),lwd=3)
lines(yrs_retro,
      SSB_mrNoLAI$rho,
      col=rgb(1,0,0,1),lwd=3)
legend("topright",
       c("with LAI","without LAI"),
       col=c(rgb(0,0,1,1),rgb(1,0,0,1)),lty=1,lwd=3)

# fbar
plot(yrs_retro,
     fbar_mrLAI$rho,
     xlab='year',ylab='mohn rho fbar',main='mohn rho fbar',type='l',col=rgb(0,0,1,1),lwd=3)
lines(yrs_retro,
      fbar_mrNoLAI$rho,
      col=rgb(1,0,0,1),lwd=3)

# recruitment
plot(yrs_retro,
     rec_mrLAI$rho,
     xlab='year',ylab='mohn rho rec',main='mohn rho rec',type='l',col=rgb(0,0,1,1),lwd=3)
lines(yrs_retro,
      rec_mrNoLAI$rho,
      col=rgb(1,0,0,1),lwd=3)