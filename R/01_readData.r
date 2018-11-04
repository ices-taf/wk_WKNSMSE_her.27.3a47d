#-------------------------------------------------------------------------------
# WKHELP
#
# Author: Niels Hintzen
#         IMARES, The Netherland
#
# Performs an MSE of North Sea Herring under different TAC scenario's
#
# Date: 05-Jun-2012
#
# Build for R2.13.2, 32bits
#-------------------------------------------------------------------------------

rm(list=ls())

library(FLCore)
library(FLSAM)

#path          <- "D:/Work/Herring MSE/NSAS/"
path          <- "W:/IMARES/Data/ICES-WG/WKHerTAC/NSAS/"
inPath        <- paste(path,"Data/",sep="")
codePath      <- paste(path,"R/",sep="")
outPath       <- paste(path,"Results/",sep="")
if(substr(R.Version()$os,1,3)== "lin"){
  path        <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",path)
  inPath      <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",inPath)
  codePath    <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",codePath)
  outPath     <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",outPath)
}

  #-------------------------------------------------------------------------------
  # 0): Read 2010 data & setup assessment
  #-------------------------------------------------------------------------------

#- Read NSH stock data + stock object fixes
setwd(path)
source(file.path(codePath,"setupAssessmentObjects.r"))
source(file.path(codePath,"setupControlObject.r"))

#- Perform the assessment
NSH.sam       <- FLSAM(NSH,NSH.tun,NSH.ctrl)
name(NSH.sam) <- "North Sea Herring"

#Update stock object
NSH           <- NSH + NSH.sam
NSH@stock     <- computeStock(NSH)

  #-------------------------------------------------------------------------------
  # 1): Save 2011 data
  #-------------------------------------------------------------------------------

save(NSH,       file=paste(outPath,"NSH.RData",     sep=""))
save(NSH.ctrl,  file=paste(outPath,"NSHctrl.RData", sep=""))
save(NSH.sam,   file=paste(outPath,"NSHsam.RData",  sep=""))
save(NSH.tun,   file=paste(outPath,"NSHtun.RData",  sep=""))

