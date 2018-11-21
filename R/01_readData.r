#-------------------------------------------------------------------------------
# WKNSMSE
#
# Author: Niels Hintzen
#         IMARES, The Netherland
#
#  MSE of North Sea Herring
#
# Date: 2018/11/18
#
# Build for R3.4.3, 64bits
#-------------------------------------------------------------------------------

rm(list=ls())

library(FLSAM)
library(FLEDA)

#path          <- "D:/Work/Herring MSE/NSAS/"
#path              <- "D:/git/wk_WKNSMSE_her.27.3a47d/R/"
path              <- "F:/WKNSMSE/wk_WKNSMSE_her.27.3a47d/R"
assessment_name   <- "NSAS_WKNSMSE2018"
try(setwd(path),silent=TRUE)

# paths to different subfolders
dataPath      <- file.path(".","data/")
outPath       <- file.path(".","results/")
scriptPath    <- file.path(".","side_scripts/")

#-------------------------------------------------------------------------------
# 0): multi-fleet assessment
#-------------------------------------------------------------------------------

source(file.path(scriptPath,"setupAssessmentObjects_mf.r"))
source(file.path(scriptPath,"setupControlObject_mf.r"))

NSH3f.sam   <- FLSAM(NSHs3,
                     NSH.tun,
                     NSH3.ctrl)

# save assessment object
save(NSHs3,
     NSH.tun,
     NSH3.ctrl,
     NSH3f.sam,
     file=file.path(outPath,paste0(assessment_name,'_mf.Rdata')))

#-------------------------------------------------------------------------------
# 0): Single fleet assessment
#-------------------------------------------------------------------------------

source(file.path(scriptPath,"setupAssessmentObjects_sf.r"))
source(file.path(scriptPath,"setupControlObject_sf.r"))

#- Perform the assessment
NSH.sam       <- FLSAM(NSH,
                       NSH.tun,
                       NSH.ctrl)
name(NSH.sam) <- assessment_name

NSH@stock.n           <- NSH.sam@stock.n[,ac(range(NSH)["minyear"]:range(NSH)["maxyear"])]
NSH@harvest           <- NSH.sam@harvest[,ac(range(NSH)["minyear"]:range(NSH)["maxyear"])]

# save assessment object
save(NSH,
     NSH.tun,
     NSH.ctrl,
     NSH.sam,
     file=file.path(outPath,paste0(assessment_name,'_sf.Rdata')))
