################################################################################
# NSH_SAM Control for Assessment
#
# $Rev: 705 $
# $Date: 2012-02-14 19:02:57 +0100 (di, 14 feb 2012) $
#
# Author: HAWG model devlopment group
#
# Sets up a control object for use by Step 04 assessments i.e. the "refined data" run
#
# Developed with:
#   - R version 2.13.0
#   - FLCore 2.4
#
# To be done:
#
# Notes: Have fun running this assessment!
#
################################################################################


#- 3-fleets

NSH3.ctrl                           <- FLSAM.control(NSHs3,NSH.tun,sumFleets=dimnames(NSHs3[["residual"]]@catch)$area)

  catchRow                          <- grep("catch",rownames(NSH3.ctrl@f.vars))
  laiRow                            <- grep("LAI",rownames(NSH3.ctrl@obs.vars))
NSH3.ctrl@states["catch A",]        <- c(-1,0:6,6)
NSH3.ctrl@states["catch BD",]       <- c(7:9,10,10,10,rep(-1,3))
NSH3.ctrl@states["catch C",]        <- c(-1,11:13,14,14,14,rep(-1,2))

NSH3.ctrl@catchabilities["HERAS",ac(1:8)]   <- c(101,102,rep(103,6))
NSH3.ctrl@catchabilities["IBTS-Q3",ac(0:5)] <- c(200:205)
NSH3.ctrl@catchabilities[laiRow,1]  <- 301

NSH3.ctrl@obs.vars["catch A",]      <- c(-1,0,1,1,1,1,2,2,2)
NSH3.ctrl@obs.vars["catch BD",]     <- c(101,102,rep(103,4),rep(-1,3))
NSH3.ctrl@obs.vars["catch C",]      <- c(-1,201,202,rep(203,4),rep(-1,2))
NSH3.ctrl@obs.vars["HERAS",ac(1:8)] <- c(301,302,rep(302,4),rep(303,2))
NSH3.ctrl@obs.vars["IBTS-Q1",ac(1)] <- 401
NSH3.ctrl@obs.vars["IBTS0",ac(0)]   <- 501
NSH3.ctrl@obs.vars["IBTS-Q3",ac(0:5)]<- c(600,601,602,602,602,602)
NSH3.ctrl@obs.vars[laiRow,1]        <- 701#:604


  idx                               <- which(rownames(NSH3.ctrl@cor.obs)=="IBTS-Q3")
NSH3.ctrl@cor.obs[idx,1:5]          <- c(rep(102,5))
NSH3.ctrl@cor.obs.Flag[idx]         <- as.factor("AR")
NSH3.ctrl@cor.F                     <- as.integer(c(2,2,2))
NSH3.ctrl@f.vars["catch A",]        <- c(-1,rep(102,5),rep(103,3))
NSH3.ctrl@f.vars["catch BD",]       <- c(202,203,rep(203,4),rep(-1,3))
NSH3.ctrl@f.vars["catch C",]        <- c(-1,301,302,rep(303,4),rep(-1,2))
NSH3.ctrl@name                      <- "North Sea herring multifleet"
NSH3.ctrl                           <- update(NSH3.ctrl)

