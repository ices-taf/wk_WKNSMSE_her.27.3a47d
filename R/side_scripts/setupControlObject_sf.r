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

NSH.ctrl                                    <- FLSAM.control(NSH,NSH.tun)

catchRow                                    <- grep("catch",rownames(NSH.ctrl@f.vars))
laiRow                                      <- grep("LAI",rownames(NSH.ctrl@obs.vars))

#Variance in N random walks (set 1st one free is usually best)
NSH.ctrl@logN.vars[]                        <- c(1,rep(2,dims(NSH)$age-1))

#All fishing mortality states are free except
#oldest ages to ensure stablity
NSH.ctrl@states[catchRow,]                  <- seq(dims(NSH)$age)
NSH.ctrl@states[catchRow,ac(7:8)]           <- 101

#Group observation variances of catches to ensure stability
NSH.ctrl@obs.vars[catchRow,]                <- c(0,0,1,1,1,1,2,2,2)
NSH.ctrl@obs.vars["HERAS",ac(1:8)]          <- c(101,rep(102,5),rep(103,2))
NSH.ctrl@obs.vars["IBTS-Q1",ac(1)]          <- 201
NSH.ctrl@obs.vars["IBTS0",ac(0)]            <- 301
NSH.ctrl@obs.vars["IBTS-Q3",ac(0:5)]        <- c(400,401,rep(402,4))
NSH.ctrl@obs.vars[laiRow,1]                 <- 501

#Catchabilities of the surveys. Set LAI all to 1 value, rest can be varied
NSH.ctrl@catchabilities["HERAS",ac(1:8)]    <- c(101,102,rep(103,6))
NSH.ctrl@catchabilities["IBTS-Q1",ac(1)]    <- c(201)
NSH.ctrl@catchabilities["IBTS-Q3",ac(0:5)]  <- c(300:305)
NSH.ctrl@catchabilities[laiRow,1]           <- 401

#Add correlation correction for Q3 survey, not for HERAS
  idx                                       <- which(rownames(NSH.ctrl@cor.obs)=="IBTS-Q3")
NSH.ctrl@cor.obs[idx,1:5]                   <- c(rep(101,5))
NSH.ctrl@cor.obs.Flag[idx]                  <- as.factor("AR")

#Variance of F random walks
NSH.ctrl@f.vars[1,]                         <- c(101,101,rep(102,4),rep(103,3))
#No correlation structure among ages in F random walks
NSH.ctrl@cor.F                              <- 0
#Finalise
NSH.ctrl@name                               <- "North Sea Herring"
NSH.ctrl                                    <- update(NSH.ctrl)

