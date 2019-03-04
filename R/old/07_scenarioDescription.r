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

#- Scenario descriptions
Fmsy          <- 0.26
FIAVBB        <- list(#- New runs for PRAC
                      opt1=list(Blim=0.8e6,Bpa=1.0e6,Btrigger=1.0e6,
                                FadultA=0.25,FadultB=0.1,FjuvA=0.05,FjuvB=0.04,
                                stabilityBreak="Blim",scen="FIAVBB"),
                      opt2=list(Blim=0.8e6,Bpa=1.0e6,Btrigger=1.1e6,
                                FadultA=0.25,FadultB=0.1,FjuvA=0.05,FjuvB=0.04,
                                stabilityBreak="Blim",scen="FIAVBB"),
                      opt3=list(Blim=0.8e6,Bpa=1.0e6,Btrigger=1.2e6,
                                FadultA=0.25,FadultB=0.1,FjuvA=0.05,FjuvB=0.04,
                                stabilityBreak="Blim",scen="FIAVBB"),
                      opt4=list(Blim=0.8e6,Bpa=1.0e6,Btrigger=1.3e6,
                                FadultA=0.25,FadultB=0.1,FjuvA=0.05,FjuvB=0.04,
                                stabilityBreak="Blim",scen="FIAVBB"),
                      opt5=list(Blim=0.8e6,Bpa=1.0e6,Btrigger=1.4e6,
                                FadultA=0.25,FadultB=0.1,FjuvA=0.05,FjuvB=0.04,
                                stabilityBreak="Blim",scen="FIAVBB"),
                      opt6=list(Blim=0.8e6,Bpa=1.0e6,Btrigger=1.4e6,
                                FadultA=0.25,FadultB=0.1,FjuvA=0.05,FjuvB=0.04,
                                stabilityBreak="Blim",scen="FIAVBB"))

LTMP          <- list(opt1=list(Blim=0.8e6,Bpa=1.3e6,Btrigger=1.5e6,
                                FadultA=0.25,FadultB=0.1,FjuvA=0.05,FjuvB=0.04,
                                stabilityBreak="Blim",scen="LTMP"))
