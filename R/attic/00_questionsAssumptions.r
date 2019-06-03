#-------------------------------------------------------------------------------
# WKHERMPII
#
# Author: Niels Hintzen
#         IMARES, The Netherland
#
# Performs an MSE of North Sea Herring under different TAC scenario's
#
# Date: 03-Oct-2011
#
# Build for R2.8.1
#-------------------------------------------------------------------------------

  #-------------------------------------------------------------------------------
  # 0): Questions to ask
  #-------------------------------------------------------------------------------

  -Natural mortality in the future period (biol@m and stocks@m)
  -recruitment in future period            (biol@n[1,projPeriod])
  -how many years in average weight and n calculation (propN, propWt)
  -which year range to take variability over for variation in F-at-age (now 2001:2011)
  -do we also want error on historic catches (as we do add error to historic survey data)
  
  
  #-------------------------------------------------------------------------------
  # 1): Coding uncertainties
  #-------------------------------------------------------------------------------

  -add noise to catches in history as well
  -not for biological processes because we don't have any understanding of the uncertainty estimates there
  -do we want to do this backwards too for the biological processes???