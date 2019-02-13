# NSAS WKNSMSE 2018
# function to scale F fleet A and B

optF_FcY      <- function(  mult,         # scalor 2x1
                            fishery,         # selectivity stored as FLQuant object. Normalized between 0 and 1.
                            iYr,
                            stock.n_sf,  # stock number single fleet
                            stock.wt_sf,
                            M,            # natural mortality
                            M.spwn,
                            F.spwn,
                            mats,
                            TACs,         # TAC FLQuant object for fleets A, B and D
                            FCProp,
                            TAC_var,
                            recruit,
                            refpoints,
                            HCR){  # proportion of F for the C fleet
  
  
  # start fun
  nFleets  <- dim(TACs)[5] # 2 fleets to have F rescaled (A, B)
  nAges    <- dim(fishery@landings.sel)[1]
  strFleet <- c('A','B','C','D')
  
  # compute new F using 1 scalor across ages for each fleet
  Ffleet <- array( 0, dim=c(nAges,nFleets)) # initialize array
  for(idxFleet in 1:nFleets){
    Ffleet[,idxFleet] <- drop(fishery@landings.sel[,iYr,strFleet[idxFleet]])*mult[idxFleet]
  }
  
  Z <-  rowSums(Ffleet) + # use single fleet F at age
        drop(M[,iYr]) # M is fleet independent and all fleet fields are the same
  
  # compute catch at age (in weight)
  catchfleet <- array( 0, dim=c(nAges,nFleets)) # initialize array for catches
  for(idxFleet in 1:nFleets){
    catchfleet[,idxFleet] <- Ffleet[,idxFleet]/Z*(1-exp(-Z))*drop(stock.n_sf[,iYr]*fishery@landings.wt[,iYr,strFleet[idxFleet]])
  }
  
  catchfleet <- colSums(catchfleet)
  
  Ftot <- apply(Ffleet,1,'sum')
  
  # calculate SSB at FcY
  ssb_FcY               <- sum( drop(stock.n_sf[,iYr]*stock.wt_sf[,iYr])*
                                exp(-Ftot*drop(F.spwn[,iYr]-M[,iYr]*M.spwn[,iYr]))*drop(mats[,iYr]))

  # calculate Ftarget and F01
  if(HCR == 'A'){
    if(ssb_FcY <= refpoints$Btrigger){
      F26tar   <- refpoints$Ftarget*ssb_FcY/refpoints$Btrigger
      F01tar   <- refpoints$F01*ssb_FcY/refpoints$Btrigger
    }else{
      F26tar   <- refpoints$Ftarget
      F01tar   <- refpoints$F01
    }
  }
  
  if(HCR == 'B'){
    if(ssb_FcY <= refpoints$Btrigger && ssb_FcY > refpoints$Blim){
      F26tar   <- refpoints$Ftarget*ssb_FcY/refpoints$Btrigger
      F01tar   <- refpoints$F01*ssb_FcY/refpoints$Btrigger
    }else if(ssb_FcY < refpoints$Blim){
      F26tar <- 0.1
      F01tar    <- 0.04
    }else{
      F26tar  <- refpoints$Ftarget
      F01tar   <- refpoints$F01
    }
  }

  # calculate F26 and F01 for current F  
  F26_bar <- mean(Ftot[3:7])# using F2-6 instead of F2+ - mean(Ftot[3:length(Ftot)])
  F01_bar  <- mean(Ftot[1:2])

  # F target for the C fleet
  Ftarget_C <- sum(rowSums(Ffleet)*drop(FCProp[,iYr]))

  TACs <-drop(TACs[,iYr]) # reduce object to the year of interest
  TAC_var <- TAC_var[iYr,]
  TACs['D'] <- TACs['D']*TAC_var['Dsplit']*TAC_var['Duptake']
  
  resA <- (F26_bar-F26tar)/F26tar
  resB <- (F01_bar-F01tar)/F01tar
  resC <- (sum(Ffleet[,3])-Ftarget_C)/Ftarget_C
  resD <- (catchfleet[4]-TACs['D'])/TACs['D']
  
  # optimize based on Ftarget F0-1, TAC C fleet in NS, TAC D fleet in NS
  #res <- sqrt((c(Ftarget,F01,Ftarget_C,TACs['D']) - c(F2up_bar, F01_bar,sum(Ffleet[,3]),catchfleet[4]))^2)
  res <- sqrt(c(resA,resB,resC,resD)^2)

  return(res)
}