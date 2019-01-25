# NSAS WKNSMSE 2018
# function to scale F fleet A and B

rescaleF_Ftargets      <- function( mult,         # scalor 2x1
                                    Fsel,         # selectivity stored as FLQuant object. Normalized between 0 and 1.
                                    iYr,
                                    Ftarget,
                                    F01,
                                    catch.wt_mf,   # catch weight at age single fleet
                                    stock.n_sf,  # stock number single fleet
                                    M,            # natural mortality
                                    TACs,         # TAC FLQuant object for fleets A, B and D
                                    FCProp,
                                    TAC_var,
                                    recruit){  # proportion of F for the C fleet
  
  
  
  # start fun
  nFleets  <- dim(TACs)[5] # 2 fleets to have F rescaled (A, B)
  nAges    <- dim(Fsel)[1]
  strFleet <- c('A','B','C','D')
  
  # compute new F using 1 scalor across ages for each fleet
  Ffleet <- array( 0, dim=c(nAges,nFleets)) # initialize array
  for(idxFleet in 1:nFleets){
    if(strFleet[idxFleet] == 'B' || strFleet[idxFleet] == 'D'){
      Ffleet[,idxFleet] <- Fsel[,iYr,'BD']*mult[idxFleet]
    }else{
      Ffleet[,idxFleet] <- Fsel[,iYr,strFleet[idxFleet]]*mult[idxFleet]
    }
  }
  
  Z <-  rowSums(Ffleet) + # use single fleet F at age
    drop(M[,iYr]) # M is fleet independent and all fleet fields are the same
  
  # propagate stock number with Z
  survivors                 <- stock.n_sf[,ac(an(iYr)-1)]*exp(-Z)
  stock.n_sf[2:nAges,iYr]   <- survivors[1:(nAges-1)]
  stock.n_sf[nAges,iYr]     <- stock.n_sf[nAges,iYr] + survivors[nAges] # add plus group
  stock.n_sf[1]             <- recruit
  
  # compute catch at age (in weight)
  catchfleet <- array( 0, dim=c(nAges,nFleets)) # initialize array for catches
  for(idxFleet in 1:nFleets){
    if(strFleet[idxFleet] == 'B' || strFleet[idxFleet] == 'D'){
      catchfleet[,idxFleet] <- Ffleet[,idxFleet]/Z*(1-exp(-Z))*drop(stock.n_sf[,iYr]*catch.wt_mf[,iYr,'BD']) # I should use catch weight here
    }else{
      catchfleet[,idxFleet] <- Ffleet[,idxFleet]/Z*(1-exp(-Z))*drop(stock.n_sf[,iYr]*catch.wt_mf[,iYr,strFleet[idxFleet]]) # I should use catch weight here
    }
  }
  
  catchfleet <- colSums(catchfleet)
  
  Ftot <- apply(Ffleet,1,'sum')
  
  F2up_bar <- mean(Ftot[3:length(Ftot)])
  F01_bar  <- mean(Ftot[1:2])
  
  # compute TAC at age for the C fleet from the proportion of F
  TAC_C_IIIa <- rowSums(Ffleet)*FCProp[iYr]/Z*(1-exp(-Z))*drop(stock.n_sf[,iYr]*catch.wt_mf[,iYr,'C'])
  TAC_C_IIIa <- sum(TAC_C_IIIa)
  
  TACs <-drop(TACs[,iYr]) # reduce object to the year of interest
  TAC_var <- TAC_var[iYr,]
  TACs['D'] <- TACs['D']*TAC_var['Dsplit']*TAC_var['Duptake']
  
  # optimize based on Ftarget F0-1, TAC C fleet in NS, TAC D fleet in NS
  res <- sqrt((c(Ftarget,F01,TAC_C_IIIa,TACs['D']) - c(F2up_bar, F01_bar,catchfleet[3],catchfleet[4]))^2)
  
  return(res)
}