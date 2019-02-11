# NSAS WKNSMSE 2018
# function to scale F fleet A and B

optF_FcY      <- function(  mult,         # scalor 2x1
                            fishery,         # selectivity stored as FLQuant object. Normalized between 0 and 1.
                            iYr,
                            Ftarget,
                            F01,
                            stock.n_sf,  # stock number single fleet
                            M,            # natural mortality
                            TACs,         # TAC FLQuant object for fleets A, B and D
                            FCProp,
                            TAC_var,
                            recruit){  # proportion of F for the C fleet
  
  
  
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
  
  F2up_bar <- mean(Ftot[3:length(Ftot)])
  F01_bar  <- mean(Ftot[1:2])
  
  # F target for the C fleet
  Ftarget_C <- sum(rowSums(Ffleet)*FCProp[iYr])
  
  # compute TAC at age for the C fleet from the proportion of F
  #TAC_C_IIIa <- rowSums(Ffleet)*FCProp[iYr]/Z*(1-exp(-Z))*drop(stock.n_sf[,iYr]*fishery@landings.wt[,iYr,'C'])
  #TAC_C_IIIa <- sum(TAC_C_IIIa)
  
  TACs <-drop(TACs[,iYr]) # reduce object to the year of interest
  TAC_var <- TAC_var[iYr,]
  TACs['D'] <- TACs['D']*TAC_var['Dsplit']*TAC_var['Duptake']
  
  
  
  A <- (1 - F2up_bar/Ftarget)^2
  B <- (1 - F01_bar/F01)^2
  C <-  (1 - sum(Ffleet[,3])/Ftarget_C)^2
  D <- (1 - catchfleet[4]/TACs['D'])^2
  
  #res <- sqrt(c(A,B,C,D))
  
  # optimize based on Ftarget F0-1, TAC C fleet in NS, TAC D fleet in NS
  res <- sqrt((c(Ftarget,F01,Ftarget_C,TACs['D']) - c(F2up_bar, F01_bar,sum(Ffleet[,3]),catchfleet[4]))^2)
  
  return(res)
}