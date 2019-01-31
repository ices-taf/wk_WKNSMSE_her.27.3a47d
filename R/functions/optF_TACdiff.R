# NSAS WKNSMSE 2018
# function that Compute differences between TAC and catches computed with F scalors
#
#  For the given year iYr and given a scalor for F for each fleet (i.e. 4 scalors), the function 
# computes the differences between the inferred catches and the TACs allocated to each fleet.
#
# Note 1: this function is specific to WKNSMSE for NSAS.
# Note 2: TAC is used for fleets A, B and D. For fleet C, the TAC is computed separately as 
# for this fleet, a proportion of F is used for forecast
# Note 3: fleet B and D are combined as they examplify the same selectivity patterns

optF_TACdiff      <- function( mult,         # scalor 4x1
                               fishery,   # catch weight at age single fleet
                               stock.n_sf,  # stock number single fleet
                               M,            # natural mortality
                               iYr,          # year of interest
                               TACs,         # TAC FLQuant object for fleets A, B and D
                               FCProp,
                               TAC_var,
                               recruit){  # proportion of F for the C fleet
  
  # start fun
  nFleets  <- dim(TACs)[5] # 4 fleets, A, B. C and D
  nAges    <- dim(stock.n_sf)[1]
  strFleet <- c('A','B','C', 'D')
  
  # compute new F using 1 scalor across ages for each fleet
  Ffleet <- array( 0, dim=c(nAges,nFleets)) # initialize array
  for(idxFleet in 1:nFleets){
    Ffleet[,idxFleet] <- drop(fishery@landings.sel[,iYr,strFleet[idxFleet]])*mult[idxFleet]
  }
  
  # compute Z using scaled Fs. Z = M+Ftot with Ftot = FA+FB+FC+FD.
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
    catchfleet[,idxFleet] <- Ffleet[,idxFleet]/Z*(1-exp(-Z))*drop(stock.n_sf[,iYr]*fishery@landings.wt[,iYr,strFleet[idxFleet]])
  }
  
  
  # sum accross the ages
  catchfleet <- colSums(catchfleet)
  
  # compute TAC at age for the C fleet from the proportion of F
  TAC_C_IIIa <- rowSums(Ffleet)*FCProp[iYr]/Z*(1-exp(-Z))*drop(stock.n_sf[,iYr]*fishery@landings.wt[,iYr,'C'])
  TAC_C_IIIa <- sum(TAC_C_IIIa)
  
  
  # Fill TAC objects with data on expected catches in NS
  TACs <-drop(TACs[,iYr]) # reduce object to the year of interest
  TAC_var <- TAC_var[iYr,]
  TACs['A'] <- TACs['A'] + TAC_var['Ctransfer']*TACs['C']
  TACs['C'] <- TAC_C_IIIa
  TACs['B'] <- TACs['B']*TAC_var['Buptake']
  TACs['D'] <- TACs['D']*TAC_var['Dsplit']*TAC_var['Duptake']
  

  res <- sqrt((drop(TACs) - catchfleet)^2)
  return(res)
}