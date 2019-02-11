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

optF_ImY      <- function(  mult,         # scalor 4x1
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
  
  #Ftot <- apply(Ffleet,1,'sum')
  #F2up_bar <- mean(Ftot[3:length(Ftot)])
  #print(F2up_bar)
  
  # compute Z using scaled Fs. Z = M+Ftot with Ftot = FA+FB+FC+FD.
  #Z <-  rowSums(Ffleet) + # use single fleet F at age
  #      drop(M[,iYr]) # M is fleet independent and all fleet fields are the same
  #print(drop(M[,iYr]))
  
  # not ideal but use F from previous year for stability, otherwise there is too much to optimize
  #Z <-  rowSums(drop(fishery[,ac(an(iYr)-1)]@landings.sel)) + # use single fleet F at age
  #      drop(M[,iYr]) # M is fleet independent and all fleet fields are the same
  
  Z <-  rowSums(Ffleet) + # use single fleet F at age
        drop(M[,iYr]) # M is fleet independent and all fleet fields are the same
  
  # compute catch at age (in weight)
  catchfleet <- array( 0, dim=c(nAges,nFleets)) # initialize array for catches
  for(idxFleet in 1:nFleets){
    #print(drop(stock.n_sf[,iYr]*fishery@landings.wt[,iYr,strFleet[idxFleet]]))
    catchfleet[,idxFleet] <- Ffleet[,idxFleet]/Z*(1-exp(-Z))*drop(stock.n_sf[,iYr]*fishery@landings.wt[,iYr,strFleet[idxFleet]])
  }
  
  #print(catchfleet)
  
  # sum accross the ages
  catchfleet <- colSums(catchfleet)
  
  # compute TAC at age for the C fleet from the proportion of F
  TAC_C_IIIa <- rowSums(Ffleet)*FCProp[iYr]/Z*(1-exp(-Z))*drop(stock.n_sf[,iYr]*fishery@landings.wt[,iYr,'C'])
  TAC_C_IIIa <- sum(TAC_C_IIIa)
  
  # F target for the C fleet in NS
  Ftarget_C <- sum(rowSums(Ffleet)*FCProp[iYr])
  
  
  # Fill TAC objects with data on expected catches in NS
  TACs <-drop(TACs[,iYr]) # reduce object to the year of interest
  #print(TACs)
  TAC_var <- TAC_var[iYr,]
  TACs['A'] <- TACs['A'] + TAC_var['Ctransfer']*TACs['C']
  TACs['C'] <- TAC_C_IIIa
  TACs['B'] <- TACs['B']*TAC_var['Buptake']
  TACs['D'] <- TACs['D']*TAC_var['Dsplit']*TAC_var['Duptake']

  res <- sqrt((c(TACs['A'],TACs['B'],Ftarget_C,TACs['D']) - c(catchfleet[1], catchfleet[2],sum(Ffleet[,3]),catchfleet[4]))^2)
  
  return(res)
}