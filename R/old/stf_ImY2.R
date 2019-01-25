stf_ImY2      <- function(NSH.sim,
                         stocks,
                         fisheryFuture,
                         TAC,
                         TAC_var,
                         FCPropIts,
                         recruit,
                         FuY){
  
  ImY <- FuY[1]
  FcY <- FuY[2]
  CtY <- FuY[3]
  
  ############# initialize stf FLStock object #############
  dsc         <- "North Sea Herring"
  nam         <- "NSAS"
  dms         <- dimnames(stocks[[1]]@m)
  dms$year    <- ac((an(ImY)-3):an(CtY))
  dms$unit    <- c("A","B","C","D")
  dms$iter    <- 1:nits
  
  nAges       <- length(dms$age)
  nits        <- length(dms$iter)
  
  # Create the stf object 
  stf         <- FLStock(name=nam,desc=dsc,FLQuant(NA,dimnames=dms))
  
  for(idxIter in 1:nits){
    for(idxFleet in dms$unit){
      # update stf object with current stock object
      stf[,,idxFleet,,,idxIter]          <- window( stocks[[idxIter]],
                                                    start=an(dms$year)[1],
                                                    end=rev(an(dms$year))[1])
      # update stf object with stock numbers from SAM object
      stf[,,idxFleet,,,idxIter]@stock.n  <- window( NSH.sim[[idxIter]]@stock.n,
                                                    start=an(dms$year)[1],
                                                    end=rev(an(dms$year))[1])
      
      # update catch.wt
      if(idxFleet == 'B' || idxFleet == 'D'){
        stf[,,idxFleet,,,idxIter]@catch.wt <- fisheryFuture[,FuY,'BD','catch.wt',,idxIter]
      }else{
        stf[,,idxFleet,,,idxIter]@catch.wt <- fisheryFuture[,FuY,idxFleet,'catch.wt',,idxIter]
      }
    }
  }
  
  # Fill slots that have no meaning for NSAS
  stf@discards.n[]          <- 0
  stf@discards[]            <- 0
  stf@discards.wt[]         <- 0
  
  
  ############################ END TEST
  
  ############# Compute F in intermediate year #############
  
  # update intermediate year with results from assessment
  
  Fscalor <- array( 0, dim=c(nFleets,nits)) # initialize array
  for(idxIter in 1:nits){
    # find F for the different fleets
    Fscalor[,idxIter] <- nls.lm(  par=runif(4), # starting point
                                  lower=rep(1e-8,4),
                                  upper=NULL,
                                  optF_TACdiff,                                  # function to optimize
                                  harvest_sf  = NSH.sim[[idxIter]]@harvest,       # single fleet FLStock object
                                  catch.wt_mf = fisheryFuture[,,,'catch.wt',,idxIter],       # catch weight at age single fleet. Using stock weight at age for now. How to get catch weight for 2018?
                                  stock.n_sf  = NSH.sim[[idxIter]]@stock.n,       # stock at age
                                  M           = stocks[[idxIter]]@m,              # natural mortality
                                  Fsel        = fisheryFuture[,,,'sel',,idxIter], # selectivity stored as FLQuant object. Normalized between 0 and 1.
                                  iYr         = ImY,                              # year of interest
                                  TACs        = TAC[,,,,,idxIter],             # TAC FLQuant object for fleets A, B and D
                                  FCProp      = FCPropIts[,idxIter],
                                  TAC_var     = TAC_var,
                                  recruit     = recruit[idxIter],
                                  nls.lm.control(ftol = (.Machine$double.eps),maxiter = 1000), # optimizer control object
                                  jac=NULL)$par
  }
  
  ############# update stf object #############
  
  # update F
  #stf@harvest[,ImY,'A'] <- t(apply(fisheryFuture[,ImY,'A','sel'],1,'*',Fscalor[1,]))
  #stf@harvest[,ImY,'B'] <- t(apply(fisheryFuture[,ImY,'BD','sel'],1,'*',Fscalor[2,]))
  #stf@harvest[,ImY,'C'] <- t(apply(fisheryFuture[,ImY,'C','sel'],1,'*',Fscalor[3,]))
  #stf@harvest[,ImY,'D'] <- t(apply(fisheryFuture[,ImY,'BD','sel'],1,'*',Fscalor[4,]))
  
  
  # compute stock.n, catch.n and landing.n
  for(idxIter in 1:nits){
    stf@harvest[,ImY,'A',,,idxIter] <- fisheryFuture[,ImY,'A','sel',,idxIter]*Fscalor[1,idxIter]
    stf@harvest[,ImY,'B',,,idxIter] <- fisheryFuture[,ImY,'BD','sel',,idxIter]*Fscalor[2,idxIter]
    stf@harvest[,ImY,'C',,,idxIter] <- fisheryFuture[,ImY,'C','sel',,idxIter]*Fscalor[3,idxIter]
    stf@harvest[,ImY,'D',,,idxIter] <- fisheryFuture[,ImY,'BD','sel',,idxIter]*Fscalor[4,idxIter]
    
    Z <-  rowSums(drop(stf@harvest[,ImY,,,,idxIter])) + # sum accross the fleets
      drop(stf[,ImY,1,,,idxIter]@m) # M is the same for all fleets in the stf object
    
    # propagate stock number with Z
    survivors                             <- stf[,ac(an(ImY)-1),1,,,idxIter]@stock.n*exp(-Z) # stock.n is the same for all fleets in the stf object
    stf@stock.n[2:nAges,ImY,,,,idxIter]   <- survivors[1:(nAges-1)]
    stf@stock.n[nAges,ImY,,,,idxIter]     <- stf[nAges,ImY,1,,,idxIter]@stock.n + survivors[nAges]
    stf@stock.n[1,ImY,,,,idxIter]         <- recruit[idxIter]
    
    for(idxFleet in 1:nFleets){  
      stf@catch.n[,ImY,idxFleet,,,idxIter]     <-  stf@stock.n[,ImY,idxFleet,,,idxIter]*
        (1-exp(-Z))*
        stf@harvest[,ImY,idxFleet,,,idxIter]/Z
      stf@catch[,ImY,idxFleet,,,idxIter]       <- computeCatch(stf[,ImY,idxFleet,,,idxIter])
      stf@landings.n[,ImY,idxFleet,,,idxIter]  <- stf@catch.n[,ImY,idxFleet,,,idxIter]
      stf@landings[,ImY,idxFleet,,,idxIter]    <- computeLandings(stf[,ImY,idxFleet,,,idxIter])
    }
  }
  
  return(stf)
  
}