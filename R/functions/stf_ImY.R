stf_ImY      <- function(stocks,
                         fishery,
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
  dms         <- dimnames(stocks@m)
  dms$year    <- ac((an(ImY)-1):an(CtY))
  dms$unit    <- c("A","B","C","D")
  
  nAges       <- length(dms$age)
  nits        <- length(dms$iter)
  
  # Create the stf object 
  stf         <- FLStock(name=nam,desc=dsc,FLQuant(NA,dimnames=dms))
  units(stf)  <- units(stocks)
  
  # initialize stf object
  stf[,,dms$unit]           <- window(stocks,
                                      start=an(dms$year)[1],
                                      end=rev(an(dms$year))[1])
  
  stf[,dms$year[1]]@harvest <- fishery@landings.sel[,dms$year[1]]
  stf[,FuY]@catch.wt <- fishery@landings.wt[,FuY]
  stf[,FuY]@m        <- stf[,ac(an(ImY)-1)]@m # copy M from terminal year to intermediate/forecast/continuation years
  stf[,FuY]@stock.wt        <- stf[,ac(an(ImY)-1)]@stock.wt # copy stock.wt from terminal year to intermediate/forecast/continuation years
  stf[,FuY]@harvest.spwn  <- stf[,ac(an(ImY)-1)]@harvest.spwn
  stf[,FuY]@m.spwn        <- stf[,ac(an(ImY)-1)]@m.spwn
  stf[,FuY]@mat           <- stf[,ac(an(ImY)-1)]@mat # copy maturity from terminal year to intermediate/forecast/continuation years
  
  # Fill slots that have no meaning for NSAS
  stf@discards.n[]          <- 0
  stf@discards[]            <- 0
  stf@discards.wt[]         <- 0
  
  # compute stock in ImY, use sf harvest in previous year
  Z <- drop(stocks@harvest[,ac(an(ImY)-1)]) + drop(stf[,ImY,1]@m)
  
  # propagate stock number with Z, only fill first slot
  survivors                             <- drop(stf@stock.n[,ac(an(ImY)-1),1])*exp(-Z) # stock.n is the same for all fleets in the stf object, taking first element
  stf@stock.n[2:nAges,ImY,1]            <- survivors[1:(nAges-1),]
  stf@stock.n[nAges,ImY,1]                <- stf[nAges,ImY,1]@stock.n + survivors[nAges,]
  stf@stock.n[1,ImY,1]         <- recruit
  
  #print(iter(stf@stock.n[,ImY,1],10))
  
  # copy stock for all fleets
  stf@stock.n[,ImY,2] <- stf@stock.n[,ImY,1]
  stf@stock.n[,ImY,3] <- stf@stock.n[,ImY,1]
  stf@stock.n[,ImY,4] <- stf@stock.n[,ImY,1]
  
  ############# Compute F in intermediate year #############
  
  #if("doParallel" %in% (.packages()))
  #  detach("package:doParallel",unload=TRUE)
  #if("foreach" %in% (.packages()))
  #  detach("package:foreach",unload=TRUE)
  #if("iterators" %in% (.packages()))
  #  detach("package:iterators",unload=TRUE)
  
  #require(doParallel)
  #ncores <- detectCores()-3
  #ncores <- ifelse(iters<ncores,nits,ncores)
  #cl <- makeCluster(ncores) #set up nodes
  #clusterEvalQ(cl,library(FLSAM))
  #clusterEvalQ(cl,library(stockassessment))
  #clusterEvalQ(cl,library(minpack.lm))
  #clusterEvalQ(cl,source('E:/wk_WKNSMSE_her.27.3a47d/R/functions/optF_ImY.R'))
  #clusterEvalQ(cl,library(stats))
  #clusterEvalQ(cl,library(FLEDA))
  #registerDoParallel(cl)
  
  
  #r<- foreach(idxIter=1:nits) %dopar% {nls.lm(  par=runif(4), # starting point
  #                                              lower=rep(1e-8,4),
  #                                              upper=NULL,
  #                                              optF_ImY,                                  # function to optimize
  #                                              fishery     = fishery[,,,,,idxIter],       # catch weight at age single fleet. Using stock weight at age for now. How to get catch weight for 2018?
  #                                              stock.n_sf  = iter(stf@stock.n[,,1],idxIter),       # stock number single fleet, taking first element, M is fleet independent
  #                                              M           = iter(stf@m[,,1],idxIter),              # natural mortality, taking first element, M is fleet independent
  #                                              iYr         = ImY,                              # year of interest
  #                                              TACs        = TAC[,,,,,idxIter],             # TAC FLQuant object for fleets A, B and D
  #                                              FCProp      = iter(FCPropIts,idxIter),
  #                                              TAC_var     = TAC_var,
  #                                              recruit     = recFuture[idxIter],
  #                                              nls.lm.control(ftol = (.Machine$double.eps),maxiter = 1000), # optimizer control object
  #                                              jac=NULL)$par}
  
  #if("doParallel" %in% (.packages()))
  #  detach("package:doParallel",unload=TRUE)
  #if("foreach" %in% (.packages()))
  #  detach("package:foreach",unload=TRUE)
  #if("iterators" %in% (.packages()))
    #detach("package:iterators",unload=TRUE)
  
  Fscalor <- array( 0, dim=c(nFleets,nits)) # initialize array
  
  #for(idxIter in 1:nits){Fscalor[,idxIter]<-r[[idxIter]]}
  
  for(idxIter in 1:nits){
    Fscalor[,idxIter]<-nls.lm(  par=runif(4), # starting point
                                lower=rep(1e-8,4),
                                upper=NULL,
                                optF_ImY,                                  # function to optimize
                                fishery     = fishery[,,,,,idxIter],       # catch weight at age single fleet. Using stock weight at age for now. How to get catch weight for 2018?
                                stock.n_sf  = iter(stf@stock.n[,,1],idxIter),       # stock number single fleet, taking first element, M is fleet independent
                                M           = iter(stf@m[,,1],idxIter),              # natural mortality, taking first element, M is fleet independent
                                iYr         = ImY,                              # year of interest
                                TACs        = TAC[,,,,,idxIter],             # TAC FLQuant object for fleets A, B and D
                                FCProp      = iter(FCPropIts,idxIter),
                                TAC_var     = TAC_var,
                                recruit     = recFuture[idxIter],
                                nls.lm.control(ftol = (.Machine$double.eps),maxiter = 1000), # optimizer control object
                                jac=NULL)$par
  }
  
  #print(Fscalor)
  
  stf@harvest[,ImY,'A'] <- t(apply(fishery@landings.sel[,ImY,'A'],1,'*',Fscalor[1,]))
  stf@harvest[,ImY,'B'] <- t(apply(fishery@landings.sel[,ImY,'B'],1,'*',Fscalor[2,]))
  stf@harvest[,ImY,'C'] <- t(apply(fishery@landings.sel[,ImY,'C'],1,'*',Fscalor[3,]))
  stf@harvest[,ImY,'D'] <- t(apply(fishery@landings.sel[,ImY,'D'],1,'*',Fscalor[4,]))
  
  harvestAll <- apply(stf@harvest[,ImY],6,'rowSums')
  
  #Z <- harvestAll + drop(stf[,ImY,1]@m) # copy M from previous year (i.e. terminal year of assessment)
  #Z <- harvestAll + drop(stf[,ImY,1]@m) # copy M from previous year (i.e. terminal year of assessment)
  
  Z <- harvestAll + drop(stf[,ImY,1]@m) # copy M from previous year (i.e. terminal year of assessment)
  
  Zmat <- replicate(4, Z)
  Zmat <- aperm(Zmat, c(1,3,2))
  
  #print(dim((1-exp(-Zmat))))
  
  stf@catch.n[,ImY] <- drop(stf@harvest[,ImY])*drop(stf@stock.n[,ImY])*drop((1-exp(-Zmat))/Zmat)
  stf@landings.n[,ImY] <- drop(stf@harvest[,ImY])*drop(stf@stock.n[,ImY])*drop((1-exp(-Zmat))/Zmat)
  
  stf@catch       <- computeCatch(stf)
  stf@landings    <- computeCatch(stf)
  stf@stock       <- computeStock(stf)
  
  return(stf)
  
}