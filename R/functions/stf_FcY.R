# NSAS MSE. Compute stf for the forecast year.

stf_FcY <- function(stf,
                    fishery,
                    TAC,
                    TAC_var,
                    FCPropIts,
                    recruit,
                    HCR,
                    referencePoints,
                    FuY){
  
  ImY <- FuY[1]
  FcY <- FuY[2]
  CtY <- FuY[3]
  
  # calculate SSB at ImY
  Ftot_fleets <- apply(stf@harvest[,ImY],c(1,2,4,5,6),'sum')
  
  ssb_ImY <- apply(drop(stf@stock.n[,ImY,1]*stf@stock.wt[,ImY,1])*
                     exp(-drop(Ftot_fleets)*drop(stf@harvest.spwn[,ImY,1]) - drop(stf@m[,ImY,1]*stf@m.spwn[,ImY,1]))*
                     drop(stf@mat[,ImY,1]),
                   2,
                   'sum')
  
  if(HCR == 'A'){
    F2plusIter <-  array(NA,dim=c(1,nits))
    F01Iter     <-  array(NA,dim=c(1,nits))
    for(idxIter in 1:nits){
      if(ssb_ImY[idxIter] <= referencePoints$Btrigger){
        F2plusIter[idxIter]  <- referencePoints$Ftarget*ssb_ImY[idxIter]/referencePoints$Btrigger
        F01Iter[idxIter]   <- referencePoints$F01*ssb_ImY[idxIter]/referencePoints$Btrigger
      }else{
        F2plusIter[idxIter]  <- referencePoints$Ftarget
        F01Iter[idxIter]   <- referencePoints$F01
      }
    }
  }
  
  if(HCR == 'B'){
    F2plusIter <-  array(NA,dim=c(1,nits))
    F01Iter     <-  array(NA,dim=c(1,nits))
    for(idxIter in 1:nits){
      if(ssb_ImY[idxIter] <= referencePoints$Btrigger && ssb_ImY[idxIter] > referencePoints$Blim){
        F2plusIter[idxIter]  <- referencePoints$Ftarget*ssb_ImY[idxIter]/referencePoints$Btrigger
        F01Iter[idxIter]   <- referencePoints$F01*ssb_ImY[idxIter]/referencePoints$Btrigger
      }else if(ssb_ImY[idxIter] < referencePoints$Blim){
        F2plusIter[idxIter] <- 0.1
        F01Iter[idxIter]    <- 0.04
      }else{
        F2plusIter[idxIter]  <- referencePoints$Ftarget
        F01Iter[idxIter]   <- referencePoints$F01
      }
    }
  }

  # use harvest from ImY
  Z <- drop(apply(stf@harvest[,ImY],c(1,2,4,5,6),'sum')) + drop(stf[,FcY,1]@m)
  
  # propagate stock number with Z, only fill first slot
  survivors                             <- drop(stf@stock.n[,ImY,1])*exp(-Z) # stock.n is the same for all fleets in the stf object, taking first element
  stf@stock.n[2:nAges,FcY,1]            <- survivors[1:(nAges-1),]
  stf@stock.n[nAges,FcY,1]              <- stf[nAges,FcY,1]@stock.n + survivors[nAges,]
  stf@stock.n[1,FcY,1]         <- recruit
  
  # copy stock for all fleets
  stf@stock.n[,FcY,2] <- stf@stock.n[,FcY,1]
  stf@stock.n[,FcY,3] <- stf@stock.n[,FcY,1]
  stf@stock.n[,FcY,4] <- stf@stock.n[,FcY,1]
  
  if("doParallel" %in% (.packages()))
    detach("package:doParallel",unload=TRUE)
  if("foreach" %in% (.packages()))
    detach("package:foreach",unload=TRUE)
  if("iterators" %in% (.packages()))
    detach("package:iterators",unload=TRUE)
  
  require(doParallel)
  ncores <- detectCores()-1
  #ncores <- ifelse(iters<ncores,nits,ncores)
  cl <- makeCluster(ncores) #set up nodes
  clusterEvalQ(cl,library(FLSAM))
  clusterEvalQ(cl,library(stockassessment))
  clusterEvalQ(cl,library(minpack.lm))
  clusterEvalQ(cl,source('E:/wk_WKNSMSE_her.27.3a47d/R/functions/optF_FcY.R'))
  clusterEvalQ(cl,library(stats))
  clusterEvalQ(cl,library(FLEDA))
  registerDoParallel(cl)
  
  r<- foreach(idxIter=1:nits) %dopar% {nls.lm(  par=runif(4), # starting point
                                                lower=rep(1e-8,4),
                                                upper=c(1,1,1,1),
                                                optF_FcY,    
                                                fishery     = fishery[,,,,,idxIter],
                                                iYr         = FcY,
                                                Ftarget     = F2plusIter[idxIter],
                                                F01         = F01Iter[idxIter],
                                                stock.n_sf  = iter(stf@stock.n[,,1],idxIter),       
                                                M           = iter(stf@m[,,1],idxIter),             
                                                TACs        = TAC[,,,,,idxIter],             
                                                FCProp      = iter(FCPropIts,idxIter),
                                                TAC_var     = TAC_var,
                                                recruit     = recFuture[idxIter],
                                                nls.lm.control(ftol = (.Machine$double.eps),maxiter = 1000), # optimizer control object
                                                jac=NULL)$par
  }
  
  if("doParallel" %in% (.packages()))
    detach("package:doParallel",unload=TRUE)
  if("foreach" %in% (.packages()))
    detach("package:foreach",unload=TRUE)
  if("iterators" %in% (.packages()))
    detach("package:iterators",unload=TRUE)
  
  Fscalor <- array( 0, dim=c(nFleets,nits)) # initialize array
  
  for(idxIter in 1:nits){Fscalor[,idxIter]<-r[[idxIter]]}
  
  stf@harvest[,FcY,'A'] <- t(apply(fishery@landings.sel[,FcY,'A'],1,'*',Fscalor[1,]))
  stf@harvest[,FcY,'B'] <- t(apply(fishery@landings.sel[,FcY,'B'],1,'*',Fscalor[2,]))
  stf@harvest[,FcY,'C'] <- t(apply(fishery@landings.sel[,FcY,'C'],1,'*',Fscalor[3,]))
  stf@harvest[,FcY,'D'] <- t(apply(fishery@landings.sel[,FcY,'D'],1,'*',Fscalor[4,]))
  
  harvestAll <- apply(stf@harvest[,FcY],6,'rowSums')
  
  Z <- harvestAll + drop(stf[,FcY,1]@m)
  
  Zmat <- replicate(4, Z)
  Zmat <- aperm(Zmat, c(1,3,2))
  
  stf@catch.n[,FcY] <- drop(stf@harvest[,FcY])*drop(stf@stock.n[,FcY])*(1-exp(-Zmat))/Zmat
  stf@landings.n[,FcY] <- drop(stf@harvest[,FcY])*drop(stf@stock.n[,FcY])*(1-exp(-Zmat))/Zmat
  
  stf@catch       <- computeCatch(stf)
  stf@landings    <- computeCatch(stf)
  stf@stock       <- computeStock(stf)
  
  return(list(stf,F2plusIter,F01Iter))
  
}