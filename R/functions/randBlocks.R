randBlocks <- function(fecYears,projPeriod,nits){
  
  nyrs <- length(projPeriod)
  
  saveBlcks_nYears     <- matrix(NA,nrow=nits,ncol=nyrs)
  saveBlcks_startYear  <- matrix(NA,nrow=nits,ncol=nyrs)
  saveBlcks_rev        <- matrix(NA,nrow=nits,ncol=nyrs)
  
  for(idxIter in 1:nits){
    saveBlcks_nYears[idxIter,] <- randNums(nyrs,1,length(fecYears)-3,nyrs)
  }
  
  # randomize the starting years
  for(idxIter in 1:nits){
    nYearsTemp <- saveBlcks_nYears[idxIter,
                                   !is.na(saveBlcks_nYears[idxIter,])] # select only non NaN values
    # loop on non NaN values
    for(idxYear in 1:length(nYearsTemp)){
      
      if(round(runif(1,0,1)) == 0){
        saveBlcks_startYear[idxIter,idxYear] <- round(runif(n = 1, 
                                                            min=min(an(fecYears)), 
                                                            max=(max(an(fecYears))-nYearsTemp[idxYear]))) # random number for the starting year
        saveBlcks_rev[idxIter,idxYear] <- 0
      }else{
        saveBlcks_startYear[idxIter,idxYear] <- round(runif(n = 1, 
                                                            min=(min(an(fecYears))+nYearsTemp[idxYear]), 
                                                            max=max(an(fecYears)))) # random number for the starting year
        saveBlcks_rev[idxIter,idxYear] <- 1
      }
    }
  }
  
  # fill in 2 arrays with start and end year for each block
  saveBlcksYrs <- array(NA,dim=c(nits,nyrs,2),dimnames=list(nits=1:nits,blocks=1:nyrs,strtstp=c("start","stop")))
  
  # simply loop on start years and number of years with a randomizer for a possible reverse order 
  for(idxIter in 1:nits){
    for(idxYear in 1:nyrs){
      if(!is.na(saveBlcks_startYear[idxIter,idxYear])){
        if(saveBlcks_rev[idxIter,idxYear] == 1){ # reverse
          tempVal <- c( saveBlcks_startYear[idxIter,idxYear],
                        saveBlcks_startYear[idxIter,idxYear] - saveBlcks_nYears[idxIter,idxYear])
          tempVal[2] <- tempVal[2]+1
        }else{ # normal order
          tempVal <- c(saveBlcks_startYear[idxIter,idxYear],
                       saveBlcks_startYear[idxIter,idxYear] + saveBlcks_nYears[idxIter,idxYear])
          tempVal[2] <- tempVal[2]-1
        }
        
        saveBlcksYrs[idxIter,idxYear,"start"] <- tempVal[1]
        saveBlcksYrs[idxIter,idxYear,"stop"]  <- tempVal[2]
      }
    }
  }
  
  yrStrngsIter <- apply(saveBlcksYrs,1,function(x){mapply(seq,from=na.omit(x[,1]),to=na.omit(x[,2]))})
  
  # if each sampled block was the same length, it creates a matrix instead of a list, gives an error below.
  # so added fix to change any matrices into lists first
  for (i in 1:nits) {
    if (!is.list(yrStrngsIter[[i]]))  {
      tmp <- list()
      for (cC in 1:ncol(yrStrngsIter[[i]])) tmp[[cC]]<- yrStrngsIter[[i]][,cC]
      yrStrngsIter[[i]] <- tmp
      rm(tmp)
    } 
  }
  
  yrStrngs    <- lapply(yrStrngsIter,function(x){do.call(c,x)})
  
  for (i in 1:nits) {
    yrStrngs[[i]] <- yrStrngs[[i]][1:length(projPeriod)]
  }
  
  return(yrStrngs)
}