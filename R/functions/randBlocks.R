randBlocks <- function(fecYears,projPeriod,nits){
  
  nyrs <- length(projPeriod)
  
  saveBlcks_nYears     <- matrix(NA,nrow=nits,ncol=nyrs)
  saveBlcks_startYear  <- matrix(NA,nrow=nits,ncol=nyrs)
  
  # generate random samples for the number of years
  for(idxIter in 1:nits){
    saveBlcks_nYears[idxIter,] <- randNums(nyrs,1,length(fecYears)-1,nyrs)
  }
  
  # randomize the starting years
  for(idxIter in 1:nits){
    nYearsTemp <- saveBlcks_nYears[idxIter,
                                   !is.na(saveBlcks_nYears[idxIter,])] # select only non NaN values
    
    # loop on non NaN values
    for(idxYear in 1:length(nYearsTemp)){
      saveBlcks_startYear[idxIter,idxYear] <- round(runif(n = 1, 
                                                          min=min(an(fecYears)), 
                                                          max=max(an(fecYears))-nYearsTemp[idxYear])) # random number for the starting year
    }
  }
  
  # fill in 2 arrays with start and end year for each block
  saveBlcksYrs <- array(NA,dim=c(nits,nyrs,2),dimnames=list(nits=1:nits,blocks=1:nyrs,strtstp=c("start","stop")))
  
  # simply loop on start years and number of years with a randomizer for a possible reverse order 
  for(idxIter in 1:nits){
    for(idxYear in 1:nyrs){
      tempVal <- c(saveBlcks_startYear[idxIter,idxYear],
                   saveBlcks_startYear[idxIter,idxYear] + saveBlcks_nYears[idxIter,idxYear])
      
      if(round(runif(1,0,1)) == 0) # randomize the inverse orderind of the years
        tempVal <- tempVal[2:1]
      
      saveBlcksYrs[idxIter,idxYear,"start"] <- tempVal[1]
      saveBlcksYrs[idxIter,idxYear,"stop"]  <- tempVal[2]
    }
  }
  
  return(saveBlcksYrs)
}