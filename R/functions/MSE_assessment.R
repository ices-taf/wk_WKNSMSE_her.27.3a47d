MSE_assessment <- function(stk0,
                           idx0,
                           sam0.ctrl,
                           escapeRuns,
                           resInit = NULL){
  
  nits <- dims(stk0)$iter
  continueRuns        <- which(!(1:dims(stk0)$iter) %in% escapeRuns)
  if (length(resInit) == 0){  
    res <- FLSAM.MSE(stk0,idx0,sam0.ctrl,return.sam=T) #first year only, we can create an FLSAMs object for all iterations similar to the final NSAS 2018 assessment too if that is needed
  }
  if (length(resInit) != 0){
    res <- FLSAM.MSE(stk0,idx0,sam0.ctrl,starting.sam=resInit,return.sam=T) #here we add the starting values for each iteration
  }
  
  trouble <- data.frame(iter = 1:nits , failure =  unlist(lapply(res , function(x) is.na(x))))
  trouble <- trouble[trouble$failure == T,]
  
  #if we only have one iteration with problems, we try to run that one again, but now without starting values
  if (dim(trouble)[1] == 1){
    tres    <-  try(FLSAM(iter(stk0,trouble$iter),
                          FLIndices(lapply(idx0 , function(x) iter(x,trouble$iter))),
                          sam0.ctrl,silent=T))
    #if that didn't work, we'll just save it as an empty FLSAM object
    if(class(tres)=="try-error"){
      res[[trouble$iter]] <- new("FLSAM")
    } else {
      res[[trouble$iter]] <- tres
    }
  }
  #if we have more iterations with trouble, we can make use of the FLSAM.MSE routine
  if (dim(trouble)[1] > 1){
    resTrouble          <- FLSAM.MSE(iter(stk0,trouble$iter),
                                     FLIndices(lapply(idx0,function(x) iter(x,trouble$iter))),
                                     sam0.ctrl,return.sam=T)
    counter <- 1
    for(ii in trouble$iter){
      if(is.na(resTrouble[[counter]])){
        res[[ii]] <- new("FLSAM")
      } else {
        res[[ii]] <-  resTrouble[[counter]]
      }
      counter <- counter + 1
    }
  }
  
  
  #Update the stock object and store which runs to escape
  stk0                      <- window(stk0,end=max(an(dimnames(stk0@harvest)$year))+1)
  for(ii in 1:nits){
    if(!is.na(res[[ii]]@harvest[1,1,drop=T])){
      iter(stk0@harvest,ii) <- res[[ii]]@harvest[,dimnames(stk0@harvest)$year]
      iter(stk0@stock.n,ii) <- res[[ii]]@stock.n[,dimnames(stk0@harvest)$year]
    } else {
      escapeRuns <- sort(unique(c(escapeRuns,ii)))
    }
  }
  
  ret <- list("stk" = stk0,"resInit" = res, "escapeRuns" = escapeRuns)
  
  return(ret)
}