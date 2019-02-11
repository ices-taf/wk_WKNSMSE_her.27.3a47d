make.arma.resid.lst <- function(arima.fit.list, age, years , recres){
  
  #    arima.fit.list= arima.fit.lst ;  age= 0 ; years= 1:10  ; recres = rec.res
  
  ##----------------------------
  ## modified by TB for 2018 MSE :
  ## simulate structured recruitment residuals on log-scale
  ## separately for each iteration of the stock
  ## in order to avoid discontinuity at the start of the simulated series, a selection of
  ## simulated series is operated for each iteration by keeping only those simulated
  ## data that start a +-10% of the last ofbserved SRdev
  ##
  ## original comment :
  ## arima.fit.list is a list of fitted arima models (each corresponding to an iteration
  ## standard deviation of the log recruitment
  ## age is a numeric scalar age at recruitment
  ## years is a numeric vector of years for the residuals
  ## nit is a numeric scalar for the number of iterations
  ## Note: only works for ARMA components here so no drif, integration
  ##----------------------------
  nyr <- length(years)
  nit <- length(arima.fit.list)
  ## container matrix
  srDev.mat <- matrix(NA, nrow = nit, ncol = nyr)
  
  for(i in 1:nit){
    arima.fit <-arima.fit.list[[i]]
    ## coefficients
    theta <- coef(arima.fit)
    ## model set up for simulations, not fully general
    model.list <- list()
    if("ar1" %in% names(theta)){
      model.list$ar <- an(theta[grep("ar", names(theta))])
    }
    if("ma1" %in% names(theta)){
      model.list$ma <- an(theta[grep("ma", names(theta))])
    }
    ## standard deviations of the innovations
    sd.innov <- sqrt(arima.fit$sigma2)
    
    ## last observed deviation
    srdev  <- iter(recres,i)
    srdev  <- tail(c(srdev),1)
    
    ## sample (until the output is close enough to the last observed)
    OK <- F
    while(!OK)
    {
      isrDev      <- an(arima.sim(n = nyr,
                                  model = model.list,
                                  sd = sd.innov))
      
      diff<-exp(isrDev[1]) - exp(srdev)
      if (abs(diff) < 0.1 )  OK <- T
    }
    
    
    srDev.mat[i, ]  <- isrDev
  }
  ## create the FLQuant
  srDev <- FLQuant(an(t(srDev.mat)),
                   dimnames = list(year = years, age = age, iter = 1:nit))
  return(srDev)
}