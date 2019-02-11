make.arma.resid <- function(arima.fit, age, years, nit){
  ##----------------------------
  ## simulate structured recruitment residuals on log-scale
  ## not used as moved to time-varying parameters
  ##
  ## arima fit is a fitted arima model
  ## standard deviation of the log recruitment
  ## age is a numeric scalar age at recruitment
  ## years is a numeric vector of years for the residuals
  ## nit is a numeric scalar for the number of iterations
  ## Note: only works for ARMA components here so no drif, integration
  ##----------------------------
  nyr <- length(years)
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
  ## container matrix
  srDev.mat <- matrix(NA, nrow = nit, ncol = nyr)
  ## sample
  for(i in 1:nit){
    srDev.mat[i, ] <- an(arima.sim(n = nyr,
                                   model = model.list,
                                   sd = sd.innov))
  }
  ## create the FLQuant
  srDev <- FLQuant(an(t(srDev.mat)),
                   dimnames = list(year = years, age = age, iter = 1:nit))
  return(srDev)
}