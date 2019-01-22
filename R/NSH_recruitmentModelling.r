library(FLCore)
library(FLSAM)
#- Load stock / MonteCarlo Stock
load("D:/Repository/ICES_HAWG/wg_HAWG/NSAS/results/NSH_HAWG2018_sf.Rdata")
load("D:\\Repository\\ICES_HAWG\\wg_HAWG\\NSAS\\results\\NSHmc.RData")

#stkMC = monteCarlo stock with 1000 iters
stkMC <- NSHmc

ny <- 20        # number of years to project - Usually 20
dy <- range(stkMC)["maxyear"]       # data year
ay <- dy                            # assessment year
iy <- ay+1                          # initial projections year (also intermediate)
fy <- iy + ny -1                    # final year


#- For 1 iteration only
#NSH.SR  <- FLSR(rec = rec(NSH)[,ac((range(NSH)["minyear"]+1): range(NSH)["maxyear"])],
#              	ssb = ssb(NSH)[,ac((range(NSH)["minyear"])  :(range(NSH)["maxyear"]-1))],
#                model='segreg')
#NSH.RI  <- FLSR(
#              	rec = rec(NSH)[,ac((range(NSH)["minyear"]+1): range(NSH)["maxyear"])],
#              	ssb = ssb(NSH)[,ac((range(NSH)["minyear"])  :(range(NSH)["maxyear"]-1))],
#              	model='ricker')


#NSH.SR <- fmle(NSH.SR)
#NSH.RI <- fmle(NSH.BH)

#- Generate empty object
stkMC.sr <- fmle(as.FLSR(stkMC,model='segreg'))

# now do the actual parameter estimation per iteration of the 1000 iterations
# 85% being ricker, 15% being segreg
itersSR <- sample(1:dims(stkMC)$iter,round(0.15*dims(stkMC)$iter),replace=F)
itersRI <- which(!(1:dims(stkMC)$iter) %in% itersSR)

for (its in itersSR)  iter(params(stkMC.sr),its)  <- params(fmle(FLSR(rec = rec(iter(stkMC,its))[,ac((range(stkMC)["minyear"]+1): range(stkMC)["maxyear"])],
                                                                    	ssb = ssb(iter(stkMC,its))[,ac((range(stkMC)["minyear"])  :(range(stkMC)["maxyear"]-1))],
                                                                      model='segreg')))
for (its in itersRI)  iter(params(stkMC.sr),its)  <- params(fmle(FLSR(rec = rec(iter(stkMC,its))[,ac((range(stkMC)["minyear"]+1): range(stkMC)["maxyear"])],
                                                                    	ssb = ssb(iter(stkMC,its))[,ac((range(stkMC)["minyear"])  :(range(stkMC)["maxyear"]-1))],
                                                                      model='ricker')))

# THIS IS MODIFIED SO THAT AN ARIMA MODEL IS FITTED FOR THE RESIDUALS OF EACH ITERATION
# AND USE TO PRODUCE THE FUTURE DEVIATIONS FOR THE CORRESPONDING ITERATION
# I want to add something that takes autocorrelation in SR relationship into account
# to do this I use an arima model
### S/R residuals - with autocorrelation
rec.res <- residuals(stkMC.sr)[,ac((range(stkMC)["minyear"]):(range(stkMC)["maxyear"]-1))]

# autoregressive model order 1
set.seed(108)
# a list with one model per iteration
arima.fit.lst <- lapply(as.list(1:dims(stkMC)$iter) ,  function(its) {arima(an(iter(rec.res,its)), order = c(1, 0, 0))})

# create autocorrelation in residuals and propagate throughout stock into the future
# from initial year of projections (iy) to last of projections (ny-1)
sr.res  <- make.arma.resid.lst(arima.fit.lst, age = 0, years = iy:(iy + ny-1) , rec.res)

save(itersSR,itersRI,sr.res,stkMC.sr,file=file.path(outPath,"SRresiduals.RData"))

#- in the OM predict recruitment following either ricker or segreg and multiply that prediction with sr.res
# e.g. yr = 2018, then ssb of 2017 produced rec in 2018
for(i in 1:dims(stkMC)$iter){
  if(i %in% itersSR)
    rec <- ifelse(c(ssb(iter(stkMC[,ac(yr-1)],i)))<=params(iter(stkMC.sr,i))["b"],
                  params(iter(stkMC.sr,i))["a"] * c(ssb(iter(stkMC[,ac(yr-1)],i))),
                  params(iter(stkMC.sr,i))["a"] * params(iter(stkMC.sr,i))["b"])
  if(i %in% itersRI)
    rec <- params(iter(stkMC.sr,i))["a"] * c(ssb(iter(stkMC[,ac(yr-1)],i))) * exp(-params(iter(stkMC.sr,i))["b"] * c(ssb(iter(stkMC[,ac(yr-1)],i))))
}
rec     <- rec * exp(sr.res[,ac(yr)])






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

