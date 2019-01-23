path              <- "D:/git/wk_WKNSMSE_her.27.3a47d/R/"

try(setwd(path),silent=TRUE)

outPath       <- file.path(".","results/")

# load workspace with 
# 1) stock and tuning objects from initial assessment (NSH and NSH.tun)
# 2) stock and tuning random sample objects (stocks and surveys)
load(file.path(outPath,'test.RData'))


#------------------------------------------------------------------------------#
# Test 1: running FLSAM the initial assessment objects. This works fine
#------------------------------------------------------------------------------#
res1 <- FLSAM(NSH,NSH.tun,NSH.ctrl)

warnings() # show warnings
assign("last.warning", NULL, envir = baseenv()) # clear warnings

#------------------------------------------------------------------------------#
# Test 2: running FLSAM with stock and surveys from iter 2 of random samples, 
# including newly generated 2018 data
#
# This gives the following error message:
# Error in solve.default(h, g) : 
# system is computationally singular: reciprocal condition number = 7.52536e-119
#------------------------------------------------------------------------------#
iYr                 <- '2018'
NSH.ctrl@range[5]   <- an(iYr)
NSH.ctrl@residuals  <- F
nits                <- length(stocks)

idxIter <- 2

NSH_test <- window(  stocks[[idxIter]],
                     start=1947,
                     end=an(iYr))

NSH.tun_test <- surveys

surveyNames <- names(NSH.tun_test) # get all the survey names
for(idxSurvey in 1:length(surveys)){
  minYearSurvey     <- min(as.numeric(colnames(NSH.tun_test[[surveyNames[idxSurvey]]]@index)))
  NSH.tun_test[[idxSurvey]] <- window( NSH.tun_test[[idxSurvey]][,,,,,idxIter],
                                       start=minYearSurvey,
                                       end=an(iYr))
}

res2 <- FLSAM(NSH_test,NSH.tun_test,NSH.ctrl)

warnings() # show warnings
assign("last.warning", NULL, envir = baseenv()) # clear warnings

#------------------------------------------------------------------------------#
# Test 3: running FLSAM with stock and surveys from iter 2 of random samples, 
# not including new 2018 data (i.e. up to 2017)
#
# No problem issued
#------------------------------------------------------------------------------#
iYr                 <- '2017'
NSH.ctrl@range[5]   <- an(iYr)
NSH.ctrl@residuals  <- F
nits                <- length(stocks)

idxIter <- 2

NSH_test <- window(  stocks[[idxIter]],
                     start=1947,
                     end=an(iYr))

NSH.tun_test <- surveys

surveyNames <- names(NSH.tun_test) # get all the survey names
for(idxSurvey in 1:length(surveys)){
  minYearSurvey     <- min(as.numeric(colnames(NSH.tun_test[[surveyNames[idxSurvey]]]@index)))
  NSH.tun_test[[idxSurvey]] <- window( NSH.tun_test[[idxSurvey]][,,,,,idxIter],
                                       start=minYearSurvey,
                                       end=an(iYr))
}

res3 <- FLSAM(NSH_test,NSH.tun_test,NSH.ctrl)

warnings() # show warnings
assign("last.warning", NULL, envir = baseenv()) # clear warnings

#------------------------------------------------------------------------------#
# Test 4: running FLSAM with stock from iter 2 of random samples and NSH.tun, 
# including new 2018 data for the stock object
#
# Warnings issued
#------------------------------------------------------------------------------#
iYr                 <- '2018'
NSH.ctrl@range[5]   <- an(iYr)
NSH.ctrl@residuals  <- F
nits                <- length(stocks)

idxIter <- 2

NSH_test <- window(  stocks[[idxIter]],
                     start=1947,
                     end=an(iYr))

NSH.tun_test <- surveys

surveyNames <- names(NSH.tun_test) # get all the survey names
for(idxSurvey in 1:length(surveys)){
  minYearSurvey     <- min(as.numeric(colnames(NSH.tun_test[[surveyNames[idxSurvey]]]@index)))
  NSH.tun_test[[idxSurvey]] <- window( NSH.tun_test[[idxSurvey]][,,,,,idxIter],
                                       start=minYearSurvey,
                                       end=an(iYr))
}

res3 <- FLSAM(NSH_test,NSH.tun,NSH.ctrl)

warnings() # show warnings
assign("last.warning", NULL, envir = baseenv()) # clear warnings
