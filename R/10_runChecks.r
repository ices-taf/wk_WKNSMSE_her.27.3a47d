iStocks <- stocks[,1:idxyrmin1]
iSurvivors <- survivors
iFishery <- tmp_fishery
iTAC <- TAC
iScenario <- scen
iFailRuns <- NULL
iHistMaxYr <- histMaxYr




#- Checks
0.15/(0.7*0.6) +0.1
0.15/0.7*((ssb-0.8e6)/1e6)+0.1


#- Biology versus stocks
ssb(stocks[,projPeriod])
ssb(stocks[,projPeriod]) / ssbb(biol,unitSums(f))[,projPeriod]

#- True Fs
quantMeans(unitSums(f[ac(0:1),projPeriod]))
quantMeans(unitSums(f[ac(2:6),projPeriod]))

#- Forecasted Fs
quantMeans(unitSums(fSTF[ac(0:1),projPeriod]))
quantMeans(unitSums(fSTF[ac(2:6),projPeriod]))

#- TAC implementation
computeCatch(fishery)[,projPeriod] / TAC[,projPeriod]

#- TAC variability
range(TAC[,ac(2013:2022),1] / TAC[,ac(2012:2021),1],na.rm=T)



HCRTAC[,projPeriod] - TAC[,projPeriod] #-noIAV
TAC[,ac(2013:2022),1] - ((TAC[,ac(2012:2021),1] + HCRTAC[,ac(2013:2022),1]) / 2) #meanTAC

idxssb  <- which(HCRSSB[,ac(2015)] > 1.5e6)
quantMeans(unitSums(fSTF[ac(2:6),ac(2015),,,,idxssb])) #FIAV

TAC[,ac(2013),1] - 0.9 * HCRTAC[,ac(2013),1]
TAC[,ac(2014),1] - ((TAC[,ac(2013),1] / 0.9 - TAC[,ac(2013),1]) + (1.1 * HCRTAC[,ac(2014),1]))
TAC[,ac(2015),1] - ((TAC[,ac(2014),1] / 1.1 - TAC[,ac(2014),1]) + (1.1 * HCRTAC[,ac(2015),1]))



## check initial single/multi-fleet NSH.sam objects against those from HAWG2018

