#-------------------------------------------------------------------------------
# WKHELP
#
# Author: Niels Hintzen
#         IMARES, The Netherland
#
# Performs an MSE of North Sea Herring under different TAC scenario's
#
# Date: 02-Sep-2012
#
# Build for R2.13.2
#-------------------------------------------------------------------------------

rm(list=ls())

library(FLSAM)
library(MASS)
library(msm)
wine <- F

path          <- "I:/WKHELP/"
inPath        <- "I:/WKHELP/data/"
codePath      <- "I:/WKHELP/R/"
outPath       <- "I:/WKHELP/Results/"
if(substr(R.Version()$os,1,3)== "lin"){
  path        <- sub("N:/","/media/n/",path)
  inPath      <- sub("N:/","/media/n/",inPath)
  codePath    <- sub("N:/","/media/n/",codePath)
  outPath     <- sub("N:/","/media/n/",outPath)
}
tmpPath <- "/home/hintz001/WKHELP/"

##-------------------------------------------------------------------------------
## Setup array to save results
##-------------------------------------------------------------------------------


lngt                                      <- length(grep("2022",dir(tmpPath)))
projPeriod                                <- 2013:2022
diags                                     <- matrix(NA,nrow=21,ncol=lngt,dimnames=list(criteria=c("Scenario","Btrigger","FA","FJ","Stability","PA","2022SSB","2022F","meanF01","meanF26","meanSSB",
                                                                                                  "High Yield A","High Yield B","2013 Yield A","2013 Yield B","meanrelTACIAV A",
                                                                                                  "meanTACIAV A","IAVrestrict up#","IAVrestrict down#","TACup","TACdown"),
                                                                                       1:lngt))



##-------------------------------------------------------------------------------
## Load results
##-------------------------------------------------------------------------------

counter <- 1
for(iScen in "FIAV"){# c("Bank","BB","F3years","FIAV","juv","LTMP","meanTAC","noIAV")){
  fls <- dir(tmpPath)
  for(iFile in fls[grep(iScen,fls)][grep("2022",fls[grep(iScen,fls)])][4:5]){
    #print(file.path(tmpPath,iFile))
    print(counter)
    load(file.path(tmpPath,iFile))

diags["Scenario",counter]   <- iScen
diags["Btrigger",counter]   <- mpPoints$Btrigger/1e6
diags["FA",counter]         <- mpPoints$FadultA
diags["FJ",counter]         <- mpPoints$FjuvA
diags["Stability",counter]  <- mpPoints$stabilityBreak


##-------------------------------------------------------------------------------
## Diagnostics on results
##-------------------------------------------------------------------------------
diags["PA",counter]                   <- length(unique(which(ssbb(biol,unitSums(f))[,ac(projPeriod)]<0.8e6,arr.ind=T)[,"dim6"]))/dims(biol)$iter  *100
diags["2022SSB",counter]              <- round(c(apply(ssbb(biol,unitSums(f))[,ac(2022)],1:5,median,na.rm=T)))
diags["meanSSB",counter]              <- round(median(c(apply(ssbb(biol,unitSums(f))[,ac(projPeriod[1]:projPeriod[length(projPeriod)-1])],3:6,mean,na.rm=T))))
diags["meanF01",counter]              <- round(median(c(apply(apply(unitSums(f)[ac(0:1),ac(projPeriod[1]:projPeriod[length(projPeriod)-1])],2:6,mean),    3:6,mean,na.rm=T))),3)
diags["meanF26",counter]              <- round(median(c(apply(apply(unitSums(f)[ac(2:6),ac(projPeriod[1]:projPeriod[length(projPeriod)-1])],2:6,mean),    3:6,mean,na.rm=T))),3)
diags["2022F",counter]                <- round(apply(apply(unitSums(f)[ac(2:6),ac(2022)],2:6,mean),1:5,median,na.rm=T),3)
diags["High Yield A",counter]         <- round(median(c(apply(computeLandings(fishery[,ac(projPeriod[1]:projPeriod[length(projPeriod)-1]),  "A"]),3:6,mean,na.rm=T))))
diags["High Yield B",counter]         <- round(median(c(apply(computeLandings(fishery[,ac(projPeriod[1]:projPeriod[length(projPeriod)-1]),  "B"]),3:6,mean,na.rm=T))))
diags["2013 Yield A",counter]         <- round(apply(computeLandings(fishery[,ac(2013),  "A"]),2,median,na.rm=T))
diags["2013 Yield B",counter]         <- round(apply(computeLandings(fishery[,ac(2013),  "B"]),2,median,na.rm=T))
diags["meanrelTACIAV A",counter]      <- round(median(c(apply(abs(TAC[,ac(projPeriod[2]:rev(projPeriod)[2]),"A"] - TAC[,ac(projPeriod[1]:rev(projPeriod)[3]),"A"]) /TAC[,ac(projPeriod[2]:rev(projPeriod)[2]),"A"] * 100,3:6,mean,na.rm=T))),3)
diags["meanTACIAV A",counter]         <- round(median(c(apply(abs(TAC[,ac(projPeriod[2]:rev(projPeriod)[2]),"A"] - TAC[,ac(projPeriod[1]:rev(projPeriod)[3]),"A"]) ,3:6,mean,na.rm=T))))

IAVUp   <- which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)-1]),"A"] == 1.15* TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-2]),"A"],arr.ind=T)
IAVDown <- which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)-1]),"A"] == 0.85* TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-2]),"A"],arr.ind=T)

  #- Average number of times the IAV rule is applied upwards or downwards
if((nrow(IAVUp)) > 0 ){
  a <- IAVUp
  diags["IAVrestrict up#",counter]    <- median(aggregate(a[,"dim2"],by=list(a[,"dim6"]),function(x){length(x)})$x)/10*100
}
if((nrow(IAVDown)) > 0 ){
  a <- IAVDown
  diags["IAVrestrict down#",counter]  <- median(aggregate(a[,"dim2"],by=list(a[,"dim6"]),function(x){length(x)})$x)/10*100
}

  #- Which TAC of the runs go up and which go down
resUp   <- which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)-1]),"A"] > TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-2]),"A"],arr.ind=T)
resDown <- which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)-1]),"A"] < TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-2]),"A"],arr.ind=T)

  #- Mean increase in TAC is TAC goes up, or mean decrease in TAC is TAC goes down
diags["TACup",counter]   <- round(mean((TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)-1]),"A"] - TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-2]),"A"])@.Data[which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)-1]),"A"] > TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-2]),"A"])]))
diags["TACdown",counter] <- round(mean((TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-2]),"A"] - TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)-1]),"A"])@.Data[which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)-1]),"A"] < TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-2]),"A"])]))

counter <- counter + 1
}
}

write.csv(t(diags),file=paste(outPath,"tables_diags.csv",sep=""),row.names=T)
