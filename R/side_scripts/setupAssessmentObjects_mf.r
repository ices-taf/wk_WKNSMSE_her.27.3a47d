require(reshape2)
# setup data path
data.source         <- file.path(".","data")    #Data source, not code or package source!!!

### ============================================================================
### Prepare stock object for assessment
### ============================================================================
#Load object
NSH                 <- readFLStock(file.path(data.source, "index.txt"),no.discards=TRUE,quiet=FALSE)

#Catch is calculated from: catch.wt * catch.n, however, the reported landings are
#normally different (due to SoP corrections). Hence we overwrite the calculate landings
# are we not using catches data then?
NSH@catch           <- NSH@landings
units(NSH)[1:17]    <- as.list(c(rep(c("tonnes","thousands","kg"),4), 
                                 rep("NA",2),"f",rep("NA",2)))

#Set object details
NSH@name                              <- "NSH_HAWG2018"
range(NSH)[c("minfbar","maxfbar")]    <- c(2,6)
NSH                                   <- setPlusGroup(NSH,NSH@range["max"])

#Historical data is only provided for ages 0-8 prior to 1960 (rather than 0-9 after)
#We therefore need to fill in the last ages by applying the following assumptions
#  - weight in the catch and weight in the stock at age 9 is the same as 
#    in period 1960-1983 (constant in both cases)
#  - catch at age reported as 8+ is split evenly between age 8 and 9
#  - natural mortality at age 9 is 0.1 (same as age 8)
#  - proportion mature at age 9 is 1.0 (same as age 8)
#  - harvest.spwn and m.spwn are the same as elsewhere
hist.yrs                      <- as.character(1947:1959)
NSH@catch.wt                  <- NSH@landings.wt #Automatic population of catch.wt introduces NAS
NSH@catch.wt["9",hist.yrs]    <- 0.271
NSH@landings.wt["9",hist.yrs] <- 0.271
NSH@catch.n["9",hist.yrs]     <- NSH@catch.n["8",hist.yrs]/2
NSH@catch.n["8",hist.yrs]     <- NSH@catch.n["9",hist.yrs]
NSH@landings.n["9",hist.yrs]  <- NSH@landings.n["8",hist.yrs]/2
NSH@landings.n["8",hist.yrs]  <- NSH@landings.n["9",hist.yrs]
NSH@stock.wt["9",hist.yrs]    <- 0.312
NSH@m["9",hist.yrs]           <- 0.1
NSH@mat["9",hist.yrs]         <- 1

#No catches of age 9 in 1977 so stock.wt does not get filled there.
#Hence, we copy the stock weight for that age from the previous year.
#Note that because we use a fixed stock.wt prior to 1983, there is no
#need to use averaging or anything fancier
NSH@stock.wt["9","1977"]              <- NSH@stock.wt["9","1976"]

#Use a running mean(y-2,y-1,y) of input wests (i.e. west_raw) to calculate west
NSH@stock.wt[,3:dim(NSH@stock.wt)[2]] <- (NSH@stock.wt[,3:(dim(NSH@stock.wt)[2]-0)] +
                                            NSH@stock.wt[,2:(dim(NSH@stock.wt)[2]-1)] +
                                            NSH@stock.wt[,1:(dim(NSH@stock.wt)[2]-2)]) / 3

### ============================================================================
### Prepare Natural Mortality estimates 
### ============================================================================
#Read in estimates from external file
M2            <- read.csv(file.path(".","data","Smoothed_span50_M_NotExtrapolated_NSASSMS2016.csv"),header=TRUE)

colnames(M2)  <- sub("X","",colnames(M2))
rownames(M2)  <- M2[,1]
M2            <- M2[,-1]# Trim off first column as it contains 'ages'
M2            <- M2[,apply(M2,2,function(x){all(is.na(x))==F})] # keep only years with data


#Extract key data from default assessment
NSHM2           <- NSH
NSHM2@m[]       <- NA
yrs             <- dimnames(NSHM2@m)$year
yrs             <- yrs[which(yrs %in% colnames(M2))]
NSHM2@m[,yrs][] <- as.matrix(M2)

#- Apply 5 year running average
extryrs         <- dimnames(NSHM2@m)$year[which(!dimnames(NSHM2@m)$year %in% yrs)]
extryrsfw       <- extryrs[which(extryrs > max(an(yrs)))]
extryrsbw       <- extryrs[which(extryrs <= max(an(yrs)))]
ages            <- dimnames(NSHM2@m)$age
extrags         <- names(which(apply(M2,1,function(x){all(is.na(x))==T})==T))
yrAver          <- 1
for(iYr in as.numeric(rev(extryrs))){
  for(iAge in ages[!ages%in%extrags]){
    if(iYr %in% extryrsbw) NSHM2@m[ac(iAge),ac(iYr)] <-
        yearMeans(NSHM2@m[ac(iAge),ac((iYr+1):(iYr+yrAver)),],na.rm=T)
    if(iYr %in% extryrsfw) NSHM2@m[ac(iAge),ac(iYr)] <-
        yearMeans(NSHM2@m[ac(iAge),ac((iYr-1):(iYr-yrAver)),],na.rm=T)
  }
}
if(length(extrags)>0){
  for(iAge in extrags)
    NSHM2@m[ac(iAge),]          <- NSHM2@m[ac(as.numeric(min(sort(extrags)))-1),]
}
#Write new M values into the original stock object
addM      <- 0.11 #M profiling based on 2018 benchmark meeting
NSH@m     <- NSHM2@m + addM


### ============================================================================
### Prepare index object for assessment
### ============================================================================
#Load and modify all numbers at age data
NSH.tun   <- readFLIndices(file.path(data.source,"fleet.txt"))
NSH.tun   <- lapply(NSH.tun,function(x) {x@type <- "number"; return(x)})
NSH.tun[["IBTS0"]]@range["plusgroup"] <- NA

# LAI index: read in raw LAI data
surveyLAI     <- read.table(file.path(data.source,"lai.txt"),stringsAsFactors=F,header=T)
ORSH          <- subset(surveyLAI,Area == "Or/Sh")
CNS           <- subset(surveyLAI,Area == "CNS")
BUN           <- subset(surveyLAI,Area == "Buchan")
SNS           <- subset(surveyLAI,Area == "SNS")

#- Put data into FLR format
formatLAI <- function(x,minYr,maxYr){
              x <- dcast(x[,c("Year","LAIUnit","L..9")],Year ~ LAIUnit,value.var="L..9")
              rownames(x) <- x$Year+1900
              missingYears <- (minYr:maxYr)[which(!minYr:maxYr %in% rownames(x))]
              if(length(missingYears)>0){
                xextra <- cbind(Year=missingYears,do.call(cbind,lapply(as.list(1:(ncol(x)-1)),function(y){return(rep(NA,length(missingYears)))})))
                rownames(xextra) <- missingYears
                colnames(xextra) <- colnames(x)
                x <- rbind(x,xextra)
              }
              x <- x[sort.int(rownames(x),index.return=T)$ix,]
              x <- x[,-grep("Year",colnames(x))]
              colnames(x) <- 0:(ncol(x)-1)
              x <- as.matrix(x)
              attr(x,"time") <- c(0.67,0.67)
             return(x)}

ORSH                <- formatLAI(ORSH,1972,range(NSH)["maxyear"]);
CNS                 <- formatLAI(CNS,1972,range(NSH)["maxyear"]);
BUN                 <- formatLAI(BUN,1972,range(NSH)["maxyear"]);
SNS                 <- formatLAI(SNS,1972,range(NSH)["maxyear"]);

FLORSH              <- FLIndex(index=FLQuant(t(ORSH),dimnames=list(age=colnames(ORSH),year=rownames(ORSH),unit="ORSH",season="all",area="unique",iter="1")))
FLCNS               <- FLIndex(index=FLQuant(t(CNS),dimnames=list(age=colnames(CNS),year=rownames(CNS),unit="CNS",season="all",area="unique",iter="1")))
FLBUN               <- FLIndex(index=FLQuant(t(BUN),dimnames=list(age=colnames(BUN),year=rownames(BUN),unit="BUN",season="all",area="unique",iter="1")))
FLSNS               <- FLIndex(index=FLQuant(t(SNS),dimnames=list(age=colnames(SNS),year=rownames(SNS),unit="SNS",season="all",area="unique",iter="1")))
range(FLORSH)[6:7]  <- range(FLCNS)[6:7] <- range(FLBUN)[6:7] <- range(FLSNS)[6:7] <- c(0.67,0.67)
name(FLORSH)        <- "LAI-ORSH"
name(FLCNS)         <- "LAI-CNS"
name(FLBUN)         <- "LAI-BUN"
name(FLSNS)         <- "LAI-SNS"
type(FLORSH)        <- type(FLCNS) <- type(FLBUN) <- type(FLSNS) <- "partial"
FLORSH@index@.Data[which(is.na(FLORSH@index))] <- -1
FLCNS@index@.Data[which(is.na(FLCNS@index))] <- -1
FLBUN@index@.Data[which(is.na(FLBUN@index))] <- -1
FLSNS@index@.Data[which(is.na(FLSNS@index))] <- -1
NSH.tun[(length(NSH.tun)+1):(length(NSH.tun)+4)] <- c(FLORSH,FLBUN,FLCNS,FLSNS)
names(NSH.tun)[(length(NSH.tun)-3):(length(NSH.tun))] <- paste("LAI",c("ORSH","BUN","CNS","SNS"),sep="-")

### ============================================================================
### Apply plusgroup to all data sets
### ============================================================================
pg <- 8

#- This function already changes the stock and landings.wts correctly
NSH <- setPlusGroup(NSH,pg)

caaF  <- read.table("./data/canum_mf.txt",stringsAsFactors=F,header=T)
caA   <- matrix(subset(caaF,Fleet=="A")[,"numbers"],nrow=length(0:9),ncol=length(1997:2017),dimnames=list(age=0:9,year=1997:2017))*1000
caB   <- matrix(subset(caaF,Fleet=="B")[,"numbers"],nrow=length(0:9),ncol=length(1997:2017),dimnames=list(age=0:9,year=1997:2017))*1000
caC   <- matrix(subset(caaF,Fleet=="C")[,"numbers"],nrow=length(0:9),ncol=length(1997:2017),dimnames=list(age=0:9,year=1997:2017))*1000
caD   <- matrix(subset(caaF,Fleet=="D")[,"numbers"],nrow=length(0:9),ncol=length(1997:2017),dimnames=list(age=0:9,year=1997:2017))*1000
cwA   <- matrix(subset(caaF,Fleet=="A")[,"weight"],nrow=length(0:9),ncol=length(1997:2017),dimnames=list(age=0:9,year=1997:2017))
cwB   <- matrix(subset(caaF,Fleet=="B")[,"weight"],nrow=length(0:9),ncol=length(1997:2017),dimnames=list(age=0:9,year=1997:2017))
cwC   <- matrix(subset(caaF,Fleet=="C")[,"weight"],nrow=length(0:9),ncol=length(1997:2017),dimnames=list(age=0:9,year=1997:2017))
cwD   <- matrix(subset(caaF,Fleet=="D")[,"weight"],nrow=length(0:9),ncol=length(1997:2017),dimnames=list(age=0:9,year=1997:2017))

#- Functions to set the plusgroups
setMatrixPlusGroupN <- function(x){
  x[nrow(x)-1,] <- colSums(x[(nrow(x)-1):(nrow(x)),],na.rm=T)
  return(x[-nrow(x),])}
setMatrixPlusGroupWt <- function(x,y){
  y[nrow(x)-1,] <- colSums(x[(nrow(x)-1):(nrow(x)),]*y[(nrow(y)-1):(nrow(y)),],na.rm=T) / colSums(x[(nrow(x)-1):(nrow(x)),],na.rm=T)
  y[is.nan(y)] <- 0
  return(y[-nrow(y),])}


NSH3f <- FLCore::expand(NSH,area=c("A","BD","C"))
NSH3f@catch.n[,ac(1997:range(NSH)["maxyear"]),,,"A"]   <- setMatrixPlusGroupN(caA)
NSH3f@catch.n[,ac(1997:range(NSH)["maxyear"]),,,"BD"]   <- setMatrixPlusGroupN(caB+caD)
NSH3f@catch.n[,ac(1997:range(NSH)["maxyear"]),,,"C"]   <- setMatrixPlusGroupN(caC)
NSH3f@catch.wt[,ac(1997:range(NSH)["maxyear"]),,,"A"]  <- setMatrixPlusGroupWt(caA,cwA)
NSH3f@catch.wt[,ac(1997:range(NSH)["maxyear"]),,,"BD"]  <- setMatrixPlusGroupWt(caB+caD,(cwB*caB + cwD*caD)/(caB+caD))
NSH3f@catch.wt[,ac(1997:range(NSH)["maxyear"]),,,"C"]  <- setMatrixPlusGroupWt(caC,cwC)
NSH3f@landings.wt <- NSH3f@catch.wt
NSH3f@catch.n[,ac(1947:1996),,,-1]    <- NA
NSH3f@catch.wt[,ac(1947:1996),,,-1]   <- NA

#We don't believe the closure catch data, so put it to NA
NSH@catch.n[,ac(1978:1979)]           <- NA
NSHsum <- NSH[,ac(1947:1996)]
NSHres3 <- NSH3f[,ac(1997:2017)]; NSHres3@landings.n[] <- NSHres3@catch.n
NSHs3        <- FLStocks(residual=NSHres3,sum=NSHsum)

NSH.tun[["HERAS"]]@index[ac(pg),]     <- quantSums(NSH.tun[["HERAS"]]@index[ac(pg:dims(NSH.tun[["HERAS"]]@index)$max),])
NSH.tun[["HERAS"]]                    <- trim(NSH.tun[["HERAS"]],age=dims(NSH.tun[["HERAS"]]@index)$min:pg)
NSH.tun[["HERAS"]]@range["plusgroup"] <- pg

#idxSCAI <- which(names(NSH.tun)=="SCAI")
#NSH.tun <- NSH.tun[-idxSCAI]

NSH.tun[["IBTS-Q3"]] <- trim(NSH.tun[["IBTS-Q3"]],age=0:5)
NSH.tun[["IBTS-Q1"]] <- trim(NSH.tun[["IBTS-Q1"]],age=1)
NSH.tun[["HERAS"]] <- trim(NSH.tun[["HERAS"]],age=1:8)

### ============================================================================
### Closure data deletion
### ============================================================================



