
rm(list=ls())

library(FLSAM)
library(MASS)
library(msm)
wine <- F

#- Set paths
  codePath        <- "W:/IMARES/Data/ICES-WG/WKHerTAC/NSAS/R/"
  dataPath        <- "W:/IMARES/Data/ICES-WG/WKHerTAC/NSAS/Data/"
  outPath         <- "W:/IMARES/Data/ICES-WG/WKHerTAC/NSAS/Results/"
  NSASPath        <- "W:/IMARES/Data/ICES-WG/WKHerTAC/NSAS/Results/"
  WBSSPath        <- "W:/IMARES/Data/ICES-WG/WKHerTAC/WBSS/Results/"
  combPathCode    <- "W:/IMARES/Data/ICES-WG/WKHerTAC/NSASWBSS/R/"
  combPathResults <- "W:/IMARES/Data/ICES-WG/WKHerTAC/NSASWBSS/Results/"
  combPathData    <- "W:/IMARES/Data/ICES-WG/WKHerTAC/NSASWBSS/Data/"
  FLMETApath      <- "W:/IMARES/Data/ICES-WG/WKHerTAC/FLMeta/"

  #- Paths for Nemo
  if((substr(R.Version()$os,1,3)== "lin" & length(dir("/media"))>0) | (substr(R.Version()$os,1,3)=="min" & length(dir("/media"))>0)){
    codePath      <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",codePath)
    dataPath      <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",dataPath)
    outPath       <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",outPath)
    NSASPath      <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",NSASPath)
    WBSSPath      <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",WBSSPath)
    combPathCode  <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",combPathCode)
    combPathResults<- sub("W:/IMARES/Data/ICES-WG/","/media/w/",combPathResults)
    combPathData  <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",combPathData)
    FLMETApath    <- sub("W:/IMARES/Data/ICES-WG/","/media/w/",FLMETApath)
  }

  #- Source the FLMeta package code
  source(file.path(FLMETApath,"FLMETA.r"))
  source(paste(combPathCode,"04_forecastScenariosNSAS.r", sep="")) #NSAS
  source(paste(combPathCode,"04_forecastScenariosWBSS.r", sep="")) #WBSS


tmpPath<- combPathResults

##-------------------------------------------------------------------------------
## Setup array to save results
##-------------------------------------------------------------------------------

# list of runs to analyse
scens   <-  c(paste("WBSS",1:6,sep=""),paste("NSAS",1:6,sep=""))

for(iScen in scens){

  load(file.path(combPathResults,paste(iScen,".RData",sep="")))
  TAC[,ac(2013),"WB",,"D"] <-   TAC[,ac(2014),"WB",,"D"]
#-------------------------------------------------------------------------------
# 9): Report figures
#-------------------------------------------------------------------------------

yrs   <- 2013:2034
cl    <- 1.1
ca    <- 1
fonts <- 1

projPeriod                                <- 2014:2034
f<-(landings.sel(fishery)) # true fishing mortality at age for both stocks


#- Test ranges
#ceiling(quantile(c(ssbb(biol[,,1],areaSums(f[,,1]))[,ac(yrs)]),probs=c(0.025,0.975))    [2]*1.25/1e6)*1e6
#ceiling(quantile(c(ssbb(biol[,,2],areaSums(f[,,2]))[,ac(yrs)]),probs=c(0.025,0.975))    [2]*1.25/1e6)*1e6
#ceiling(quantile(c(unitSums(computeLandings(fishery)[,ac(yrs),,,"A"])),probs=c(0.025,0.975))[2]*1.25/1e5)*1e5
#ceiling(quantile(c(unitSums(computeLandings(fishery)[,ac(yrs),,,"B"])),probs=c(0.025,0.975))[2]*1.25/1e4)*1e4
#ceiling(quantile(c(unitSums(computeLandings(fishery)[,ac(yrs),,,"C"])),probs=c(0.025,0.975))[2]*1.25/1e5)*1e5
#ceiling(quantile(c(unitSums(computeLandings(fishery)[,ac(yrs),,,"D"])),probs=c(0.025,0.975))[2]*1.25/1e4)*1e4
#ceiling(quantile(c(unitSums(computeLandings(fishery)[,ac(yrs),,,"F"])),probs=c(0.025,0.975))[2]*1.25/1e4)*1e4
#
#ceiling(quantile(c(quantMeans(areaSums(f[ac(2:6),ac(yrs),1]))),probs=c(0.025,0.975))    [2]*1.25*10) /10
#ceiling(quantile(c(quantMeans(areaSums(f[ac(3:6),ac(yrs),2]))),probs=c(0.025,0.975))    [2]*1.25*10) /10



yrangeSSBNSAS <- c(0,4e6);
yrangeSSBWBSS <- c(0,4e5)
yrangeLanA    <- c(0,7e5);
yrangeLanB    <- c(0,40000)
yrangeLanC    <- c(0,1e5);
yrangeLanD    <- c(0,10000)
yrangeLanF    <- c(0,50000)
yrangeFNSAS   <- c(0,0.5)
yrangeFWBSS   <- c(0,0.5)

pdf(file=paste(tmpPath,paste(iScen,"_overview",sep=""),".pdf" ,sep=""))
#pdf(file=paste(tmpPath,paste(iScen,"_","0.3150000Blim",sep=""),".pdf" ,sep=""))
#pdf(file=paste(tmpPath,paste(iScen,"_","0.251500000",sep=""),".pdf" ,sep=""))

par(mfrow=c(3,1),oma=c(6,6,2,6),mar=c(1,0.5,0,0.5))

  #---------
  #- Landings
  #---------
rLandfTot <- apply(areaSums(computeLandings(fishery[,ac(yrs)]))@.Data,1:5,quantile,probs=c(0.025,0.5,0.975),na.rm=T)
#rLands    <- apply(computeLandings(stock[,ac(yrs)])@.Data,1:5,quantile,probs=c(0.025,0.5,0.975),na.rm=T)
#rTAC      <- apply(TAC[,ac(yrs)]@.Data,1:5,quantile,probs=c(0.025,0.5,0.975),na.rm=T)

xrange  <- range(yrs)
yrange1 <- c(0,max(pretty(max(rLandfTot[,,,1,,]))))
yrange2 <- c(0,max(pretty(max(rLandfTot[,,,2,,]))))

plot(y=rLandfTot["50%",,ac(yrs),1,,],x=yrs,xlim=xrange,ylim=yrange1,xlab="",ylab="",
     type="l",lwd=1,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n")
axis(2,las=1,at=pretty(yrange1),labels=pretty(yrange1)/1000,cex=ca,font=fonts)
mtext(text=expression(paste("Landings (",10^3," tonnes)",sep="")),side=2,at=(yrange1[2]-yrange1[1])/2+yrange1[1],outer=F,cex=cl,line=4,font=fonts)
grid(); box()
lines(y=rLandfTot["50%",,ac(yrs),1,,],x=yrs,lty=1,lwd=2)
lines(y=rLandfTot["2.5%",,ac(yrs),1,,],x=yrs,lty=3,lwd=1)
lines(y=rLandfTot["97.5%",,ac(yrs),1,,],x=yrs,lty=3,lwd=1)
par(new=T)
plot(y=rLandfTot["50%",,ac(yrs),2,,],x=yrs,xlim=xrange,ylim=yrange2,xlab="",ylab="",
     type="l",lwd=1,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n",col="red")
axis(4,las=1,at=pretty(yrange2),labels=pretty(yrange2)/1000,cex=ca,font=fonts,col="red",col.ticks="red",col.axis="red")
lines(y=rLandfTot["50%",,ac(yrs),2,,],x=yrs,lty=1,lwd=2,col="red")
lines(y=rLandfTot["2.5%",,ac(yrs),2,,],x=yrs,lty=3,lwd=1,col="red")
lines(y=rLandfTot["97.5%",,ac(yrs),2,,],x=yrs,lty=3,lwd=1,col="red")
legend("bottomleft",legend=c("Area IV","Area IIIa & 22-24"),col=c("black","red"),lwd=2,lty=1,box.lty=0)
text(x=xrange[1],y=yrange2[2],pos=1,labels="(A)",font=fonts,cex=cl)

  #------------------
  # True F
  #------------------
rFpNSAS   <- apply(apply(areaSums(f)[ac(2:6),,1],2:6,mean,na.rm=T)@.Data,1:5,quantile,probs=c(0.025,0.5,0.975),na.rm=T)
rFpWBSS   <- apply(apply(areaSums(f)[ac(3:6),,2],2:6,mean,na.rm=T)@.Data,1:5,quantile,probs=c(0.025,0.5,0.975),na.rm=T)
#rFs   <- apply(fbar(stocks)@.Data,1:5,quantile,probs=c(0.025,0.5,0.975),na.rm=T)
#rFstf <- apply(apply(unitSums(fSTF)[ac(2:6),],2:6,mean,na.rm=T)@.Data,1:5,quantile,probs=c(0.025,0.5,0.975),na.rm=T)


xrange  <- range(yrs)
plot(rFpNSAS["50%",,ac(yrs),,,]~yrs,xlim=xrange,ylim=yrangeFNSAS,xlab="",ylab="",
     type="l",lwd=1,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n")
abline(h=0.26,col="blue",lty=2,lwd=2)
text(x=yrs[1]+0.5,y=0.26,labels="Fmsy NSAS",las=1,col="blue",font=fonts,pos=1)
axis(2,las=1,at=pretty(yrangeFNSAS),labels=pretty(yrangeFNSAS),cex=ca,font=fonts)
mtext(text=expression(paste(F[2-6]," (",year^-1,")",sep="")),side=2,at=(yrangeFNSAS[2]-yrangeFNSAS[1])/2+yrangeFNSAS[1],outer=F,cex=cl,line=4,font=fonts)
grid(); box()
lines(rFpNSAS["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2)
lines(rFpNSAS["2.5%",,ac(yrs),,,]~yrs,lty=3,lwd=1)
lines(rFpNSAS["97.5%",,ac(yrs),,,]~yrs,lty=3,lwd=1)
par(new=T)
plot(rFpWBSS["50%",,ac(yrs),,,]~yrs,xlim=xrange,ylim=yrangeFWBSS,xlab="",ylab="",
     type="l",lwd=1,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n",col="red")
abline(h=0.28,col="darkgreen",lty=2,lwd=2)
text(x=yrs[11]-0.5,y=0.28,labels="Fmsy WBSS",las=1,col="darkgreen",font=fonts,pos=3)
lines(rFpWBSS["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2,col="red")
lines(rFpWBSS["2.5%",,ac(yrs),,,]~yrs,lty=3,lwd=1,col="red")
lines(rFpWBSS["97.5%",,ac(yrs),,,]~yrs,lty=3,lwd=1,col="red")
axis(4,las=1,at=pretty(yrangeFWBSS),labels=pretty(yrangeFWBSS),cex=ca,font=fonts,col="red",col.ticks="red",col.axis="red")
mtext(text=expression(paste(F[3-6]," (",year^-1,")",sep="")),side=4,at=(yrangeFWBSS[2]-yrangeFWBSS[1])/2+yrangeFWBSS[1],outer=F,cex=cl,line=4,font=fonts,col="red")
text(x=xrange[1],y=yrangeFWBSS[2],pos=1,labels="(B)",font=fonts,cex=cl)


  #------------------
  # True SSB
  #------------------
rSSBpNSAS <- apply(ssbb(biol[,,1],areaSums(f[,,1]))@.Data,1:5,quantile,probs=c(0.025,0.5,0.975),na.rm=T)
rSSBpWBSS <- apply(ssbb(biol[,,2],areaSums(f[,,2]))@.Data,1:5,quantile,probs=c(0.025,0.5,0.975),na.rm=T)


xrange  <- range(yrs)

plot(rSSBpNSAS["50%",,ac(yrs),,,]~yrs,xlim=xrange,ylim=yrangeSSBNSAS,xlab="",ylab="",
     type="l",lwd=1,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n")
axis(2,las=1,at=pretty(yrangeSSBNSAS),labels=pretty(yrangeSSBNSAS)/1000,cex=ca,font=fonts)
axis(1,las=1,at=pretty(xrange),labels=pretty(xrange),cex=ca,font=fonts)
mtext(text=expression(paste("SSB (",10^3," tonnes)",sep="")),side=2,at=(yrangeSSBNSAS[2]-yrangeSSBNSAS[1])/2+yrangeSSBNSAS[1],outer=F,cex=cl,line=4,font=fonts)
abline(h=0.8e6,col="blue",lwd=2,lty=2);
text(x=yrs[1]+0.5,y=0.8e6,labels="Blim NSAS",las=1,col="blue",font=fonts,pos=1)
grid(); box()
lines(rSSBpNSAS["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2)
lines(rSSBpNSAS["2.5%",,ac(yrs),,,]~yrs,lty=3,lwd=1)
lines(rSSBpNSAS["97.5%",,ac(yrs),,,]~yrs,lty=3,lwd=1)
  #-Reference level Blim & Bpa
par(new=T)
plot(rSSBpWBSS["50%",,ac(yrs),,,]~yrs,xlim=xrange,ylim=yrangeSSBWBSS,xlab="",ylab="",
     type="l",lwd=1,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n",col="red")
grid(); box()
axis(4,las=1,at=pretty(yrangeSSBWBSS),labels=pretty(yrangeSSBWBSS)/1000,cex=ca,font=fonts,col="red",col.ticks="red",col.axis="red")
lines(rSSBpWBSS["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2,col="red")
lines(rSSBpWBSS["2.5%",,ac(yrs),,,]~yrs,lty=3,lwd=1,col="red")
lines(rSSBpWBSS["97.5%",,ac(yrs),,,]~yrs,lty=3,lwd=1,col="red")
abline(h=90000,col="darkgreen",lwd=2,lty=2)
text(x=yrs[11]-0.5,y=90000,labels="Blim WBSS",las=1,col="darkgreen",font=fonts,pos=3)
text(x=xrange[1],y=yrangeSSBWBSS[2],pos=1,labels="(C)",font=fonts,cex=cl)

  #- Labels x-axis
mtext(text=expression(Years),side=1,at=0.5,outer=T,cex=cl,line=2.5,font=fonts)
#mtext(at=c(0.25,0.75),text=c("Blim","Btrigger"),outer=T,line=0.2,cex=1.1,font=fonts,side=3)
#mtext(at=c(0.25,0.75),text=c("Fmsy = 0.2","Fmsy = 0.2"),outer=T,line=0.2,cex=1.1,font=fonts,side=3)
dev.off()
#savePlot(paste(tmpPath,paste(iScen,"_",fls22scen[idxFiles[1]],sep=""),".png" ,sep=""),type="png")

}

for(iScen in scens){

  load(file.path(combPathResults,paste(iScen,".RData",sep="")))


pdf(file=paste(tmpPath,paste(iScen,"_TACS",sep=""),".pdf" ,sep=""))
par(mfrow=c(3,2),oma=c(6,6,2,6),mar=c(1,2.5,0,2.5))
rLandf    <- apply(unitSums(computeLandings(fishery[,ac(yrs)]))@.Data,1:5,quantile,probs=c(0.025,0.5,0.975),na.rm=T)
TAC[,ac(yrs),1,,c("C","D")] <- 0
rTACf     <- apply(TAC[,ac(yrs)]@.Data,       1:5,quantile,probs=c(0.025,0.5,0.975),na.rm=T)
rOutf     <- apply(outtakeTAC[,ac(yrs[-1])]@.Data,1:5,quantile,probs=c(0.025,0.5,0.975),na.rm=T)
rTACf[is.na(rTACf)==T] <- 0
rLandf[is.na(rLandf)==T] <- 0
rOutf[is.na(rOutf)==T] <- 0

  #---------
  #- Landings by fleet & species & TAC by fleet
  #---------
xrange  <- range(yrs)
for(iFsh in c("A","B","C","D","F")){
  yrange2 <- c(0,max(pretty(max(rLandf["50%",,,,,iFsh]))))
  plot(y=rLandf["50%",,ac(yrs),,,iFsh],x=yrs,xlim=xrange,ylim=yrange2,xlab="",ylab="",
       type="l",lwd=1,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n")
  lines(y=rTACf["50%",,ac(yrs),1,,iFsh],    x=yrs,    lty=2)
  lines(y=rOutf["50%",,ac(yrs[-1]),1,,iFsh],x=yrs[-1],lty=3)
  axis(2,las=1,at=pretty(yrange2),labels=pretty(yrange2)/1000,cex=ca,font=fonts)
  if(iFsh %in% c("A","C","F"))
    mtext(text=expression(paste("Landings (",10^3," tonnes)",sep="")),side=2,at=(yrange2[2]-yrange2[1])/2+yrange2[1],outer=F,cex=cl,line=4,font=fonts)
  grid(); box()
  par(new=T)
  yrange2 <- c(0,max(pretty(max(c(rOutf["50%",,ac(yrs[-1]),2,,iFsh],rTACf["50%",,ac(yrs),2,,iFsh])))))
  plot(y=rTACf["50%",,ac(yrs),2,,iFsh],x=yrs,xlim=xrange,ylim=yrange2,xlab="",ylab="",
       type="l",lwd=1,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n",col="red",lty=2)
  lines(y=rOutf["50%",,ac(yrs[-1]),2,,iFsh],x=yrs[-1],lty=3,col="red")
  axis(4,las=1,at=pretty(yrange2),labels=pretty(yrange2)/1000,cex=ca,font=fonts,col="red",col.ticks="red",col.axis="red")
  text(x=xrange[1],y=yrange2[2],pos=1,labels=paste("",iFsh,"",sep=""),font=fonts,cex=cl)
  mtext(text=expression(Years),side=1,at=0.5,outer=T,cex=cl,line=2.5,font=fonts)
  if(iFsh %in% c("D","F"))
    axis(1,las=1,at=pretty(xrange),labels=pretty(xrange),cex=ca,font=fonts)
}
plot.new()
legend("center",legend=c("Landings","Advised catch NSAS","Advised catch WBSS","Outtake NSAS","Outtake WBSS"),col=c(rep("black",2),"red","black","red"),lwd=1,lty=c(1,2,2,3,3),box.lty=0)
dev.off()

}
