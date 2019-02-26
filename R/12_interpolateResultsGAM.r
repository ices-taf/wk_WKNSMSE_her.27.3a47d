library(mgcv)
library(lattice)

#- Read in data
risk  <- as.matrix(read.csv("D:/Downloads/hertest.csv",header=T,stringsAsFactor=F))
colnames(risk)[-1] <- seq(1e6,2e6,0.2e6)
rownames(risk) <- risk[,1]
risk  <- risk[,-1]

#- Turn data into data.frame
risk            <- as.data.frame(as.table(risk))
colnames(risk)  <- c("Ftarget","Btrigger","Risk")
risk$Ftarget    <- as.numeric(as.character(risk$Ftarget))
risk$Btrigger   <- as.numeric(as.character(risk$Btrigger))

#- Predict data
predData <- expand.grid(Ftarget = seq(0.15,0.35,0.01),Btrigger=seq(0.8e6,2.4e6,0.1e6))
predData$predRisk <- predict.gam(gam(Risk ~  s(Ftarget,Btrigger),data=risk),newdata=predData)
predData        <- predData[which(predData$predRisk >= 0.05),]

#- Plot the data
levelplot(predRisk ~  Ftarget * Btrigger, data=predData,main="Risk")