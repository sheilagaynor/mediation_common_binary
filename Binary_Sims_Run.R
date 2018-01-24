Args <- commandArgs(trailingOnly = TRUE)
lowBd <-as.numeric(Args[[2]]); upBd <-as.numeric(Args[[3]]); n <-as.numeric(Args[[4]])

#Necessary functions, packages
library(boot);library(MASS);library(NORMT3)

source("BinSims_Func.R")

#Start loop across prevalence values
out <- matrix(NA, nrow=2, ncol=11)
for(forInd in lowBd:upBd) {
  set.seed(33)  
  const <-  seq(-2.75,1,0.005)[forInd]
  simMatrix <- matrix(NA, nrow=1000, ncol=8)
  #Store y prevalence
  prev <- inv.logit(const + 0.4*0.4 + .5*inv.logit(0.1 + 0.5*0.4 + 0.4*0.3*0.4) + 0.25*(0.3*0.4)) 
  for (i in 1:1000){ 
    #Simulate Data
    a <- rnorm(n, 0.4, 0.75);  con <- rnorm(n,0.3*a, 0.75) #Exposure, Covariate
    m <- rbinom(n, 1, inv.logit(0.1 + 0.5*a + 0.4*con)); y <- rbinom(n, 1, inv.logit(const + .4*a + .5*m  + .25*con))  #Mediator, Outcome
    cCon<- median(con)
    #Fit regression models
    logitYModel <- glm(y ~ a + m + con, na.action=na.omit, family=binomial(link="logit"))
    logitMModel <- glm(m ~ a + con, na.action=na.omit, family=binomial(link="logit"))
    logitYModelCoef <- coef(logitYModel); logitMModelCoef <- coef(logitMModel)
    #Get proposed estimator
    simMatrix[i,1] <-  indirectEffect(t0=logitYModelCoef[1],t1=logitYModelCoef[2],t2=logitYModelCoef[3],t4=logitYModelCoef[4],
                                      c=cCon,b0=logitMModelCoef[1],b1=logitMModelCoef[2],b2=logitMModelCoef[3])
    simMatrix[i,2] <-    directEffect(t0=logitYModelCoef[1],t1=logitYModelCoef[2],t2=logitYModelCoef[3],t4=logitYModelCoef[4],
                                      c=cCon,b0=logitMModelCoef[1],b1=logitMModelCoef[2],b2=logitMModelCoef[3])
    #Get delta estimator
    simMatrix[i,3] <-  indirectDelta(varMLogitModel=vcov(logitMModel),varYLogitModel=vcov(logitYModel),
                                     t0=logitYModelCoef[1],t1=logitYModelCoef[2],t2=logitYModelCoef[3],t4=logitYModelCoef[4],
                                      c=cCon,b0=logitMModelCoef[1],b1=logitMModelCoef[2],b2=logitMModelCoef[3])
    simMatrix[i,4] <-    directDelta(varMLogitModel=vcov(logitMModel),varYLogitModel=vcov(logitYModel),
                                     t0=logitYModelCoef[1],t1=logitYModelCoef[2],t2=logitYModelCoef[3],t4=logitYModelCoef[4],
                                      c=cCon,b0=logitMModelCoef[1],b1=logitMModelCoef[2],b2=logitMModelCoef[3])
    #Run bootstrap
    IEBoot<- rep(NA,500) ; DEBoot<- rep(NA,500)
    for (j in 1:500){ 
      #Get bootstrap sample
      tempID <- sample(1:n, replace = TRUE)
      aBoot <- a[tempID]; conBoot <- con[tempID]; cConBoot <- median(conBoot)
      mBoot <- m[tempID]; yBoot <- y[tempID]
      #Fit regression models
      logitYModelBoot <- glm(yBoot ~ aBoot + mBoot + conBoot, na.action=na.omit, family=binomial(link="logit"))
      logitMModelBoot <- glm(mBoot ~ aBoot +  conBoot, na.action=na.omit, family=binomial(link="logit"))
      logitYModelCoefBoot <- coef(logitYModelBoot); logitMModelCoefBoot <- coef(logitMModelBoot)
      #Get proposed estimator
      IEBoot[j] <-  indirectEffect(t0=logitYModelCoefBoot[1],t1=logitYModelCoefBoot[2],t2=logitYModelCoefBoot[3],t4=logitYModelCoefBoot[4],
                                        c=cConBoot,b0=logitMModelCoefBoot[1],b1=logitMModelCoefBoot[2],b2=logitMModelCoefBoot[3])
      DEBoot[j] <-    directEffect(t0=logitYModelCoefBoot[1],t1=logitYModelCoefBoot[2],t2=logitYModelCoefBoot[3],t4=logitYModelCoefBoot[4],
                                        c=cConBoot,b0=logitMModelCoefBoot[1],b1=logitMModelCoefBoot[2],b2=logitMModelCoefBoot[3])
    }
    simMatrix[i,5]   <- sd(IEBoot, na.rm = TRUE);     simMatrix[i,6]  <- sd(DEBoot, na.rm = TRUE)
    simMatrix[i,7]   <- mad(IEBoot, na.rm = TRUE);    simMatrix[i,8]  <- mad(DEBoot, na.rm = TRUE)
  }
  out <- rbind( out, c(prev, sd(simMatrix[,1], na.rm = TRUE), sd(simMatrix[,2], na.rm = TRUE),
                       colMeans(simMatrix[,c(1:8)], na.rm = TRUE)) )
}

write.table(out[3:dim(out)[1],], file=paste("binSimn", n, "_up", upBd, ".txt", sep=""),
            sep="\t", row.names=FALSE, quote = FALSE)

