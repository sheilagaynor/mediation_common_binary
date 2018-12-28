Args <- commandArgs(trailingOnly = TRUE)
lowBd <-as.numeric(Args[[2]]); upBd <-as.numeric(Args[[3]]); n <-as.numeric(Args[[4]])

#Necessary functions, packages
library(boot);library(MASS);library(NORMT3)

source("ContSims_Func.R")

#Start loop across prevalence values
out <- matrix(NA, nrow=2, ncol=17)
for(forInd in lowBd:upBd) { 
  set.seed(33)  
  const <-  seq(-5,1.5,0.005)[forInd]
  simMatrix <- matrix(NA, nrow=1000, ncol=12)
  trueT0 <- const;  trueT1 <- 0.4;  trueT2 <- 0.5;  trueT4 <- 0.25
  trueB0 <- 0.1; trueB1 <- 0.5; trueB2 <- 0.4; trueSig <- 0.75
  trueA <- 0.4; trueC <- 0.3*trueA
  #Store y prevalence
  prevVec <- rep(NA,1000) 
  #Approximation method to get truth
  approxIE <- exp(  logit(integrate(ie1Function,-Inf, Inf,t0=trueT0,t1=trueT1,t2=trueT2,
                                    t4=trueT4,c=trueC,sigma=trueSig,b0=trueB0,b1=trueB1,b2=trueB2)$value) -
                      logit(integrate(ie2Function,-Inf, Inf,t0=trueT0,t1=trueT1,t2=trueT2,
                                      t4=trueT4,c=trueC,sigma=trueSig,b0=trueB0,b1=trueB1,b2=trueB2)$value) )
  approxDE <- exp( logit(integrate(de1Function,-Inf, Inf,t0=trueT0,t1=trueT1,t2=trueT2,
                                   t4=trueT4,c=trueC,sigma=trueSig,b0=trueB0,b1=trueB1,b2=trueB2)$value) -
                     logit(integrate(de2Function,-Inf, Inf,t0=trueT0,t1=trueT1,t2=trueT2,
                                     t4=trueT4,c=trueC,sigma=trueSig,b0=trueB0,b1=trueB1,b2=trueB2)$value) )
  for (i in 1:1000){ 
    #Simulate Data
    a <- rnorm(n, trueA, trueSig);  con <- rnorm(n,0.3*trueA, trueSig) #Exposure, Covariate
    m <- rnorm(n, trueB0 + trueB1*a + trueB2*con, trueSig); 
    y <- rbinom(n, 1, inv.logit(trueT0 + trueT1*a + trueT2*m  + trueT4*con))  #Mediator, Outcome
    cCon<- median(con)
    prevVec[i] <- mean(y)
    #Fit regression models
    logitModel <- glm(y ~ a + m  + con, na.action=na.omit, family=binomial(link="logit"))
    probitModel <- glm(y ~ a + m  + con, na.action=na.omit, family=binomial(link="probit"))
    sEst <- median(coef(probitModel)/coef(logitModel), na.rm = TRUE)
    linearModel <- lm(m ~ a + con, na.action=na.omit)
    logitModelCoef <- coef(logitModel); linearModelCoef <- coef(linearModel)
    #Get proposed estimator
    simMatrix[i,1] <-  indirectEffectProposed(s=sEst,t0=logitModelCoef[1],t1=logitModelCoef[2],t2=logitModelCoef[3],t4=logitModelCoef[4],
                                              c=cCon,sigma=sqrt(mse(linearModel)),b0=linearModelCoef[1],b1=linearModelCoef[2],b2=linearModelCoef[3])
    simMatrix[i,2] <-    directEffectProposed(s= sEst,t0=logitModelCoef[1],t1=logitModelCoef[2],t2=logitModelCoef[3],t4=logitModelCoef[4],
                                              c=cCon,sigma=sqrt(mse(linearModel)),b0=linearModelCoef[1],b1=linearModelCoef[2],b2=linearModelCoef[3])
    #Get delta estimator
    simMatrix[i,3] <-  indirectDelta(varLinModel=vcov(linearModel),varLogitModel=vcov(logitModel),
                                     s= sEst,t0=logitModelCoef[1],t1=logitModelCoef[2],t2=logitModelCoef[3],t4=logitModelCoef[4],
                                     c=cCon,sigma=sqrt(mse(linearModel)),b0=linearModelCoef[1],b1=linearModelCoef[2],b2=linearModelCoef[3])
    simMatrix[i,4] <-    directDelta(varLinModel=vcov(linearModel),varLogitModel=vcov(logitModel),
                                     s= sEst,t0=logitModelCoef[1],t1=logitModelCoef[2],t2=logitModelCoef[3],t4=logitModelCoef[4],
                                     c=cCon,sigma=sqrt(mse(linearModel)),b0=linearModelCoef[1],b1=linearModelCoef[2],b2=linearModelCoef[3])
    #See coverage of delta CI
    simMatrix[i,5] <-  as.numeric((Re(simMatrix[i,1] - 1.96*simMatrix[i,3]) <= approxIE ) & (Re(simMatrix[i,1] + 1.96*simMatrix[i,3]) >= approxIE ))
    simMatrix[i,6]  <- as.numeric((Re(simMatrix[i,2] - 1.96*simMatrix[i,4]) <= approxDE ) & (Re(simMatrix[i,2] + 1.96*simMatrix[i,4]) >= approxDE ))
    #Run bootstrap
    IEBoot<- rep(NA,500) ; DEBoot<- rep(NA,500)
    for (j in 1:500){
      #Get bootstrap sample
      tempID <- sample(1:n, replace = TRUE)
      aBoot <- a[tempID]; conBoot <- con[tempID]; cConBoot <- median(conBoot)
      mBoot <- m[tempID]; yBoot <- y[tempID]
      #Fit regression models
      logitModelBoot <- glm(yBoot ~ aBoot + mBoot + conBoot, na.action=na.omit, family=binomial(link="logit"))
      probitModelBoot <- glm(yBoot ~ aBoot + mBoot + conBoot, na.action=na.omit, family=binomial(link="probit"))
      sEstBoot <- median(coef(probitModelBoot)/coef(logitModelBoot), na.rm = TRUE)
      linearModelBoot <- lm(mBoot ~ aBoot + conBoot, na.action=na.omit)
      logitModelCoefBoot <- coef(logitModelBoot); linearModelCoefBoot <- coef(linearModelBoot)
      #Get proposed estimator
      IEBoot[j] <-  indirectEffectProposed(s= sEstBoot,t0=logitModelCoefBoot[1],t1=logitModelCoefBoot[2],t2=logitModelCoefBoot[3],t4=logitModelCoefBoot[4],
                                           c=cConBoot,sigma=sqrt(mse(linearModelBoot)),b0=linearModelCoefBoot[1],b1=linearModelCoefBoot[2],b2=linearModelCoefBoot[3])
      DEBoot[j] <-    directEffectProposed(s= sEstBoot,t0=logitModelCoefBoot[1],t1=logitModelCoefBoot[2],t2=logitModelCoefBoot[3],t4=logitModelCoefBoot[4],
                                           c=cConBoot,sigma=sqrt(mse(linearModelBoot)),b0=linearModelCoefBoot[1],b1=linearModelCoefBoot[2],b2=linearModelCoefBoot[3])
    }
    simMatrix[i,7]   <- sd(IEBoot, na.rm = TRUE);     simMatrix[i,8]  <- sd(DEBoot, na.rm = TRUE)
    simMatrix[i,9]   <- mad(IEBoot, na.rm = TRUE);    simMatrix[i,10]  <- mad(DEBoot, na.rm = TRUE)
    simMatrix[i,11] <- as.numeric((quantile(IEBoot,0.025, na.rm = TRUE) <= approxIE ) & (quantile(IEBoot,0.975, na.rm = TRUE) >= approxIE ))
    simMatrix[i,12] <- as.numeric((quantile(DEBoot,0.025, na.rm = TRUE) <= approxDE ) & (quantile(DEBoot,0.975, na.rm = TRUE) >= approxDE ))
    
  }
  prev <- mean(prevVec)
  out <- rbind( out, c(prev, approxIE, approxDE, sd(simMatrix[,1], na.rm = TRUE), sd(simMatrix[,2], na.rm = TRUE),
                       colMeans(simMatrix[,c(1:12)], na.rm = TRUE)) )
}

write.table(out[3:dim(out)[1],], file=paste("contSimn", n, "_up", upBd, ".txt", sep=""),
            sep="\t", row.names=FALSE, quote = FALSE)

