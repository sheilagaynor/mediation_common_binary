methData <- read.delim("~/Desktop/Xihong/Final_Application/NAS_Data.txt", stringsAsFactors=FALSE)

library(boot);library(MASS);library(NORMT3)
indirectEffectProposed <-function(s,t0,t1,t2,t4,t5,t6,t7,t8,t9,c1,c2,c3,c4,c5,c6,sigma,b0,b1,b2,b3,b4,b5,b6,b7){ 
  return( oddCalc(pnorm( (s*t0 + s*t1 + s*t4*c1 + s*t5*c2 + s*t6*c3 + s*t7*c4 + s*t8*c5 + s*t9*c6 + s*t2*(b0 + b1 + b2*c1 + b3*c2 + b4*c3 + b5*c4 + b6*c5 + b7*c6)) / (sqrt(1 + (s*t2*sigma)^2)))) / 
          oddCalc(pnorm( (s*t0 + s*t1 + s*t4*c1 + s*t5*c2 + s*t6*c3 + s*t7*c4 + s*t8*c5 + s*t9*c6 + s*t2*(b0 +      b2*c1 + b3*c2 + b4*c3 + b5*c4 + b6*c5 + b7*c6)) / (sqrt(1 + (s*t2*sigma)^2))  ))) }
directEffectProposed   <-function(s,t0,t1,t2,t4,t5,t6,t7,t8,t9,c1,c2,c3,c4,c5,c6,sigma,b0,b1,b2,b3,b4,b5,b6,b7){ 
  return( oddCalc(pnorm( (s*t0 + s*t1 + s*t4*c1 + s*t5*c2 + s*t6*c3 + s*t7*c4 + s*t8*c5 + s*t9*c6 + s*t2*(b0 + b1 + b2*c1 + b3*c2 + b4*c3 + b5*c4 + b6*c5 + b7*c6)) / (sqrt(1 + (s*t2*sigma)^2)))) / 
          oddCalc(pnorm( (s*t0 +        s*t4*c1 + s*t5*c2 + s*t6*c3 + s*t7*c4 + s*t8*c5 + s*t9*c6 + s*t2*(b0 +      b2*c1 + b3*c2 + b4*c3 + b5*c4 + b6*c5 + b7*c6)) / (sqrt(1 + (s*t2*sigma)^2))  ))) }
DEBoot <- rep(NA, 1000); IEBoot <- rep(NA, 1000) 
ageMedian <- median(as.numeric(as.character(methData$age)), na.rm=T); C1Median <- median(as.numeric(as.character(methData$C1)), na.rm=T)
C2Median <- median(as.numeric(as.character(methData$C2)), na.rm=T); C3Median <- median(as.numeric(as.character(methData$C3)), na.rm=T)
C4Median <- median(as.numeric(as.character(methData$C4)), na.rm=T); C5Median <- median(as.numeric(as.character(methData$C5)), na.rm=T)
mse <- function(sm){ mean(sm$residuals^2) }; erfinv <-  function(x) { qnorm((1 + x) /2) / sqrt(2) }; oddCalc<-function(p){ return((p/(1-p))) }

#For Proposed method
effectMatrix <- matrix(NA, nrow=62, ncol=8)
rowCount <- 0
#Obtain results 
for (i in 26:87) {
  MEDIATOR <- names(methData)[i]
  rowCount <- 1 + rowCount
  #####################
  ###Point Estimates###
  #####################
  temp <- data.frame(methData[1:675, c(98, which(names(methData)==MEDIATOR), 16, 6, 92:96)])
  nameExp <- names(methData)[i]
  effectMatrix[rowCount,8] <- nameExp
  #Fit two regression models
  logitModel <- (glm((obsInd) ~ as.numeric(smoke) + as.numeric(get(MEDIATOR)) + 
                       as.numeric(age) + as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + as.numeric(C5), 
                       na.action=na.omit, family="binomial", data=temp))
  probitModel <- (glm((obsInd) ~ as.numeric(smoke) + as.numeric(get(MEDIATOR)) + 
                       as.numeric(age) + as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + as.numeric(C5), 
                     na.action=na.omit, family=binomial(link="probit"), data=temp))
  sEst <- median(coef(probitModel)/coef(logitModel), na.rm = TRUE)
  logitModelCoef <- coef(logitModel)
  linearModel <- lm( as.numeric(get(MEDIATOR)) ~ as.numeric(smoke) + as.numeric(age) +
                       as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + as.numeric(C5) , data=temp)
  logitModelCoef <- coef(logitModel); linearModelCoef <- coef(linearModel)
  #Obtain effect estimates
  effectMatrix[rowCount,1] <-  indirectEffectProposed(s=sEst,t0=logitModelCoef[1],t1=logitModelCoef[2],t2=logitModelCoef[3],t4=logitModelCoef[4],
                                            t5=logitModelCoef[5],t6=logitModelCoef[6],t7=logitModelCoef[7],t8=logitModelCoef[8],t9=logitModelCoef[9],
                                            c1=ageMedian,c2=C1Median,c3=C2Median,c4=C3Median,c5=C4Median,c6=C5Median,
                                            sigma=sqrt(mse(linearModel)),b0=linearModelCoef[1],b1=linearModelCoef[2],b2=linearModelCoef[3],
                                            b3=linearModelCoef[4],b4=linearModelCoef[5],b5=linearModelCoef[6],b6=linearModelCoef[7],b7=linearModelCoef[8])
  effectMatrix[rowCount,2] <-    directEffectProposed(s=sEst,t0=logitModelCoef[1],t1=logitModelCoef[2],t2=logitModelCoef[3],t4=logitModelCoef[4],
                                            t5=logitModelCoef[5],t6=logitModelCoef[6],t7=logitModelCoef[7],t8=logitModelCoef[8],t9=logitModelCoef[9],
                                            c1=ageMedian,c2=C1Median,c3=C2Median,c4=C3Median,c5=C4Median,c6=C5Median,
                                            sigma=sqrt(mse(linearModel)),b0=linearModelCoef[1],b1=linearModelCoef[2],b2=linearModelCoef[3],
                                            b3=linearModelCoef[4],b4=linearModelCoef[5],b5=linearModelCoef[6],b6=linearModelCoef[7],b7=linearModelCoef[8])
  #Get testing p-value
  effectMatrix[rowCount,3] <- max(summary(logitModel)$coefficients[3,4], summary(linearModel)$coefficients[2,4])

  
  #####################
  ##Bootstrapped MAD###
  #####################
  for (j in 1:500){
    tempID <- sample(1:675, replace = TRUE)
    tempBoot <- data.frame(methData[tempID, c(98, which(names(methData)==MEDIATOR), 16, 6, 92:96)])
    names(tempBoot) <- names(temp)
    #Get bootstrap sample expectations
    ageMedian2 <- median(as.numeric(as.character(methData$age[tempID])), na.rm=T); C1Median2 <- median(as.numeric(as.character(methData$C1[tempID])), na.rm=T)
    C2Median2 <- median(as.numeric(as.character(methData$C2[tempID])), na.rm=T); C3Median2 <- median(as.numeric(as.character(methData$C3[tempID])), na.rm=T)
    C4Median2 <- median(as.numeric(as.character(methData$C4[tempID])), na.rm=T); C5Median2 <- median(as.numeric(as.character(methData$C5[tempID])), na.rm=T)
    #Fit regression models
    logitModel2 <- (glm((obsInd) ~ as.numeric(smoke) + as.numeric(get(MEDIATOR)) +
                          as.numeric(age) + as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + as.numeric(C5),
                        na.action=na.omit, family="binomial", data=tempBoot))
    probitModel2 <- (glm((obsInd) ~ as.numeric(smoke) + as.numeric(get(MEDIATOR)) +
                          as.numeric(age) + as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + as.numeric(C5),
                        na.action=na.omit, family=binomial(link="probit"), data=tempBoot))
    sEst2 <- median(coef(probitModel2)/coef(logitModel2), na.rm = TRUE)
    linearModel2 <- lm( as.numeric(get(MEDIATOR)) ~ as.numeric(smoke) + as.numeric(age) +
                          as.numeric(C1) + as.numeric(C2) + as.numeric(C3) + as.numeric(C4) + as.numeric(C5) , data=tempBoot)
    logitModelCoef2 <- coef(logitModel2); linearModelCoef2 <- coef(linearModel2)
    #Store regression output for use in estimator
    IEBoot[j] <-  indirectEffectProposed(s=sEst2,t0=logitModelCoef2[1],t1=logitModelCoef2[2],t2=logitModelCoef2[3],t4=logitModelCoef2[4],
                                              t5=logitModelCoef2[5],t6=logitModelCoef2[6],t7=logitModelCoef2[7],t8=logitModelCoef2[8],t9=logitModelCoef2[9],
                                              c1=ageMedian2,c2=C1Median2,c3=C2Median2,c4=C3Median2,c5=C4Median2,c6=C5Median2,
                                              sigma=sqrt(mse(linearModel)),b0=linearModelCoef2[1],b1=linearModelCoef2[2],b2=linearModelCoef2[3],
                                              b3=linearModelCoef2[4],b4=linearModelCoef2[5],b5=linearModelCoef2[6],b6=linearModelCoef2[7],b7=linearModelCoef2[8])
    DEBoot[j] <-    directEffectProposed(s=sEst2,t0=logitModelCoef2[1],t1=logitModelCoef2[2],t2=logitModelCoef2[3],t4=logitModelCoef2[4],
                                              t5=logitModelCoef2[5],t6=logitModelCoef2[6],t7=logitModelCoef2[7],t8=logitModelCoef2[8],t9=logitModelCoef2[9],
                                              c1=ageMedian2,c2=C1Median2,c3=C2Median2,c4=C3Median2,c5=C4Median2,c6=C5Median2,
                                              sigma=sqrt(mse(linearModel)),b0=linearModelCoef2[1],b1=linearModelCoef2[2],b2=linearModelCoef2[3],
                                              b3=linearModelCoef2[4],b4=linearModelCoef2[5],b5=linearModelCoef2[6],b6=linearModelCoef2[7],b7=linearModelCoef2[8])
 }
  effectMatrix[rowCount,4] <- quantile(IEBoot,0.025, na.rm = TRUE);  effectMatrix[rowCount,5] <- quantile(IEBoot,0.975, na.rm = TRUE)
  effectMatrix[rowCount,6] <- quantile(DEBoot,0.025, na.rm = TRUE);  effectMatrix[rowCount,7] <- quantile(DEBoot,0.975, na.rm = TRUE)
}
effects <- data.frame(effectMatrix)
names(effects) <- c("IE_Est", "DE_Est", "PValue", "IE_Low", "IE_Up", "DE_Low", "DE_Up", "Site")

write.table(effects, "~/Desktop/Xihong/Final_Application/NAS_Results.txt", sep="\t", row.names=FALSE, quote = FALSE)
