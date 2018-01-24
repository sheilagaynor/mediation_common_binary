#Necessary functions, packages
library(boot);library(MASS);library(NORMT3)

source("BinSims_Func.R")

#Start loop across prevalence values
out <- matrix(NA, nrow=length(1:751), ncol=3)
for(i in 1:751) {
  set.seed(33)  
  const <-  seq(-2.75,1,0.005)[i]
    #Prevalence
    out[i,1] <- inv.logit(const + 0.4*0.4 + .5*inv.logit(0.1 + 0.5*0.4 + 0.4*0.3*0.4) + 0.25*(0.3*0.4)) 
    #Estimators
    out[i,2] <-  indirectEffect(t0=const,t1=0.4,t2=0.5,t4=0.25,
                                      c=.12,b0=0.1,b1=0.5,b2=0.4)
    out[i,3] <-    directEffect(t0=const,t1=0.4,t2=0.5,t4=0.25,
                                      c=.12,b0=0.1,b1=0.5,b2=0.4)
  }

write.table(out, file="binSimTruth.txt",
            sep="\t", row.names=FALSE, quote = FALSE)

