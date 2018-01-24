###Mediation method
###Applied normative aging study
library(sas7bdat)

#Read in methylation, phenotype data
methPhen <- read.sas7bdat("/Users/sheilagaynor/Desktop/Xihong/Paper_NAS_Application/sg_04oct16.sas7bdat")
methPhen2 <- read.sas7bdat("/Users/sheilagaynor/Desktop/Xihong/Paper_NAS_Application/sg_03mar17.sas7bdat")
NAS450k_samplesmetadata_2014.06.12 <- read.csv("~/Desktop/Xihong/Paper_NAS_Application/NAS450k_samplesmetadata_2014-06-12.csv")

#Get the methylation IDs, remove any duplicates and take most recent methylation
idListMeth <- NAS450k_samplesmetadata_2014.06.12[,1:3]
idListMeth <- idListMeth[order(idListMeth$yr_vis, decreasing = TRUE),]
idListMeth <- idListMeth[which(duplicated(idListMeth$id,fromLast = FALSE)==FALSE),]

#Obtain years exposure/outcome, to ensure correct ordering
methPhen$obstr <- methPhen$fev1 / methPhen$fvc
idListExp <- methPhen[,c(1,133,132,214,3)] #smoking and years
idListOut <- methPhen[,c(1,41, 58, 39, 54, 62, 63, 60, 56, 214, 40, 215)] #potential outcomes

#Loop through to check which dates and observations should be taken 
fullData <- matrix(NA, nrow=675, ncol=16)
for (idIndex in 1:dim(idListMeth)[1]){
  fullData[idIndex, 1] <- as.character(idListMeth[idIndex,3]) #Save sample name (450k)
  fullData[idIndex, 2] <- idListMeth[idIndex,1] #Save id 
  fullData[idIndex, 3] <- idListMeth[idIndex,2] #Save methylation year
  #Obtain all outcomes, exposures under the same id
  #Obtain only outcomes after methyl, exposures before methyl
  outcomes <- idListOut[which(idListOut$id == idListMeth[idIndex,1]),]
  outcomes <- outcomes[which(outcomes$year >= idListMeth$yr_vis[idIndex]),]
  exposures <- idListExp[which(idListExp$id == idListMeth[idIndex,1]),]
  exposures <- exposures[which(exposures$year <= idListMeth$yr_vis[idIndex]),]
  #For those with exposure data, make sure it precedes the methylation mediator
  if ( dim(exposures)[1] > 0){
  fullData[idIndex, 4] <- exposures$smk[which.max(exposures$year)]
  fullData[idIndex, 5] <- exposures$yr_smk[which.max(exposures$year)]
  fullData[idIndex, 6] <- exposures$age[which.max(exposures$year)]
  fullData[idIndex, 7] <- max(outcomes$stroke, na.rm = TRUE)
  fullData[idIndex, 8] <- max(outcomes$emphys, na.rm = TRUE)
  fullData[idIndex, 9] <- max(outcomes$chd, na.rm = TRUE)
  fullData[idIndex, 10] <- max(outcomes$asthma, na.rm = TRUE)
  fullData[idIndex, 11] <- max(outcomes$cough, na.rm = TRUE)
  fullData[idIndex, 12] <- max(outcomes$coughcr, na.rm = TRUE)
  fullData[idIndex, 13] <- max(outcomes$hayfev, na.rm = TRUE)
  fullData[idIndex, 14] <- max(outcomes$bronch, na.rm = TRUE)
  outcomes <- outcomes[is.na(outcomes$obstr) == FALSE, ]
  if (dim(outcomes)[1] > 0){
  fullData[idIndex, 15] <- outcomes$obstr[which.max(outcomes$year)]
  }
  }}
#Create indicator from fev1/fvc ratio
fullData[,16] <- ifelse(fullData[,15] < .7, 1, 0)
methPhenotypes <- data.frame(fullData)
names(methPhenotypes) <- c("samplename_450k", "id", "yr_vis", "smk", "yr_smk", "age", "stroke",
                           "emphys", "chdEvent", "asthma", "cough", "coughcr", "hayfev", "bronch", "obstr", "obsInd" )

#Read in methylation data
load("/Users/sheilagaynor/Desktop/Xihong/Paper_NAS_Application/datdbbmiq_beta_2015-07-08.Rdata")
methSites <- c("cg00310412","cg01692968", "cg01731783", "cg01899089", "cg01901332", "cg01940273", "cg02451831",
               "cg02657160", "cg03329539", "cg03547355", "cg03636183", "cg03991871", "cg04885881", "cg05284742", "cg05575921",
               "cg05951221", "cg06060868", "cg06126421", "cg06644428", "cg07123182", "cg08709672", "cg09935388", "cg11207515",
               "cg11231349", "cg11314684", "cg11660018", "cg12075928", "cg12803068", "cg12806681", "cg12876356", "cg13193840",
               "cg13976502", "cg14580211", "cg14753356", "cg14817490", "cg15342087", "cg18316974", "cg19572487", "cg19859270",
               "cg20295214", "cg21121843", "cg21161138", "cg21566642", "cg21611682", "cg21913886", "cg22132788", "cg22851561",
               "cg23079012", "cg23161492", "cg23576855", "cg23771366", "cg23916896", "cg24090911", "cg24859433", "cg24996979",
               "cg25189904", "cg25648203", "cg25949550", "cg26271591", "cg26703534", "cg26963277", "cg27241845")
methMatrix <- mdata[which(row.names(mdata) %in% methSites),]
rm(mdata)
methFrame <- data.frame(t(methMatrix))
rm(methMatrix)

#Merge methylation and phenotype
NAS450k_samplesmetadata_2014.06.12 <- read.csv("~/Desktop/Xihong/Paper_NAS_Application/NAS450k_samplesmetadata_2014-06-12.csv")
methFrame$samplename_450k <- row.names(methFrame)
meth <- merge(NAS450k_samplesmetadata_2014.06.12,methFrame,by="samplename_450k")
meth2 <- merge(meth, methPhen2, by="samplename_450k")
methData <- merge(methPhenotypes,meth2,by="samplename_450k")
methData$smoke <- ifelse(methData$smk== 1, 0, 1)
names(methData)[92:96] <- c("C1","C2", "C3", "C4", "C5") 


write.table(methData, "~/Desktop/Xihong/Final_Application/NAS_Data.txt", sep="\t", row.names=FALSE, quote = FALSE)

