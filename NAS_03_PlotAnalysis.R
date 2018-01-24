NAS_Results <- read.delim("~/Desktop/Xihong/Final_Application/NAS_Results.txt", stringsAsFactors=FALSE)

#First p-value plot
library(ggplot2)
setwd("~/Desktop/Xihong/Final_Application/")
logPval <- -log(as.numeric(NAS_Results$PValue),10)
plotMat <- data.frame(cbind(NAS_Results$Site, (logPval), (1:62)))
names(plotMat) <- c("Site", "logP", "Ind")
plotMat$logP <- as.character(plotMat$logP)
plotMat$logP <- as.numeric(plotMat$logP)

plot1 <- ggplot(plotMat, aes(x=Site, y=logP)) + 
  geom_point(colour="#7570b3") +
  geom_hline(size=0.6,aes(yintercept=-log(0.05/62)),colour="dimgray") +
  labs(x="CpG Site", y="-Log P-Value") +
  theme(axis.text.x=element_text(angle = 90, hjust = 0, size=7),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
ggsave('~/Desktop/Xihong/Final_Application/Figure8.eps', plot=plot1,
       device="eps", width = 8, height = 3, units = 'in', dpi = 800)



#Second CI plot
NAS_Results$ID1 <- seq(1,124,2)
NAS_Results$ID2 <- seq(1.4,124.4,2)
plotData <- NAS_Results

plot2 <- ggplot(plotData, aes(x = ID1, y = DE_Est)) +
  geom_point(aes(x = plotData$ID2, y = plotData$IE_Est, color="#1b9e77")) +
  geom_point(aes(x = plotData$ID1, y = plotData$DE_Est, color="#d95f02")) +
  geom_errorbar(aes(ymax = plotData$IE_Up,ymin = plotData$IE_Low,x=ID2),
                width = .425, colour=("#d95f02")) +
  geom_errorbar(aes(ymax = plotData$DE_Up,ymin = plotData$DE_Low,x=ID1),
                width = .425, colour=("#1b9e77")) +
  theme(legend.position = "bottom") +
  scale_color_manual(name = "Effect", labels = c("NIE", "NDE"), 
                     values= c("#d95f02", "#1b9e77")) +
  scale_x_continuous(breaks=c(seq(1,124,2)) , labels=plotData$Site) +
  theme(axis.text.x=element_text(angle = 90, hjust = 0, size=7),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "bottom") +
  geom_hline(yintercept = 1,colour="dimgray") +
  xlab("CpG Site") + ylab("Effect CI") 


ggsave('~/Desktop/Xihong/Final_Application/Figure7.eps', plot=plot2,
       device="eps", width = 8, height = 4.5, units = 'in', dpi = 800)






