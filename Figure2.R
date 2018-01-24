library(ggplot2)

#1000 samples; set location of this simulation output as working directory
setwd("~/Cont1000/")
table_out <- list.files("~/Cont1000/")
figData1000 <- read.delim(table_out[1], stringsAsFactors=FALSE)
for(file in table_out[2:length(table_out)]){
  figData1000 <- rbind(figData1000, read.delim(file, stringsAsFactors=FALSE))
}
figData1000 <- data.frame(sapply(figData1000, Re))
names(figData1000) <- c("prev", "trueIE", "trueDE", "sdIE", "sdDE", "propIE", "propDE",
                    "deltaIE", "deltaDE", "covIE_delta", "covDE_delta", "bootIE_se",
                    "bootDE_se", "bootIE_mad", "bootDE_mad", "covIE_boot", "covDE_boot")
figData1000 <- figData1000[order(figData1000$prev),]
figData1000$IEPercBias <- 100*(figData1000$propIE-figData1000$trueIE)/figData1000$trueIE
figData1000$DEPercBias <- 100*(figData1000$propDE-figData1000$trueDE)/figData1000$trueDE
figData1000$IEBias <- figData1000$propIE - figData1000$trueIE
figData1000$DEBias <- figData1000$propDE - figData1000$trueDE

#500 samples
setwd("~/Cont500/")
table_out <- list.files("~/Cont500/")
figData500 <- read.delim(table_out[1], stringsAsFactors=FALSE)
for(file in table_out[2:length(table_out)]){
  figData500 <- rbind(figData500, read.delim(file, stringsAsFactors=FALSE))
}
figData500 <- data.frame(sapply(figData500, Re))
names(figData500) <- c("prev", "trueIE", "trueDE", "sdIE", "sdDE", "propIE", "propDE",
                        "deltaIE", "deltaDE", "covIE_delta", "covDE_delta", "bootIE_se",
                        "bootDE_se", "bootIE_mad", "bootDE_mad", "covIE_boot", "covDE_boot")
figData500 <- figData500[order(figData500$prev),]
figData500$IEPercBias <- 100*(figData500$propIE-figData500$trueIE)/figData500$trueIE
figData500$DEPercBias <- 100*(figData500$propDE-figData500$trueDE)/figData500$trueDE
figData500$IEBias <- figData500$propIE - figData500$trueIE
figData500$DEBias <- figData500$propDE - figData500$trueDE
figData1000 <- figData1000[which(figData1000$prev %in% figData500$prev),]


library(ggplot2);library(reshape);library(grid)
library(gridExtra);library(lattice)

plotData <- melt(data.frame(figData500$prev, figData500$IEBias, figData1000$IEBias), id.vars="figData500.prev")
names(plotData)[2] <- "Sample"
g1 <- ggplot(plotData, aes(x = figData500.prev, y = value, colour = (Sample))) + 
  geom_smooth(se = FALSE) +  
  labs(list(title = "A", x = "Prevalence", y = "Absolute bias")) +
  scale_color_manual(values=c("#46ACC8", "#B40F20"), 
                     labels=c("n=500", "n=1000")) +
  theme(panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = -0.25, face = 'bold')) +
  xlim(0.1,0.75) +  ylim(-.1,.1) 


plotData <- melt(data.frame(figData500$prev, figData500$DEBias, figData1000$DEBias), id.vars="figData500.prev")
names(plotData)[2] <- "Sample"
g3 <- ggplot(plotData, aes(x = figData500.prev, y = value, colour = (Sample))) + 
  geom_smooth(se = FALSE) +  
  labs(list(title = "B", x = "Prevalence", y = "Absolute bias")) +
  scale_color_manual(values=c("#46ACC8", "#B40F20"), 
                     labels=c("n=500", "n=1000")) +
  theme(panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = -0.25,face = 'bold')) +
  xlim(0.1,0.75) +  ylim(-.1,.1) 


plotData <- melt(data.frame(figData500$prev, figData500$IEPercBias, figData1000$IEPercBias), id.vars="figData500.prev")
names(plotData)[2] <- "Sample"
g2 <- ggplot(plotData, aes(x = figData500.prev, y = value, colour = (Sample))) + 
  geom_smooth(se = FALSE) +  
  labs(list(title = "C", x = "Prevalence", y = "Percent bias")) +
  scale_color_manual(values=c("#46ACC8", "#B40F20"), 
                     labels=c("n=500", "n=1000")) +
  theme(panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = -0.16,face = 'bold')) +
  xlim(0.1,0.75)  +  ylim(-2,4) 


plotData <- melt(data.frame(figData500$prev, figData500$DEPercBias, figData1000$DEPercBias), id.vars="figData500.prev")
names(plotData)[2] <- "Sample"
g4 <- ggplot(plotData, aes(x = figData500.prev, y = value, colour = (Sample))) + 
  geom_smooth(se = FALSE) +  
  labs(list(title = "D", x = "Prevalence", y = "Percent bias")) +
  scale_color_manual(values=c("#46ACC8", "#B40F20"), 
                     labels=c("n = 500", "n = 1000"),
                     name="") +
  theme(panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = -0.16,face = 'bold'),
        legend.position="bottom",
        legend.text=element_text(size=11)) +
  xlim(0.1,0.75) +  ylim(-2,4) 

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend<-g_legend(g4)


ggsave('~/Figure2.eps', plot=grid.arrange(arrangeGrob( g1 + theme(legend.position="none"),
                                                         g3 + theme(legend.position="none"),
                                                         g2 + theme(legend.position="none"),
                                                         g4 + theme(legend.position="none"),
                                                         nrow=2),mylegend, nrow=2, heights=c(20,2)),
       device="eps", width = 6, height = 6, units = 'in', dpi = 800)

