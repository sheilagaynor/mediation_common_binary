library(ggplot2)
#1000 samples; set working directory to where this output is
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


###########
#Error estimate plots
plotData <- melt(data.frame(figData500$prev, figData500$deltaIE, figData500$bootIE_se,
                            figData500$bootIE_mad, figData500$sdIE), id.vars="figData500.prev")
names(plotData)[2] <- "Sample"
g1 <- ggplot(plotData, aes(x = figData500.prev, y = value, colour = (Sample))) + 
  geom_smooth(se = FALSE) +  
  labs(list(title = "A", x = "Prevalence", y = "NIE OR error")) +
  scale_color_manual(values=c("#d7191c","#fdae61","#abdda4","#2b83ba")) +
  xlim(0.1,0.75) +  ylim(0.035,0.125) +
  theme(panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = -0.24, face = 'bold')) 

plotData <- melt(data.frame(figData500$prev, figData500$deltaDE, figData500$bootDE_se,
                            figData500$bootDE_mad, figData500$sdDE), id.vars="figData500.prev")
names(plotData)[2] <- "Sample"
g3 <- ggplot(plotData, aes(x = figData500.prev, y = value, colour = Sample)) + 
  geom_smooth(se = FALSE) +  
  labs(list(title = "C", x = "Prevalence", y = "NDE OR error", 
            colour = "")) +
  scale_color_manual(values=c("#d7191c","#fdae61","#abdda4","#2b83ba"), 
                     labels=c("Delta Method", "Bootstrap SD", "Bootstrap MAD", "Empirical")) +
  xlim(0.1,0.75) +  ylim(0.125,0.275) + 
  theme(panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = -0.215, face = 'bold'),
        legend.position="bottom")

plotData <- melt(data.frame(figData1000$prev, figData1000$deltaIE, figData1000$bootIE_se,
                            figData1000$bootIE_mad, figData1000$sdIE), id.vars="figData1000.prev")
names(plotData)[2] <- "Sample"
g2 <- ggplot(plotData, aes(x = figData1000.prev, y = value, colour = (Sample))) + 
  geom_smooth(se = FALSE) +  
  labs(list(title = "B", x = "Prevalence", y = "NIE OR error")) +
  scale_color_manual(values=c("#d7191c","#fdae61","#abdda4","#2b83ba")) +
  xlim(0.1,0.75) +  ylim(0.035,0.125) +
  theme(panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = -.255, face = 'bold')) 

plotData <- melt(data.frame(figData1000$prev, figData1000$deltaDE, figData1000$bootDE_se,
                            figData1000$bootDE_mad, figData1000$sdDE), id.vars="figData1000.prev")
names(plotData)[2] <- "Sample"
g4 <- ggplot(plotData, aes(x = figData1000.prev, y = value, colour = Sample)) + 
  geom_smooth(se = FALSE) +  
  labs(list(title = "D", x = "Prevalence", y = "NDE OR error", 
            colour = "")) +
  scale_color_manual(values=c("#d7191c","#fdae61","#abdda4","#2b83ba"), 
                     labels=c("Delta Method", "Bootstrap SD", "Bootstrap MAD", "Empirical")) +
  xlim(0.1,0.75) +  ylim(0.125,0.275) + 
  theme(panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="bottom",
        legend.text=element_text(size=11),
        plot.title = element_text(hjust = -.235, face = 'bold'))


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend<-g_legend(g4)


ggsave('~/Figure3.eps', plot=grid.arrange(arrangeGrob( g1 + theme(legend.position="none"),
                                    g2 + theme(legend.position="none"),
                                    g3 + theme(legend.position="none"),
                                    g4 + theme(legend.position="none"),
                                    nrow=2),mylegend, nrow=2, heights=c(20,2)),
       device="eps", width = 6, height = 6, units = 'in', dpi = 800)


