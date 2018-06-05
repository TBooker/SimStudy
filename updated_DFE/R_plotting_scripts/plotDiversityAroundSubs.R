rm(list=ls())
library(ggplot2)
cbPalette <- c("#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")

dan<-read.csv('project/2.simulations/diversity_around_substitutions/DansData/analysedData.csv')
dan <- dan[ which(dan$stat == 'ratio'),]
dan <- dan[ which(dan$Substituted == 'TRUE'),]
dan$label<-dan$type
dan$label<- factor(dan$label, levels = c('0-fold','4-fold'),
                   labels = c('M. m. castaneus: 0-fold','M. m. castaneus: 4-fold'))

simD_neAnc <- read.csv('project/2.simulations/updated_DFE/LongRuns/substitutionPlots/Ne-anc.sfs.subs.analysed')
simD_neAnc$dfe<- 'Model B'
simD_fulluSFS <- read.csv('project/2.simulations/updated_DFE/LongRuns/substitutionPlots/full_usfs.sfs.analysed')
simD_fulluSFS$dfe<- 'Model A'
com <- rbind( simD_fulluSFS, simD_neAnc )


com$label<- factor(com$label, levels = c('Non-Synonymous','Synonymous'),
                  labels = c('Simulation: Non-Synonymous','Simulation: Synonymous'))
# Saved the plot as 11.96 (Width) x 4.28 (Height)
cairo_pdf('~/project/simulation_paper_stuff/plots/SubstitutionPlot.pdf' , width = 8.47  , height = 6.15) 
ggplot(data = com, aes(x = pos/1000, y = actual , fill = label, col = label))+
  geom_line()+
  geom_ribbon(aes(ymax=upper, ymin = lower), alpha = 0.8)+
  facet_grid(dfe~.)+
  ylab(expression(pi))+
  xlab('Distance from Substituted Site (Kbp)')+
  geom_line(data = dan, aes(x = mid, y= value*0.15), lwd = 0.8)+
  scale_fill_manual('Class of Substitution',values= cbPalette)+
  scale_color_manual('Class of Substitution',values= cbPalette)+
  theme_bw()+
  theme(
    text=element_text(family="Trebuchet MS"),
    axis.title.x = element_text(size=14,angle=0),
    axis.title.y = element_text(size=20,angle=0,vjust=0.5),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    legend.text = element_text(size =13, face = 'italic'),
    legend.title = element_text(size =15, face = 'bold'),
    strip.text.y = element_text(size = 15)
  )
dev.off()
# Saved the plot as 11.96 (Width) x 4.28 (Height) 

#simD_fulluSFS <- read.csv('project/2.simulations/diversity_around_substitutions/Kfam_0.01/divCorrect/N1000_x4000_divCorrect.sites.csv')
#plot(simD_fulluSFS$pos/1000 , log(simD_fulluSFS$sites/1e6), col = simD_fulluSFS$label)
#points(dan$sites/1e6~dan$mid, col = 'red')
#abline(h = 10)
