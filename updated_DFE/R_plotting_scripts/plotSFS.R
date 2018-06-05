rm(list=ls())
### Plotting scripts to get uSFS plots for simulation paper

sfs <- na.omit(read.csv('~/project/simulation_paper_stuff/site_frequency_spectrum_data.csv'))
str(sfs)
library(ggplot2)
sfs$Class.1 <- factor( sfs$Class.1 , levels = c('Neutral','Uncorrected', 'Corrected'),
                  labels = c('Neutral','Uncorrected', 'Corrected'))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cairo_pdf('TIFFplots/uSFS_comparison.pdf', height = 8, width = 9)

ggplot(data = sfs, aes( x = Alleles, y = Proportion, fill = Class.1, col= Class.1))+
  geom_bar( stat = 'identity', position = 'dodge', alpha = 0.7, width = 0.7)+
    facet_grid( Group~. )+
    scale_y_sqrt('Proportion of Sites')+
#      scale_y_log10()+
      scale_fill_manual('', values = cbPalette)+
      scale_color_manual(values = cbPalette, guide=FALSE)+
      scale_x_continuous(breaks= seq(1,19))+
      xlab('Number of Derived Allele Copies')+
      theme_bw()+
      theme(
        text=element_text(family="Trebuchet MS"),
        axis.title.x = element_text(size=16,angle=0),
        axis.title.y = element_text(size=17,angle=90,vjust = 0.5),
        axis.text.x = element_text(size=14,angle=0),
        axis.text.y = element_text(size=14,angle=0),
        strip.text.y = element_text(size = 17, face = 'bold'),
        legend.text = element_text(size =13),
        legend.key = element_blank())
dev.off()


sfs2 <- sfs[sfs$Group == '0-fold',]

sfs2 <- sfs2[sfs2$Class != '0-fold_corrected',]

sfs2$Class.1 <- factor( sfs2$Class.1 , levels = c('Neutral','Uncorrected'),
                       labels = c('4-fold','0-fold'))
ggplot(data = sfs2, aes( x = Alleles, y = Proportion, fill = Class.1, col= Class.1))+
  geom_bar( stat = 'identity', position = 'dodge', alpha = 0.7, width = 0.7)+
  #    scale_y_sqrt()+
  #      scale_y_log10()+
  scale_fill_manual('', values = cbPalette)+
  scale_color_manual(values = cbPalette, guide=FALSE)+
  scale_x_continuous(breaks= seq(1,19))+
  xlab('Number of Derived Allele Copies')+
  theme_bw()+
  theme(
#    text=element_text(family="Trebuchet MS"),
    axis.title.x = element_text(size=14,angle=0),
    axis.title.y = element_text(size=17,angle=90,vjust = 0.5),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    strip.text.y = element_text(size = 17, face = 'bold'),
    legend.text = element_text(size =13),
    legend.key = element_blank())


rm(list= ls())

sfs <- read.csv('~/project/simulation_paper_stuff/neutral_fit_v_observed.csv')
cbPalette <- c("white","maroon")

library(reshape2)
sfs1<-melt(sfs, id = c('Class', 'Alleles'))
sfs1<- sfs1[sfs1$Alleles>0,]
sfs1<- sfs1[sfs1$Alleles<20,]

sfs1$variable <- factor( sfs1$variable , levels = c('Expected','Observed.SFS'),
                       labels = c('Expected','Observed'))

cairo_pdf('TIFFplots/Neutral_uSFS_fit.pdf', height = 6, width = 9)

ggplot(data = sfs1, aes(x= Alleles, y= value, fill = variable))+
  geom_bar( stat = 'identity', position = 'dodge', width = 0.7, col = 'black')+
  facet_grid(Class~.)+
  scale_y_sqrt('Proportion of sites',breaks = c(0.001,0.003,0.006, 0.009))+
  scale_fill_manual('uSFS', values = cbPalette)+
  scale_x_continuous(breaks = seq(1,19))+
  xlab('Number of Derived Allele Copies')+
  theme_bw()+
  theme(
    text=element_text(family="Trebuchet MS"),
    axis.title.x = element_text(size=14,angle=0),
    axis.title.y = element_text(size=17,angle=90, vjust = 0.5),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    strip.text.y = element_text(size = 17, face = 'bold'),
    legend.text = element_text(size =13),
    legend.title = element_text(size =15),
    legend.key = element_blank())

dev.off()

rm(list=ls())
### Plotting scripts to get uSFS plots for simulation paper

sfs <- na.omit(read.csv('~/project/simulation_paper_stuff/DFE_fits_NewNw.csv'))
sfs1<-melt(sfs, id = c('Class', 'Method', 'Alleles'))

sfs1<- sfs1[sfs1$Alleles<20,]


sfs1<- sfs1[sfs1$Alleles>0,]

sfs1$Class <- factor( sfs1$Class , levels = c('0-fold','UTR', 'CNE'),
                       labels = c('0-fold','UTR', 'CNE'))
sfs1$Method <- factor( sfs1$Method , levels = c( 'Observed', 'Full uSFS','Ne-anc'),
                       labels = c( 'Observed', 'Expected: Model A','Expected: Model B'))

cbPalette <- c( "#1b9e77", '#d95f02', '#7570b3')

cairo_pdf('TIFFplots/selected_fit.pdf', height = 8, width = 9)

ggplot(data = sfs1, aes( x = Alleles, y = value, fill = Method))+
  geom_bar( stat = 'identity', position = 'dodge', alpha = 0.6, width = 0.9,col='black')+
  facet_grid( Class~. )+
  scale_y_sqrt('Proportion of Sites' , breaks = c(0.001, 0.003, 0.006))+
  #      scale_y_log10()+
  scale_fill_manual('', values = cbPalette)+
  scale_color_manual(values = cbPalette, guide=FALSE)+
  scale_x_continuous(breaks= seq(1,19))+
  xlab('Number of Derived Allele Copies')+
  theme_bw()+
  theme(
    text=element_text(family="Trebuchet MS"),
    axis.title.x = element_text(size=14,angle=0),
    axis.title.y = element_text(size=17,angle=90,vjust = 0.5),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    strip.text.y = element_text(size = 17, face = 'bold'),
    legend.text = element_text(size =13),
    legend.key = element_blank())

dev.off()
