rm(list = ls())
library(ggplot2)
# Now plotting additional lots for supplementary section 

expo <- read.csv('project/2.simulations/exponential_DFE/melted_data.csv')


vnames <-list(
  'pi' = expression(pi),
  'pi2' = expression(pi[ratio])
)

vlabeller <- function(variable,value){
  return(vnames[value])
}



ggplot(data = expo, aes(x=mid/1000, y=lower))+
  geom_ribbon( aes( ymax = upper, ymin = lower), fill= "#CC79A7", alpha = 0.8)+
  facet_grid(stat~., scales = 'free',  labeller = vlabeller)+
  theme_bw()+
  xlab("Distance to Exon (Kb)")+
  ylab('')+
  theme(
    text=element_text(family="Trebuchet MS"),
    axis.title.x = element_text(size=14,angle=0),
    axis.title.y = element_text(size=14,angle=0),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
  strip.text.y = element_text(size=15, face = 'bold'),
        legend.text = element_text(size =13)
  )






##########################################################
# Plot the Kim and Stephan expectation
##########################################################
rm(list=ls())

library(ggplot2)

cbPalette <- c("#CC79A7", "#E69F00", "red")

#plot(exons$pi_red~exons$dist, xlim=c(-100000,-10000),ylim=c(0.92,1.03))

bgs<-read.csv("project/2.simulations/analysis_for_RMS/Ne-anc/summary_BGS.NeAnc.Exons.csv")
bgs<-bgs[order(bgs$mid),]
bgs$lab<-"Background Selection"
bgs$red<-(bgs$pi/0.0083)

sel<-read.csv("project/2.simulations/analysis_for_RMS/Ne-anc/summary_SWEEP.NeAnc.Exons.csv")
sel<-sel[order(sel$mid),]
sel$lab<-'Selective Sweeps'
sel$red<-sel$pi/0.0083

both<-read.csv("project/2.simulations/analysis_for_RMS/Ne-anc/summary_BGS+SWEEP.NeAnc.Exons.csv")
both<-both[order(both$mid),]
both$red<-1/(bgs$red^-1 + (1-(sel$red)))
both$lab<-'Kim and Stephan'

all1 <- rbind(bgs,sel,both)

all1$lab <- factor(all1$lab, levels = c('Background Selection','Selective Sweeps', 'Kim and Stephan'),
                  labels = c('Background Selection','Selective Sweeps', 'Kim and Stephan'))

both$lab2<-'null'

ggplot(data= all1, aes(x = mid/1000))+
  geom_ribbon(data = both, aes(ymax = pi_upper/0.0083, ymin = pi_lower/0.0083, group = lab2 , fill = 'test 1'), alpha = 0.8, lwd = 0.8)+
  geom_line( aes(y = red, group = lab, colour = lab), lwd = 0.75 )+
  scale_y_continuous(limits = c(0.65,1.0))+
  ylab(expression(pi/theta))+
  xlab('Distance to Exon (Kb)')+
  scale_colour_manual('', values = cbPalette)+
  scale_fill_manual('', values = c('#56B4E9'), guide = FALSE)+
  theme_bw()+
  theme(
    text=element_text(family="Trebuchet MS"),
    axis.title.x = element_text(size=14,angle=0),
    axis.title.y = element_text(size=17,angle=0),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    legend.text = element_text(size =13),
    legend.key = element_blank())


#####################################

rm(list=ls())

hall<- read.csv('project/2.simulations/halligan_2013_DFE/halligan2013.sfs.stats.csv')
hall$lab <- 'Halligan et al. 2013'
str(hall)
bgs <- read.csv('project/2.simulations/updated_DFE/Exons/analysis/summary_BGS.full_usfs.csv')
bgs$lab <- 'Current Study'
#comb <- rbind(hall, bgs)
str(bgs)
library(ggplot2)



pallete = c("#CC79A7",'#009E73')
cairo_pdf('~/project/simulation_paper_stuff/plots/redux/Halligan_dDFE.pdf', height = 4, width = 8)

ggplot(data = hall, aes(x= position/1000, y = pi))+
  geom_ribbon( aes( ymax = pi_up, ymin = pi_low, fill = lab), alpha = 0.7)+
  geom_ribbon(data = bgs, aes( x = mid/1000, ymax = pi_upper, ymin = pi_lower, fill = lab), alpha = 0.7)+
  scale_fill_manual('dDFE Estimate', values = pallete)+
  xlab('Distance from Exon (Kb)')+
  theme_bw()+
  ylab(expression(pi))+
  theme(
    text=element_text(family="Trebuchet MS"),
    axis.title.x = element_text(size=14,angle=0),
    axis.title.y = element_text(size=17,angle=0, vjust = 0.5),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    legend.text = element_text(size =13),
    legend.key = element_blank())

dev.off()
