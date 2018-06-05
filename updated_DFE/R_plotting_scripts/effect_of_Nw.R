rm(list=ls())
library(ggplot2)
library(extrafont)
cbPalette <- c("#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")

exons<-read.csv('~/project/1.stats_from_exons/intervals/analysis/strategy_2/TEST')
exons$pi_red <- exons$pi/mean(exons[which(abs(exons$dist) > 75000 ),]$pi)
exons$lab<-'M. m. castaneus'
#plot(exons$pi_red~exons$dist, xlim=c(-100000,-10000),ylim=c(0.92,1.03))

bgs<-read.csv("~/project/2.simulations/updated_DFE/Exons/Ne-anc/BGS/BGS.Ne-anc.Exons.grand.csv")
bgs<-bgs[order(bgs$mid),]
bgs$pi_red <- bgs$pi/mean(bgs[which(abs(bgs$mid) > 75000 ),]$pi)
bgs$lab<-"Deleterious Mutations"
bgs$red<-(bgs$pi/0.0083)

sel<-read.csv("~/project/2.simulations/updated_DFE/Exons/Ne-anc/SWEEP/SWEEP.Ne-anc.Exons.grand.csv")
sel<-sel[order(sel$mid),]
sel$pi_red <- sel$pi/mean(sel[which(abs(sel$mid) > 75000 ),]$pi)
sel$lab<-"Advantageous Mutations"
sel$red<-(sel$pi/0.0083)


both<-read.csv("~/project/2.simulations/updated_DFE/Exons/Ne-anc/BGS+SWEEP/BGS+SWEEP.grand.csv")
both<-both[order(both$mid),]
both$pi_red <- both$pi/mean(both[which(abs(both$mid) > 75000 ),]$pi)
both$lab<-"Both"
both$red<-(both$pi/0.0083)

all1 <- rbind(bgs,sel,both)

# 
# all <- rbind(all1,all2)
# all$lab <- factor(all$lab, levels = c('Deleterious Mutations','Advantageous Mutations', 'Both'),
#                   labels = c('Deleterious Mutations','Advantageous Mutations', 'Both'))
# 
# 
# vnames <-list(
#   'Deleterious Mutations' = 'Deleterious Mutations',
#   'Advantageous Mutations' = 'Advantageous Mutations',
#   'Both' ='Both',
#   'Full uSFS' = 'Model A',
#   'uSFS Excl. Divergence' = 'Model B'
#   #  'Ne-anc' = bquote(N[e-anc])
# )
# 
# vlabeller <- function(variable,value){
#   return(vnames[value])
# }
# 

########################
## Now I have all the data in a usable format
########################
cairo_pdf('~/project/2.simulations/updated_DFE/Nw_comparison.pdf', height =6 , width =12)

ggplot(data=all1,aes(x=mid/1000,y=pi_red, col = lab))+
  geom_line()+
  geom_line(data=exons,aes(x=dist/1000,y=pi_red), col='black')+
  scale_y_continuous(limits=c(0.75,1.1))+
  scale_colour_manual('', values = cbPalette)+
  xlab("Distance from Exon (Kb)")+
  ylab(expression(pi / pi[Ref]))+
#  facet_grid(dfe~.)+
  #  facet_grid(dfe~. , labeller = vlabeller)+
  theme_bw()+
  theme(
    text=element_text(family="Trebuchet MS"),
    axis.title.x = element_text(size=14,angle=0),
    axis.title.y = element_text(size=20,angle=0,vjust=0.5),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    legend.text = element_text(size =13),
    strip.text.y = element_text(size = 15)
  )

dev.off()

rm(list=ls())
bgs<-read.csv("~/project/2.simulations/updated_DFE/Exons/Ne-anc/BGS/BGS.Ne-anc.Exons.grand.csv")
bgs$red<-(bgs$pi/0.0083)
bgs$lab <- 'BGS'
sel<-read.csv("~/project/2.simulations/updated_DFE/Exons/Ne-anc/SWEEP/SWEEP.Ne-anc.Exons.grand.csv")
sel$red<-(sel$pi/0.0083)
sel$lab <- 'SSW'
both<-read.csv("~/project/2.simulations/updated_DFE/Exons/Ne-anc/BGS+SWEEP/BGS+SWEEP.grand.csv")
both$red<-(both$pi/0.0083)
both$lab <- 'BOTH'
all<-rbind(bgs,sel,both)
ggplot(data=all, aes(x=mid,y=red,col=lab))+
  geom_line()
