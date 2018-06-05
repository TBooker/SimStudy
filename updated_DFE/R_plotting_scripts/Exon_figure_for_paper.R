rm(list=ls())
library(ggplot2)
library(extrafont)
cbPalette <- c("#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")

exons<-read.csv('~/project/1.stats_from_exons/intervals/analysis/strategy_2/TEST')

exons$pi_red <- exons$pi/mean(exons[which(abs(exons$dist) > 75000 ),]$pi)

#plot(exons$pi_red~exons$dist, xlim=c(-100000,-10000),ylim=c(0.92,1.03))

bgs<-read.csv("~/project/2.simulations/updated_DFE/Exons/analysis/summary_BGS.Ne-anc.csv")
bgs<-bgs[order(bgs$mid),]
bgs$pi_red <- bgs$pi/mean(bgs[which(abs(bgs$mid) > 75000 ),]$pi)
bgs$pi_red_up <- bgs$pi_up/mean(bgs[which(abs(bgs$mid) > 75000 ),]$pi)
bgs$pi_red_low <- bgs$pi_low/mean(bgs[which(abs(bgs$mid) > 75000 ),]$pi)
bgs$lab<-"Deleterious Mutations"
bgs$red<-(bgs$pi/0.0083)
sqrt(sum((bgs$pi_red - exons$pi_red)^2)/length(exons$pi_red))
sqrt(sum((bgs$pi_red - exons$pi_red)^2))

sel<-read.csv("~/project/2.simulations/updated_DFE/Exons/analysis/summary_SWEEP.Ne-anc.csv")
sel<-sel[order(sel$mid),]
sel$pi_red <- sel$pi/mean(sel[which(abs(sel$mid) > 75000 ),]$pi)
sel$pi_red_up <- sel$pi_up/mean(sel[which(abs(sel$mid) > 75000 ),]$pi)
sel$pi_red_low <- sel$pi_low/mean(sel[which(abs(sel$mid) > 75000 ),]$pi)
sel$lab<-"Advantageous Mutations"
sel$red<-sel$pi/0.0083

sqrt(sum((sel$pi_red - exons$pi_red)^2)/length(exons$pi_red))
sqrt(sum((sel$pi_red - exons$pi_red)^2))

both<-read.csv("~/project/2.simulations/updated_DFE/Exons/analysis/summary_BGS+SWEEP.Ne-anc.csv")
both<-both[order(both$mid),]
both$pi_red <- both$pi/mean(both[which(abs(both$mid) > 60000 ),]$pi)
both$pi_red_up <- both$pi_up/mean(both[which(abs(both$mid) > 60000 ),]$pi)
both$pi_red_low <- both$pi_low/mean(both[which(abs(both$mid) > 60000 ),]$pi)
both$red<-1/(bgs$red^-1 + (1-(sel$red)))
both$lab<-"Both"



all1 <- rbind(bgs,sel,both)
all1$dfe <- 'uSFS Excl. Divergence'

bgs<-read.csv("~/project/2.simulations/updated_DFE/Exons/analysis/summary_BGS.full_usfs.csv")
bgs<-bgs[order(bgs$mid),]
bgs$pi_red <- bgs$pi/mean(bgs[which(abs(bgs$mid) > 75000 ),]$pi)
bgs$pi_red_up <- bgs$pi_up/mean(bgs[which(abs(bgs$mid) > 75000 ),]$pi)
bgs$pi_red_low <- bgs$pi_low/mean(bgs[which(abs(bgs$mid) > 75000 ),]$pi)
bgs$lab<-"Deleterious Mutations"
bgs$red<-(bgs$pi/0.0083)
sqrt(sum((bgs$pi_red - exons$pi_red)^2))


sel<-read.csv("~/project/2.simulations/updated_DFE/Exons/analysis/summary_SWEEP.full_usfs.csv")
sel<-sel[order(sel$mid),]
sel$pi_red <- sel$pi/mean(sel[which(abs(sel$mid) > 75000 ),]$pi)
sel$pi_red_up <- sel$pi_up/mean(sel[which(abs(sel$mid) > 75000 ),]$pi)
sel$pi_red_low <- sel$pi_low/mean(sel[which(abs(sel$mid) > 75000 ),]$pi)
sel$lab<-"Advantageous Mutations"
sel$red<-sel$pi/0.0083
sqrt(sum((sel$pi_red - exons$pi_red)^2))

both<-read.csv("~/project/2.simulations/updated_DFE/Exons/analysis/summary_BGS+SWEEP.full_usfs.csv")
both<-both[order(both$mid),]
both$pi_red <- both$pi/mean(both[which(abs(both$mid) > 75000 ),]$pi)
both$pi_red_up <- both$pi_up/mean(both[which(abs(both$mid) > 75000 ),]$pi)
both$pi_red_low <- both$pi_low/mean(both[which(abs(both$mid) > 75000 ),]$pi)
both$red<-1/(bgs$red^-1 + (1-(sel$red)))
both$lab<-"Both"

sqrt(sum((both$pi_red - exons$pi_red)^2))

all2 <- rbind(bgs,sel,both)
all2$dfe <- 'Full uSFS'



all <- rbind(all1,all2)
all$lab <- factor(all$lab, levels = c('Deleterious Mutations','Advantageous Mutations', 'Both'),
                  labels = c('Deleterious Mutations','Advantageous Mutations', 'Both'))


vnames <-list(
  'Deleterious Mutations' = 'Deleterious Mutations',
  'Advantageous Mutations' = 'Advantageous Mutations',
  'Both' ='Both',
  'Full uSFS' = 'Model A',
  'uSFS Excl. Divergence' = 'Model B'
  #  'Ne-anc' = bquote(N[e-anc])
)

vlabeller <- function(variable,value){
  return(vnames[value])
}


########################
## Now I have all the data in a usable format
########################
## Scaled pi
cairo_pdf('~/project/simulation_paper_stuff/plots/ExonPiScaled.pdf', height = 8, width = 10)

ggplot(data=all,aes(x=mid/1000,y=pi_red))+
  # geom_line(col='black',lty=3)+
  geom_ribbon(aes(ymin=pi_red_low, ymax=pi_red_up, fill=lab ) ,alpha=0.8)+
  geom_line(data=exons,aes(x=dist/1000,y=pi_red))+
  geom_line(aes(y=1.0),lty=2)+
  scale_y_continuous(limits=c(0.75,1.1))+
  scale_fill_manual('', values = cbPalette)+
  xlab("Distance from Exon (Kb)")+
  ylab(expression(pi / pi[Ref]))+
  facet_grid(dfe~. , labeller = vlabeller)+
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
## Unscaled pi
cairo_pdf('~/project/simulation_paper_stuff/plots/ExonPiUnscaled.pdf', height = 8, width = 10)

ggplot(data=all,aes(x=mid/1000,y=pi))+
  geom_ribbon(aes(ymin=pi_lower, ymax=pi_upper, fill=lab ) ,alpha=0.8)+
  geom_line(data=exons,aes(x=dist/1000,y=pi))+
  geom_line(aes(y=0.0083), lty = 2)+
  scale_fill_manual('', values = cbPalette)+
  xlab("Distance from Exon (Kb)")+
  ylab(expression(pi))+
  facet_grid(dfe~. , labeller = vlabeller)+
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

cairo_pdf('~/project/simulation_paper_stuff/plots/ExonPi2.pdf', height = 8, width = 10)

ggplot(data=all, aes(x=mid/1000,y=pi2))+
  # geom_line(col='black',lty=3)+
  geom_ribbon(aes(ymin=pi2_lower, ymax=pi2_upper, fill=lab ) ,alpha=0.8)+
  geom_line(data=exons, aes(x=dist/1000, y=pi2))+
  scale_x_continuous(limits = c(-50,50))+
  scale_fill_manual('', values = cbPalette)+
  xlab("Distance from Exon (Kb)")+
  ylab(expression(pi[ratio]))+
  facet_grid(dfe~. , labeller = vlabeller, scales = 'free')+
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


#cairo_pdf('~/project/simulation_paper_stuff/plots/ExonXsi.pdf', height = 8, width = 8)

ggplot(data=all, aes(x=mid/1000,y=tajima))+
  # geom_line(col='black',lty=3)+
  geom_ribbon(aes(ymin=tajima_lower, ymax=tajima_upper, fill=lab ) ,alpha=0.8)+
  geom_line(data=exons, aes(x=dist/1000, y=tajima))+
  scale_x_continuous(limits = c(-50,50))+
  scale_fill_manual('', values = cbPalette)+
  xlab("Distance from Exon (Kb)")+
  ylab("Tajima's D")+
  facet_grid(dfe~. , labeller = vlabeller, scales = 'free')+
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
#dev.off()


cairo_pdf('~/project/simulation_paper_stuff/plots/ExonXsi.pdf', height = 8, width = 10)
ggplot(data=all, aes(x=mid/1000,y=xsi))+
  # geom_line(col='black',lty=3)+
  geom_ribbon(aes(ymin=xsi_lower+0.02, ymax=xsi_upper+0.02, fill=lab ) ,alpha=0.8)+
  geom_line(data=exons, aes(x=dist/1000, y=xsi))+
  scale_x_continuous(limits = c(-50,50))+
  scale_y_continuous(limits = c(0.725,0.79))+
  scale_fill_manual('', values = cbPalette)+
  xlab("Distance from Exon (Kb)")+
  ylab(expression(Xi))+
  facet_grid(dfe~. , labeller = vlabeller)+
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

############################## Stop




## Plot unscaled Pi
ggplot(data=all,aes(x=mid/1000,y=pi))+
  # geom_line(col='black',lty=3)+
  geom_ribbon(aes(ymin=pi_lower, ymax=pi_upper, fill=lab) ,alpha=0.8)+
  geom_line(data=exons,aes(x=dist/1000,y=pi))+
  #  geom_line(aes(y=1.0),lty=2)+
  scale_y_continuous(limits=c(0.006,0.009))+
  scale_x_continuous(limits=c(-100,100))+
  xlab("Distance from Exon (Kb)")+
  ylab(expression(pi))+
  facet_grid(dfe~lab)+
  theme_bw()+
  theme(
    axis.title.y = element_text(size=16,angle=0)
  )




###########################################
#    Variable recombination rates
#
###########################################
vaR<-read.csv("~/project/2.simulations/all_elements/all/eddieRuns/analysis2/BGS+SWEEP.take2.stats.csv")
vaR$lab<-'Variable'
mean(vaR[which(abs(vaR$position) > 60000 ),]$pi)/0.0083

conR<-read.csv("~/project/2.simulations/all_elements/all/eddieRuns/analysis2/constant_rho.stats.csv")
conR$lab<-'Constant'
mean(conR[which(abs(conR$position) > 60000 ),]$pi)/0.0083

rho<-rbind(vaR,conR)

cbPalette <- c( "#0072B2","#CC79A7", "#D55E00")

ggplot(data=rho,aes(x=position/1000,y=pi))+
  geom_ribbon(aes(ymin=pi_low, ymax=pi_up, fill=lab),alpha=0.8)+
  scale_fill_manual('Recombination Rates',values=c("orangered4",'blue'))+
  #  scale_y_continuous(limits=c(0.8,1.1))+
  scale_x_continuous(limits=c(-100,100))+
  xlab("Distance from Exon (Kbp)")+
  ylab(expression(pi))+
  scale_fill_manual('Recombination Rates', values = cbPalette)+
  #  facet_grid(lab~.)+
  theme_bw()+
  theme(
    text=element_text(family="Trebuchet MS"),
    axis.title.x = element_text(size=14,angle=0),
    axis.title.y = element_text(size=20,angle=0),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    legend.text = element_text(size =13),
    strip.text.y = element_text(size = 15)
  )
