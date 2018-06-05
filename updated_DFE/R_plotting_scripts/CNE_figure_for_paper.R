rm(list=ls())
library(ggplot2)

cbPalette <- c("#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")

cne<-read.csv('~/project/1.stats_from_exons/intervals/analysis/strategy_2/cne/subtracted_beds/CNE_nCpG_unfolded.stats.4.csv')
cne$pi_k <- cne$pi/cne$t2
cne$pi_red <- cne$pi_k/mean(cne[which(abs(cne$dist) > 4000 ),]$pi_k)
cne$pi_scale <- cne$pi_red*0.0083

bgs<-read.csv("~/project/2.simulations/updated_DFE/CNEs/analysis/summary_BGS.Ne-anc.csv")
bgs$sites <- NA
bgs<-bgs[order(bgs$mid),]
bgs<- bgs[which(abs(bgs$mid) < max(cne$dist)+1), ]
bgs$pi_red <- bgs$pi/mean(bgs[which(abs(bgs$mid) > 4000 ),]$pi)
bgs$pi_red_up <- bgs$pi_up/mean(bgs[which(abs(bgs$mid) > 4000 ),]$pi)
bgs$pi_red_low <- bgs$pi_low/mean(bgs[which(abs(bgs$mid) > 4000 ),]$pi)
bgs$lab<-'Deleterious Mutations'
sqrt(sum((bgs$pi_red - cne$pi_red)^2)/length(cne$pi_red))

sel<-read.csv("~/project/2.simulations/updated_DFE/CNEs/analysis/summary_SWEEP.Ne-anc.csv")
sel<-sel[order(sel$mid),]
sel<- sel[which(abs(sel$mid) < max(cne$dist)+1), ]
sel$pi_red <- sel$pi/mean(sel[which(abs(sel$mid) > 4000 ),]$pi)
sel$pi_red_up <- sel$pi_up/mean(sel[which(abs(sel$mid) > 4000 ),]$pi)
sel$pi_red_low <- sel$pi_low/mean(sel[which(abs(sel$mid) > 4000 ),]$pi)
sel$lab<-'Advantageous Mutations'
sqrt(sum((sel$pi_red - cne$pi_red)^2)/length(cne$pi_red))
str(sel)

both<-read.csv("~/project/2.simulations/updated_DFE/CNEs/analysis/summary_BGS+SWEEP.Ne-anc.csv")
both$sites <- NA
both<-both[order(both$mid),]
both<- both[which(abs(both$mid) < max(cne$dist)+1), ]
both$pi_red <- both$pi/mean(both[which(abs(both$mid) > 4000 ),]$pi)
both$pi_red_up <- both$pi_up/mean(both[which(abs(both$mid) > 4000 ),]$pi)
both$pi_red_low <- both$pi_low/mean(both[which(abs(both$mid) > 4000 ),]$pi)
both$lab<-'Both'
sqrt(sum((both$pi_red - cne$pi_red)^2)/length(cne$pi_red))


all1<-rbind(bgs,sel,both)
all1$dfe <- 'Ne-anc'

bgs<-read.csv("~/project/2.simulations/updated_DFE/CNEs/analysis/summary_BGS.full_usfs.csv")
bgs<-bgs[order(bgs$mid),]
bgs<- bgs[which(abs(bgs$mid) < max(cne$dist)+1), ]
bgs$pi_red <- bgs$pi/mean(bgs[which(abs(bgs$mid) > 4000 ),]$pi)
bgs$pi_red_up <- bgs$pi_up/mean(bgs[which(abs(bgs$mid) > 4000 ),]$pi)
bgs$pi_red_low <- bgs$pi_low/mean(bgs[which(abs(bgs$mid) > 4000 ),]$pi)
bgs$lab<-'Deleterious Mutations'
sqrt(sum((bgs$pi_red - cne$pi_red)^2)/length(cne$pi_red))

sel<-read.csv("~/project/2.simulations/updated_DFE/CNEs/analysis/summary_SWEEP.full_usfs.csv")
sel<-sel[order(sel$mid),]
sel<- sel[which(abs(sel$mid) < max(cne$dist)+1), ]
sel$pi_red <- sel$pi/mean(sel[which(abs(sel$mid) > 4000 ),]$pi)
sel$pi_red_up <- sel$pi_up/mean(sel[which(abs(sel$mid) > 4000 ),]$pi)
sel$pi_red_low <- sel$pi_low/mean(sel[which(abs(sel$mid) > 4000 ),]$pi)
sel$lab<-'Advantageous Mutations'
sqrt(sum((sel$pi_red - cne$pi_red)^2)/length(cne$pi_red))

both<-read.csv("~/project/2.simulations/updated_DFE/CNEs/analysis/summary_BGS+SWEEP.full_usfs.csv")
both<-both[order(both$mid),]
both<- both[which(abs(both$mid) < max(cne$dist)+1), ]
both$pi_red <- both$pi/mean(both[which(abs(both$mid) > 4000 ),]$pi)
both$pi_red_up <- both$pi_up/mean(both[which(abs(both$mid) > 4000 ),]$pi)
both$pi_red_low <- both$pi_low/mean(both[which(abs(both$mid) > 4000 ),]$pi)
both$lab<-'Both'
sqrt(sum((both$pi_red - cne$pi_red)^2)/length(cne$pi_red))

all2<-rbind(bgs,sel,both)
all2$dfe <- 'FulluSFS'

all<-rbind(all1,all2)

all$lab <- factor(all$lab, levels = c('Deleterious Mutations','Advantageous Mutations', 'Both'),
                  labels = c('Deleterious Mutations','Advantageous Mutations', 'Both'))


vnames <-list(
  'Deleterious Mutations' = 'Deleterious Mutations',
  'Advantageous Mutations' = 'Advantageous Mutations',
  'Both' ='Both',
  'FulluSFS' = 'Model A',
  'Ne-anc' = 'Model B'
  #  'Ne-anc' = bquote(N[e-anc])
)
vlabeller <- function(variable,value){
  return(vnames[value])
}

########################
## Now I have all the data in a usable format
########################
## Scaled pi

cairo_pdf('~/project/simulation_paper_stuff/plots/redux/CNEPiScaled.pdf', height = 8, width = 9)
ggplot(data=all,aes(x=mid/1000,y=pi_red))+
  # geom_line(col='black',lty=3)+
  geom_ribbon(aes(ymin=pi_red_low, ymax=pi_red_up, fill=lab ) ,alpha=0.8)+
  geom_line(data=cne,aes(x=dist/1000,y=pi_red))+
  geom_line(aes(y=1.0),lty=2)+
  scale_y_continuous(limits=c(0.83,1.05),position = 'left')+
  scale_fill_manual('', values = cbPalette)+
  xlab("Distance from CNE (Kb)")+
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
cairo_pdf('~/project/simulation_paper_stuff/plots/redux/CNEPiUnscaled.pdf', height = 8, width = 9)
ggplot(data=all,aes(x=mid/1000, y=pi))+
  geom_line(aes(y = 0.0083) ,col='black',lty=2)+
  geom_ribbon(aes(ymin=pi_lower, ymax=pi_upper, fill=lab ) ,alpha=0.8)+
  geom_line(data=cne, aes(x=dist/1000, y=pi_scale))+
  scale_fill_manual('', values = cbPalette)+
  xlab("Distance from CNE (Kb)")+
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

cairo_pdf('~/project/simulation_paper_stuff/plots/redux/CNEPi2.pdf', height = 8, width = 9)
ggplot(data=all, aes(x=mid/1000,y=pi2))+
  # geom_line(col='black',lty=3)+
  geom_ribbon(aes(ymin=pi2_lower, ymax=pi2_upper, fill=lab ) ,alpha=0.8)+
  geom_line(data=cne, aes(x=dist/1000, y=pi2))+
  # scale_y_continuous(limits = c(0.75,1.1))+
  scale_fill_manual('', values = cbPalette)+
  xlab("Distance from CNE (Kb)")+
  ylab(expression(pi[ratio]))+
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

ggplot(data=all, aes(x=mid/1000,y=tajima))+
  # geom_line(col='black',lty=3)+
  geom_ribbon(aes(ymin=tajima_lower, ymax=tajima_upper, fill=lab ) ,alpha=0.8)+
  geom_line(data=cne, aes(x=dist/1000, y=tajima+0.46))+
  # scale_y_continuous(limits = c(0.75,1.1))+
  scale_fill_manual('', values = cbPalette)+
  xlab("Distance from CNE (Kb)")+
  ylab(expression("Tajima's D"))+
  facet_grid(dfe~. , labeller = vlabeller)+
  theme_bw()+
  theme(
    text=element_text(family="Trebuchet MS"),
    axis.title.x = element_text(size=14,angle=0),
    axis.title.y = element_text(size=20,angle=90,vjust=0.5, face = 'italic'),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    legend.text = element_text(size =13),
    strip.text.y = element_text(size = 15)
  )


ggplot(data=all, aes(x=mid/1000,y=xsi))+
  # geom_line(col='black',lty=3)+
  geom_ribbon(aes(ymin=xsi_lower, ymax=xsi_upper, fill=lab ) ,alpha=0.8)+
  geom_line(data=cne, aes(x=dist/1000, y=xsi))+
  # scale_y_continuous(limits = c(0.75,1.1))+
  scale_fill_manual('', values = cbPalette)+
  xlab("Distance from CNE (Kb)")+
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
############# STOP















