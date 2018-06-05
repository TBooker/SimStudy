
rm(list=ls())
library(ggplot2)


cne<-read.csv('~/project/1.stats_from_exons/intervals/analysis/strategy_2/cne/subtracted_beds/CNE_nCpG_unfolded.stats.4.csv')
cne$pi_k <- cne$pi/cne$t2
cne$pi_red <- cne$pi_k/mean(cne[which(abs(cne$dist) > 4000 ),]$pi_k)
cne$pi_scale <- cne$pi_red*0.0083
cne$lab<-'CNEs'

bgs<-read.csv("~/project/2.simulations/updated_DFE/CNEs/analysis/summary_BGS.Ne-anc.csv")
bgs$sites <- NA
bgs<-bgs[order(bgs$mid),]
bgs<- bgs[which(abs(bgs$mid) < max(cne$dist)+1), ]
bgs$pi_red <- bgs$pi/mean(bgs[which(abs(bgs$mid) > 4000 ),]$pi)
bgs$pi_red_up <- bgs$pi_up/mean(bgs[which(abs(bgs$mid) > 4000 ),]$pi)
bgs$pi_red_low <- bgs$pi_low/mean(bgs[which(abs(bgs$mid) > 4000 ),]$pi)
bgs$lab<-'Deleterious Mutations'

sel<-read.csv("~/project/2.simulations/updated_DFE/CNEs/analysis/summary_SWEEP.Ne-anc.csv")
sel<-sel[order(sel$mid),]
sel<- sel[which(abs(sel$mid) < max(cne$dist)+1), ]
sel$pi_red <- sel$pi/mean(sel[which(abs(sel$mid) > 4000 ),]$pi)
sel$pi_red_up <- sel$pi_up/mean(sel[which(abs(sel$mid) > 4000 ),]$pi)
sel$pi_red_low <- sel$pi_low/mean(sel[which(abs(sel$mid) > 4000 ),]$pi)
sel$lab<-'Advantageous Mutations'


both<-read.csv("~/project/2.simulations/updated_DFE/CNEs/analysis/summary_BGS+SWEEP.Ne-anc.csv")
both$sites <- NA
both<-both[order(both$mid),]
both<- both[which(abs(both$mid) < max(cne$dist)+1), ]
both$pi_red <- both$pi/mean(both[which(abs(both$mid) > 4000 ),]$pi)
both$pi_red_up <- both$pi_up/mean(both[which(abs(both$mid) > 4000 ),]$pi)
both$pi_red_low <- both$pi_low/mean(both[which(abs(both$mid) > 4000 ),]$pi)
both$lab<-'Both'


all1<-rbind(bgs,sel,both)
all1$dfe <- 'Model B'

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
all2$dfe <- 'Model A'

allCNE<-rbind(all1,all2)


rm(list=setdiff(ls(), c("allCNE","cne")))

cbPalette <- c("#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")

exons<-read.csv('~/project/1.stats_from_exons/intervals/analysis/strategy_2/TEST')

exons$pi_red <- exons$pi/mean(exons[which(abs(exons$dist) > 75000 ),]$pi)

exons$lab <- 'Exons'
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
all1$dfe <- 'Model B'

bgs<-read.csv("~/project/2.simulations/updated_DFE/Exons/analysis/summary_BGS.full_usfs.csv")
bgs<-bgs[order(bgs$mid),]
bgs$pi_red <- bgs$pi/mean(bgs[which(abs(bgs$mid) > 75000 ),]$pi)
bgs$pi_red_up <- bgs$pi_up/mean(bgs[which(abs(bgs$mid) > 75000 ),]$pi)
bgs$pi_red_low <- bgs$pi_low/mean(bgs[which(abs(bgs$mid) > 75000 ),]$pi)
bgs$lab<-"Deleterious Mutations"
bgs$red<-(bgs$pi/0.0083)
sqrt(sum((bgs$pi_red - exons$pi_red)^2))


sel<-read.csv("~/project/2.simulations/updated_DFE/Exons/analysis/summary_SWEEP.Ne-anc.csv")
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
all2$dfe <- 'Model A'



allExons <- rbind(all1,all2)
allExons$lab <- factor(all$lab, levels = c('Deleterious Mutations','Advantageous Mutations', 'Both'),
                  labels = c('Deleterious Mutations','Advantageous Mutations', 'Both'))

#pdf('test',width=12, height=8)
####################################
####################################
##              Tajima's D
####################################
####################################

allExons$label <- NULL
ExonsTajima <- data.frame('lab' = allExons$lab, 'distance' = allExons$mid, 'tajima' = allExons$tajima,'tajima_lower' = allExons$tajima_lower,'tajima_upper' = allExons$tajima_upper, element = 'Exons', dfe = allExons$dfe)
CNETajima <- data.frame('lab' = allCNE$lab,'distance' = allCNE$mid, 'tajima' = allCNE$tajima,'tajima_lower' = allCNE$tajima_lower,'tajima_upper' = allCNE$tajima_upper, element = 'CNEs', dfe = allCNE$dfe)
combined <- rbind(ExonsTajima, CNETajima)

combined$lab <- factor(combined$lab, levels = c('Deleterious Mutations','Advantageous Mutations', 'Both'),
                       labels = c('Deleterious Mutations','Advantageous Mutations', 'Both'))

empiric <- data.frame('dist' = c(exons$dist, cne$dist),'tajima' = c(exons$tajima, cne$tajima),'element' = c(exons$lab, cne$lab))

vnames <-list(
  'Deleterious Mutations' = 'Deleterious Mutations',
  'Advantageous Mutations' = 'Advantageous Mutations',
  'Both' ='Both',
  'FulluSFS' = 'Model A',
  'Ne-anc' = 'Model B',
  'Exons' = 'Exons',
  'CNEs' = 'CNEs'
  )

vlabeller <- function(variable,value){
  return(vnames[value])
}
#pdf('test',width=12, height=8)
cairo_pdf('~/project/simulation_paper_stuff/plots/redux/Tajima_ExonCNE.pdf', width = 9.75, height = 6)
ggplot(data=combined, aes(x=distance/1000,y=tajima))+
  geom_ribbon(aes(ymin=tajima_lower, ymax=tajima_upper, fill=lab ) ,alpha=0.8)+
  scale_fill_manual('', values = cbPalette)+
  xlab("Distance from Element (Kb)")+
  ylab("Tajima's D")+
  geom_line(data = empiric, aes(x=dist/1000, y=tajima))+
  facet_grid(dfe~element , scales = 'free_x')+
  theme_bw()+
  theme(
    text=element_text(family="Trebuchet MS"),
    axis.title.x = element_text(size=14,angle=0),
    axis.title.y = element_text(size=20,vjust=0.5, face= 'italic'),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    legend.text = element_text(size =13),
    strip.text.y = element_text(size = 15),
    strip.text.x = element_text(size = 15)
  )
dev.off()
####################################
##              Xsi
####################################
####################################


Exonsxsi <- data.frame('lab' = allExons$lab, 'distance' = allExons$mid, 'xsi' = allExons$xsi,'xsi_lower' = allExons$xsi_lower,'xsi_upper' = allExons$xsi_upper, element = 'Exons', dfe = allExons$dfe)
CNExsi <- data.frame('lab' = allCNE$lab,'distance' = allCNE$mid, 'xsi' = allCNE$xsi,'xsi_lower' = allCNE$xsi_lower,'xsi_upper' = allCNE$xsi_upper, element = 'CNEs', dfe = allCNE$dfe)

combined <- rbind(Exonsxsi, CNExsi)

combined$lab <- factor(combined$lab, levels = c('Deleterious Mutations','Advantageous Mutations', 'Both'),
                       labels = c('Deleterious Mutations','Advantageous Mutations', 'Both'))

empiric <- data.frame('dist' = c(exons$dist, cne$dist),'xsi' = c(exons$xsi, cne$xsi),'element' = c(exons$lab, cne$lab))

vnames <-list(
  'Deleterious Mutations' = 'Deleterious Mutations',
  'Advantageous Mutations' = 'Advantageous Mutations',
  'Both' ='Both',
  'FulluSFS' = 'Full uSFS',
  'Ne-anc' = 'uSFS Excl. Divergence',
  'Exons' = 'Exons',
  'CNEs' = 'CNEs'
)

vlabeller <- function(variable,value){
  return(vnames[value])
}


ggplot(data=combined, aes(x=distance/1000,y=xsi))+
  geom_ribbon(aes(ymin=xsi_lower+0.02, ymax=xsi_upper+0.02, fill=lab ) ,alpha=0.8)+
  scale_fill_manual('', values = cbPalette)+
  xlab("Distance from Element (Kb)")+
  ylab(expression(Xi))+
  geom_line(data = empiric, aes(x=dist/1000, y=xsi))+
  facet_grid(dfe~element , scales = 'free_x')+
  theme_bw()+
  theme(
 #   text=element_text(family="Trebuchet MS"),
    axis.title.x = element_text(size=14,angle=0),
    axis.title.y = element_text(size=20,angle=0,vjust=0.5, face= 'italic'),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    legend.text = element_text(size =13),
    strip.text.y = element_text(size = 15),
    strip.text.x = element_text(size = 15)
  )

####################################
####################################
##              Pi2
####################################
####################################

Exonspi2 <- data.frame('lab' = allExons$lab, 'distance' = allExons$mid, 'pi2' = allExons$pi2,'pi2_lower' = allExons$pi2_lower,'pi2_upper' = allExons$pi2_upper, element = 'Exons', dfe = allExons$dfe)
CNEpi2 <- data.frame('lab' = allCNE$lab,'distance' = allCNE$mid, 'pi2' = allCNE$pi2,'pi2_lower' = allCNE$pi2_lower,'pi2_upper' = allCNE$pi2_upper, element = 'CNEs', dfe = allCNE$dfe)

combined <- rbind(Exonspi2, CNEpi2)

combined$lab <- factor(combined$lab, levels = c('Deleterious Mutations','Advantageous Mutations', 'Both'),
                       labels = c('Deleterious Mutations','Advantageous Mutations', 'Both'))

empiric <- data.frame('dist' = c(exons$dist, cne$dist),'pi2' = c(exons$pi2, cne$pi2),'element' = c(exons$lab, cne$lab))



ggplot(data=combined, aes(x=distance/1000,y=pi2))+
  geom_ribbon(aes(ymin=pi2_lower+0.02, ymax=pi2_upper+0.02, fill=lab ) ,alpha=0.8)+
  scale_fill_manual('', values = cbPalette)+
  xlab("Distance from Element (Kb)")+
  ylab(expression(pi[ratio]))+
  geom_line(data = empiric, aes(x=dist/1000, y=pi2))+
  facet_grid(dfe~element , scales = 'free_x')+
  theme_bw()+
  theme(
 #   text=element_text(family="Trebuchet MS"),
    axis.title.x = element_text(size=14,angle=0),
    axis.title.y = element_text(size=20,angle=0,vjust=0.5, face= 'italic'),
    axis.text.x = element_text(size=12,angle=0),
    axis.text.y = element_text(size=12,angle=0),
    legend.text = element_text(size =13),
    strip.text.y = element_text(size = 15),
    strip.text.x = element_text(size = 15)
  )

#dev.off()




