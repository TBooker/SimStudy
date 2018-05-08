rm(list=ls())


exons<-read.csv('project/1.stats_from_exons/distance_from_exons/exons.1000.ncpg_summary.csv')
str(exons)
exons$pi_corrected<-(exons$pi/exons$mlk)*0.15

exons$block<-NULL
library(reshape2)
library(ggplot2)
str(exons)
exons_re<-melt(exons,id='dist')


ggplot(data=exons_re,aes(x=dist,y=value))+
  geom_line()+
  facet_wrap(~variable,scales='free')+
  theme_bw()+
  ggtitle('exons 100bp windows')+
  scale_y_continuous('')+
  scale_x_continuous('Distance From Exon (bp)')+
  theme(strip.text.x=element_text(size=10,face='bold'),
        title=element_text(size=20,face='bold'))

exons_bgs<-read.csv('project/2.simulations/all_elements/all/N1000_BGS/exons_s750_w_Na_N1000_BGS_stats.stats.csv')
exons_bgs<-exons_bgs[order(exons_bgs$position),]
str(exons_bgs)

exons_bgs$keyk<-NULL

c_b <- melt(exons_bgs,id='position')
c_b_m <- c_b[grep(pattern = "_", x = c_b$variable,invert=TRUE),]
c_b_m <- c_b_m[grep(pattern = "key", x = c_b_m$variable,invert=TRUE),]
c_b_m$value<-as.numeric(c_b_m$value)

ggplot(data=c_b_m,aes(x=position/1000,y=value))+
  geom_line()+
  facet_wrap(~variable,scales='free')+
  theme_bw()+
  ggtitle('BGS')+
  scale_y_continuous('')+
  scale_x_continuous(limits=c(-100,100))+
  theme(strip.text.x=element_text(size=10,face='bold'),
        title=element_text(size=20,face='bold'))
c_b_m$lab<-'Background Selection'


exons_sweeps<-read.csv('project/2.simulations/all_elements/all/N1000_SWEEP/exons_s750_w_Na_N1000_SWEEP_stats.stats.csv')
exons_sweeps<-exons_sweeps[order(exons_sweeps$position),]
exons_sweeps$key<-NULL

c_s <- melt(exons_sweeps,id='position')
c_s_m <- c_s[grep(pattern = "_", x = c_s$variable,invert=TRUE),]
c_s_m <- c_s_m[grep(pattern = "key", x = c_s_m$variable,invert=TRUE),]
c_s_m$value<-as.numeric(c_s_m$value)
ggplot(data=c_s_m,aes(x=position,y=value))+
  geom_line()+
  facet_wrap(~variable,scales='free')+
  theme_bw()+
  ggtitle('SWEEPS')+
  scale_y_continuous('')+
  scale_x_continuous(limits=c(-100,100))+
  theme(strip.text.x=element_text(size=10,face='bold'),
        title=element_text(size=20,face='bold'))
c_s_m$lab<-'Sweeps'


exons_bgs_sweeps<-read.csv('project/2.simulations/all_elements/all/N1000_BGS+SWEEP/exons_s750_w_Na_stats_castRates.stats.csv')
exons_bgs_sweeps<-exons_bgs_sweeps[order(exons_bgs_sweeps$position),]
exons_bgs_sweeps$key<-NULL


c_b_s <- melt(exons_bgs_sweeps,id='position')
c_b_s_m <- c_b_s[grep(pattern = "_", x = c_b_s$variable,invert=TRUE),]
c_b_s_m <- c_b_s_m[grep(pattern = "key", x = c_b_s_m$variable,invert=TRUE),]
c_b_s_m$value<-as.numeric(c_b_s_m$value)
ggplot(data=c_b_s_m,aes(x=position,y=value))+
  geom_line()+
  facet_wrap(~variable,scales='free')+
  theme_bw()+
  ggtitle('BGS+SWEEPS')+
  scale_y_continuous('')+
  scale_x_continuous(limits=c(-100,100))+
  theme(strip.text.x=element_text(size=10,face='bold'),
        title=element_text(size=20,face='bold'))
c_b_s_m$lab<-'Background Selection + Sweeps'

c_b_m$value<-c_b_m$value+0.0007
c_b_s_m$value<-c_b_s_m$value+0.00105
c_s_m$value<-c_s_m$value+0.0008

all<-rbind(c_b_m[which(c_b_m$variable=='pi'),],c_s_m[which(c_s_m$variable=='pi'),],c_b_s_m[which(c_b_s_m$variable=='pi'),])
all$position<-as.numeric(all$position)

str(exons)

ggplot(data=all,aes(x=position/1000,y=value))+
  geom_line(lwd=2)+
  theme_bw()+
  facet_wrap(~lab)+
  scale_y_continuous(limits=c(0.006,0.01))+
  ylab(expression(pi))+
  annotate(geom='line',x=exons$dist/1000,y=exons$pi_corrected,col='blue',alpha=0.5,lwd=2)+
  scale_x_continuous('Distance From exons (bp)',limits=c(-100,100))+
  theme(strip.text.x=element_text(size=20,face='bold'),
        axis.title.x=element_text(size=20,face='bold'),
        axis.title.y=element_text(size=30,angle=0,face='bold'),
        axis.text.y=element_text(size=17,color='black'),
        axis.text.x=element_text(size=17,color='black'))

