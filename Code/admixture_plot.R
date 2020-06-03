library(pophelper)
library(grid)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

setwd('./admixture_out/')
tbl<-read.table('Merge.9.Q')
colnames(tbl)<-paste0('Cluster',seq(1:9))
pop<-read.table('Merge.pop',col.names='pop')
id<-read.table('../merged_plink/Merge.fam',header=F)%>%select(id=V2)
tbl<-cbind(id,pop,tbl)


par(mar=c(1.5,4,2.5,2),cex.lab=0.75,cex.axis=0.6)
barplot(t(as.matrix(tbl)),
        col=brewer.pal(9,"Set1"),zlab="Anc.Proportions",
        border=NA,space=0)




Qfile<-list.files(path='./admixture_out',pattern="\\.Q",full.names=T)
q<-read.table(Qfile)
colnames(q)<-paste0('Cluster',seq(1:9))
pop<-read.table('./merged_plink/Merge.pop')%>%rename(pop=V1)
id<-read.table('./merged_plink/Merge.fam')%>%select(id=V2)
q<-cbind(id,pop,q)
m<-melt(q,id.vars=c("id","pop"),value.name='value')
asw<-m%>%filter(pop=='-')
arrange(asw,desc(Cluster1))
ggplot(asw,aes(id,value,fill=variable))+
  geom_bar(position='stack',stat='identity')

p<-plotQ(q_pop,
         exportplot=F,returnplot=T,
         sortind='all',
         ordergrp=T,
         showlegend=T)
p


View(q)
barplot(t(as.matrix(q)),col=rainbow(9),xlab="Individual #",ylab="Ancestry",border=NA)

