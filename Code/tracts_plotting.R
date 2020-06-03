#!/usr/bin/env Rscript
#setwd('~/OneDrive/job_env/Data')
populations<-c('ASW','GBR','CEU','TSI','IBS','GWD','ESN','MSL','YRI','LWK')
ASW<-c('ASW')
European<-c('GBR','CEU','TSI','IBS')
African<-c('GWD','ESN','MSL','YRI','LWK')

library(unikn)
eur_col<-"#FBD154"
asw_col<-"#B978AA"
afr_col<-"#2776BB"
pal_superpops<-newpal(col=c(eur_col,asw_col,afr_col),names=c("EUR","ASW","AFR"))

eur_pal<-newpal(col=c("#f5db73","#e05eb7"))%>%seecol(n=4)%>%newpal(names=European)
afr_pal<-newpal(col=c("#85e69d","#3e419e"))%>%seecol(n=5)%>%newpal(names=African)
afr_other_pal<-seecol(c(afr_pal,"#c2c2c2"))%>%newpal(names=c(African,"Other"))
pal_pops<-seecol(c(eur_pal,afr_pal))

write.table(pal_superpops,'pal_superpops.txt',quote=F,row.names=F,sep="\n")

library(plyr)
library(dplyr)
library(ggpubr)
library(qtl)
library(reshape2)
library(ggplot2)
library(tidyr)
library(lattice)
library(data.table)
library(stringr)
library(chromoMap)
library(htmlwidgets)
library(webshot)
library(scales)
library(gridExtra)
library(taRifx)
library(scatterpie)
library(rlist)
library(gtable)

moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", 
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}


#Slave trade data
slaves<-read.csv("slave_trade_data.csv",header=T)%>%select(-Totals)%>%filter(Year!="Totals",Year!="1860")
slaves<-slaves%>%gather(pop,n,GWD,MSL,ESN_YRI,LWK,Other)
slaves$nk<-slaves$n/1000
slaves<-slaves%>%group_by(pop)%>%mutate(cum=cumsum(nk))
slaves$Year<-as.numeric(as.character(slaves$Year))
order.slaves<-data.frame(pop=c("Other","ESN_YRI","GWD","MSL","LWK"),order=c(1:5),
                         col=c("Other","ESN","GWD","MSL","LWK"))
slaves<-left_join(slaves,order.slaves,by="pop")%>%as.data.frame()

slaves.line<-ggplot(slaves,(aes(x=Year,y=cum,fill=reorder(col,order))))+
  geom_area(position="stack")+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(breaks=seq(from=1660,to=1820,by=20),
                     expand=c(0,0))+
  scale_fill_manual(values=afr_other_pal,
                    labels=c("Other African \nsource populations \n(West-Central African)",
                             "The Bight of Benin, \n The Bight of Biafra \n& The Gold Coast\n(includes ESN & YRI)",
                             "Senegambia \n(includes GWD)",
                             "Sierra Leone \n& The Windward Coast \n(includes MSL)",
                             "South-East Africa \n(includes LWK)"))+
  labs(fill="Source region:",
       y="Cumulative number of captives disembarked (thousands)",
       x="Year")+
  theme_classic()+
  theme(legend.position="top",
        axis.title=element_text(size=15),
        axis.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=11))
##ggsave(slaves.line,file="~/OneDrive/fyp_docs/Graphs_final/slaves_time.png")


percents.slaves<-slaves%>%select(pop,cum)%>%group_by(pop)%>%filter(cum==max(cum))%>%unique()
percents.slaves$percents<-(percents.slaves$cum/sum(percents.slaves$cum))*100
order.slaves<-data.frame(pop=c("Other","ESN_YRI","GWD","MSL","LWK"),order=c(1:5),
                         col=c("Other","ESN","GWD","MSL","LWK"))
percents.slaves<-left_join(percents.slaves,order.slaves,by="pop")
percents.slaves$col<-factor(percents.slaves$col,levels =percents.slaves$col[order(percents.slaves$order,decreasing=TRUE)])

slaves.per.pie<-ggplot(percents.slaves,aes(x=factor(1),y=cum,fill=factor(col)))+
  coord_polar('y',start=0,direction=1)+
  geom_col(color="white",show.legend=F)+
  scale_fill_manual(values=rev(afr_other_pal))+
  scale_color_manual(values=c("white"))+
  labs(fill="Population")+
  theme_void()
##ggsave(slaves.per.pie,file="~/OneDrive/fyp_docs/Graphs_final/slaves_percent_pie.png")

#With LWK
afr_pop_sum<-pop_sum%>%filter(ancestral_superpop=="African")
afr_pop_sum<-afr_pop_sum%>%mutate(total=sum(percent))
afr_pop_sum<-afr_pop_sum%>%mutate(percent=(percent/total)*100)
afr_pop_sum<-afr_pop_sum%>%select(pop=ancestral_pop,rf.per=percent)
ESN_YRI<-data.frame(pop=1,rf.per=1)
ESN_YRI$rf.per<-afr_pop_sum%>%filter(pop=="ESN"|pop=="YRI")%>%summarise(rf.per=sum(rf.per))
ESN_YRI$pop<-"ESN"
afr_pop_sum<-afr_pop_sum%>%filter(pop!="ESN"&pop!="YRI")
afr_pop_sum<-rbind(afr_pop_sum,ESN_YRI)

history<-percents.slaves%>%ungroup()%>%select(pop=col,hist.per=percents)
history_vs_rfmix<-left_join(afr_pop_sum,history,by="pop")

ESN_YRI<-long%>%filter(ancestral_pop=="ESN"|ancestral_pop=="YRI")%>%
  ddply(.(ID),numcolwise(sum))%>%mutate(ancestral_pop="ESN")%>%
  select(ID,ancestral_pop,value)
afr_pop_inds<-long%>%filter(ancestral_pop!="ESN"|ancestral_pop!="YRI")%>%
  select(ID,ancestral_pop,value)
afr_pop_inds<-rbind(ESN_YRI,afr_pop_inds)
afr_pop_inds<-ddply(afr_pop_inds,.(ancestral_pop,ID),numcolwise(mean))

afr_pop_inds<-afr_pop_inds%>%mutate(rf.per=value*100)%>%select(pop=ancestral_pop,rf.per)%>%filter(pop%in%African)

history_vs_rfmix<-left_join(afr_pop_inds,history,by="pop")%>%filter(pop!="YRI")

wLWK_r2<-summary(lm(history_vs_rfmix$hist.per~history_vs_rfmix$rf.per))$r.squared

w_lwk<-ggplot(history_vs_rfmix,aes(x=rf.per,y=hist.per))+
  geom_smooth(method=lm,se=T,color="black")+
  geom_jitter(aes(color=pop),width=0,height=0.1,size=3)+
  labs(x="Average global estimates of African ancestry proportions \nin modern ASW individuals",
       y="Historical African ancestry contributions during the \nTransatlantic Slave Trade")+
  scale_color_manual(values=afr_pal,labels=c("ESN & YRI","GWD","LWK","MSL"))+
  labs(color="Population")+
  theme_classic()+
  scale_x_continuous(expand=c(0.01,0))+
  scale_y_continuous(expand=c(0.01,0))+
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20)
        )

##ggsave(w_lwk,file="~/OneDrive/fyp_docs/Graphs_final/hist_vs_rfmix.png")

#Without LWK
history_vs_rfmix_noLWK<-history_vs_rfmix%>%filter(pop!="LWK")

w.oLWK_r2<-summary(lm(history_vs_rfmix_noLWK$hist.per~history_vs_rfmix_noLWK$rf.per))$r.squared

w.o_lwk<-ggplot(history_vs_rfmix_noLWK,aes(x=rf.per,y=hist.per))+
  geom_smooth(method=lm,se=T,color="black")+
  geom_jitter(aes(color=pop),width=0,height=0.2,size=3)+
  labs(x="Average global estimates of African ancestry proportions\nin modern ASW individuals",
       y="Historical African ancestry contributions during the\nTransatlantic Slave Trade")+
  scale_color_manual(values=afr_pal,labels=c("ESN & YRI","GWD","MSL"))+
  labs(color="Population")+
  scale_x_continuous(expand=c(0.01,0))+
  scale_y_continuous(expand=c(0.01,0))+
  theme_classic()+
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=20),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20))


##ggsave(w.o_lwk,file="~/OneDrive/fyp_docs/Graphs_final/hist_vs_rfmix.noLWK.png")

rf_vs_hist<-grid.arrange(w_lwk,w.o_lwk,nrow=2)


setwd('./rfmixout')
#RFMIX output
chromosomes<-c(1:22)
for (chr in chromosomes) {
  file<-paste0('rfmix_output_chr',chr,'.fb.tsv')
  df.name<-paste('rfmixout',chr,sep="")
  assign(df.name,read.table(file=file,header=T))
}
df.names <- ls(pattern='^rfmixout\\d+')
rfmixout<-'row.names<-'(do.call(rbind,mget(df.names)),NULL)

long<-melt(setDT(rfmixout),id.vars=colnames(rfmixout)[1:5],variable.name="INDIVIDUAL") %>%
  separate(col='INDIVIDUAL',into=c('ID','hap','ancestral_pop'),remove=F,extra='warn')
long$hap<-long$hap%>%str_remove('hap')
long<-long%>%mutate(ancestral_superpop=ifelse(ancestral_pop%in%African,'African',
                                              ifelse(ancestral_pop%in%European,'European',NA)))
long<-long%>%select(ID,physical_position,genetic_position,chromosome,hap,ancestral_pop,ancestral_superpop,value)
long_for_sum<-long%>%select(ID,physical_position,genetic_position,chromosome,hap,ancestral_superpop,value)
summed<-ddply(long_for_sum,.(ID,physical_position,genetic_position,chromosome,hap,ancestral_superpop),numcolwise(sum))
summed<-summed%>%unite(ID_pos,ID,physical_position,sep="_",remove=F)

#Load chromosome coordinates file (2 sets, 1 for each haplotype)
coords<-read.table('../chromomap_chromosome_coords.txt',
                   sep='\t',
                   header=F)
coordsA<-coords%>%mutate(V1=paste(coords$V1,'A',sep=""))
coordsB<-coords%>%mutate(V1=paste(coords$V1,'B',sep=""))
write.table(coordsA,
            file='chromomap_chromosome_coordsA.txt',
            sep='\t',
            quote=F,
            row.names=F,
            col.names=F)
write.table(coordsB,
            file='chromomap_chromosome_coordsB.txt',
            sep='\t',
            quote=F,
            row.names=F,
            col.names=F)

#Assign ancestry based on >.5 ancestry proportion (categorical)
one.ind.test<-data.frame(ID=long$ID%>%sample(1))
chromomap.df<-right_join(x=summed,y=one.ind.test,by='ID')
chromomap.df<-chromomap.df%>%
  filter(ancestral_superpop=='African')%>%
  select(ID_pos,chromosome,pos_start=physical_position,pos_end=physical_position,value,hap)
chromomap.df$putative_pop<-ifelse(chromomap.df$value<0.5,'European',
                                  ifelse(chromomap.df$value>0.5,'African',NA))
chromomap.df<-chromomap.df%>%select(-value)
chromomap.df1<-chromomap.df%>%filter(hap==1)%>%
  select(-hap)%>%mutate(chromosome=paste(chromosome,'A',sep=''))
chromomap.df2<-chromomap.df%>%filter(hap==2)%>%
  select(-hap)

write.table(chromomap.df1,
            file='chromomap_annotation1.txt',
            sep='\t',
            quote=F,
            row.names=F,
            col.names=F)
write.table(chromomap.df2,
            file='chromomap_annotation2.txt',
            sep='\t',
            quote=F,
            row.names=F,
            col.names=F)

plot<-chromoMap(c('chromomap_chromosome_coordsA.txt','../chromomap_chromosome_coords.txt'),
          c('chromomap_annotation1.txt','chromomap_annotation2.txt'),
          ploidy = 2,
          data_based_color_map = T,
          data_type = 'categorical',
          legend=c(T,T),
          lg_x=100,
          lg_y=250,
          data_colors=list(c(afr_col,eur_col),c(afr_col,eur_col)),
          chr_text=c(F,T),
          chr_color = 'lightgrey')
plot
#Save the widget as .html format
#setwd('../../Graphs')
#saveWidget(plot,file='ancestry_proportions_chromomap.html')
#Save as non-interactive .png
#webshot(url="ancestry_proportions_chromomap.html" , file="ancestry_proportions_chromomap.png", delay = 0.2)

one.ind.test
global.one.ind<-data.frame(1)
global.one.ind$afr<-chromomap.df%>%select(putative_pop)%>%filter(putative_pop=='African')%>%nrow()
global.one.ind$eur<-chromomap.df%>%select(putative_pop)%>%filter(putative_pop=='European')%>%nrow()
global.one.ind<-t(global.one.ind)%>%as.data.frame()
global.one.ind<-data.frame(anc=c("African","European"),
           percent=c(77.5,22.4))

pop_sum%>%
  ggplot(aes(x=factor(1),y=percent,fill=factor(ancestral_pop)))+
  geom_col(aes(colour=ancestral_superpop))+
  geom_text(aes(label=paste0(round(percent,0),'%'),colour=ancestral_superpop),
            position=position_stack(vjust=0.5),size=5,show.legend=F)+
  coord_polar('y',start=0,direction=-1)+
  theme()+
  scale_fill_manual(values=rev(pal_superpops))+
  scale_colour_manual(values=c("#e3eeff","#6e4466"))+
  labs(fill="Ancestral population")+
  theme_void()


ind<-ggplot(global.one.ind,aes(x=factor(1),y=percent,fill=factor(anc)))+
  geom_col(aes(fill=anc))+
  scale_fill_manual(values=c(afr_col,eur_col))+
  geom_text(aes(label=paste0(round(percent,0),'%')),
            position=position_stack(vjust=0.5),size=10,show.legend=F)+
  coord_polar('y',start=0,direction=1)+
  theme()+
  labs(fill="Ancestral population")+
  theme_void()
##ggsave(ind,file="~/Desktop/ind.test.pie.png")

#Do plots for highest and lowest individuals
##Highest
high_afr<-rfmix_props%>%filter(African==max(African))%>%select(ID)
one.ind.test<-data.frame(ID=high_afr[[1]])
chromomap.df<-right_join(x=summed,y=one.ind.test,by='ID')
chromomap.df<-chromomap.df%>%
  filter(ancestral_superpop=='African')%>%
  select(ID_pos,chromosome,pos_start=physical_position,pos_end=physical_position,value,hap)
chromomap.df$putative_pop<-ifelse(chromomap.df$value<0.5,'European',
                                  ifelse(chromomap.df$value>0.5,'African',NA))
chromomap.df<-chromomap.df%>%select(-value)
chromomap.df1<-chromomap.df%>%filter(hap==1)%>%
  select(-hap)%>%mutate(chromosome=paste(chromosome,'A',sep=''))
chromomap.df2<-chromomap.df%>%filter(hap==2)%>%
  select(-hap)

write.table(chromomap.df1,
            file='chromomap_annotation1.txt',
            sep='\t',
            quote=F,
            row.names=F,
            col.names=F)
write.table(chromomap.df2,
            file='chromomap_annotation2.txt',
            sep='\t',
            quote=F,
            row.names=F,
            col.names=F)

plot<-chromoMap(c('chromomap_chromosome_coordsA.txt','../chromomap_chromosome_coords.txt'),
                c('chromomap_annotation1.txt','chromomap_annotation2.txt'),
                ploidy = 2,
                data_based_color_map = T,
                data_type = 'categorical',
                legend=c(T,T),
                lg_x=100,
                lg_y=250,
                data_colors=list(c(afr_col,eur_col),c(afr_col,eur_col)),
                chr_text=c(F,T),
                chr_color = 'lightgrey')
plot
#Save the widget as .html format
setwd('../../Graphs')
#saveWidget(plot,file='ancestry_proportions_chromomap.html')
#Save as non-interactive .png
#webshot(url="ancestry_proportions_chromomap.html" , file="ancestry_proportions_chromomap.png", delay = 0.2)

one.ind.test
global.one.ind<-data.frame(1)
global.one.ind$afr<-chromomap.df%>%select(putative_pop)%>%filter(putative_pop=='African')%>%nrow()
global.one.ind$eur<-chromomap.df%>%select(putative_pop)%>%filter(putative_pop=='European')%>%nrow()
global.one.ind<-t(global.one.ind)%>%as.data.frame()
global.one.ind<-data.frame(anc=c("African","European"),
                           percent=c(98.1,1.9))


ind<-ggplot(global.one.ind,aes(x=factor(1),y=percent,fill=factor(anc)))+
  geom_col(aes(fill=anc))+
  scale_fill_manual(values=c(afr_col,eur_col))+
  coord_polar('y',start=0,direction=1)+
  theme()+
  theme_void()
##ggsave(ind,file="~/Desktop/ind.test.pie.png")


#Lowest
high_eur<-rfmix_props%>%filter(European==max(European))%>%select(ID)
one.ind.test<-data.frame(ID=high_eur[[1]])
chromomap.df<-right_join(x=summed,y=one.ind.test,by='ID')
chromomap.df<-chromomap.df%>%
  filter(ancestral_superpop=='African')%>%
  select(ID_pos,chromosome,pos_start=physical_position,pos_end=physical_position,value,hap)
chromomap.df$putative_pop<-ifelse(chromomap.df$value<0.5,'European',
                                  ifelse(chromomap.df$value>0.5,'African',NA))
chromomap.df<-chromomap.df%>%select(-value)
chromomap.df1<-chromomap.df%>%filter(hap==1)%>%
  select(-hap)%>%mutate(chromosome=paste(chromosome,'A',sep=''))
chromomap.df2<-chromomap.df%>%filter(hap==2)%>%
  select(-hap)

write.table(chromomap.df1,
            file='chromomap_annotation1.txt',
            sep='\t',
            quote=F,
            row.names=F,
            col.names=F)
write.table(chromomap.df2,
            file='chromomap_annotation2.txt',
            sep='\t',
            quote=F,
            row.names=F,
            col.names=F)

plot<-chromoMap(c('chromomap_chromosome_coordsA.txt','../chromomap_chromosome_coords.txt'),
                c('chromomap_annotation1.txt','chromomap_annotation2.txt'),
                ploidy = 2,
                data_based_color_map = T,
                data_type = 'categorical',
                legend=c(T,T),
                lg_x=100,
                lg_y=250,
                data_colors=list(c(eur_col,afr_col),c(eur_col,afr_col)),
                chr_text=c(F,T),
                chr_color = 'lightgrey')
plot
#Save the widget as .html format
#setwd('../../Graphs')
#saveWidget(plot,file='ancestry_proportions_chromomap.html')
#Save as non-interactive .png
#webshot(url="ancestry_proportions_chromomap.html" , file="ancestry_proportions_chromomap.png", delay = 0.2)

one.ind.test
global.one.ind<-data.frame(1)
global.one.ind$afr<-chromomap.df%>%select(putative_pop)%>%filter(putative_pop=='African')%>%nrow()
global.one.ind$eur<-chromomap.df%>%select(putative_pop)%>%filter(putative_pop=='European')%>%nrow()
global.one.ind<-t(global.one.ind)%>%as.data.frame()
global.one.ind<-data.frame(anc=c("African","European"),
                           percent=c(0.8,99.2))


ind<-ggplot(global.one.ind,aes(x=factor(1),y=percent,fill=factor(anc)))+
  geom_col(aes(fill=anc))+
  scale_fill_manual(values=c(afr_col,eur_col))+
  coord_polar('y',start=0,direction=1)+
  theme()+
  theme_void()

#Assign ancestry based on continuous ancestry proportions
chromomap.df<-right_join(x=summed,y=one.ind.test,by='ID')
chromomap.df<-chromomap.df%>%
  filter(ancestral_superpop=='African')%>%
  select(ID_pos,chromosome,pos_start=physical_position,pos_end=physical_position,value,hap)
chromomap.df1<-chromomap.df%>%filter(hap==1)%>%
  select(-hap)%>%mutate(chromosome=paste(chromosome,'A',sep=''))
chromomap.df2<-chromomap.df%>%filter(hap==2)%>%
  select(-hap)%>%mutate(chromosome=paste(chromosome,'B',sep=''))

write.table(chromomap.df1,
            file='chromomap_annotation1.txt',
            sep='\t',
            quote=F,
            row.names=F,
            col.names=F)
write.table(chromomap.df2,
            file='chromomap_annotation2.txt',
            sep='\t',
            quote=F,
            row.names=F,
            col.names=F)

plot<-chromoMap(c('chromomap_chromosome_coordsA.txt','chromomap_chromosome_coordsB.txt'),
                c('chromomap_annotation1.txt','chromomap_annotation2.txt'),
                ploidy = 2,
                data_based_color_map = T,
                data_type='numeric',
                chr_color='lightgrey',
                data_colors=list(c(eur_col,afr_col),c(eur_col,afr_col)),
                legend=c(T,F),
                chr_text=c(T,T),
                lg_x=100,
                lg_y=250)
plot
#Save the widget as .html format
setwd('../../Graphs')
#saveWidget(plot,file='ancestry_proportions_chromomap.html')
#Save as non-interactive .png
#webshot(url="ancestry_proportions_chromomap.html" , file="../../Graphs/ancestry_proportions_chromomap.png", delay = 0.2)


#Calculating global ancestry proportions
#Population-level
library(plyr)
library(dplyr)
pop_sum<-ddply(long,.(ancestral_pop,ancestral_superpop),numcolwise(sum))%>%select(-c(physical_position,genetic_position,chromosome))
pop_sum$percent<-pop_sum$value/sum(pop_sum$value)*100
pop_sum$ancestral_pop<-factor(pop_sum$ancestral_pop,
                        levels = pop_sum$ancestral_pop[order(pop_sum$percent,decreasing=TRUE)])

pop_pie<-pop_sum%>%
  ggplot(aes(x=factor(1),y=percent,fill=factor(ancestral_pop)))+
  coord_polar('y',start=0,direction=-1)+
  geom_col(color="white")+
  geom_text(aes(label=paste0(round(percent,0),'%'),color=ancestral_superpop),
            position=position_stack(vjust=0.5),size=5,show.legend=F)+
  scale_fill_manual(values=rev(pal_pops))+
  scale_color_manual(values=c("white","black"))+
  labs(fill="Ancestral \npopulation")+
  theme_void()
pop_pie
##ggsave(pop_pie,file='~/Desktop/pop_pie.png')


#Superpopulation-level
superpop_bar.df<-long%>%select(ID,value,ancestral_superpop)
superpop_bar.df<-ddply(superpop_bar.df,.(ID,ancestral_superpop),numcolwise(sum))%>%
  mutate(African=ifelse(ancestral_superpop=='African',value/(value+lead(value)),NA))%>%
  filter('African'==ancestral_superpop)%>%
  mutate(European=1-African)%>%
  select(ID,African,European)
rfmix_props<-superpop_bar.df
superpop_bar.df<-melt(superpop_bar.df,id.vars="ID")

superpop_bar<-ggplot(superpop_bar.df,aes(x=value,y=variable))+
  stat_boxplot(geom='errorbar',width=0.5)+
  geom_boxplot(width=0.6)+
  geom_jitter(aes(col=variable),size=2,alpha=0.9,width=0,height=0.15)+
  theme_classic()+
  labs(x='Ancestry proportion')+
  theme(legend.position="none",
        axis.title.y = element_blank(),
        axis.title.x=element_text(size=15),
        axis.text=element_text(colour="black",size=13))+
  scale_color_manual(values=c(afr_col,eur_col))+
  coord_equal(ratio=0.3)
superpop_bar

superpop_sum<-ddply(summed,.(ancestral_superpop),numcolwise(sum))
superpop_sum$percent<-superpop_sum$value/sum(superpop_sum$value)*100
superpop_sum$ancestral_pop<-factor(superpop_sum$ancestral_pop,
                              levels = superpop_sum$ancestral_pop[order(superpop_sum$percent,decreasing=TRUE)])

superpop_pie<-superpop_sum%>%
  ggplot(aes(x=factor(1),y=percent,fill=factor(ancestral_superpop)))+
  geom_col()+
  geom_text(aes(label=paste0(round(percent,0),'%')),
            position=position_stack(vjust=0.5),
            color=c("white","black"),size=6)+
  coord_polar('y',start=0,direction=-1)+
  theme_void()+
  theme(legend.position = "none")+
  scale_fill_manual(values=c(afr_col,eur_col))
superpop_pie

grid.arrange(superpop_bar,superpop_pie,nrow=1,widths=c(1.5,1))
superpop_props_plot<-arrangeGrob(superpop_bar,superpop_pie,nrow=1,widths=c(1.5,1))
##ggsave(superpop_props_plot,file='~/Desktop/superpops_props_plot.png')

#RFMIXout - global ancestry proportions
seecol(pal_superpops)

chromosomes<-c(1:22)
Qrfmixout<-list()
for (chr in chromosomes) {
  file<-paste0('rfmix_output_chr',chr,'.rfmix.Q')
  df<-read.table(file=file,header=F)
  df$chr<-chr
  colnames(df)<-c("ID","CEU","ESN","GBR","GWD","IBS","LWK","MSL","TSI","YRI","chr")
  Qrfmixout[[chr]]<-df
}
Qrfmixout<-'row.names<-'(do.call(rbind,Qrfmixout),NULL)%>%as.data.frame()
Qrfmixout<-Qrfmixout%>%group_by(ID)%>%mutate(African=sum(ESN,GWD,LWK,MSL,YRI)/22,European=sum(CEU,GBR,IBS,TSI)/22)
Qrfmixout<-Qrfmixout%>%gather(CEU:YRI,key="pop",value="prop",factor_key=T)%>%as.data.frame()

keys <- colnames(Qrfmixout)[!grepl('prop',colnames(Qrfmixout))]
Qrfmixout<- as.data.table(Qrfmixout)
dt<-Qrfmixout[,list(mean=mean(prop)),keys]

dt<-dt%>%
  arrange(factor(pop, levels = populations), desc(mean))

xorder.df<-data.frame(pop=populations[2:10],order=1:9)
dt<-left_join(dt,order.df,by="pop")


rfmix_global<-ggplot(data=dt,aes(x=reorder(ID, -African),y=mean*100,fill=reorder(pop,order)))+
  geom_col(width=1)+
  scale_fill_manual(values=seecol(pal_pops))+
  theme_classic()+
  theme(axis.text.x=element_blank())+
  scale_y_continuous(expand=c(0,0))+
  labs(fill="Ancestral \npopulation",
       x="ASW Individuals",
       y="Proportion of assigned tracts (%)")
##ggsave(rfmix_global,file="~/Desktop/rfmix_global.png")

#Tracts distribution
tracts<-summed%>%filter(ancestral_superpop=='African')%>%
  select(-ID_pos,-ancestral_superpop)%>%
  mutate(assigned_pop=ifelse(value<0.5,'European',
                             ifelse(value>0.5,'African',NA)))
tracts<-tracts%>%arrange(ID,chromosome,hap,physical_position)%>%select(-value)
tracts<-tracts%>%unite(id_chr_hap,ID,chromosome,hap,remove=F)
tracts<-tracts%>%mutate(changes=ifelse(id_chr_hap!=lag(id_chr_hap,default='start'),'start',
                                       ifelse(id_chr_hap!=lead(id_chr_hap,default='end'),'end',
                                              ifelse(assigned_pop!=lag(assigned_pop),'tract_intersection',NA))))
tracts<-tracts%>%filter(!is.na(changes))
tracts<-tracts%>%mutate(length.bp=ifelse(changes!='end',lead(physical_position)-physical_position,NA),
                        length.cM=ifelse(changes!='end',lead(genetic_position)-genetic_position,NA))%>%
  filter(!is.na(length.bp))

tracts_hist<-ggplot(data=tracts,aes(x=length.bp/1000000,fill=assigned_pop))+
  geom_histogram(binwidth=3,alpha=0.8)+
  theme_classic()+
  labs(x="Length (Mbp)",
       y="Number of blocks",
       fill="Ancestral \nsuperpopulation")+
  theme(legend.position=c(0.9,0.8),
        legend.text = element_text(size=12),
        legend.title=element_text(size=14),
        axis.text = element_text(size=17),
        axis.title=element_text(size=17))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=c(afr_col,eur_col))+
  expand_limits(x = c(0, 255))
tracts_hist
##ggsave(tracts_hist,file="~/Desktop/tracts_length_hist.png")

#Calculating blocks per individual versus global European ancestry proportion
#global_anc_per_individual<-summed%>%select(ID,ancestral_superpop,value)%>%unite(ID_superpop,ID,ancestral_superpop,remove=F)
#global_anc_per_individual<-ddply(global_anc_per_individual,.(ID_superpop,ID,ancestral_superpop),numcolwise(mean))

global_anc_per_individual<-Qrfmixout%>%select(-c(chr,pop,prop))%>%unique()%>%gather(superpop,value,African,European)%>%
  unite(ID_superpop,ID,superpop,remove=F)

blocks_per_individual<-tracts%>%unite(ID_superpop,ID,assigned_pop,remove=F)
blocks_per_individual<-data.frame(table(blocks_per_individual$ID_superpop))%>%select(ID_superpop=Var1,blocks=Freq)

time<-left_join(global_anc_per_individual,blocks_per_individual,by="ID_superpop")

#Modelling number of blocks at different times since admixture
L.cM<-3435 #whole genome length in cM
Time<-c(rep(seq(1,5,by=1),each=101))
value<-c(seq(0,1,by=0.01))
model<-data.frame(Time,value,L.cM)
model$blocks<-(2*2*0.01*Time*L.cM*value*(1-value))+2*22
model<-model%>%mutate(label=ifelse(value==0.5,Time,NA))

tracts_model_plot<-ggplot(NULL,aes(x=value,y=blocks))+
  geom_point(data=model,size=0.2)+
  geom_text(data=model,aes(label=label),nudge_y=8)+
  geom_point(data=time,aes(col=superpop))+
  theme_classic()+
  theme(legend.position=c(0.9,0.8),
        legend.text = element_text(size=12),
        legend.title=element_text(size=14),
        axis.text = element_text(size=12),
        axis.title=element_text(size=14))+
  labs(x='Genome-wide ancestry',y="Number of blocks",colour="Ancestral \nsuperpopulation")+
  scale_color_manual(values=c(afr_col,eur_col))

grid.arrange(tracts_hist,tracts_model_plot,nrow=2)
tracts_hist.model<-arrangeGrob(tracts_hist,tracts_model_plot,nrow=2)
##ggsave(tracts_hist.model,file="~/Desktop/tracts_hist_model.png")

#Creating TRACTS input file [chrom][begin][end][assignment][cmBegin][cmEnd]
#Reading in genetic maps for all chromosomes and combining
setwd('..')
tracts_input<-tracts%>%
  mutate(hap=LETTERS[as.numeric(hap)])%>%
  unite(id_hap,ID,hap,remove=F)%>%
  select(id_hap,chrom=chromosome,begin=physical_position,length.bp,cmBegin=genetic_position,length.cM,assignment=assigned_pop)%>%
  mutate(end=begin+length.bp,cmEnd=cmBegin+length.cM)%>%
  select(id_hap,chrom,begin,end,assignment,cmBegin,cmEnd)
tracts_input_list<-tracts_input%>%group_split(id_hap,keep=F)
names(tracts_input_list)<-unique(tracts_input$id_hap)

for( i in 1:length(tracts_input_list)){
  write.table(tracts_input_list[[i]],
              paste0('2pops/tracts_in/',names(tracts_input_list)[i],'.bed'),
              quote=F,sep='\t',row.names=F)
}


#Producing TRACTS output graphic
bins<-read.table('./2pops/out_bins',header=F,sep='\t')
dat<-read.table('./2pops/out_dat',header=F,sep='\t')
mig<-read.table('./2pops/out_mig',header=F,sep='\t')%>%select(EUR=V1,AFR=V2) #migration matrix (sizes/locations of pie charts)
pred<-read.table('./2pops/out_pred',header=F,sep='\t') #predicted counts in each bin by the model
pars<-read.table('./2pops/out_pars',header=F,sep='\t')%>%select(pp=V1) #time of the pulse
bins2<-read.table('./2pops/out2_bins',header=F,sep='\t')
dat2<-read.table('./2pops/out2_dat',header=F,sep='\t')
mig2<-read.table('./2pops/out2_mig',header=F,sep='\t')%>%select(EUR=V1,AFR=V2)
pred2<-read.table('./2pops/out2_pred',header=F,sep='\t')
pars2<-read.table('./2pops/out2_pars',header=F,sep='\t')%>%select(pp=V1,px=V2,eur_px=V3) 

#pp.py
##param - 0.03305922571234303
##likelihoods found:  [-1614.8278961725737, -326.0848869818305]
#pp_px.py
##params - 0.07684188369069068	0.039908476998231815	0.2006785935865797
#likelihoods found: -330.029
#likelihood_ratio=log(-330/-1614)

#pp graphing
mig.matrix<-mig%>%
  mutate(GA=1:nrow(mig))%>%
  arrange(-GA)%>%
  mutate(yr=2014-(GA*25),order=1:nrow(mig))
mig.matrix<-mig.matrix%>%gather(anc,contrib,EUR,AFR,factor_key = T)
mig.matrix<-mig.matrix%>%group_by(anc)%>%mutate(cum=cumsum(contrib))%>%
  group_by(yr)%>%mutate(percent=cum/sum(cum))%>%group_by(order)
pie1<-mig%>%mutate(GA=1:nrow(mig))%>%
  arrange(-GA)%>%
  mutate(yr=2014-GA*25,pos=1,
         order=1:nrow(mig),
         contrib=as.numeric(0.5*(EUR+AFR)/(sum(EUR)+sum(AFR))))%>%
  group_by(GA)
pp<-ggplot()+
  geom_bar(data=mig.matrix,stat='identity',width=1,
           aes(x=order-0.5,y=percent*2,fill=anc))+
  scale_x_continuous(breaks=c(0:nrow(mig)),
                     labels=c(unique(mig.matrix$yr),"2014"))+
  geom_scatterpie(data=pie1,
                  aes(x=order-1,y=pos+1.5,r=contrib),
                  cols=c("EUR","AFR"),color=NA)+
  coord_fixed()+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.x=element_line(),
        axis.text.x=element_text())+
  #annotate("text",x=6,y=2.6,label="Single-pulse model",hjust="right",fontface=2,size=5.5)+
  #annotate("text",x=6,y=2.3,label="L = -1615",hjust="right",size=5)+
  theme(axis.title.x=element_text(vjust=-0.35,hjust=0.52),
        legend.position='none')+
  labs(fill="Ancestral \npopulation",
       x="Year")+
  scale_fill_manual(values=c(eur_col,afr_col))
pp

#pp_px graphing
mig.matrix2<-mig2%>%
  mutate(GA=1:nrow(mig2))%>%
  arrange(-GA)%>%
  mutate(yr=2014-GA*25,order=1:nrow(mig2))%>%
  gather(anc,contrib,EUR,AFR,factor_key = T)%>%
  group_by(anc)%>%
  mutate(cum=cumsum(contrib))%>%
  group_by(yr)%>%
  mutate(percent=cum/sum(cum))%>%
  group_by(order)
pie2<-mig2%>%mutate(GA=1:nrow(mig2))%>%
  arrange(-GA)%>%
  mutate(yr=2014-GA*25,pos=1,
         order=1:nrow(mig2),
         contrib=as.numeric((EUR+AFR)/(sum(EUR)+sum(AFR))))%>%
  group_by(GA)

pp_px<-ggplot()+
  geom_bar(data=mig.matrix2,stat='identity',width=1,
           aes(x=order+0.5,y=percent*3.2,fill=anc))+
  scale_x_continuous(breaks=seq(nrow(pie2)+1),
                     labels=c(unique(mig.matrix2$yr),"2014"))+
  geom_scatterpie(data=pie2,
                  aes(x=order,y=pos+3.,r=contrib),
                  cols=c("EUR","AFR"),color=NA)+
  coord_fixed()+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.x=element_line(),
        axis.text.x=element_text())+
  #annotate("text",x=10,y=4.5,label="Two-pulse model",hjust="right",fontface=2,size=5.5)+
  #annotate("text",x=10,y=3.7,label="L = -326",hjust="right",size=5)+
  theme(axis.title.x=element_text(),
        legend.position="bottom")+
  labs(fill="Ancestral population",x="Year")+
  scale_fill_manual(values=c(eur_col,afr_col))
pp_px

##ggsave(pp,file="~/Desktop/pp_graph.png")
##ggsave(pp_px,file="~/Desktop/pp_px_graph.png")


grid.arrange(pp,pp_px,nrow=2,ncol=1,heights=c(1,4))
tracts_models<-arrangeGrob(pp,pp_px,nrow=2,ncol=1,heights=c(1,4))
##ggsave(tracts_models,file="~/Desktop/tracts_models.png")

#CV error across different K values for the unsupervised ADMIXTURE run
cv<-read.table("admixture_unsupervised.Kx/runs1/K_search.CV.txt")
colnames(cv)<-c("K","CV")
min<-cv %>%
  group_by(K) %>%
  summarise(CV = min(CV))
min$color<-1
cv<-left_join(cv,min,by="CV")
cv$color=as.character(cv$color)

cv_plot<-ggplot(data=cv,aes(x=K,y=CV))+
  geom_line(size=2.5)+
  geom_point(size=5,aes(color=color))+
  scale_color_manual(values=c("red","black"))+
  scale_x_continuous(breaks=seq(1,12,1))+
  scale_y_continuous(breaks=seq(0.5,0.8,0.05))+
  theme_classic()+
  theme(axis.text=element_text(size=20),
      axis.title=element_text(size=22),
      legend.position = "none")+
  labs(y="Cross-validation value",
       x="Number of ancestral population (K) passed to ADMIXTURE")
#ggsave(cv_plot,file="~/OneDrive/fyp_docs/Graphs_final/cv_plot.png",width=10,height=8)

#Comparing RFMIX and ADMIXTURE global ancestry estimates
rfmix_props<-long%>%select(ID,value,ancestral_superpop)
rfmix_props<-ddply(rfmix_props,.(ID,ancestral_superpop),numcolwise(sum))%>%
  mutate(African=ifelse(ancestral_superpop=='African',value/(value+lead(value)),NA))%>%
  filter('African'==ancestral_superpop)%>%
  mutate(European=1-African)%>%
  select(ID,European,African)
rfmix_props<-rfmix_props%>%mutate(programme="RFMix")

admix.q<-read.table('./admixture_unsupervised.K2/Merge.2.Q',header=F,row.names=NULL)
fam<-read.table('./admixture_unsupervised.K2/Merge.fam',sep='\t',header=F)%>%select(V2)
pop<-read.table('./admixture_unsupervised.K2/Merge.pop',sep='\t',header=F)
admix_props<-cbind(fam,pop,admix.q)
colnames(admix_props)<-c('ID','pop','European','African')
admix_props<-admix_props%>%mutate(programme='ADMIXTURE')%>%filter("-"==pop)%>%select(-pop)

compare<-rbind(admix_props,rfmix_props)
compare<-gather(compare,population,value,European:African,factor_key=T)
compare<-spread(compare,programme,value)

compare_plot<-ggplot(data=compare,aes(x=ADMIXTURE,y=RFMix,col=population))+
  geom_point()+
  theme_classic()+
  theme(title = element_blank(),
        legend.position=c(0.32,0.902),
        legend.text=element_text(size=18),
        axis.title= element_text(size=20),
        axis.text=element_text(size=18))+
  labs(x="ADMIXTURE estimates of global ancestry",y="RFMix estimates of global ancestry")+
  stat_cor(aes(label = ..rr.label..),show.legend = F,size=6)+
  scale_color_manual(values=c(eur_col,afr_col))+
  coord_fixed()
compare_plot
##ggsave(compare_plot,file="~/Desktop/compare_rf_vs_admix.png")

#Assortative mating

#Do a permutation test for genome-wide AMI to test for baseline ancestry-assortative mating and calculate the AMI values for the gene sets under random mating
library(splitstackshape)

all<-long%>%select(chr=chromosome,pos=physical_position)%>%
  unite(chr.pos,chr,pos,sep=".",remove=F)%>%
  sample_n(1000)
df<-tracts_input%>%unite(id_hap_chr,id_hap,chrom,remove=F)

shuffle<-data.frame(id_hap_chr=df$id_hap_chr)
shuffle.counts<-shuffle%>%table()%>%as.data.frame()
colnames(shuffle.counts)<-c("id_hap_chr","count")
shuffle<-left_join(unique(shuffle),shuffle.counts,by="id_hap_chr")
shuffle<-shuffle%>%separate(col=id_hap_chr,into=c("id","hap","chr"),sep="\\_",remove=F)%>%
  select(id_hap_chr,chr,count)
shuffle$id_hap_chr<-ave(shuffle$id_hap_chr,shuffle$chr,FUN = sample)
shuffle<-expandRows(shuffle,"count",count.is.col = TRUE,drop=T)
shuffle$s.id_hap_chr<-shuffle$id_hap_chr

df<-cbind(shuffle$s.id_hap_chr,df)%>%select(id_hap_chr="shuffle$s.id_hap_chr",chrom,begin,end,assignment,cmBegin,cmEnd)
df<-df%>%separate(col=id_hap_chr,into=c("id","hap","chr"),sep="_")%>%
  select(-chr)%>%
  unite(id_hap,id,hap,remove=T)

shuffled_tracts_list<-list()
for (i in 1:nrow(all)) {
  l<-filter(df,chrom==all[i,]$chr & begin<all[i,]$pos & end>all[i,]$pos)
  l<-l%>%mutate(chr.pos=all[i,1])
  shuffled_tracts_list[[i]]<-l
}
names(shuffled_tracts_list)<-as.character(all$chr.pos)
shuffled_tracts_list<-lapply(shuffled_tracts_list,separate,col=id_hap,into=c("id","hap"),sep="\\_")
shuffled_tracts_list<-lapply(shuffled_tracts_list,select,id,hap,assignment,chrom,chr.pos)
shuffled_tracts_list<-lapply(shuffled_tracts_list,spread,hap,assignment)
shuffled_tracts_list<-shuffled_tracts_list[sapply(shuffled_tracts_list, nrow)>0]
shuffled_tracts_list<-lapply(shuffled_tracts_list,mutate,genotype=ifelse(A==B,A,'het'))
shuffled_tracts_list<-lapply(shuffled_tracts_list,filter,!is.na(genotype))


shuffled_hom_vs_het<-list()
for (i in 1:length(shuffled_tracts_list)) {
  chr.pos<-names(shuffled_tracts_list)[i]
  chr<-strsplit(chr.pos, "[.]")[[1]][1]
  o_het<-shuffled_tracts_list[[i]]%>%filter(genotype=="het")%>%nrow()
  n<-shuffled_tracts_list[[i]]%>%nrow()
  o_hom_afr<-shuffled_tracts_list[[i]]%>%filter(genotype=="African")%>%nrow()
  o_hom_eur<-shuffled_tracts_list[[i]]%>%filter(genotype=="European")%>%nrow()
  df<-data.frame(chr.pos=names(shuffled_tracts_list)[i],chr,
                 o_hom_afr,o_hom_eur,n=n,o_het=o_het)
  shuffled_hom_vs_het[[i]]<-df
}
shuffled_hom_vs_het<-`row.names<-`(do.call(rbind,shuffled_hom_vs_het), NULL)
shuffled_hom_vs_het<-shuffled_hom_vs_het%>%filter(n==61)%>%
  mutate(o_hom=n-o_het)

ancestral_fractions<-superpop_sum%>%select(anc=ancestral_superpop,fraction=percent)%>%
  mutate(fraction= fraction / 100)%>%
  spread(anc,fraction)
shuffled_hom_vs_het<-cbind(shuffled_hom_vs_het,ancestral_fractions)
shuffled_hom_vs_het<-shuffled_hom_vs_het%>%mutate(e_het=2*African*European*n,
                                        e_hom=((African*African+European*European)*n))
shuffled_hom_vs_het<-shuffled_hom_vs_het%>%mutate(ami=log((o_hom/e_hom)/(o_het/e_het)))%>%
  arrange(ami)
shuffled_hom_vs_het<-shuffled_hom_vs_het%>%mutate(hom_afr_per=o_hom_afr/(o_hom_afr+o_hom_eur),
                                        hom_eur_per=o_hom_eur/(o_hom_afr+o_hom_eur))

shuffled_pooled<-data.frame(ami=mean(shuffled_hom_vs_het$ami),
                       stdev=sd(shuffled_hom_vs_het$ami),
                       samples=mean(shuffled_hom_vs_het$n))%>%
  mutate(CI.95=(1.96*stdev)/sqrt(samples))

ggplot(shuffled_hom_vs_het,aes(x=chr.pos,y=ami))+
  geom_hline(aes(yintercept=0),linetype="dashed")+
  geom_point()+
  coord_flip()+
  theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank())+
  labs(y="Assortative Mating Index")

#Calculate the average genome-wide assortative mating index (based on broad ancestral homozygosity)
#all<-long%>%select(chr=chromosome,pos=physical_position)%>%
#  unite(chr.pos,chr,pos,sep=".",remove=F)%>%
#  sample_n(1000)

df<-tracts_input
all_tracts_list<-list()
for (i in 1:nrow(all)) {
  l<-filter(df,chrom==all[i,]$chr & begin<all[i,]$pos & end>all[i,]$pos)
  l<-l%>%mutate(chr.pos=all[i,1])
  #assign(paste0(all[i,1],"_matches"),l)
  all_tracts_list[[i]]<-l
}
names(all_tracts_list)<-as.character(all$chr.pos)
all_tracts_list<-lapply(all_tracts_list,separate,col=id_hap,into=c("id","hap"),sep="\\_")
all_tracts_list<-lapply(all_tracts_list,select,id,hap,assignment,chr.pos)
all_tracts_list<-lapply(all_tracts_list,spread,hap,assignment)
all_tracts_list<-all_tracts_list[sapply(all_tracts_list, nrow)>0]
all_tracts_list<-lapply(all_tracts_list,mutate,genotype=ifelse(A==B,A,'het'))
all_tracts_list<-lapply(all_tracts_list,filter,!is.na(genotype))

all_hom_vs_het<-list()
for (i in 1:length(all_tracts_list)) {
  chr.pos<-names(all_tracts_list)[i]
  o_het<-all_tracts_list[[i]]%>%filter(genotype=="het")%>%nrow()
  n<-all_tracts_list[[i]]%>%nrow()
  o_hom_afr<-all_tracts_list[[i]]%>%filter(genotype=="African")%>%nrow()
  o_hom_eur<-all_tracts_list[[i]]%>%filter(genotype=="European")%>%nrow()
  df<-data.frame(chr.pos=names(all_tracts_list)[i],
                 o_hom_afr,o_hom_eur,n=n,o_het=o_het)
  all_hom_vs_het[[i]]<-df
}
all_hom_vs_het<-`row.names<-`(do.call(rbind,all_hom_vs_het), NULL)
all_hom_vs_het<-all_hom_vs_het%>%filter(n==61)%>%
  mutate(o_hom=n-o_het)

ancestral_fractions<-superpop_sum%>%select(anc=ancestral_superpop,fraction=percent)%>%
  mutate(fraction= fraction / 100)%>%
  spread(anc,fraction)
all_hom_vs_het<-cbind(all_hom_vs_het,ancestral_fractions)
all_hom_vs_het<-all_hom_vs_het%>%mutate(e_het=2*African*European*n,
                                        e_hom=((African*African+European*European)*n))
all_hom_vs_het<-all_hom_vs_het%>%mutate(ami=log((o_hom/e_hom)/(o_het/e_het)))%>%
  arrange(ami)
all_hom_vs_het<-all_hom_vs_het%>%mutate(hom_afr_per=o_hom_afr/(o_hom_afr+o_hom_eur),
                                        hom_eur_per=o_hom_eur/(o_hom_afr+o_hom_eur))
all_hom_vs_het<-all_hom_vs_het%>%mutate(e_hom_afr=(African*African*n),
                                e_hom_eur=(European*European*n),
                                ah_afr=(o_hom_afr-e_hom_afr)/e_hom_afr,
                                ah_eur=(o_hom_eur-e_hom_eur)/e_hom_eur)

all_ah<-all_hom_vs_het%>%summarise(ah_afr=sum(ah_afr),
                                            ah_eur=sum(ah_eur))
all_ah<-all_ah%>%gather(ah_afr,ah_eur,key="anc",value="ah")
all_ah<-all_ah%>%mutate(Ancestry=ifelse(anc=="ah_afr","African","European"))
all_ah_plot<-ggplot(all_ah,aes(y=ah,fill=Ancestry))+
  geom_col(position="dodge",aes(x=1))+
  geom_vline(aes(xintercept=0))+
  labs(y="Ancestry homozygosity")+
  scale_fill_manual(values=c(afr_col,eur_col))+
  scale_x_discrete(limits=c("Genome-wide"))+
  theme_classic()+
  theme(legend.position=c(0.25,0.85),
        axis.text.x=element_text(size=15),
        axis.title.x=element_blank(),
        legend.title=element_text(size=15),
        legend.text=element_text(size=13),
        axis.text.y=element_text(size=15),
        axis.title.y=element_text(size=15))
ggsave(all_ah_plot,file="~/Desktop/all_ancestry_homozygosity_plot.jpeg")

all_ah$phen<-"Genome-wide"
ah$set<-"genes"
all_ah$set<-"genome"

library(lattice)

ah<-rbind(all_ah,ah)
ggplot(ah,aes(y=ah,x=phen,fill=Ancestry))+
  geom_col(position="dodge")+
  facet_grid(~set,
             scales="free_y")


all_ah_plot<-ggplot(NULL,aes(y=ah,x=phen,fill=Ancestry))+
  geom_col(data=ah,position="dodge")+
  geom_col(data=all_ah,position="dodge")
  geom_vline(aes(xintercept=0))+
  labs(y="Ancestry homozygosity")+
  scale_fill_manual(values=c(afr_col,eur_col))+
  scale_x_discrete(limits=c("Genome-wide"))+
  theme_classic()+
  theme(legend.position=c(0.8,0.8),
        axis.text.x=element_text(size=15),
        axis.title.x=element_blank(),
        legend.title=element_text(size=15),
        legend.text=element_text(size=13),
        axis.text.y=element_text(size=15),
        axis.title.y=element_text(size=15))

all_pooled<-data.frame(ami=mean(all_hom_vs_het$ami),
                       stdev=sd(all_hom_vs_het$ami),
                       samples=mean(all_hom_vs_het$n))%>%
  mutate(CI.95=(1.96*stdev)/sqrt(samples))


ggplot(all_hom_vs_het,aes(x=chr.pos,y=ami))+
  geom_hline(aes(yintercept=0),linetype="dashed")+
  geom_point()+
  coord_flip()+
  theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank())+
  labs(y="Assortative Mating Index")

#Do AMI test for genes of interest
sample_info<-read.table('20130606_g1k.ped',sep='\t',header=T)
paired_ASW<-sample_info%>%filter('ASW'==Population,c('father','mother')==Relationship)%>%arrange(Family.ID)%>%
  select(Family.ID,ID=Individual.ID,Relationship)
paired_ASW<-left_join(global_anc_per_individual,paired_ASW,by="ID")%>%filter(!is.na(Family.ID))
spouses<-c("2368","2424","2431","2434","2436") #Families with both mother and father in dataset
paired_ASW<-paired_ASW%>%filter(Family.ID%in%spouses,superpop=="European")

couples<-paired_ASW%>%select(-ID_superpop,-ID)%>%spread(key="Relationship",value)
ggplot(couples,aes(x=mother*100,y=father*100))+
  geom_point(size=5)+
  stat_cor(aes(label = ..r.label..), label.x = 75, size=15)+
  geom_smooth(method='lm',se=F,size=2,color="darkgrey")+
  theme_classic()+
  scale_x_continuous(expand=c(0.05,0))+
  scale_y_continuous(expand=c(0.05,0))+
  labs(y="Father's global African ancestry (%)",
       x="Mother's global African ancestry (%)")+
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25))

paired_ASW%>%spread(key="Relationship",ID)
paired_ASW%>%select(-ID_superpop)%>%spread(key="superpop",value)


#Alleles (pigmentation - polygenic phenotype) and their significant trait-associated SNPs.  
geneset<-list(MC1R<-data.frame(gene="MC1R",rsid="rs1805007",phen="skin"),
              OCA2<-data.frame(gene="OCA2",rsid="rs12913832",phen="skin"),
              TYR<-data.frame(gene="TYR",rsid="rs1042602",phen="skin"),
              TYRP1<-data.frame(gene="TYRP1",rsid="rs1408799",phen="skin"),
              DCT<-data.frame(gene="DCT",rsid="rs1407995",phen="skin"),
              SLC45A2<-data.frame(gene="SLC45A2",rsid="rs16891982",phen="skin"),
              SLC24A5<-data.frame(gene="SLC24A5",rsid="rs1426654",phen="skin"),
              KITLG<-data.frame(gene="KITLG",rsid="rs642742",phen="skin"),
              IRF4<-data.frame(gene="IRF4",rsid="rs12203592",phen="skin"),
              ASIP<-data.frame(gene="ASIP",rsid="rs1015362",phen="skin"),
              TCHH<-data.frame(gene="TCHH",rsid="rs11803731",phen="hair"),
              PADI3<-data.frame(gene="PADI3",rsid="rs11203346",phen="hair"),
              MATP<-data.frame(gene="MATP",rsid="rs16891982",phen="eye"),
              HERC2<-data.frame(gene="HERC2",rsid="rs1129038",phen="eye"))
geneset<-rbindlist(geneset)
write.table(geneset$rsid,'assortative_mating_analysis/geneset_rsids.txt',
            quote=F,row.names = F,col.names = F) #(run block in run.sh to extract matches)

#Count the number of homozygous and heterozygous individuals at these positions
matched<-read.table('assortative_mating_analysis/matched.txt')
colnames(matched)<-c('chr','pos','rsid','eur_af','afr_af','ref','alt')
matched<-left_join(geneset,matched,by='rsid')
matched<-remove.factors(matched)

df<-tracts_input
gene_tracts_list<-list()
for (i in 1:nrow(matched)) {
  l<-filter(df,chrom==matched[i,]$chr & begin<matched[i,]$pos & end>matched[i,]$pos)
  l<-l%>%mutate(gene=matched[i,1],phen=matched[i,3])
  #assign(paste0(matched[i,1],"_matches"),l)
  gene_tracts_list[[i]]<-l
}
names(gene_tracts_list)<-as.character(matched$gene)
gene_tracts_list<-lapply(gene_tracts_list,separate,col=id_hap,into=c("id","hap"),sep="\\_")
gene_tracts_list<-lapply(gene_tracts_list,select,id,hap,assignment,gene,phen)
gene_tracts_list<-lapply(gene_tracts_list,spread,hap,assignment)
gene_tracts_list<-lapply(gene_tracts_list,mutate,genotype=ifelse(A==B,A,
                                                                 ifelse(A!=B,'het',NA)))
hom_vs_het<-list()
for (i in 1:nrow(matched)) {
  gene<-names(gene_tracts_list)[i]
  o_het<-gene_tracts_list[[i]]%>%filter(genotype=="het")%>%nrow()
  n<-gene_tracts_list[[i]]%>%nrow()
  o_hom_afr<-gene_tracts_list[[i]]%>%filter(genotype=="African")%>%nrow()
  o_hom_eur<-gene_tracts_list[[i]]%>%filter(genotype=="European")%>%nrow()
  df<-data.frame(gene=names(gene_tracts_list)[i],phen=gene_tracts_list[[i]][[3]][[1]],
                 o_hom_afr,o_hom_eur,n=n,o_het=o_het)
  hom_vs_het[[i]]<-df
}
hom_vs_het<-`row.names<-`(do.call(rbind,hom_vs_het), NULL)
hom_vs_het<-hom_vs_het%>%mutate(o_hom=n-o_het)

ancestral_fractions<-superpop_sum%>%select(anc=ancestral_superpop,fraction=percent)%>%
  mutate(fraction= fraction / 100)%>%
  spread(anc,fraction)
hom_vs_het<-cbind(hom_vs_het,ancestral_fractions)
hom_vs_het<-hom_vs_het%>%mutate(e_het=2*African*European*n,
                                e_hom=((African*African+European*European)*n),
                                e_hom_afr=(African*African*n),
                                e_hom_eur=(European*European*n))
hom_vs_het<-hom_vs_het%>%mutate(ami=log((o_hom/e_hom)/(o_het/e_het)))%>%
  arrange(ami)
hom_vs_het<-hom_vs_het%>%mutate(hom_afr_per=o_hom_afr/(o_hom_afr+o_hom_eur),
                                hom_eur_per=o_hom_eur/(o_hom_afr+o_hom_eur))
hom_vs_het<-hom_vs_het%>%mutate(ah_afr=(o_hom_afr-e_hom_afr)/e_hom_afr,
                                ah_eur=(o_hom_eur-e_hom_eur)/e_hom_eur)

ah<-hom_vs_het%>%group_by(phen)%>%summarise(ah_afr=sum(ah_afr),
                                            ah_eur=sum(ah_eur))
ah<-ah%>%gather(ah_afr,ah_eur,key="anc",value="ah")
ah<-ah%>%mutate(Ancestry=ifelse(anc=="ah_afr","African","European"))
ah_plot<-ggplot(ah,aes(x=ah,y=phen,fill=Ancestry))+
  geom_col(position="dodge")+
  geom_vline(aes(xintercept=0))+
  labs(y="Trait",
       x="Ancestry homozygosity")+
  scale_fill_manual(values=c(afr_col,eur_col))+
  scale_y_discrete(limits=c("eye","hair","skin"))+
  theme_classic()+
  theme(axis.line.y =element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        legend.position=c(0.8,0.15),
        axis.text.x=element_text(size=20),
        axis.title.x=element_text(size=20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=18))
ggsave(ah_plot,file="~/Desktop/ancestry_homozygosity_plot.jpeg")

hom_vs_het$set<-hom_vs_het$phen
pooled<-data.frame(ami=mean(hom_vs_het$ami),
                   stdev=sd(hom_vs_het$ami),
                   samples=mean(hom_vs_het$n))%>%
  mutate(CI.95=(1.96*stdev)/sqrt(samples))

ggplot(hom_vs_het,aes(y=ami))+
  geom_hline(aes(yintercept=0),linetype="dashed")+
  geom_hline(aes(yintercept=0.216,linetype="dashed"))+
  geom_boxplot(aes(x=phen))+
  geom_point(aes(x=phen))+
  coord_flip()+
  theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank())+
  labs(y="Assortative Mating Index")

#Do AMI test for genes of interest using shuffled random mating genomes
df<-tracts_input%>%unite(id_hap_chr,id_hap,chrom,remove=F)

shuffle<-data.frame(id_hap_chr=df$id_hap_chr)
shuffle.counts<-shuffle%>%table()%>%as.data.frame()
colnames(shuffle.counts)<-c("id_hap_chr","count")
shuffle<-left_join(unique(shuffle),shuffle.counts,by="id_hap_chr")
shuffle<-shuffle%>%separate(col=id_hap_chr,into=c("id","hap","chr"),sep="\\_",remove=F)%>%
  select(id_hap_chr,chr,count)
shuffle$id_hap_chr<-ave(shuffle$id_hap_chr,shuffle$chr,FUN = sample)
shuffle<-expandRows(shuffle,"count",count.is.col = TRUE,drop=T)
shuffle$s.id_hap_chr<-shuffle$id_hap_chr

df<-cbind(shuffle$s.id_hap_chr,df)%>%select(id_hap_chr="shuffle$s.id_hap_chr",chrom,begin,end,assignment,cmBegin,cmEnd)
df<-df%>%separate(col=id_hap_chr,into=c("id","hap","chr"),sep="_")%>%
  select(-chr)%>%
  unite(id_hap,id,hap,remove=T)

shuffled_gene_tracts_list<-list()
for (i in 1:nrow(matched)) {
  l<-filter(df,chrom==matched[i,]$chr & begin<matched[i,]$pos & end>matched[i,]$pos)
  l<-l%>%mutate(gene=matched[i,1],phen=matched[i,3])
  shuffled_gene_tracts_list[[i]]<-l
}
names(shuffled_gene_tracts_list)<-as.character(matched$gene)
shuffled_gene_tracts_list<-lapply(shuffled_gene_tracts_list,separate,col=id_hap,into=c("id","hap"),sep="\\_")
shuffled_gene_tracts_list<-lapply(shuffled_gene_tracts_list,select,id,hap,assignment,gene,phen)
shuffled_gene_tracts_list<-lapply(shuffled_gene_tracts_list,spread,hap,assignment)
shuffled_gene_tracts_list<-lapply(shuffled_gene_tracts_list,mutate,genotype=ifelse(A==B,A,
                                                                 ifelse(A!=B,'het',NA)))
shuffled_gene_hom_vs_het<-list()
for (i in 1:nrow(matched)) {
  gene<-names(shuffled_gene_tracts_list)[i]
  o_het<-shuffled_gene_tracts_list[[i]]%>%filter(genotype=="het")%>%nrow()
  n<-shuffled_gene_tracts_list[[i]]%>%nrow()
  o_hom_afr<-shuffled_gene_tracts_list[[i]]%>%filter(genotype=="African")%>%nrow()
  o_hom_eur<-shuffled_gene_tracts_list[[i]]%>%filter(genotype=="European")%>%nrow()
  df<-data.frame(gene=names(shuffled_gene_tracts_list)[i],phen=shuffled_gene_tracts_list[[i]][[3]][[1]],
                 o_hom_afr,o_hom_eur,n=n,o_het=o_het)
  shuffled_gene_hom_vs_het[[i]]<-df
}
shuffled_gene_hom_vs_het<-`row.names<-`(do.call(rbind,shuffled_gene_hom_vs_het), NULL)
shuffled_gene_hom_vs_het<-shuffled_gene_hom_vs_het%>%mutate(o_hom=n-o_het)

ancestral_fractions<-superpop_sum%>%select(anc=ancestral_superpop,fraction=percent)%>%
  mutate(fraction= fraction / 100)%>%
  spread(anc,fraction)
shuffled_gene_hom_vs_het<-cbind(shuffled_gene_hom_vs_het,ancestral_fractions)
shuffled_gene_hom_vs_het<-shuffled_gene_hom_vs_het%>%mutate(e_het=2*African*European*n,
                                e_hom=((African*African+European*European)*n))
shuffled_gene_hom_vs_het<-shuffled_gene_hom_vs_het%>%mutate(ami=log((o_hom/e_hom)/(o_het/e_het)))%>%
  arrange(ami)
shuffled_gene_hom_vs_het<-shuffled_gene_hom_vs_het%>%mutate(hom_afr_per=o_hom_afr/(o_hom_afr+o_hom_eur),
                                hom_eur_per=o_hom_eur/(o_hom_afr+o_hom_eur))
shuffled_gene_hom_vs_het$set<-shuffled_gene_hom_vs_het$phen
shuffled_gene_pooled<-data.frame(ami=mean(hom_vs_het$ami),
                   stdev=sd(hom_vs_het$ami),
                   samples=mean(hom_vs_het$n))%>%
  mutate(CI.95=(1.96*stdev)/sqrt(samples))


shuffled_skin.x<-shuffled_gene_hom_vs_het%>%filter(phen=="skin")%>%select(ami)
shuffled_hair.x<-shuffled_gene_hom_vs_het%>%filter(phen=="hair")%>%select(ami)
shuffled_eye.x<-shuffled_gene_hom_vs_het%>%filter(phen=="eye")%>%select(ami)

wilcox.test(y=shuffled_skin.x[,1],x=skin.x[,1],alternative="two.sided")
wilcox.test(y=shuffled_hair.x[,1],x=hair.x[,1],alternative="two.sided")
wilcox.test(y=shuffled_eye.x[,1],x=eye.x[,1])

df1<-shuffled_gene_hom_vs_het%>%mutate(svg="shuffle")
df2<-hom_vs_het%>%mutate(svg="gene")%>%select(-c(ah_afr,ah_eur,e_hom_afr,e_hom_eur))
genevshuff.df<-rbind(df1,df2)

gene_v_shuffle<-ggplot(data=genevshuff.df)+
  geom_vline(aes(xintercept=0))+
  geom_density(aes(x=ami,fill=svg,alpha=0.5))+
  geom_jitter(aes(x=ami,y=0.5,color=svg),width=0,height=0.1)+
  facet_grid(rows=vars(set),
             scales="free_y")+
  theme_classic()+
  scale_fill_manual(values=c("lightpink","lightblue"))+
  scale_color_manual(values=c(asw_col,afr_col))+
  scale_y_continuous(expand=c(0,0))+
  theme(legend.position = "none",
        axis.line.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.spacing = unit(1,"cm"),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_text(size=20),
        axis.text.x=element_text(size=20))+
  labs(x="Assortative Mating Index (AMI)")
##ggsave(gene_v_shuffle,file="~/Desktop/gene_vs_shuffle_plot.png")

genevshuff<-ggplot(NULL,group=set)+
  geom_vline(aes(xintercept = 0),size=0.5)+
  geom_jitter(data=shuffled_gene_hom_vs_het,aes(y=set,x=ami),width=0,height=0.1,colour=afr_col)+
  geom_jitter(data=hom_vs_het,aes(y=set,x=ami),width=0,height=0.1,colour=asw_col)+
  theme_classic()+
  facet_grid(~phen)+
  theme(axis.line.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.text=element_text(colour="black",size=18),
        axis.title.x=element_text(size=20))+
  scale_y_discrete(labels=c("Genome\nwide","Shuffled"))+
  labs(x="Assortative Mating Index (AMI)")

genevshuff.density<-ggplot(NULL,aes(alpha=0.5))+
  geom_vline(aes(xintercept = 0),size=0.5)+
  geom_density(data=shuffled_gene_hom_vs_het,aes(x=ami),fill="lightblue")+
  geom_density(data=hom_vs_het,aes(x=ami),fill="lightpink")+
  theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(colour="black",size=18),
        axis.text.y=element_blank(),
        axis.title.x=element_text(size=20),
        legend.position = "none")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  labs(x="Assortative Mating Index (AMI)")



#Do a permutation test for ancestry-based assortative mating (permutation of random sites to generate null distributions of gene set AMI values expected given the observed genome-wide levels of ancestry-based assortative mating)
#14 genes were used in total (2 for eye, 2 for hair, 2 for skin)
gene_sets_summary<-matched%>%select(phen)%>%table()%>%as.data.frame()
permutations_n2<-vector("numeric",1000)
for (i in 1:1000){
  set<-sample_n(all_hom_vs_het,2)
  ami<-mean(set$ami)
  permutations_n2[i]<-ami
}
mean_n2<-mean(permutations_n2)
permutations_n10<-vector("numeric",1000)
for (i in 1:1000){
  set<-sample_n(all_hom_vs_het,2)
  ami<-mean(set$ami)
  permutations_n10[i]<-ami
}
mean_n10<-mean(permutations_n10)

hair.perm<-data.frame(ami=permutations_n2,size=2,phen="hair")
eye.perm<-data.frame(ami=permutations_n2,size=2,phen="eye")
skin.perm<-data.frame(ami=permutations_n10,size=10,phen="skin")
perm.genesets<-rbind(hair.perm,eye.perm,skin.perm)

#Graph permuted versus calculated AMI for genesets
library(lattice)
ggplot(hom_vs_het,aes(x=ami))+
  geom_vline(aes(xintercept=0))+
  geom_density(fill="lightpink",alpha=0.5)+
  facet_grid(rows=vars(phen),scales="free_y")+
  geom_density(data=perm.genesets,aes(x=ami),fill="lightblue",alpha=0.5)+
  geom_jitter(aes(y=0.2),color=asw_col,height=0.1,width=0)+
  theme_classic()+
  theme(axis.title.y=element_blank(),
        axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=12))+
  scale_y_continuous(expand=c(0,0))

perm_hair.y<-perm.genesets%>%filter(phen=="hair")%>%select(ami)
perm_eye.y<-perm.genesets%>%filter(phen=="eye")%>%select(ami)
perm_skin.y<-perm.genesets%>%filter(phen=="skin")%>%select(ami)

wilcox.test(perm_hair.y[,1],hair.x[,1])
wilcox.test(perm_eye.y[,1],hair.x[,1])
wilcox.test(perm_skin.y[,1],hair.x[,1])


#Test of the significance of AMIs of the ancestry cue genes
detach(package:plyr)
cue.geneset<-hom_vs_het%>%group_by(phen)%>%summarize(stdev=sd(ami),
                                                     ami=mean(ami),
                                             samples=mean(n),
                                             hom_afr=mean(o_hom_afr),
                                             hom_eur=mean(o_hom_eur),
                                             African=mean(African),
                                             European=mean(European))%>%
  rename(set="phen")

all.geneset<-all_hom_vs_het%>%summarise(stdev=sd(ami),
                                        ami=mean(ami),
                                        samples=mean(n),
                                        hom_afr=mean(o_hom_afr),
                                        hom_eur=mean(o_hom_eur),
                                        African=mean(African),
                                        European=mean(European))%>%
  mutate(set="all")

shuffled.geneset<-shuffled_hom_vs_het%>%summarise(stdev=sd(ami,na.rm=T),
                                                  ami=mean(ami),
                                                  samples=mean(n),
                                                  hom_afr=mean(o_hom_afr),
                                                  hom_eur=mean(o_hom_eur),
                                                  African=mean(African),
                                                  European=mean(European))%>%
  mutate(set="shuffle")
set<-rbind(cue.geneset,all.geneset,shuffled.geneset)%>%
  mutate(CI=qnorm(0.95)*(stdev/sqrt(samples)),
         lower=ami-CI,
         upper=ami+CI)

#Statistical testing

amis<-data.frame(set=hom_vs_het$phen,
           ami=hom_vs_het$ami)
amis.all<-all_hom_vs_het%>%select(ami)%>%mutate(set='all')
amis<-rbind(amis.all,amis)
amis$set<-factor(amis$set)

TukeyHSD(aov(ami~set,amis))


genes.x<-hom_vs_het%>%select(ami)
skin.x<-hom_vs_het%>%filter(phen=='skin')%>%select(ami)
hair.x<-hom_vs_het%>%filter(phen=='hair')%>%select(ami)
eye.x<-hom_vs_het%>%filter(phen=='eye')%>%select(ami)
all.y<-all_hom_vs_het$ami

wilcox.test(genes.x[,1],all.y) 
genes.p<-wilcox.test(genes.x[,1],all.y)$p.value
#Welch Two Sample t-test

#data:  genes.x and all.p
#t = 0.59951, df = 13.746, p-value = 0.5586
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.1054133  0.1870101
#sample estimates:
#  mean of x mean of y 
#0.2563503 0.2155518 

t.test(skin.x,all.y)
skin.p<-t.test(skin.x,all.y)$p.value
skin.t<-t.test(skin.x,all.y)$statistic
#Welch Two Sample t-test

#data:  skin.x and all.p
#t = 0.36888, df = 9.2557, p-value = 0.7205
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.1803494  0.2509800
#sample estimates:
#  mean of x mean of y 
#0.2508671 0.2155518 

t.test(hair.x,all.y)
hair.p<-t.test(hair.x,all.y)$p.value
hair.t<-t.test(hair.x,all.y)$statistic
#Welch Two Sample t-test

#data:  hair.x and all.p
#t = 1.7484, df = 1.1532, p-value = 0.3051
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.3267044  0.4772657
#sample estimates:
#  mean of x mean of y 
#0.2908325 0.2155518 

t.test(eye.x,all.y)
eye.p<-t.test(eye.x,all.y)$p.value
eye.t<-t.test(eye.x,all.y)$statistic
#Welch Two Sample t-test
#
#data:  eye.x and all.p
#t = 2.9871, df = 579, p-value = 0.002936
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.01155254 0.05591128
#sample estimates:
#  mean of x mean of y 
#0.2492838 0.2155518 

all.y<-all_hom_vs_het$ami
shuffled.y<-shuffled_hom_vs_het$ami
t.test(all.y,shuffled.y)
all.p<-t.test(all.y,shuffled.y)$p.value
all.t<-t.test(all.y,shuffled.y)$statistic


set<-set%>%mutate(p=c(skin.p,hair.p,eye.p,NA,NA),
                  t=c(skin.t,hair.t,eye.t,NA,NA))%>%
  as.data.frame()
set<-left_join(set,gene_sets_summary,by=c("set"="."))
set<-set%>%rename(n=Freq)%>%
  mutate(labels=c("Skin","Hair","Eyes",NA,NA))
set$order<-c(1:5)
set$n<-c(10,2,2,NA,NA)
set$set<-factor(set$set, levels = set$set[order(set$order)])


#Plotting
hom_vs_het$set<-hom_vs_het$phen
hom_vs_het$order<-ifelse(hom_vs_het$set=="skin",1,ifelse(hom_vs_het$set=="hair",2,3))
all_hom_vs_het$set<-factor("all")
all_hom_vs_het$order<-factor(4)
shuffled_hom_vs_het$set<-factor("shuffle")
shuffled_hom_vs_het$order<-factor(5)
set$set <- factor(set$set, as.character(set$set))

r.hom_vs_het<-rbind(select(hom_vs_het,ami,set,order),
      select(all_hom_vs_het,ami,set,order),
      select(shuffled_hom_vs_het,ami,set,order))
r.hom_vs_het$set<-factor(r.hom_vs_het$set)     


ami_plot<-ggplot(NULL,aes(x=ami,y=reorder(set,desc(order))),group=set)+
  geom_vline(aes(xintercept = 0,linetype="zero"),size=1.2)+
  geom_vline(data=all_hom_vs_het,aes(xintercept=mean(ami),linetype="mean"),alpha=0.5)+
  geom_vline(data=shuffled_hom_vs_het,aes(xintercept=mean(ami),linetype="mean"),alpha=0.2)+
  geom_jitter(data=r.hom_vs_het,height=0.1,aes(alpha=0.8),color=afr_col)+
  geom_point(data=set,aes(x=ami),size=7)+
  geom_errorbar(data=subset(set,set!='all'&set!='shuffle'),aes(xmin=lower,xmax=upper,y=set),width=0.4,size=0.5)+
  geom_text(data=set,size=6,aes(x=0.6,label=ifelse(is.na(p),"",paste0("P = ",round(p,3),", N = ",n))),hjust="left")+
  theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.text=element_text(colour="black",size=18),
        axis.title.x=element_text(size=20),
        legend.position = "none")+
  scale_y_discrete(labels=c("Shuffled","Genome\nwide","Eye\ncolour","Hair","Skin\npigment"))+
  labs(x="Assortative Mating Index (AMI)")+
  scale_linetype_manual(values=c(zero="solid",mean="dashed"))  
##ggsave(ami_plot,file="~/Desktop/ami_plot.png")


allvshuff<-ggplot(NULL)+
  geom_vline(aes(xintercept = 0),size=0.5)+
  geom_jitter(data=shuffled_hom_vs_het,aes(y=set,x=ami),width=0.035,height=0.3,colour=afr_col)+
  geom_jitter(data=all_hom_vs_het,aes(y=set,x=ami),width=0.035,height=0.3,colour=asw_col)+
  theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.text=element_text(colour="black",size=18),
        axis.title.x=element_text(size=20),
        legend.position = "none")+
  scale_y_discrete(labels=c("Genome\nwide","Shuffled"))+
  labs(x="Assortative Mating Index (AMI)")

allvshuff.density<-ggplot(NULL,aes(alpha=0.5))+
  geom_vline(aes(xintercept = 0),size=0.5)+
  geom_density(data=shuffled_hom_vs_het,aes(x=ami),fill="lightblue")+
  geom_density(data=all_hom_vs_het,aes(x=ami),fill="lightpink")+
  theme(axis.line.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.text=element_text(colour="black",size=18),
        axis.title.x=element_text(size=20),
        legend.position = "none")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  labs(x="Assortative Mating Index (AMI)")

grid.arrange(allvshuff,allvshuff.density)
gen_wide_AMI<-arrangeGrob(allvshuff,allvshuff.density)
##ggsave(gen_wide_AMI,file="~/Desktop/AMI_permutation.png")

ami_dist_plot<-ggplot(NULL,aes(x=ami,y=reorder(factor(set),order),group=set))+
  geom_vline(aes(xintercept = 0,linetype="zero"),size=1.2)+
  geom_vline(data=all_hom_vs_het,aes(xintercept=mean(ami),linetype="mean"))+
  geom_point(data=subset(set,set!='all'&set!='shuffle'),aes(x=ami),size=7)+
  geom_errorbar(data=subset(set,set!='all'),aes(xmin=lower,xmax=upper),width=0.4,size=0.5)+
  geom_jitter(data=hom_vs_het,aes(x=ami),col="darkgrey",width=0,height=0.1)+
  geom_text(data=set,size=6,aes(x=0.6,label=ifelse(is.na(p),"",paste0("P = ",round(p,3),", N = ",n))),hjust="left")+
  geom_jitter(data=shuffled_hom_vs_het,aes(x=ami,y=reorder(set,order)),width=0.4,size=0.5)+
  geom_jitter(data=all_hom_vs_het,aes(x=ami),width=0.035,height=0.3)+
  theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.text=element_text(colour="black",size=18),
        axis.title.x=element_text(size=20),
        legend.position = "none")+
  scale_y_discrete(labels=c("Shuffled","Genome\nwide","Eye\ncolour","Hair","Skin\npigment"))+
  labs(x="Assortative Mating Index (AMI)")+
  scale_linetype_manual(values=c(zero="solid",mean="dashed"))
ami_dist_plot
##ggsave(ami_dist_plot,file="~/Desktop/ami_sets.png")


anc_homozyg<-set%>%gather(anc,hom_n,hom_afr,hom_eur)%>%
  select(set,anc,hom_n)
anc_homozyg_plot<-ggplot(anc_homozyg,aes(x=set,y=hom_n,fill=anc))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
  geom_bar(stat='identity',position='dodge')+
  scale_fill_manual(values=c(afr_col,eur_col))+
  theme_classic()+
  labs(y="Average ancestral homozygosity counts within each set")+
  scale_y_continuous(expand=c(0,0))+
  coord_flip()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank())

grid.arrange(ami_dist_plot,anc_homozyg_plot,nrow=1)

#Chromopainter input file write
setwd('chromopainter_in/individuals/')

options(scipen=999)

physical.to.genetic<-rfmixout%>%
  select(chromosome,physical_position,genetic_position)%>%
  unite(chrom.pos,chromosome,physical_position,sep=".",remove=F) #positions with genetic position info

#Haplotype_infile
files<-list.files()

individuals<-list()
haplotypes<-list()
for (i in 1:length(files)) {
  file<-files[i]
  ind<-gsub(file,pattern='.txt$',replacement="")
  df<-read.table(file=file,header=F,sep='\t')
  colnames(df)<-c('rsid','chr','pos','haps')
  df<-df%>%unite(chr.pos,chr,pos,sep='.',remove=F)
  df<-df%>%separate(col=haps,into=c(paste0(ind,"_A"),paste0(ind,"_B")))
  df<-df%>%filter(chr.pos%in%physical.to.genetic$chrom.pos)
  individuals[[i]]<-df
  names(individuals)[i]<-paste(ind)
  haplotypes[[i]]<-df[,c(5,6)]
}
haplotypes.bound<-list.cbind(haplotypes)
all_na<-function(x) any(!is.na(x))
haplotypes.bound<-haplotypes.bound%>%select_if(all_na)

line1<-ncol(haplotypes.bound)
line3<-c('P',individuals[[1]][['pos']])%>%as.character()
line2<-length(line3)
haplotypes<-t(haplotypes.bound)%>%as.data.frame()
haplotypes<-haplotypes%>%unite(col='hap',sep="")
haplotypes<-haplotypes[,1]
haplotype_infile<-c(line1,
                    line2,
                    paste(line3,collapse=" "),
                    haplotypes)

haplotype_infile<-list(line1,line2,line3,haplotypes)
lapply(haplotype_infile,cat,file='../haplotype_infile',
       sep=" ",append=T,row.names=F,quote=F)

#Recom_rate_infile
recom<-data.frame(start.pos=individuals[[1]][['pos']],chrom=individuals[[1]][['chr']])%>%
  unite(chrom.pos,chrom,start.pos,sep=".",remove=F)
physical.to.genetic.match<-left_join(recom,physical.to.genetic,by="chrom.pos")%>%
  select(start.pos,physical_position,genetic_position)
physical.to.genetic.match<-physical.to.genetic.match%>%
  mutate(recom.rate.perbp=ifelse(physical_position>lead(physical_position),-9,
                                 genetic_position/(lead(physical_position)-physical_position)))
recom_rate_infile<-physical.to.genetic.match%>%
  select(start.pos,recom.rate.perbp)
write.table(recom_rate_infile,file="../recom_rate_infile",sep=" ",quote=F,row.names=F)

#label_infie
sample_info<-read.table('../../20130606_g1k.ped',sep='\t',header=T)%>%
  select(id=Individual.ID,Population)
samples<-data.frame(id.hap=names(haplotypes.bound))%>%
  separate(col=id.hap,into=c('id','hap'))%>%
  select(id)%>%unique()
label_infile<-left_join(samples,sample_info,by="id")%>%
  mutate(include=1)
write.table(label_infile,file="../label_infile",sep=" ",quote=F,row.names=F,col.names=F)

#population_list_infile
population_list_infile<-data.frame(Population=label_infile$Population)%>%
  mutate(d_r=ifelse(Population%in%European,'D',
                    ifelse(Population%in%African,'D',
                           ifelse(Population%in%ASW,'R',NA))))%>%
  unique()
write.table(population_list_infile,file="../population_list_infile",sep=" ",quote=F,row.names=F,col.names=F)





####GRAVEYARD


