#!/usr/bin/env Rscript
library(dplyr,warn.conflicts=F)
library(tidyr,warn.conflicts=F)

populations<-c('ASW','GBR','CEU','TSI','IBS','GWD','ESN','MSL','YRI','LWK')
ASW<-c('ASW')
European<-c('GBR','CEU','TSI','IBS')
African<-c('GWD','ESN','MSL','YRI','LWK')

sample_info<-read.csv("20130606_g1k.ped",sep='\t',header=T)%>%select(Individual.ID,Population)
write.table(sample_info,file='pca/pca.clst.txt')

sample_info<-data.frame(lapply(sample_info,function(x) {
  gsub("ASW","-",x)
}))

fam<-read.table('./merged_plink/Merge.fam')
inds<-fam%>%select(Individual.ID=V2)
pop_file<-left_join(inds,sample_info,by='Individual.ID')
pop_file.K9<-pop_file%>%select(Population)

write.table(pop_file.K9,file='./merged_plink/Merge9.pop',quote=F,row.names=F,col.names=F)


pop_file.K2<-pop_file.K9%>%mutate(superpop=ifelse(Population%in%African,'AFR',
                                                  ifelse(Population%in%European,'EUR',"-")))%>%
  select(superpop)

write.table(pop_file.K2,file='./merged_plink/Merge2.pop',quote=F,row.names=F,col.names=F)



