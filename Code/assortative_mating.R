setwd('~/job_env/Data')
library(dplyr)

sample_info<-read.table('20130606_g1k.ped',sep='\t',header=T)
paired_ASW<-sample_info%>%filter('ASW'==Population,c('father','mother')==Relationship)%>%arrange(Family.ID)

#Alleles (pigmentation - polygenic phenotype) and their significant trait-associated SNPs.  
MC1R<-data.frame(gene="MC1R",rsid=c("rs1805007","rs1805005","rs2228479"))
OC2<-data.frame(gene="OC2",rsid=c("rs12913832"))
TYR<-data.frame(gene="TYR",rsid=c("rs1042602","rs1126809"))
TYRP1<-data.frame(gene="TYRP1",rsid=c("rs1408799"))
DCT<-data.frame(gene="DCT",rsid=c("rs1407995","rs2031526"))
SLC45A2<-data.frame(gene="SLC45A2",rsid=c("rs16891982","rs26722"))
SLC24A5<-data.frame(gene="SLC24A5",rsid=c("rs1426654"))
KITLG<-data.frame(gene="KITLG",rsid=c("rs642742","rs12821256"))
IRF4<-data.frame(gene="IRF4",rsid=c("rs12203592"))
geneset<-rbind(MC1R,OC2,TYR,TYRP1,DCT,SLC45A2,SLC24A5,KITLG,IRF4)

positions<-read.csv('assortative_mating_analysis/mart_export.txt')
positions<-positions%>%select(rsid=Variant.name,chromosome=Chromosome.scaffold.name,pos=Chromosome.scaffold.position.start..bp.)

#Need to filter to genes that pass the ancestry genotype threshold and to prune for LD.
matched<-read.table('assortative_mating_analysis/matched.txt')
colnames(matched)<-c('chr','rsid','eur_af','afr_af','ref','alt')
matched<-left_join(matched,geneset,by='rsid')



