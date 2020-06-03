setwd('chromopainter_in/individuals/')
library(rlist)

options(scipen = 999)

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
  individuals[[i]]<-df
  names(individuals)[i]<-paste(ind)
  haplotypes[[i]]<-df[,c(5,6)]
}
haplotypes.bound<-list.cbind(haplotypes)

all_na<-function(x) any(!is.na(x))
haplotypes.bound<-haplotypes.bound%>%select_if(all_na)

line1<-ncol(haplotypes.bound)
line3<-c('P',individuals[[1]][['pos']])
line2<-length(line3)
haplotypes<-t(haplotypes.bound)

haplotype_infile<-list(line1,line2,line3,haplotypes)
lapply(haplotype_infile,cat,"\n",file='../haplotype_infile.txt',append=T)  


#Recom_rate_infile
recom<-data.frame(start.pos=individuals[[1]][['pos']])

