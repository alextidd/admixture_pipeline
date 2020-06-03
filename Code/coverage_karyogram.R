#!/usr/bin/env Rscript
library(ggplot2,warn.conflicts = F)
library(dplyr,warn.conflicts = F)
library(tidyr,warn.conflicts = F)
library(lattice,warn.conflicts = F)
library(data.table,warn.conflicts = F)
library(stringr,warn.conflicts = F)
library(chromoMap,warn.conflicts = F)
library(plotly,warn.conflicts = F)
library(webshot,warn.conflicts = F)
numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
} 

#Read in the position information and generate the annotation file
chromosomes<-dir(path='./coverage/',pattern="^chr")%>%numextract()
chromosomes<-gsub("\\.","",chromosomes)%>%as.numeric()
for (chr in chromosomes) {
  file<-paste0('./coverage/chr',chr,'.txt')
  df.name<-paste('chr',chr,sep="")
  assign(df.name,read.table(file=file)%>%
           mutate('chr'=chr,start=V1,stop=V1)%>%
           mutate(POS=paste('chr',chr,'.',V1,sep=""))%>%
           select(POS,chr,start,stop))
}
df.names <- ls(pattern='^chr\\d+')
coverage<-'row.names<-'(do.call(rbind,mget(df.names)),NULL)
write.table(coverage,file='./coverage/coverage.txt',col.names = F,row.names=F,sep='\t')

print(paste('The number of AIM sites matched across the genome is ',nrow(coverage),sep=""))

#Create a chromoMap point annotation karyogram of covered sites
plot<-chromoMap('chromomap_chromosome_coords.txt','./coverage/coverage.txt',
          chr_color='#B4B4B4',
          anno_col = '#95FE67',
          segment_annotation = T,
          canvas_width=1400)
#Green means AIM within 2.5*10^6 bp
plot

#Save the widget as .html format
htmlwidgets::saveWidget(as_widget(plot),'../Graphs/coverage_plot.html',selfcontained = F)
#Save as non-interactive .png
webshot::install_phantomjs()
webshot("../Graphs/coverage_plot.html" , "../Graphs/coverage_plot.pdf", delay = 0.2)


