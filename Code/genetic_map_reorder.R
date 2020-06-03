#!/usr/bin/env Rscript
setwd('./1000_Genomes_Phase3')
library(dplyr)
library(plyr)
library(stringr)

filelist<-list.files(pattern=".*.txt")
datalist<-lapply(filelist, function(x)read.table(x, header=T))
names(datalist)<-list.files(pattern=".*.txt")
remove <- c("genetic_map_chr", "_combined_b37.txt")
chrs<-str_remove_all(filelist, paste(remove, collapse = "|"))
datalist<-mapply(cbind,datalist,'chr'=chrs,SIMPLIFY=F)
columns<-c('chr','position','Genetic_Map.cM.')
result<-lapply(datalist,'[', ,columns)

dir.create('./reordered_rfmixv2')
lapply(names(result),
       function(x,result) write.table(result[[x]],paste("./reordered_rfmixv2/",x, sep = ""),
                                      col.names=F,row.names=F,sep='\t',quote=F),
       result)