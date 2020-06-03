#!/usr/bin/env Rscript
library(dplyr)

filelist<-c(sprintf("%s_reference.sample",seq(from=1,to=22)),
            sprintf("%s_query.sample",seq(from=1,to=22)))
datalist<-lapply(filelist, function(x)read.table(x, header=T))
names(datalist)<-filelist
datalist<-lapply(datalist,transform,father=0,mother=0,plink_pheno=0)
columns<-c('ID_1','ID_2','missing','father','mother','sex','plink_pheno')
reordered<-lapply(datalist,'[', ,columns)

lapply(names(reordered),function(df) write.table(reordered[[df]],file=paste0(df),quote=F,row.names=F,col.names=T))

