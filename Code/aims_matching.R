#!/usr/bin/env Rscript
suppressMessages(library(plyr))
suppressMessages(library(dplyr))

#Reading in the file names for the rsIDs of interest and the chromosomes (within a for loop)
args<-commandArgs(trailingOnly=TRUE)
data_rsIDs_filename<-args[1]
chr<-args[2]
data_rsIDs<-read.table(data_rsIDs_filename,header=T)%>%as.data.frame()
colnames(data_rsIDs)<-"SNP_ID"

#Reading in the AIM panel to be used for filtering sites
aims_rsIDs<-read.table('./AIMs/aims_rsIDs.txt',header=T)

#Joining to find rsIDs shared by the data and the AIMs list
matched_rsIDs<-inner_join(data_rsIDs,aims_rsIDs,by="SNP_ID")
write.table(matched_rsIDs,
            file=paste("./AIMs/chr",chr,"_all_matched_aims.txt",sep=""),
            row.names=F,col.names=F,
            quote=F,sep='\t')
