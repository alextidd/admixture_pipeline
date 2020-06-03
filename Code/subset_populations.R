#!/usr/bin/env Rscript
library(dplyr,warn.conflicts=F)

args<-commandArgs(trailingOnly=TRUE)
number_of_individuals<-tail(args,n=1)
number_of_individuals<-number_of_individuals[1]
populations<-head(args,-1)

print(populations)
sample_info<-read.table("20130606_g1k.ped",sep="\t",header=T)
sample_info_subset<-sample_info%>%subset(Population%in%populations)

print(number_of_individuals)

sample_info_subset<-if(number_of_individuals!='ALL'){
  sample_info_subset%>%
    group_by(Population)%>%
    sample_n(size = as.numeric(number_of_individuals))%>%
    as.data.frame()
} else {
  sample_info_subset  
  }

paste0("The number of sample IDs in the subset is ",nrow(sample_info_subset))%>%print()
ind_pop_panel<-sample_info_subset%>%select(Individual.ID)

pops_prefix<-populations%>%paste(collapse='_',sep='')

write.table(ind_pop_panel,
            file = paste0(pops_prefix,"_ids.txt"),
            sep = "\t",
            row.names = F,
            col.names=F,quote=F)

