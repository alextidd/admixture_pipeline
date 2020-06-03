#!/usr/bin/env Rscript
ref_populations<-commandArgs(trailingOnly=TRUE)

print(ref_populations)
sample_info<-read.csv("20130606_sample_info.csv")
library(dplyr,warn.conflicts=F)

ref_sample_info_subset<-sample_info%>%subset(Population%in%ref_populations)
ref_ind_pop_panel<-ref_sample_info_subset%>%select(Sample,Population)

pops_prefix<-ref_populations%>%paste(collapse='_',sep='')

write.table(ref_ind_pop_panel,
            file = paste0("filtered_",pops_prefix,"_reference_sample_map.txt"),
            sep = "\t",
            row.names = F,
            col.names=F,quote=F)
