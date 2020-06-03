#!/usr/bin/env Rscript
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
library(readxl)
read_excel_allsheets<-function(file,tibble=F){
  sheets <- readxl::excel_sheets(file)
  x <- lapply(sheets, function(X) readxl::read_excel(file, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

#STUDY 1 - Load in panel excel file from Reich et al., 2015 as a list and then rbind
file<-"reich2015aims.xls"
sheets<-read_excel_allsheets('./AIMs/reich2015aims.xls')
aims_rsIDs1<-ldply(sheets,rbind)%>%select("SNP_ID")

#STUDY 2 - Load in panel csv from Tandon, Patterson & Reich, 2015
aims_rsIDs2<-read.csv('./AIMs/reich2007aims.csv')%>%select("SNP_ID")

#Combine together all AIMs from both studies
aims_rsIDs<-rbind(aims_rsIDs1,aims_rsIDs2)%>%unique()

#Write combined panel file
write.table(aims_rsIDs,
            file='./AIMs/aims_rsIDs.txt',
            col.names=T,row.names=F,
            quote=F)
