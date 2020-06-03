#Graphics for poster
library(ggplot2)
library(dplyr)
library(tidyr)
library(lattice)
library(data.table)
setwd('~/job_env/Data')

rfmixout<-read.table(file='RFMIX_output.fb.tsv',header=T)
#hist(rfmixout$physical_position,breaks=1000)
rfmixout$physical_position_old<-rfmixout$physical_position
rfmixout$physical_position<-(rfmixout$physical_position)/50818468

moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", 
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}

rfmixout<-rfmixout[moveme(names(rfmixout),"physical_position_old first")]

long<-melt(setDT(rfmixout),id.vars=colnames(rfmixout)[1:5],variable.name="INDIVIDUAL") %>%
  separate(col='INDIVIDUAL',into=c('ID','hap','ancestral_pop'),remove=F,extra='warn')
long<-long%>%rename('ancestry_proportion'=value)
long<-long%>%unite(ID_pos_pop,ID,physical_position,ancestral_pop,sep="_",remove=F)

head(long)
long.test.plot.list<-long$ID%>%
  as.data.frame()%>%
  distinct()%>%
  top_n(1)
colnames(long.test.plot.list)<-'ID'
long.test.plot<-right_join(x=long,y=long.test.plot.list,by='ID')
long.test.plot%>%distinct(ID)%>%nrow()

p<-ggplot(long.test.plot,aes(x=physical_position,y=ancestry_proportion)) +
  geom_area(aes(col=ancestral_pop))+
  facet_grid(rows=vars(hap))
p
  
long.dt<-long%>%
  select(ID_pos_pop,ancestry_proportion)%>%
  data.table()
long.dt<-long.dt[, lapply(.SD,sum), by = ID_pos_pop]

long.diploid<-long.dt%>%
  separate(col='ID_pos_pop',into=c('ID','physical_position','population'),sep="_",remove=F,extra='warn')%>%
  as.data.frame()

long.diploid.small.list<-long.diploid$ID%>%
  as.data.frame()%>%
  distinct()%>%
  top_n(1)
colnames(long.diploid.small.list)<-'ID'
long.diploid.small<-right_join(x=long.diploid,y=long.diploid.small.list,by='ID')


p<-ggplot(long.diploid.small,aes(x=physical_position,y=ancestry_proportion,group=population)) +
  geom_line(aes(color=population))+
  facet_grid(rows=vars(ID))
p

p<-ggplot(long.diploid.small,aes(x=physical_position,y=ancestry_proportion)) +
  geom_histogram(aes(fill=population))+
  facet_grid(rows=vars(ID))
p

p<-ggplot(rfmixout_NA20294,
          aes(x=physical_position,y=CEU))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0),breaks=seq(0,2,by=1)) +
  theme_classic()+
  labs(title="Individual NA20294",x=expression(Position~on~Chromosome~22~(x10^7)),y="Ancestry proportions")+
  theme(panel.border=element_blank(),
        plot.title=element_text(hjust=0.5,size=25,face="bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=22))+
  geom_area(fill="#f0b387")+
  geom_ribbon(aes(ymin=CEU,ymax=2),fill="#8facda")
setwd('../Graphs/')
png(filename='./plot_NA20294.png',width=900,height=300)
print(p)
dev.off()

#NA19835
rfmixout_NA19835<-rfmixout%>%
  select(physical_position,
         NA19835...hap1...CEU,
         NA19835...hap1...YRI,
         NA19835...hap2...CEU,
         NA19835...hap2...YRI)%>%
  rename(CEU1=NA19835...hap1...CEU,
         YRI1=NA19835...hap1...YRI,
         CEU2=NA19835...hap2...CEU,
         YRI2=NA19835...hap2...YRI)

rfmixout_NA19835$CEU<-rowSums(rfmixout_NA19835[,c("CEU1","CEU2")])
rfmixout_NA19835$YRI<-rowSums(rfmixout_NA19835[,c("YRI1","YRI2")])

p<-ggplot(rfmixout_NA19835,
          aes(x=physical_position,y=CEU,fill=))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0),breaks=seq(0,2,by=1)) +
  theme_classic()+
  labs(title="Individual NA19835",x=expression(Position~on~Chromosome~22~(x10^7)),y="Ancestry proportions")+
  theme(panel.border=element_blank(),
        plot.title=element_text(hjust=0.5,size=25,face="bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=22))+
  geom_area(fill="#f0b387")+
  geom_ribbon(aes(ymin=CEU,ymax=2),fill="#8facda")
setwd('~/OneDrive/FYP/Graphs/')
png(filename='./plot_NA19835.png',width=900,height=300)
print(p)
dev.off()

#NA19904
rfmixout_NA19904<-rfmixout%>%
  select(physical_position,
         NA19904...hap1...CEU,
         NA19904...hap1...YRI,
         NA19904...hap2...CEU,
         NA19904...hap2...YRI)%>%
  rename(CEU1=NA19904...hap1...CEU,
         YRI1=NA19904...hap1...YRI,
         CEU2=NA19904...hap2...CEU,
         YRI2=NA19904...hap2...YRI)

rfmixout_NA19904$CEU<-rowSums(rfmixout_NA19904[,c("CEU1","CEU2")])
rfmixout_NA19904$YRI<-rowSums(rfmixout_NA19904[,c("YRI1","YRI2")])

p<-ggplot(rfmixout_NA19904,
          aes(x=physical_position,y=CEU,fill=))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0),breaks=seq(0,2,by=1)) +
  theme_classic()+
  labs(title="Individual NA19904",x=expression(Position~on~Chromosome~22~(x10^7)),y="Ancestry proportions")+
  theme(panel.border=element_blank(),
        plot.title=element_text(hjust=0.5,size=25,face="bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=22))+
  geom_area(fill="#f0b387")+
  geom_ribbon(aes(ymin=CEU,ymax=2),fill="#8facda")
setwd('~/OneDrive/FYP/Graphs/')
png(filename='./plot_NA19904.png',width=900,height=300)
print(p)
dev.off()



NA19835<-gather(data=rfmixout,
       key='source',
       value=NA19835,CEU:YRI,
       factor_key=TRUE)
View(NA19835)
