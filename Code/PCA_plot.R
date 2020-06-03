#!/usr/bin/env Rscript
library(dplyr,warn.conflicts = F)
library(R.devices,warn.conflicts=F)
library(PCAviz)
library(cowplot)
library(ggplot2)

library(unikn)
eur_col<-"#FBD154"
asw_col<-"#CC79A7"
afr_col<-"#2776BB"
pal_pops<-newpal(col=c(eur_col,asw_col,afr_col),names=c("EUR","ASW","AFR"))
pal_superpops<-seecol(pal_pops,n=9)%>%newpal(names=c("GBR","CEU","TSI","IBS","MSL","ESN","LWK","GWD","YRI"))

#Read in the pedigree file
ped<-read.table('20130606_g1k.ped',sep='\t',header=T)%>%select(Family.ID,Individual.ID,Population)

##Read in individual coordinates on PCs and eigenvalues
PCA<-read.table('pca/pca.eigs')
nPCs<-20
names(PCA)<-c("ID",paste("PC",(1:nPCs),sep=""),"case.control")
PCA<-PCA%>%select(-case.control)
eig.val<-sqrt(unlist(read.table('pca/pca.eval'))[1:nPCs])
sum.eig<-sum(unlist(read.table('pca/pca.eval')))

#Read in SNP weightings matrix
snpeigs<-read.table(paste('pca/pca.snpeigs'))
names(snpeigs)<-c("ID","chr","pos",paste("PC",(1:nPCs),sep=""))
snpeigs$chr<-factor(snpeigs$chr)
rownames(snpeigs)<-snpeigs$ID
snpeigs<-snpeigs[,-1]

#Smartpca pushes the plink family and individual IDs together so we need to extract out the IDs
tmp<-unlist(sapply(as.character(PCA$ID),strsplit,":"))
ids<-tmp[seq(2,length(tmp),by=2)]
PCA$ID<-ids

#Read in the group/cluster labels
clst<-ped
clst_unord<-PCA%>%select(Individual.ID=ID)
clst_unord<-left_join(clst_unord,clst,by='Individual.ID')
clst_unord<-clst_unord%>%select(Population)
PCA<-cbind(PCA,clst_unord)
names(PCA)[ncol(PCA)]<-'Population'

#Create superpop annotations
ASW<-c('ASW')
European<-c('GBR','CEU','TSI','IBS')
African<-c('GWD','ESN','MSL','YRI','LWK')
PCA<-PCA%>%mutate(superpop=ifelse(Population%in%African,'African',
                                              ifelse(Population%in%European,'European',
                                                      ifelse(Population%in%ASW,'ASW',NA))))

#Calculate median positions for ancestral populations
median1<-PCA%>%group_by(superpop)%>%summarise(median(PC1))
corrected<-median1[1,2]-median1[3,2]


#Build the PCAviz object
hgdp<-pcaviz(dat=PCA,sdev=eig.val,
             var=sum.eig,rotation=snpeigs)
hgdp<-pcaviz_abbreviate_var(hgdp,'Population')

#Make PCA plots - grid
geom.point.summary.params<-list(
  shape=16,stroke=1,size=5,
  show.legend=F)
plot1<-plot(hgdp,
            coords=paste0("PC",c(1,2)),color='Population',
            geom.point.summary.params=geom.point.summary.params,
            scale.pc.axes=0.6)
plot2<-plot(hgdp,
            coords=paste0("PC",c(2,3)),color='Population',
            geom.point.summary.params=geom.point.summary.params,
            scale.pc.axes=1)
plot3<-plot(hgdp,
            coords=paste0("PC",c(3,4)),color='Population',
            geom.point.summary.params=geom.point.summary.params,
            scale.pc.axes=0.6)
plot4<-plot(hgdp,
            coords=paste0("PC",c(4,5)),color='Population',
            geom.point.summary.params=geom.point.summary.params,
            scale.pc.axes=0.6)
plot5<-plot(hgdp,
            coords=paste0("PC",c(5,6)),color='Population',
            geom.point.summary.params=geom.point.summary.params,
            scale.pc.axes=0.6)

pca_grid<-plot_grid(plot1+theme(legend.position="none"),
                    plot2+theme(legend.position="none"),
                    plot3+theme(legend.position="none"),
                    plot4+theme(legend.position="none"),
                    nrow=2,align='vh',axis="tlr")
pca_grid
plot_legend<-get_legend(plot1)
ggsave(pca_grid,file='~/Desktop/pca_grid.png')
ggsave(plot_legend,file='~/Desktop/pca_grid_legend.png')

#Make PCA plots - violin
violin<-pcaviz_violin(hgdp,pc.dim=paste0('PC',c(1:3)),
              plot.grid.params=list(nrow=3))
violin



ggsave(violin,file="~/Desktop/pca_violin.png")

#Make PCA plots - loadings
for (i in 1:5) {
  plotname <- paste0("plot",i)
  plot<-pcaviz_loadingsplot(hgdp,
                            pc.dim=paste0("PC",i),
                            min.rank=0.8,gap=200,color="chr",
                            geom.point.params=list(show.legend=F)
                            )+
    xlab("SNPs") +
    ylab(paste0("PC",i," loading"))
  assign(plotname,plot)
}
#Grep common legend
plot<-pcaviz_loadingsplot(hgdp,
                          pc.dim=paste0("PC",i),
                          min.rank=0.8,gap=200,color="chr"
                          )+
  guides(color=guide_legend(nrow=1,byrow=T))+
  theme(legend.position="bottom",
        legend.justification="center",
        legend.spacing.x=unit(0.103,"cm"),
        legend.title=element_blank())
plot_legend<-get_legend(plot)
#Plot loadings
prow<-plot_grid(plot1,
                plot2,
                plot3,
                plot4,
                plot5,
                NULL,
                ncol=1,align='vh',
                rel_heights=c(1,1,1,1,1,0.3))

pca_loadings<-prow+draw_grob(plot_legend,x=0.525,y=-0.46,width=0,height=1)

ggsave(pca_loadings,file="~/Desktop/pca_loadings.png")


#PC1 for husbands and wives
PCA.paired<-paired_ASW%>%select(ID,Family.ID,Relationship)
PCA.paired<-left_join(PCA.paired,PCA,by="ID")%>%select(Family.ID,Relationship,PC1)
PCA.paired<-PCA.paired%>%spread(key="Relationship",PC1)
ggplot(PCA.paired,aes(x=mother,y=father))+
  geom_point(size=5)+
  stat_cor(size=15)+
  geom_smooth(method='lm',se=F,size=2,color="darkgrey")+
  theme_classic()+
  scale_x_continuous(expand=c(0.05,0))+
  scale_y_continuous(expand=c(0.05,0))+
  labs(y="Father's PC1",
       x="Mother's PC1")+
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25))

ggplot(couples,aes(x=mother*100,y=father*100))+
  geom_point(size=5)+
  stat_cor(aes(label = ..r.label..), label.x = 75, size=15)+
  geom_smooth(method='lm',se=F,size=2,color="darkgrey")+
  theme_classic()+
  scale_x_continuous(expand=c(0.05,0))+
  scale_y_continuous(expand=c(0.05,0))+
  labs(y="Father's global African ancestry (%)",
       x="Mother's global African ancestry (%)")+
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25))

#Build the PCAviz object
paired_plot<-pcaviz(dat=PCA.paired,sdev=eig.val,
             var=sum.eig,rotation=snpeigs)
paired_plot<-pcaviz_abbreviate_var(paired_plot,'Population')

#Make PCA plots - grid
geom.point.summary.params<-list(
  shape=16,stroke=1,size=5,
  show.legend=F)
plot1<-plot(paired_plot,
            coords=paste0("PC",c(1,2)),color='Population',
            geom.point.summary.params=geom.point.summary.params,
            scale.pc.axes=0.6)