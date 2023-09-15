#Fig. 4a
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(indicspecies)
library(superheat)
library(readxl)
#VSURF
meta <- fread("input/meta_pwy.tbl")
gapseq <- read_excel("input/pathway150.xlsx", sheet = "t")
tax<-read_excel("input/pathway150.xlsx", sheet = "anno")
tax_1<-tax[,c(1,25)]
rownames(gapseq) <- gapseq$ID
gapseq<-gapseq[-72,] #move HQB-92
gapseq_1 <- gapseq[,apply(gapseq, 2, function(x){(length(unique(x))>1)})]
#--------------------------------------------------------------------------------------------
#F23
phyno.f23<-tax[,c(10)]
rownames(phyno.f23)<-tax$ID
##k-means
library(factoextra)
library(cluster)
a<-fviz_nbclust(phyno.f23, kmeans, method = "wss", linecolor  = "black",nboot=1000)+
  geom_vline(xintercept = 3, linetype = 2,colour = 'lightgrey')  #k=3
a
b<-fviz_nbclust(phyno.f23, kmeans, method = "gap_stat",linecolor  = "black") #k=3
b
c<-fviz_nbclust(phyno.f23, kmeans, method = "silhouette",
                linecolor  = "black",nboot=1000) #k=3
c

library(ggpubr)
abc<-ggarrange(a,b,c,ncol=3,labels = c('a','b','c'))
abc
#ggsave('abc.pdf',abc,width = 10,height = 5)
##gap_stat <- clusGap(phyno,
##                    FUN = kmeans,
##                   nstart = 100,
##                    K.max = 10,
##                    B = 100)
##fviz_gap_stat(gap_stat)   #k=6
#
## kmeans-clust
set.seed(666)
#
km <- kmeans(phyno.f23, centers = 3, nstart = 1000)
km
##plot-cluster
# plot(silhouette(km$cluster, daisy(phyno.f23)),  
#      border=NA,
#      col=c("#EE2617FF", "#F2A241FF", "#558934FF"#,'#0E54B6FF'
#      ),
#      main = 'Silhouette width of each cluster',
#      xlab='Silhouette width',
#      #asp = 0.01
# )
# 
# aggregate(phyno.f23, by=list(cluster=km$cluster), mean)
# #
# final_data <- cbind(phyno.f23, cluster = km$cluster)
# 
# final_data<-final_data[order(final_data$cluster),]
# final_data<-mutate(final_data,ID=rownames(final_data))
# final_data_1<-final_data[,c(3,2)]
# ISA_input<- merge(final_data_1,gapseq_1, by = "ID")
# ISA_pathway <- ISA_input[,3:ncol(ISA_input)]
# ISA_groups <- ISA_input$cluster
# 
# ISA_output.f23 <- multipatt(x = ISA_pathway,
#                             cluster = ISA_groups,
#                             func = "IndVal.g",
#                             #restcomb = c(1:25,262), #262 is just the grouping of amo genomes together
#                             max.order = 2,
#                             control = how(nperm=9999),
#                             print.perm = TRUE)
# ##FDR
# ISA_output_tab.f23 <- data.table(path = rownames(ISA_output.f23$sign), ISA_output.f23$sign)
# ISA_output_tab.f23$FDR <- p.adjust(ISA_output_tab.f23$p.value, method = "fdr")
# ### Merge Annotations into ISA output table
# # Subfam summaries + remove columns
# meta_1<-meta[,c(1,2,4)]
# ISA_output_tab.f23 <- merge(ISA_output_tab.f23, meta_1, by.x = "path",by.y = 'id')
# ### Filter ISA output to correct FDR and stat 
# ISA_output_tab_filt.f23 <- subset(ISA_output_tab.f23, FDR <= 0.05 & stat >= 0.4)
# #writexl::write_xlsx(ISA_output_tab_filt.f23,"ISA_output_tab.f23.xlsx")
# ISA_output_tab_final.f23 <- subset(ISA_output_tab_filt.f23, index %in% c(1,2,5,6)) 
# #92 pathways
# ISA_output_tab_final.f23<-ISA_output_tab_final.f23[order(ISA_output_tab_final.f23$index),]
# #### Draw Heatmap of Significant KOs
# ### Set up data
# foo <- which(colnames(ISA_input) %in% c("ID", ISA_output_tab_final.f23$path))
# bar <- ISA_input[,foo]
# nex <- as.matrix(bar[,2:ncol(bar)],rownames = bar$ID)
# rownames(nex)<-bar$ID
# t_nex <- t(nex)%>%as.data.frame()
# #order at groups
# t_nex2<-mutate(t_nex,ID=rownames(t_nex))
# idx = match(ISA_output_tab_final.f23$path,t_nex2$ID)
# t_nex3<-t_nex2[idx,]
# t_nex3<-t_nex3[,-75]
# 
# ##order
# ISA_input$cluster<-factor(ISA_input$cluster,levels=c(1,3,2)) 
# ISA_output_tab_final.f23$index<-factor(ISA_output_tab_final.f23$index,levels=c(1,3,2,5,6))   
# final_data2<-final_data[order(final_data$ID),]
# #Fig. S34ab
# superheat(t_nex3,
#           pretty.order.cols = T,
#           pretty.order.rows = T,
#           heat.na.col = "white",
#           heat.pal = c("#E8F5E9FF", "#2E7D32FF"),
#           membership.rows = ISA_output_tab_final.f23$index,
#           membership.cols = ISA_input$cluster,
#           dist.method = "euc",
#           #linkage.method = "ward.D2",
#           scale = F,
#           bottom.label.size = 0.1,
#           bottom.label.text.size = 3,
#           #force.left.label = TRUE,
#           #left.label = "variable",
#           #bottom.label = "variable",
#           #left.label.text.size = 1,
#           left.label.size = 0.05,
#           #left.label.col = "white",
#           grid.hline.col = "black",
#           grid.vline.col = "black",
#           # bottom.label.col = c("#E68098FF","#CA4A68FF","#F3BFCBFF"),
#           legend = FALSE,
#           #smooth.heat = TRUE
#           yt.plot.type = "bar",
#           yt=final_data2$mean_5,
#           yt.axis.name = "F23 growth rate",
#           yt.cluster.col = c("#CA4A68FF","#E68098FF","#F3BFCBFF"),
#           #bottom.label.text.angle = 90
# )
########################################################################################
#select 'beneficial' and 'detrimental' groups   
final_data_sub<-subset(final_data,cluster != 3)  #remove cluster 3
final_data_sub<-final_data_sub[-which(rownames(final_data_sub)%in%c('HQB-88','HQB-440')),]
# 55 strains
#writexl::write_xlsx(final_data_sub,"group.xlsx")
#final_data_sub<-final_data_sub[order(-final_data_sub$mean_5),]
final_data_sub<-mutate(final_data_sub,ID=rownames(final_data_sub))
final_data_sub_1<-final_data_sub[,c(3,2)]


ISA_input<- merge(final_data_sub_1,gapseq_1, by = "ID")
ISA_pathway <- ISA_input[,3:ncol(ISA_input)]
ISA_groups <- ISA_input$cluster

ISA_output.f23 <- multipatt(x = ISA_pathway,
                            cluster = ISA_groups,
                            func = "IndVal.g",
                            #restcomb = c(1:25,262), #262 is just the grouping of amo genomes together
                            max.order = 2,
                            control = how(nperm=9999),
                            print.perm = TRUE)
##FDR
ISA_output_tab.f23 <- data.table(path = rownames(ISA_output.f23$sign), ISA_output.f23$sign)
ISA_output_tab.f23$FDR <- p.adjust(ISA_output_tab.f23$p.value, method = "fdr")
### Merge Annotations into ISA output table
# Subfam summaries + remove columns
meta_1<-meta[,c(1,2,4)]
ISA_output_tab.f23 <- merge(ISA_output_tab.f23, meta_1, by.x = "path",by.y = 'id')
### Filter ISA output to correct FDR and stat 
ISA_output_tab_filt.f23 <- subset(ISA_output_tab.f23, FDR <= 0.05 & stat >= 0.4)
#writexl::write_xlsx(ISA_output_tab_filt.f23,"ISA_output_tab_filt.f23.sub.xlsx")
ISA_output_tab_final.f23 <- ISA_output_tab_filt.f23 
#126 pathwasys
ISA_output_tab_final.f23<-ISA_output_tab_final.f23[order(ISA_output_tab_final.f23$index),]
#### Draw Heatmap of Significant KOs
### Set up data
foo <- which(colnames(ISA_input) %in% c("ID", ISA_output_tab_final.f23$path))
bar <- ISA_input[,foo]
nex <- as.matrix(bar[,2:ncol(bar)],rownames = bar$ID)
rownames(nex)<-bar$ID
#idx1 = match(final_data_sub_1$ID,rownames(nex))
#nex1<-nex[idx1,]
t_nex <- t(nex)%>%as.data.frame()
#order at groups
t_nex2<-mutate(t_nex,ID=rownames(t_nex))
idx = match(ISA_output_tab_final.f23$path,t_nex2$ID)
t_nex3<-t_nex2[idx,]
t_nex3<-t_nex3[,-56]

##order
ISA_input$cluster<-factor(ISA_input$cluster,levels=c(1,2)) #fast to slow
ISA_output_tab_final.f23$index<-factor(ISA_output_tab_final.f23$index,levels=c(1,2))   
final_data_sub2<-final_data_sub[order(final_data_sub$ID),]


#Fig. 34ab
superheat(t_nex3,
          pretty.order.cols = T,
          pretty.order.rows = T,
          heat.na.col = "white",
          heat.pal = c("white", "#B6DBFFFF"),
          membership.rows = ISA_output_tab_final.f23$index,
          membership.cols = ISA_input$cluster,
          dist.method = "euc",
          #linkage.method = "ward.D2",
          scale = F,
          bottom.label.size = 0.1,
          bottom.label.text.size = 3,
          force.left.label = TRUE,
          #left.label = "variable",
          #bottom.label = "variable",
          #left.label.text.size = 0.9,
          left.label.size = 0.05,
          #left.label.col = "white",
          grid.hline.col = "black",
          grid.vline.col = "black",
          bottom.label.col = c("#D6D6D6FF","#ADADADFF"),
          left.label.col = c("#D6D6D6FF","#ADADADFF"),
          legend = FALSE,
          #smooth.heat = TRUE
          yt.plot.type = "bar",
          yt=final_data_sub2$mean_5,
          yt.axis.name = "F23 growth rate",
          yt.cluster.col = c("grey"),
          #bottom.label.text.angle = 90
)

#merge pathways and growth data
merge.tax<-merge(final_data_sub,tax, by = "ID")
View(merge.tax)
merge.tax.path<-merge(merge.tax,gapseq,by = "ID")
write.csv(merge.tax.path,'merge.tax.path.csv')

#class-plot
#clear up in excel
library(ggpubr)
library(ggplot2)
library(ggforce)
library(readxl)
sum<-read_xlsx('./input/sub.cluster-vsurf.xlsx',sheet = 'sum')
ko_cluster1<-
  ggplot()+
  geom_arc_bar(data=sum,
               stat = "pie",
               aes(x0=0,y0=0,r0=0.6,r=2,
                   amount=`cluster1%`,fill=Class,
                   #explode=c(0.05,-0.05,0.05,0.05,0.05,-0.05,-0.2,0.05)
               )
               
  )+theme_void()+
  scale_fill_manual(values = c('Amine Degradation'='#006DDBFF','Amino Acid Metabolism'='#B6DBFFFF',
                               'Aromatic Compounds Metabolism'='#490092FF','Carbohydrates Metabolism'='#FFB6DBFF',
                               'Cell Structure Biosynthesis'='#009292FF','Cofactor Metabolism'='#D4419EFF',
                               'Energy Metabolism'='#6DB6FFFF','Lipid Metabolism'='#920000FF',
                               'Metabolic Regulators'='#24FF24FF','Nucleotide Metabolism'='#DB6D00FF',
                               'Secondary Metabolite Metabolism'='#FFFF6DFF','Others'='#9E9E9E'))+
  theme(legend.key.size = unit(0.1, "inches"),
        legend.title = element_text(face = "bold", size=9))
ko_cluster2<-
  ggplot()+
  geom_arc_bar(data=sum,
               stat = "pie",
               aes(x0=0,y0=0,r0=0.6,r=2,
                   amount=`cluster2%`,fill=Class,
                   #explode=c(0.05,-0.05,0.05,0.05,0.05,-0.05,-0.2,0.05)
               )
               
  )+theme_void()+
  scale_fill_manual(values = c('Amine Degradation'='#006DDBFF','Amino Acid Metabolism'='#B6DBFFFF',
                               'Aromatic Compounds Metabolism'='#490092FF','Carbohydrates Metabolism'='#FFB6DBFF',
                               'Cell Structure Biosynthesis'='#009292FF','Cofactor Metabolism'='#D4419EFF',
                               'Energy Metabolism'='#6DB6FFFF','Lipid Metabolism'='#920000FF',
                               'Metabolic Regulators'='#24FF24FF','Nucleotide Metabolism'='#DB6D00FF',
                               'Secondary Metabolite Metabolism'='#FFFF6DFF','Others'='#9E9E9E'))+
  theme(legend.key.size = unit(0.1, "inches"),
        legend.title = element_text(face = "bold", size=9))
class<-ggarrange(ko_cluster1,ko_cluster2,ncol=1,common.legend = T,legend = 'right')
class
ggsave('class-pathway.pdf',class,width = 5,height = 6)

#Fig. 4a
#VSURF
library(VSURF)
merge<-merge(final_data_sub,gapseq_1, by = "ID")
set.seed(666)
rf.path<-VSURF(x=merge[,-c(1,2,3)], y=merge$mean_5,  parallel=T)
gapseq_2<-merge[,-c(1,2,3)]
rownames(gapseq_2)<-merge$ID
rf.path.imp <- data.table(nr=rf.path$imp.mean.dec.ind, id=colnames(gapseq_2)[rf.path$imp.mean.dec.ind], imp.mean=rf.path$imp.mean.dec, imp.sd=rf.path$imp.sd.dec) 
rf.path.imp[,name:=meta$name[match(id, meta$id)]]

gapseq_3<- t(gapseq_2)
gapseq_3<- as.data.frame(gapseq_3)
gapseq_3$id<- rownames(gapseq_3)

merge.path.cluster2<-merge(rf.path.imp,gapseq_3,by='id',sort=F) %>% t
merge.path.cluster2<-data.frame(names = row.names(merge.path.cluster2), merge.path.cluster2)
write.csv(rf.path.imp,'rf.path.imp.csv')
write.csv(merge.path.cluster2,'merge.path.cluster2.csv')

#clear up in excel
#plot
library(pheatmap)
library(ggpubr)


library(readxl)
dt <- read_excel("./input/sub.cluster-vsurf.xlsx", sheet = "plot.top30.new")
#remove first col
dt1<-dt[,-c(1)]
rownames(dt1)=dt$`Pathway name`


anno_col<-read_excel("./input/sub.cluster-vsurf.xlsx", sheet = "anno")

anno_col1 <- anno_col[,c(26,24,#6,
                         11
)]
rownames(anno_col1)=anno_col$ID
anno_col1<-as.data.frame(anno_col1)
colnames(anno_col1)<-c('Phylum / Class','Gram strain',#'Log2(N2)',
                       'F23'
)

anno_colors = list('Phylum / Class'=c('Actinobacteria'='#EEE685','Bacteroidetes'='#4F94CD',
                                      'Firmicutes'='#EE7600','Alphaproteobacteria'='#00CD00',
                                      'Gammaproteobacteria'='#90EE90'),
                   'Gram strain' = c('Negative'='#7E8BD6FF','Positive'='#E84A5FFF'),
                   #'Log2(N2)'=c("#004D40FF","#E0F2F1FF"),
                   'F23'=c('#E3F2FDFF','#0D47A1FF')
)

#plot  Fig. 4a
pheatmap(dt1, cluster_cols = T, cluster_rows = F, 
         clustering_distance_cols = "euclidean",clustering_method='average',
         color= c("#F9FBE7FF","#604A76FF"),legend = F,
         fontsize = 5,fontsize_row = 3, fontsize_col = 3.2,angle_col = 45,
         cutree_cols = 2,
         scale="none",
         # border_color= 'white',
         border=F,
         #annotation_row = anno_row,
         annotation_col = anno_col1,
         annotation_colors = anno_colors,
         #gaps_row=c(4,30,35,59),
         annotation_names_row=F,annotation_legend=T,treeheight_col=10,
         cellwidth = 4, cellheight = 3,
         #filename='sig.path.t30.pdf',width=10,height=8
)

#Fig. S18  all HQiome Vsurf
#VSURF-ALL74
phyno.f23<-tax[,c(1,10)]
rownames(phyno.f23)<-tax$ID
merge74<-merge(phyno.f23,gapseq_1,by='ID')
merge74.path<-merge74[,-c(1,2)]
rownames(merge74.path)<-merge74$ID
rf.74<-VSURF(x=merge74.path, y=merge74$mean_5,  parallel=T)
gapseq_1.1<-gapseq_1[,-1]
rownames(gapseq_1.1)<-gapseq_1$ID
rf.74.imp <- data.table(nr=rf.74$imp.mean.dec.ind, id=colnames(gapseq_1.1)[rf.74$imp.mean.dec.ind], imp.mean=rf.74$imp.mean.dec, imp.sd=rf.74$imp.sd.dec) 
rf.74.imp[,name:=meta$name[match(id, meta$id)]]

gapseq_3.2<- t(gapseq_1.1)
gapseq_3.2<- as.data.frame(gapseq_3.2)
gapseq_3.2$id<- rownames(gapseq_3.2)

merge.path.74<-merge(rf.74.imp,gapseq_3.2,by='id',sort=F) %>% t
merge.path.74<-data.frame(names = row.names(merge.path.74), merge.path.74)
write.csv(rf.74.imp,'rf.74.imp.csv')
write.csv(merge.path.74,'merge.path.74.csv')

#wilcoxon test

library(psych)
sp.b150 <- corr.test(merge74$mean_5,merge74.path, method = 'spearman',adjust='BH')
sp.b150.imp <- data.table(id=colnames(sp.b150$r), r=t(sp.b150$r),absr=abs(t(sp.b150$r)),t=t(sp.b150$t), p=t(sp.b150$p))
sp.b150.imp[,name:=meta$name[match(id, meta$id)]]

library("xlsx")
write.xlsx(sp.b150.imp,'Lmarina_74pathway_sig_test.xlsx')

#PLOT Fig. S18 
#plot
library(pheatmap)
library(ggpubr)

#plot-all
library(readxl)
dt <- read_excel("input/74plot.xlsx", sheet = "pathway")

dt1<-dt[,-c(1)]
rownames(dt1)=dt$Name


anno_col<-read_excel("./input/74plot.xlsx", sheet = "anno")

anno_col1 <- anno_col[,c(10, 23, 25)]
rownames(anno_col1)=anno_col$ID
anno_col1<-as.data.frame(anno_col1)
colnames(anno_col1)<-c('Growth','Gram strain','Phylum / Class')

anno_colors = list('Growth'=c('#E3F2FDFF','#0D47A1FF'),
                   'Gram strain' = c('Negative'='#7E8BD6FF','Positive'='#E84A5FFF'),
                   'Phylum / Class'=c('Actinobacteria'='#EEE685','Bacteroidetes'='#4F94CD',
                                      'Firmicutes'='#EE7600','Alphaproteobacteria'='#00CD00',
                                      'Gammaproteobacteria'='#90EE90'))

#plot
pheatmap(dt1, cluster_cols = T, cluster_rows = F, 
         clustering_distance_cols = "euclidean",clustering_method='average',
         color= c("#F9FBE7FF","#604A76FF"),legend = F,
         fontsize = 5,fontsize_row = 3, fontsize_col = 3.2,angle_col = 45,
         cutree_cols = 2,
         scale="none",
         # border_color= 'white',
         border=F,
         #annotation_row = anno_row,
         annotation_col = anno_col1,
         annotation_colors = anno_colors,
         #gaps_row=c(4,30,35,59),
         annotation_names_row=F,annotation_legend=T,treeheight_col=10,
         cellwidth = 4, cellheight = 3,
         #filename='all_F23_sig.path.t30.pdf',width=10,height=8
)

#========================================================================================
#Fig. 5a
#N2 -VSURF
library(VSURF)
phyno.n2<-tax[,c(1,4)]
rownames(phyno.n2)<-tax$ID
merge74<-merge(phyno.n2,gapseq_1,by='ID')
merge74.path<-merge74[,-c(1,2)]
rownames(merge74.path)<-merge74$ID
rf.74<-VSURF(x=merge74.path, y=merge74$N2,  parallel=T)
gapseq_1.1<-gapseq_1[,-1]
rownames(gapseq_1.1)<-gapseq_1$ID
rf.74.imp <- data.table(nr=rf.74$imp.mean.dec.ind, id=colnames(gapseq_1.1)[rf.74$imp.mean.dec.ind], imp.mean=rf.74$imp.mean.dec, imp.sd=rf.74$imp.sd.dec) 
rf.74.imp[,name:=meta$name[match(id, meta$id)]]

gapseq_3.2<- t(gapseq_1.1)
gapseq_3.2<- as.data.frame(gapseq_3.2)
gapseq_3.2$id<- rownames(gapseq_3.2)

merge.path.74<-merge(rf.74.imp,gapseq_3.2,by='id',sort=F) %>% t
merge.path.74<-data.frame(names = row.names(merge.path.74), merge.path.74)
write.csv(rf.74.imp,'rf.n2.74.imp.csv')
write.csv(merge.path.74,'merge.n2.path.74.csv')

#GDP-mannose biosynthesis
cor.test(merge74.path$`|PWY-1781|`, merge74$N2, 
         method = "spearman")

#wilcoxon test
library(psych)
sp.b150 <- corr.test(merge74$N2,merge74.path, method = 'spearman',adjust='BH')
sp.b150.imp <- data.table(id=colnames(sp.b150$r), r=t(sp.b150$r),absr=abs(t(sp.b150$r)),t=t(sp.b150$t), p=t(sp.b150$p))
sp.b150.imp[,name:=meta$name[match(id, meta$id)]]

library("xlsx")
write.xlsx(sp.b150.imp,'CELEGANS_pathway_sig_test.xlsx')

#clear up in excel
#plot Fig. 5a
library(pheatmap)
library(ggpubr)

#plot-all
library(readxl)
dt <- read_excel("./input/N2_pathway_plot.xlsx", sheet = "path")

dt1<-dt[,-c(1)]
rownames(dt1)=dt$Name

anno_col<-read_excel("./input/N2_pathway_plot.xlsx", sheet = "anno")

anno_col1 <- anno_col[,c(4, 23, 25)]
rownames(anno_col1)=anno_col$ID
anno_col1<-as.data.frame(anno_col1)
colnames(anno_col1)<-c('Growth','Gram strain','Phylum / Class')

anno_colors = list('Growth'=c("#004D40FF","#E0F2F1FF"),
                   'Gram strain' = c('Negative'='#7E8BD6FF','Positive'='#E84A5FFF'),
                   'Phylum / Class'=c('Actinobacteria'='#EEE685','Bacteroidetes'='#4F94CD',
                                      'Firmicutes'='#EE7600','Alphaproteobacteria'='#00CD00',
                                      'Gammaproteobacteria'='#90EE90'))

#plot
pheatmap(dt1, cluster_cols = T, cluster_rows = F, 
         clustering_distance_cols = "euclidean",clustering_method='average',
         color= c("#F9FBE7FF","#604A76FF"),legend = F,
         fontsize = 5,fontsize_row = 3, fontsize_col = 3.2,angle_col = 45,
         cutree_cols = 4,
         scale="none",
         # border_color= 'white',
         border=F,
         #annotation_row = anno_row,
         annotation_col = anno_col1,
         annotation_colors = anno_colors,
         #gaps_row=c(4,30,35,59),
         annotation_names_row=F,annotation_legend=T,treeheight_col=10,
         cellwidth = 4, cellheight = 3,
         #filename='N2_sig.path.t30.pdf',width=10,height=8
)

#=====================================================================


##################################################################################
#Fig. 4b and Fig. 5b
#lifespan
library(readxl)
data <- read_excel("input/lifespan.xlsx",sheet = 'F23_plot')
#shapiro.test and  bartlett.test
shapiro.test(data$lifespan)  #  p-value = 0.9065
bartlett.test(lifespan ~ group, data = data)   # p-value = 0.1604

#growth-f23 
shapiro.test(data$d5)  #p-value = 0.0001439   
bartlett.test(d5 ~ group2, data = data)   #p-value < 2.2e-16   

#test
library(rstatix)
#ANOVA 
data %>% anova_test(lifespan ~ group)  # p=4.48e-06
#data %>% anova_test(d5 ~ group)  

#t test
lifespan_t<-data %>%     #
  t_test(lifespan ~ group, paired = F) 


library(agricolae)
#Turkey para
model2 <- aov(lifespan~group, data = data)
#out2 <- LSD.test(model2, "group", p.adj = "BH" )  #LSD
out2 <- HSD.test(model2, "group", alpha = 0.05 )  #HSD turkey
print(out2$groups)
plot(out2)
#add sig
stat_lifespan=out2$groups
data$p.lifespan=stat_lifespan[as.character(data$group),]$groups

#stat_growth=out1$groups
#data$stat2=stat_growth[as.character(data$group),]$groups


#growth  non-para
#data %>% kruskal_test(lifespan ~ group)  # P=0.0021 
data %>% kruskal_test(d5 ~ group2)   # P=0.000177 

##Kruskal Wallis test and multiple comparison
p.growth<-with(data,kruskal(d5,group,group=TRUE, p.adj='BH',main='Growth'))
#add sig
stat_growth=p.growth$groups
data$p.growth=stat_growth[as.character(data$group),]$groups

#============================================================================================================

#genus test
shapiro.test(data$lifespan)  #  p-value = 0.9065  
bartlett.test(lifespan ~ genus, data = data)   # p-value = 0.6421   

#growth-f23 
shapiro.test(data$d5)  #p-value = 0.0001439    
bartlett.test(d5 ~ genus, data = data)   #p-value < 0.05823    

#t检验---lifespan 
data %>%     
  t_test(lifespan ~ genus, paired = F) #0.0265

#  ---groth     wilcox dunntt
#p.lifespan_genus <- data %>% 
#  wilcox_test(lifespan ~ genus, p.adjust.method = "fdr")

data %>% 
  wilcox_test(d5 ~ genus, p.adjust.method = "fdr")   # p=0.00112

#dunn_test  
#p.lifespan <- data %>% 
#  dunn_test(lifespan ~ group, p.adjust.method = "fdr")

data %>% 
  dunn_test(d5 ~ genus, p.adjust.method = "fdr")      # p=0.00109

#kruskal---growth    p=6e-04
p.growth_2<-with(data,kruskal(d5,genus,group=F, p.adj='BH',main='Growth'))
p.growth_2
#vibrio and shewanella have signicant differences at genus

#=====================================================================================
library(ggplot2)
library(ggpubr)
#p.d5<-data %>% 
#  tukey_hsd(d5 ~ group, paired = T) %>% add_y_position()

#p.lifespan<- data %>%  
#  tukey_hsd(lifespan ~ group, paired = T)  %>% add_y_position(step.increase = 0.48)

#Fig. S15a
growth<-ggboxplot(
  data, x = "group", y = "d5",
  color = "genus", palette = c('#009292FF','#F9791EFF')
  , #ylim = c(0,0.8),
  ylab='Growth rate (%)',xlab = '',size=1
)+
  #stat_pvalue_manual(p.d5, label = "p.adj.signif")+
  rotate_x_text(45,size=10)+geom_text(aes(label = p.growth),size=3)
growth

#Fig. S15b
lifespan<-ggboxplot(
  data, x = "group", y = "lifespan",
  color = "genus", palette = c('#009292FF','#F9791EFF'), 
  ylim = c(5,15),
  ylab='Mean lifespan (day)',xlab = '',size=1
)+
  #stat_pvalue_manual(p.lifespan, label = "p.adj.signif",hide.ns=T)+
  rotate_x_text(45,size=10)+geom_text(aes(label = p.lifespan),size=3)
lifespan

merge<-ggarrange(growth, lifespan, common.legend = T, labels = c('a','b'))
merge
ggsave('merge3.pdf',merge,height = 4,width = 9)

#pathway corralation

#survival curve
library(tidyverse)
library(ggsurvfit)
library(gghighlight)
library(tidycmprsk)
library(paletteer) 
data <- read_excel("input/lifespan.xlsx",sheet = 'f23_ggsurvfit')
#Fig. S15c
surcical1<-survfit2(Surv(time, status) ~ species, data = data) %>% 
  ggsurvfit(size = 1)+
  add_censor_mark()+
  #add_confidence_interval() +
  #add_risktable(theme=theme_test()+
  #                theme(axis.title = element_blank(),
  #                      axis.text.x = element_blank(),
  #                      axis.ticks.x=element_blank(),
  #                      axis.text.y = element_text(color="black",size=10)),
  #              combine_groups=F)+
  add_quantile(color ="grey80",size=0.8,linetype =5)+theme_classic()+
  labs(y = "Percentage Survival",x='Day')+
  #add_pvalue(caption = "Log-rank {p.value}",location  = "annotation", x = 20)+
  guides(colour =guide_legend(ncol=2))+
  scale_color_manual(values = c(paletteer_d("ggthemes::hc_darkunica"),paletteer_d("basetheme::brutal")))+
  scale_fill_manual(values = c(paletteer_d("ggthemes::hc_darkunica"),paletteer_d("basetheme::brutal")))+
  scale_x_continuous(breaks = seq(0, 30, by = 1))
surcical1
ggsave('sur1.pdf',surcical1,width = 7,height = 3)

#shewanella vs. vibrio
# REMOVED
# data2<-subset(data,genus!='OP50')
# surcical2<-survfit2(Surv(time, status) ~ genus, data = data2) %>% 
#   ggsurvfit(size = 1)+
#   add_censor_mark()+
#   #add_confidence_interval() +
#   #add_risktable(theme=theme_test()+
#   #                theme(axis.title = element_blank(),
#   #                      axis.text.x = element_blank(),
#   #                      axis.ticks.x=element_blank(),
#   #                      axis.text.y = element_text(color="black",size=10)),
#   #              combine_groups=F)+
#   add_quantile(color ="grey80",size=0.8,linetype =5)+theme_classic()+
#   labs(y = "Percentage Survival",x='Day')+
#   add_pvalue(caption = "Log-rank {p.value}",location  = "annotation", x = 15)+  #添加P值---log rank
#   scale_color_manual(values = c('#009292FF','#F9791EFF'))+
#   scale_fill_manual(values = c('#009292FF','#F9791EFF'))+
#   scale_x_continuous(breaks = seq(0, 30, by = 1))
# surcical2
# ggsave('sur2.pdf',surcical2,width = 6,height = 3)
# 
# library(survival)
# survdiff(Surv(time, status) ~ species, data = data)
# p <- survfit2(Surv(time, status) ~ genus, data = data2) %>%
#   ggsurvfit()
# 
# #pairwise p values
# library(survminer)
# pair<-pairwise_survdiff(Surv(time, status) ~ species,
#                         data = data)
# pair.frame<-as.data.frame(pair$p.value)
# pair.frame.star<-symnum(pair$p.value, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
#                         symbols = c("****", "***", "**", "*", "+", " "),
#                         abbr.colnames = FALSE, na = "")
# write.csv(pair.frame,'pair.frame.csv')
# write.csv(pair.frame.star,'pair.frame.star.csv')

###===============================================================================================
#correlation  Lifespan and pathways for Fig. 4b
library(readxl)
library(VSURF)
library(dplyr) 
data <- read_excel("input/lifespan.xlsx",sheet = 'F23_plot')
data2<-data %>%
  group_by(group) %>%
  summarise_all(list(mean))
data2<-data2[,c(1,3,4)]

library(data.table)
meta <- fread("input/meta_pwy.tbl")

gapseq <- read_excel("input/pathway150.xlsx", sheet = "t")
tax<-read_excel("input/pathway150.xlsx", sheet = "anno")
tax_1<-tax[,c(1,25)]

idx = match(data2$group,gapseq$ID)
gapseq_1<-gapseq[idx,]
rownames(gapseq_1) <- gapseq_1$ID

gapseq_2 <- gapseq_1[,apply(gapseq_1, 2, function(x){(length(unique(x))>1)})]
set.seed(666)
rf.path<-VSURF(x=gapseq_2[,-1], y=data2$lifespan,  parallel=T)
rownames(gapseq_2)<-gapseq_2$ID
rf.path.imp <- data.table(nr=rf.path$imp.mean.dec.ind, id=colnames(gapseq_2)[rf.path$imp.mean.dec.ind], imp.mean=rf.path$imp.mean.dec, imp.sd=rf.path$imp.sd.dec) 
rf.path.imp[,name:=meta$name[match(id, meta$id)]]
write.csv(rf.path.imp,'rf.path.imp.csv')

#P test

cor.test(gapseq_2$`|RHAMCAT-PWY|`, data2$lifespan, 
         method = "spearman")

#wilcoxon
library(psych)
gapseq_2<-gapseq_2[,-1]
cor.test <- corr.test(data2$lifespan,gapseq_2, method = 'spearman',adjust='BH')
cor.test.imp <- data.table(id=colnames(cor.test$r), r=t(cor.test$r),absr=abs(t(cor.test$r)),t=t(cor.test$t), p=t(cor.test$p))
cor.test.imp[,name:=meta$name[match(id, meta$id)]]

library("xlsx")
write.xlsx(cor.test.imp,'lifespan_pathway_sig_test.xlsx')

#PWY-4261  p-value = 0.0191
#PWY-5901 p-value = 0.3237
#PWY4 p-value = 0.08253
#PWY-6938-1 p-value = 0.7183
#RHAMCAT-PWYp-value = 0.2451

#growth
set.seed(666)
rf.path.growth<-VSURF(x=gapseq_2[,-1], y=data2$d5,  parallel=T)
rf.path.growth.imp <- data.table(nr=rf.path.growth$imp.mean.dec.ind, id=colnames(gapseq_2)[rf.path.growth$imp.mean.dec.ind], imp.mean=rf.path.growth$imp.mean.dec, imp.sd=rf.path.growth$imp.sd.dec) 
rf.path.growth.imp[,name:=meta$name[match(id, meta$id)]]
write.csv(rf.path.growth.imp,'rf.growth.imp.csv')

##==================================================================================================
#plot-heatmap
#Fig. 4b
library(readxl)
library(pheatmap)
library(ggpubr)
data_p<-read_excel('input/cor_results.xlsx', sheet='f23.plot')
lifespan_p<- data_p[,c(1,2)]
pathway_p<-data_p[,c(1,3,4)]
lifespan_p<-lifespan_p[,-1]
pathway_p<-pathway_p[,-1]
rownames(lifespan_p)=data_p$ID
rownames(pathway_p)=data_p$ID
lifespan_p<-as.matrix(lifespan_p)
pathway_p<-as.matrix(pathway_p)
#lifespan
pheatmap(lifespan_p, cluster_cols = F, cluster_rows = F, 
         color= c("#E4F1E1FF", '#63A6A0FF','#0D585FFF'),legend = T,
         fontsize = 8,fontsize_row = 8, fontsize_col = 8,#angle_col = 45,
         scale="none",
         # border_color= 'white',
         border=F,border_color = "grey",
         cellwidth = 10, cellheight = 10,
         #filename='sig.path.t30.pdf',width=10,height=8
)
#pathway
pheatmap(pathway_p, cluster_cols = F, cluster_rows = F, 
         color= c("#F9FBE7FF","#604A76FF"),legend = F,
         fontsize = 8,fontsize_row = 8, fontsize_col = 8,#angle_col = 45,
         scale="none",
         # border_color= 'white',
         border=F,border_color = "grey",
         cellwidth = 10, cellheight = 10,
         #filename='lifespan_pathway.pdf',width=10,height=8
)


##==================================================================================================
#Fig. S34 
#vibrio and shewanella ON LIFESPAN
#genus
genus_id<-match(data2$group,data$group)
genus<-data[genus_id,]
genus<-genus[,c(2,5)]
data3<-merge(data2,genus,by = 'group')

ISA_input<- merge(genus,gapseq_2, by.x = "group",by.y='ID')
ISA_pathway <- ISA_input[,3:ncol(ISA_input)]
ISA_groups <- ISA_input$genus

library(indicspecies)
ISA_output <- multipatt(x = ISA_pathway,
                        cluster = ISA_groups,
                        func = "IndVal.g",
                        #restcomb = c(1:25,262), #262 is just the grouping of amo genomes together
                        max.order = 3,
                        control = how(nperm=9999),
                        print.perm = TRUE)

#FDR
ISA_output_tab <- data.table(path = rownames(ISA_output$sign), ISA_output$sign)
ISA_output_tab$FDR <- p.adjust(ISA_output_tab$p.value, method = "fdr")
## Merge Annotations into ISA output table
# Subfam summaries + remove columns
meta_1<-meta[,c(1,2,4)]
ISA_output_tab <- merge(ISA_output_tab, meta_1, by.x = "path",by.y = 'id')
## Filter ISA output to correct FDR and stat 
ISA_output_tab_filt <- subset(ISA_output_tab, FDR <= 0.05 & stat >= 0.4)
write.csv(ISA_output_tab_filt,'correlation/ISA_output_tab_filt.csv')

## Draw Heatmap
## Set up data
foo <- which(colnames(ISA_input) %in% c("group", ISA_output_tab_filt$path))
bar <- ISA_input[,foo]
nex <- as.matrix(bar[,2:ncol(bar)],rownames = bar$group)
rownames(nex)<-bar$group
t_nex <- t(nex)%>%as.data.frame()
#order by groups
t_nex2<-mutate(t_nex,ID=rownames(t_nex))
idx = match(ISA_output_tab_filt$path,t_nex2$ID)
t_nex3<-t_nex2[idx,]
t_nex3<-t_nex3[,-18]

ISA_input$genus<-factor(ISA_input$genus,levels=c('Shewanella','Vibrio'))
library(superheat)

#Fig. s34cd
superheat(t_nex3,
          #pretty.order.cols = T,
          #pretty.order.rows = T,
          heat.na.col = "white",
          heat.pal = c("white", "#2A363BFF"),##B6DBFFFF
          membership.rows = ISA_output_tab_filt$index,
          membership.cols = ISA_input$genus,
          dist.method = "euc",
          #linkage.method = "ward.D2",
          scale = F,
          bottom.label.size = 0.1,
          bottom.label.text.size = 3,
          force.left.label = TRUE,
          #left.label = "variable",
          #left.label.text.size = 1,
          left.label.size = 0.05,
          left.label.col = c('#009292FF','#F9791EFF'),
          grid.hline.col = "white",
          grid.vline.col = "white",
          bottom.label.col = c('#5DA5DAFF','#009292FF','#F9791EFF'),
          legend = FALSE,
          #smooth.heat = TRUE
)

#class-plot
#clear up in excel
library(ggpubr)
library(ggplot2)
library(ggforce)
library(readxl)
sum<-read_xlsx('input/cor_results.xlsx',sheet = 'f23.Sum')
shewanella<-
  ggplot()+
  geom_arc_bar(data=sum,
               stat = "pie",
               aes(x0=0,y0=0,r0=0.6,r=2,
                   amount=`cluster1%`,fill=Class,
                   #explode=c(0.05,-0.05,0.05,0.05,0.05,-0.05,-0.2,0.05)
               )
               
  )+theme_void()+
  scale_fill_manual(values = c('Amino Acid Metabolism'='#B6DBFFFF',
                               'Aromatic Compounds Metabolism'='#490092FF',
                               'Carbohydrates Metabolism'='#FFB6DBFF',
                               'Cell Structure Biosynthesis'='#009292FF','Cofactor Metabolism'='#D4419EFF',
                               'Detoxification'='#8FA87AFF',
                               'Energy Metabolism'='#6DB6FFFF','Lipid Metabolism'='#920000FF',
                               'Nucleotide Metabolism'='#DB6D00FF','Olefins Metabolism'='#D070B9FF',
                               'Phosphorus Compounds Metabolism'='#BF616AFF','Polyamine Metabolism'='#E7D202FF',
                               'Secondary Metabolite Metabolism'='#FFFF6DFF','Sulfur Metabolism'='#66CDAAFF',
                               'Others'='#9E9E9E'))+
  theme(legend.key.size = unit(0.1, "inches"),
        legend.title = element_text(face = "bold", size=9))
shewanella

Vibrio<-
  ggplot()+
  geom_arc_bar(data=sum,
               stat = "pie",
               aes(x0=0,y0=0,r0=0.6,r=2,
                   amount=`cluster2%`,fill=Class,
                   #explode=c(0.05,-0.05,0.05,0.05,0.05,-0.05,-0.2,0.05)
               )
               
  )+theme_void()+
  scale_fill_manual(values = c('Amino Acid Metabolism'='#B6DBFFFF',
                               'Aromatic Compounds Metabolism'='#490092FF',
                               'Carbohydrates Metabolism'='#FFB6DBFF',
                               'Cell Structure Metabolism'='#009292FF','Cofactor Metabolism'='#D4419EFF',
                               'Detoxification'='#8FA87AFF',
                               'Energy Metabolism'='#6DB6FFFF','Lipid Metabolism'='#920000FF',
                               'Nucleotide Metabolism'='#DB6D00FF','Olefins Metabolism'='#D070B9FF',
                               'Phosphorus Compounds Metabolism'='#BF616AFF','Polyamine Metabolism'='#E7D202FF',
                               'Secondary Metabolite Metabolism'='#FFFF6DFF','Sulfur Metabolism'='#66CDAAFF',
                               'Others'='#9E9E9E'))+
  theme(legend.key.size = unit(0.1, "inches"),
        legend.title = element_text(face = "bold", size=9))
Vibrio
class<-ggarrange(shewanella,Vibrio,ncol=2,common.legend = T,legend = 'right')
class
ggsave('class-pathway.pdf',class,width = 5,height = 1.5)


#################################################################################################
#Fig. S35
#N2-plot
library(readxl)
library(rstatix)
library(agricolae)
data_N2 <- read_excel("input/lifespan.xlsx",sheet = 'N2_plot')
#lifespan
 
#lifespan----
shapiro.test(data_N2$lifespan)  # p=0.03428,  
bartlett.test(lifespan ~ group, data = data_N2)   # p=0.009233, 

#growth-f23 
shapiro.test(data_N2$egg)  # p=2.898e-10,   
bartlett.test(egg ~ group, data = data_N2)   #p=2.2e-16, 

data_N2 %>% kruskal_test(lifespan ~ group)  # P=0.00088 
data_N2 %>% kruskal_test(egg ~ group)   # P=0.000106 

##Kruskal Wallis test and multiple comparison
p.lifespan<-with(data_N2,kruskal(lifespan,group,group=T, p.adj='BH'))
p.growth<-with(data_N2,kruskal(egg,group,group=TRUE, p.adj='BH'))
#add sig
stat_lifespan=p.lifespan$groups
stat_growth=p.growth$groups
data_N2$p.lifespan=stat_lifespan[as.character(data_N2$group),]$groups
data_N2$p.growth=stat_growth[as.character(data_N2$group),]$groups
#============================================================================================================

#at genus
shapiro.test(data_N2$lifespan)  #  p-value = 0.03428   
bartlett.test(lifespan ~ genus, data = data_N2)   # p-value = 0.08648  

#growth-f23 ---------不正太,方差齐
shapiro.test(data_N2$egg)  #p-value = 2.898e-10 
bartlett.test(egg ~ genus, data = data_N2)   #p-value=0.4491 

data_N2 %>% 
  wilcox_test(egg ~ genus, p.adjust.method = "fdr")   # p=0.274


data_N2 %>% 
  wilcox_test(lifespan ~ genus, p.adjust.method = "fdr")   # p=0.581

#dunn_test  
data_N2 %>%  
  dunn_test(lifespan ~ genus, p.adjust.method = "fdr")   # p=0.571

data_N2 %>% 
  dunn_test(egg ~ genus, p.adjust.method = "fdr")      # p=0.270

#kruskal ---growth    p=0.2739   survival  p=0.5765
p.growth_2<-with(data_N2,kruskal(lifespan,genus,group=F, p.adj='BH',main='Growth'))
p.growth_2


#vibrio and shewanella do not have difference at genus level
#=====================================================================================
library(ggplot2)
library(ggpubr)
#Fig. S35a
#p.d5<-data %>% 
#  tukey_hsd(d5 ~ group, paired = T) %>% add_y_position()
#
#p.lifespan<- data %>%  
#  tukey_hsd(lifespan ~ group, paired = T)  %>% add_y_position(step.increase = 0.48)
growth<-ggboxplot(
  data_N2, x = "group", y = "egg",
  color = "genus", palette = c('#009292FF','#F9791EFF')
  , #ylim = c(0,0.8),
  ylab='Egg laying time (h)',xlab = '',size=1
)+
  #stat_pvalue_manual(p.d5, label = "p.adj.signif")+
  rotate_x_text(45,size=10)+geom_text(aes(label = p.growth),size=3)
growth
#Fig. S35b
lifespan<-ggboxplot(
  data_N2, x = "group", y = "lifespan",
  color = "genus", palette = c('#009292FF','#F9791EFF'), 
  #ylim = c(5,25),
  ylab='Mean lifespan (day)',xlab = '',size=1
)+
  #stat_pvalue_manual(p.lifespan, label = "p.adj.signif",hide.ns=T)+
  rotate_x_text(45,size=10)+geom_text(aes(label = p.lifespan),size=3)
lifespan

merge<-ggarrange(growth, lifespan, common.legend = T, labels = c('a','b'))
merge
ggsave('merge-N2.pdf',merge,height = 4,width = 9)

#=================================================================================
#N2  survival curve
library(tidyverse)
library(ggsurvfit)
library(gghighlight)
library(tidycmprsk)
library(paletteer) 
data <- read_excel("input/lifespan.xlsx",sheet = 'N2_ggsurvfit')

#Supplementary Fig. 29c
surcical1<-survfit2(Surv(time, status) ~ species, data = data) %>% 
  ggsurvfit(size = 1)+
  add_censor_mark()+
  #add_confidence_interval() +
  #add_risktable(theme=theme_test()+
  #                theme(axis.title = element_blank(),
  #                      axis.text.x = element_blank(),
  #                      axis.ticks.x=element_blank(),
  #                      axis.text.y = element_text(color="black",size=10)),
  #              combine_groups=F)+
  add_quantile(color ="grey80",size=0.8,linetype =5)+theme_classic()+
  labs(y = "Percentage Survival",x='Day')+
  #add_pvalue(caption = "Log-rank {p.value}",location  = "annotation", x = 20)+
  guides(colour =guide_legend(ncol=2))+
  scale_color_manual(values = c(paletteer_d("ggthemes::hc_darkunica"),paletteer_d("basetheme::brutal")))+
  scale_fill_manual(values = c(paletteer_d("ggthemes::hc_darkunica"),paletteer_d("basetheme::brutal")))+
  scale_x_continuous(breaks = seq(0, 38, by = 1))
surcical1

ggsave('sur_N2.pdf',surcical1,width = 9,height = 4)
#shewanella vs vibrio
#Supplementary Fig. 29d
data2<-subset(data,genus!='OP50')
surcical2<-survfit2(Surv(time, status) ~ genus, data = data2) %>% 
  ggsurvfit(size = 1)+
  add_censor_mark()+
  #add_confidence_interval() +
  #add_risktable(theme=theme_test()+
  #                theme(axis.title = element_blank(),
  #                      axis.text.x = element_blank(),
  #                      axis.ticks.x=element_blank(),
  #                      axis.text.y = element_text(color="black",size=10)),
  #              combine_groups=F)+
  add_quantile(color ="grey80",size=0.8,linetype =5)+theme_classic()+
  labs(y = "Percentage Survival",x='Day')+
  add_pvalue(caption = "Log-rank {p.value}",location  = "annotation", x = 15)+  #添加P值---log rank
  scale_color_manual(values = c('#009292FF','#F9791EFF'))+
  scale_fill_manual(values = c('#009292FF','#F9791EFF'))+
  scale_x_continuous(breaks = seq(0, 38, by = 1))
surcical2
ggsave('sur_N2-2.pdf',surcical2,width = 9,height = 4)

#log rank p value
library(survival)
survdiff(Surv(time, status) ~ genus, data = data2)
p <- survfit2(Surv(time, status) ~ genus, data = data2) %>%
  ggsurvfit()

#pairwise p values
library(survminer)
pair<-pairwise_survdiff(Surv(time, status) ~ species,
                        data = data)
pair.frame<-as.data.frame(pair$p.value)
pair.frame.star<-symnum(pair$p.value, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
                        symbols = c("****", "***", "**", "*", "+", " "),
                        abbr.colnames = FALSE, na = "")
write.csv(pair.frame,'pair.frame.N2.csv')
write.csv(pair.frame.star,'pair.frame.star.N2.csv')

###===============================================================================================
#Fig. 5b
#N2-pathways correlation----VSURF
library(VSURF)
library(dplyr) 
data2<-data_N2 %>%
  group_by(group) %>%
  summarise_all(list(mean))
data2<-data2[,c(1,3,4)]

library(data.table)
meta <- fread("input/meta_pwy.tbl")

gapseq <- read_excel("input/pathway150.xlsx", sheet = "t")
tax<-read_excel("input/pathway150.xlsx", sheet = "anno")
tax_1<-tax[,c(1,25)]

idx = match(data2$group,gapseq$ID)
gapseq_1<-gapseq[idx,]
rownames(gapseq_1) <- gapseq_1$ID


gapseq_2 <- gapseq_1[,apply(gapseq_1, 2, function(x){(length(unique(x))>1)})]

set.seed(666)
rf.path<-VSURF(x=gapseq_2[,-1], y=data2$lifespan,  parallel=T)
rownames(gapseq_2)<-gapseq_2$ID
rf.path.imp <- data.table(nr=rf.path$imp.mean.dec.ind, id=colnames(gapseq_2)[rf.path$imp.mean.dec.ind], imp.mean=rf.path$imp.mean.dec, imp.sd=rf.path$imp.sd.dec) 
rf.path.imp[,name:=meta$name[match(id, meta$id)]]
write.csv(rf.path.imp,'N2.rf.lifespan.imp.csv')

#p test
library(psych)
gapseq_2<-gapseq_2[,-1]
cor.test <- corr.test(data2$lifespan,gapseq_2, method = 'spearman',adjust='BH')
cor.test.imp <- data.table(id=colnames(cor.test$r), r=t(cor.test$r),absr=abs(t(cor.test$r)),t=t(cor.test$t), p=t(cor.test$p))
cor.test.imp[,name:=meta$name[match(id, meta$id)]]

library("xlsx")
write.xlsx(cor.test.imp,'CELEGANS_lifespan_pathway_sig_test.xlsx')
#PWY-6608 p-value = 0.92
#PWY-8200 p-value = 0.1038
#PWY-5078 p-value = 0.6967

#growth
set.seed(666)
rf.path.growth<-VSURF(x=gapseq_2[,-1], y=data2$egg,  parallel=T)
rf.path.growth.imp <- data.table(nr=rf.path.growth$imp.mean.dec.ind, id=colnames(gapseq_2)[rf.path.growth$imp.mean.dec.ind], imp.mean=rf.path.growth$imp.mean.dec, imp.sd=rf.path.growth$imp.sd.dec) 
rf.path.growth.imp[,name:=meta$name[match(id, meta$id)]]
write.csv(rf.path.growth.imp,'N2.rf.growth.imp.csv')

##########################################################################################################
#Fig. 4c-f
library(readxl)
library(ggpubr)
library(tidyverse)
library(ggplot2)
library(ggalt)
library(paletteer)

HQB10 <-read_excel('input/Supplementary.xlsx', sheet='bmc_hqb10')


HQB10 %>% 
  rowwise() %>% 
  mutate(mean_value=mean(c(rep1,rep2,rep3)),
         std_error=plotrix::std.error(c(rep1,rep2,rep3))) %>% 
  ungroup() -> new.HQB10

new.HQB10$strain<-factor(new.HQB10$strain,levels=c('HQB-10','DMSO','CoQ','Heme',
                                                   'Acetyl-CoA', 'Pyruvate', 'Acetaldehyde',
                                                   'Oleate'))

#plot--hqb10
fig4c<- ggplot(data=new.HQB10,aes(x=Day,y=mean_value, color=strain))+
  geom_errorbar(aes(ymin=mean_value-std_error,
                    ymax=mean_value+std_error),linetype="solid",
                width=0.1,linewidth=0.5)+
  geom_point(size=1)+
  geom_line(linewidth =0.5)+
  #geom_xspline(spline_shape = -0.5)+
  #scale_linetype_manual(values=c("solid", "dashed",'dotted'))+
  geom_vline(xintercept=5, linewidth= 0.5,color="lightgrey",linetype = "dashed")+
  scale_color_manual(values = c('black',paletteer_d("basetheme::minimal")),
                     #labels=c("Vector","mRHIM2","mRHIM1",expression(paste("mZ",alpha,"2")),"WT"),
                     #breaks = c("A","F","E","D","B"),
                     name=NULL)+
  scale_x_continuous(limits = c(2.9,10.1),
                     expand = expansion(mult=c(0,0)),
                     breaks = seq(10)
  )+
  scale_y_continuous(limits = c(-0.5,60),
                     expand = expansion(mult=c(0,0)),
                     #breaks = seq(),
  )+
  coord_cartesian(clip = "off")+
  theme_classic()+
  theme(legend.text = element_text(hjust=0),
        # legend.position = c(0.1,0.4)
        axis.text = element_text(colour = "black")
  )+
  labs(x="Day",
       y=expression(paste("L4 proportion of L. marina (%)")))
fig4c
ggsave('fig4c.pdf', fig4c, width = 5,height = 2.2)



# hqb138
HQB138 <-read_excel('input/Supplementary.xlsx', sheet='bmc_hqb138')


HQB138 %>% 
  rowwise() %>% 
  mutate(mean_value=mean(c(rep1,rep2,rep3)),
         std_error=plotrix::std.error(c(rep1,rep2,rep3))) %>% 
  ungroup() -> new.HQB138

new.HQB138$strain<-factor(new.HQB138$strain,levels=c('HQB-138','DMSO','CoQ','Heme',
                                                     'Acetyl-CoA', 'Pyruvate', 'Acetaldehyde',
                                                     'Oleate'))

#plot--138
fig4d<- ggplot(data=new.HQB138,aes(x=Day,y=mean_value, color=strain))+
  geom_errorbar(aes(ymin=mean_value-std_error,
                    ymax=mean_value+std_error),linetype="solid",
                width=0.1,linewidth=0.5)+
  geom_point(size=1)+
  geom_line(linewidth =0.5)+
  #geom_xspline(spline_shape = -0.5)+
  #scale_linetype_manual(values=c("solid", "dashed",'dotted'))+
  geom_vline(xintercept=5, linewidth= 0.5,color="lightgrey",linetype = "dashed")+
  scale_color_manual(values = c('black', paletteer_d("basetheme::void")),
                     #labels=c("Vector","mRHIM2","mRHIM1",expression(paste("mZ",alpha,"2")),"WT"),
                     #breaks = c("A","F","E","D","B"),
                     name=NULL)+
  scale_x_continuous(limits = c(2.9,10.1),
                     expand = expansion(mult=c(0,0)),
                     breaks = seq(10)
  )+
  scale_y_continuous(limits = c(-0.5,80),
                     expand = expansion(mult=c(0,0)),
                     #breaks = seq(),
  )+
  coord_cartesian(clip = "off")+
  theme_classic()+
  theme(legend.text = element_text(hjust=0),
        # legend.position = c(0.1,0.4)
        axis.text = element_text(colour = "black")
  )+
  labs(x="Day",
       y=expression(paste("L4 proportion of L. marina (%)")))
fig4d
ggsave('fig4d.pdf', fig4d, width = 5,height = 2.2)

#fig 4e PLP
p5p <-read_excel('input/Supplementary.xlsx', sheet='bmc_p5p')

p5p %>% 
  rowwise() %>% 
  mutate(mean_value=mean(c(rep1,rep2,rep3)),
         std_error=plotrix::std.error(c(rep1,rep2,rep3))) %>% 
  ungroup() -> new.p5p

new.p5p$strain<-factor(new.p5p$strain,levels=c('OP50','OP50_P5P','HQB-5','HQB-5_P5P'))

#plot--pLp
fig4e<- ggplot(data=new.p5p,aes(x=Day,y=mean_value, linetype=type, col=color))+
  geom_errorbar(aes(ymin=mean_value-std_error,
                    ymax=mean_value+std_error),linetype="solid",
                width=0.1,linewidth=0.5)+
  geom_point(size=1)+
  geom_line(linewidth =0.5)+
  #geom_xspline(spline_shape = -0.5)+
  scale_linetype_manual(values=c("solid", "dashed",'dotted'))+
  geom_vline(xintercept=5, linewidth= 0.5,color="lightgrey",linetype = "dashed")+
  scale_color_manual(values = c('#4AC6AEFF','#E27069FF'),
                     #labels=c("Vector","mRHIM2","mRHIM1",expression(paste("mZ",alpha,"2")),"WT"),
                     #breaks = c("A","F","E","D","B"),
                     name=NULL)+
  scale_x_continuous(#limits = c(3,10),
    expand = expansion(mult=c(0,0)),
    breaks = seq(10)
  )+
  scale_y_continuous(limits = c(-0.5,80),
                     expand = expansion(mult=c(0,0)),
                     #breaks = seq(),
  )+
  coord_cartesian(clip = "off")+
  theme_classic()+
  theme(legend.text = element_text(hjust=0),
        # legend.position = c(0.1,0.4)
        axis.text = element_text(colour = "black")
  )+
  labs(x="Day",
       y=expression(paste("L4 proportion of L. marina (%)")))
fig4e
ggsave('fig4e.pdf', fig4e, width = 3.5,height = 2.2)

#p value
library(rstatix)

data.P5P<-read_excel('input/Supplementary.xlsx', sheet='p5p_pvalue')
data.P5P$Group<-factor(data.P5P$Group,levels=c('Control','P5P'))

#t test
op50.ttest<-data.P5P[1:6,]
H5.ttest<-data.P5P[6:12,]
op50.ttest%>%t_test(DAY5 ~ Group)
#DAY3 0.00288   **
#DAY4 0.00196     **
#DAY5 0.114
H5.ttest%>%t_test(DAY5 ~ Group)
#DAY3 0.391
#DAY4 0.0171  *
#DAY5 0.055

#########################################################################################
#each day
data <-read_excel('input/Supplementary.xlsx', sheet = "P_value")
data2=as.data.frame(lapply(data,as.numeric))
data2[,1]=data[,1]
data2[,11]=data[,11]

HQB_10<- data2 %>% filter(group %in% 'HQB-10')
HQB_138<- data2 %>% filter(group %in% 'HQB-138')

#p-value
library(DescTools)
set.seed(666)
#HQB-10
#day4
DunnettTest(x=HQB_10$d4, g=HQB_10$strain, control = 'HQB-10')
#$`HQB-10`
#                          diff     lwr.ci    upr.ci   pval     
#Acetaldehyde-HQB-10  0.1190476 -6.5214637  6.759559 1.0000      
#CoA-HQB-10           6.3095238 -0.3309875 12.950035 0.0647 .    
#Heme-HQB-10          7.7380952  1.0975839 14.378607 0.0207 *    
#Oleate-HQB-10       -3.2142857 -9.8547971  3.426226 0.5372      
#Pyruvate-HQB-10      0.1190476 -6.5214637  6.759559 1.0000      
#   

#day5
DunnettTest(x=HQB_10$d5, g=HQB_10$strain, control = 'HQB-10')
#$`HQB-10`
#                         diff     lwr.ci    upr.ci    pval    
#Acetaldehyde-HQB-10 -1.071429  -8.748962  6.606105 0.99352    
#CoA-HQB-10          17.976190  10.298657 25.653724 3.9e-05 ***
#Heme-HQB-10         16.547619   8.870085 24.225153 0.00011 ***
#Oleate-HQB-10       -9.642857 -17.320391 -1.965324 0.01283 *  
#Pyruvate-HQB-10      3.690476  -3.987057 11.368010 0.54337    
#

#day6
DunnettTest(x=HQB_10$d6, g=HQB_10$strain, control = 'HQB-10')
# $`HQB-10`
#                            diff      lwr.ci    upr.ci   pval    
# Acetaldehyde-HQB-10  -5.0000000 -18.3961782  8.396178 0.7446    
# CoA-HQB-10           19.2857143   5.8895361 32.681892 0.0047 ** 
# Heme-HQB-10          14.0476190   0.6514408 27.443797 0.0387 *  
# Oleate-HQB-10       -19.2857143 -32.6818925 -5.889536 0.0046 ** 
# Pyruvate-HQB-10      -0.7142857 -14.1104639 12.681892 0.9999    


#day7
DunnettTest(x=HQB_10$d7, g=HQB_10$strain, control = 'HQB-10')
# $`HQB-10`
#                           diff     lwr.ci     upr.ci    pval    
# Acetaldehyde-HQB-10  -9.761905 -24.241748   4.717938 0.25105    
# CoA-HQB-10           19.285714   4.805871  33.765557 0.00843 **  
# Heme-HQB-10          13.095238  -1.384605  27.575081 0.08216 .  
# Oleate-HQB-10       -26.428571 -40.908415 -11.948728 0.00061 ***
# Pyruvate-HQB-10      -5.952381 -20.432224   8.527462 0.67431    

#day8
DunnettTest(x=HQB_10$d8, g=HQB_10$strain, control = 'HQB-10')
# $`HQB-10`
# diff      lwr.ci      upr.ci    pval    
# Acetaldehyde-HQB-10 -15.595238 -29.6921087  -1.4983675  0.0284 *  
# CoA-HQB-10           14.404762   0.3078913  28.5016325  0.0446 *  
# Heme-HQB-10           9.166667  -4.9302040  23.2635373  0.2794    
# Oleate-HQB-10       -34.642857 -48.7397278 -20.5459865 4.1e-05 ***
# Pyruvate-HQB-10     -13.214286 -27.3111563   0.8825849  0.0694 .  


#day9
DunnettTest(x=HQB_10$d9, g=HQB_10$strain, control = 'HQB-10')
# $`HQB-10`
#                           diff      lwr.ci     upr.ci    pval    
# Acetaldehyde-HQB-10 -17.261905 -30.3609704  -4.162839  0.0092 ** 
# CoA-HQB-10           13.690476   0.5914105  26.789542  0.0393 *  
# Heme-HQB-10           9.880952  -3.2181133  22.980018  0.1734    
# Oleate-HQB-10       -38.214286 -51.3133514 -25.115220 7.7e-06 ***
# Pyruvate-HQB-10     -14.404762 -27.5038276  -1.305696  0.0296 *  
# 

#day10
DunnettTest(x=HQB_10$d10, g=HQB_10$strain, control = 'HQB-10')
# $`HQB-10`
#                           diff     lwr.ci     upr.ci    pval    
# Acetaldehyde-HQB-10 -19.404762 -31.571656  -7.237868  0.0021 ** 
# CoA-HQB-10           10.595238  -1.571656  22.762132  0.0974 .  
# Heme-HQB-10           7.738095  -4.428799  19.904989  0.2972    
# Oleate-HQB-10       -41.785714 -53.952608 -29.618820 9.1e-08 ***
# Pyruvate-HQB-10     -16.071429 -28.238323  -3.904534  0.0091 ** 

#HQB-138
#day4
DunnettTest(x=HQB_138$d4, g=HQB_138$strain, control = 'HQB-138')
#$`HQB-138`
#                          diff     lwr.ci     upr.ci   pval     
#Acetaldehyde-HQB-138 -2.738095 -10.109424  4.6332331 0.7479       
#CoA-HQB-138           5.833333  -1.537995 13.2046617 0.1449       
#Heme-HQB-138          9.642857   2.271529 17.0141855 0.0098 **    
#Oleate-HQB-138       -7.500000 -14.871328 -0.1286717 0.0458 *     
#Pyruvate-HQB-138      4.404762  -2.966566 11.7760903 0.3496       
#   

#day5
DunnettTest(x=HQB_138$d5, g=HQB_138$strain, control = 'HQB-138')
#$`HQB-138`
#                           diff     lwr.ci   upr.ci    pval    
#Acetaldehyde-HQB-138  16.071429   4.245980 27.89688 0.00726 ** 
#CoA-HQB-138           19.880952   8.055504 31.70640 0.00148 ** 
#Heme-HQB-138          15.595238   3.769789 27.42069 0.00920 ** 
#Oleate-HQB-138       -21.071429 -32.896877 -9.24598 0.00079 ***
#Pyruvate-HQB-138       8.928571  -2.896877 20.75402 0.17254    
#

#day6
DunnettTest(x=HQB_138$d6, g=HQB_138$strain, control = 'HQB-138')
# $`HQB-138`
#                            diff      lwr.ci     upr.ci    pval    
# Acetaldehyde-HQB-138   7.023810  -5.8186771  19.866296  0.4275    
# CoA-HQB-138           11.785714  -1.0567723  24.628201  0.0768 .  
# Heme-HQB-138          12.261905  -0.5805818  25.104391  0.0632 .  
# Oleate-HQB-138       -41.071429 -53.9139152 -28.228942 2.3e-07 ***
# Pyruvate-HQB-138      -4.404762 -17.2472485   8.437725  0.7980    


#day7
DunnettTest(x=HQB_138$d7, g=HQB_138$strain, control = 'HQB-138')
# $`HQB-138`
#                            diff     lwr.ci    upr.ci    pval    
# Acetaldehyde-HQB-138   6.666667  -8.461652  21.79499  0.6182    
# CoA-HQB-138           11.904762  -3.223557  27.03308  0.1482    
# Heme-HQB-138          16.190476   1.062157  31.31879  0.0343 *  
# Oleate-HQB-138       -45.714286 -60.842604 -30.58597 5.4e-06 ***
# Pyruvate-HQB-138      -2.857143 -17.985462  12.27118  0.9760    

#day8
DunnettTest(x=HQB_138$d8, g=HQB_138$strain, control = 'HQB-138')
# $`HQB-138`
#                             diff     lwr.ci    upr.ci    pval    
# Acetaldehyde-HQB-138   6.9047619  -7.941718  21.75124  0.5724    
# CoA-HQB-138           13.0952381  -1.751241  27.94172  0.0921 .  
# Heme-HQB-138          17.8571429   3.010663  32.70362  0.0168 *  
# Oleate-HQB-138       -47.8571429 -62.703622 -33.01066 2.1e-06 ***
# Pyruvate-HQB-138       0.2380952 -14.608384  15.08457  1.0000    


#day9
DunnettTest(x=HQB_138$d9, g=HQB_138$strain, control = 'HQB-138')
# $`HQB-138`
#                            diff      lwr.ci    upr.ci    pval    
# Acetaldehyde-HQB-138   5.714286 -10.1411307  21.56970  0.7677    
# CoA-HQB-138           11.904762  -3.9506545  27.76018  0.1765    
# Heme-HQB-138          16.666667   0.8112503  32.52208  0.0380 *  
# Oleate-HQB-138       -51.428571 -67.2839878 -35.57316 1.9e-07 ***
# Pyruvate-HQB-138      -1.428571 -17.2839878  14.42684  0.9992    
# 

#day10
DunnettTest(x=HQB_138$d10, g=HQB_138$strain, control = 'HQB-138')
# $`HQB-138`
#                            diff      lwr.ci    upr.ci   pval    
# Acetaldehyde-HQB-138   5.714286 -10.2697771  21.69835 0.7728    
# CoA-HQB-138           10.000000  -5.9840628  25.98406 0.3107    
# Heme-HQB-138          16.190476   0.2064134  32.17454 0.0466 *  
# Oleate-HQB-138       -54.285714 -70.2697771 -38.30165  5e-08 ***
# Pyruvate-HQB-138      -1.428571 -17.4126343  14.55549 0.9992    


#CoQ vs DMSO
#t test
library(rstatix)
coq <-read_excel('input/Supplementary.xlsx', sheet = "CoQ_DMSO")
coq_hqb10<-coq %>% filter(group %in% 'HQB-10')
coq_hqb138<-coq %>% filter(group %in% 'HQB-138')
#HQB-10
#day4
coq_hqb10 %>% t_test(d4 ~ strain)  #0.379
#day5
coq_hqb10 %>% t_test(d5 ~ strain)  #0.293
#day6
coq_hqb10 %>% t_test(d6 ~ strain)  #.273
#day7
coq_hqb10 %>% t_test(d7 ~ strain)  #0.388
#day8
coq_hqb10 %>% t_test(d8 ~ strain)  #.382
#day9
coq_hqb10 %>% t_test(d9 ~ strain)  #0.56
#day10
coq_hqb10 %>% t_test(d10 ~ strain)  #0.833

#HQB-138
#day4
coq_hqb138 %>% t_test(d4 ~ strain)  #0.0576
#day5
coq_hqb138 %>% t_test(d5 ~ strain)  #0.169
#day6
coq_hqb138 %>% t_test(d6 ~ strain)  #0.0778
#day7
coq_hqb138 %>% t_test(d7 ~ strain)  #0.00313
#day8
coq_hqb138 %>% t_test(d8 ~ strain)  #0.0125
#day9
coq_hqb138 %>% t_test(d9 ~ strain)  #0.0633
#day10
coq_hqb138 %>% t_test(d10 ~ strain)  #0.137



######################################################################################
##PLP   
#p value
library(rstatix)
data.P5P<-read_excel('input/Supplementary.xlsx', sheet='p5p_pvalue')
data.P5P$Group<-factor(data.P5P$Group,levels=c('Control','P5P'))

#t test
op50.ttest<-data.P5P[1:6,]
H5.ttest<-data.P5P[6:12,]
#op50
#day3
op50.ttest%>%t_test(DAY3 ~ Group)   #0.00288
#day4
op50.ttest%>%t_test(DAY4 ~ Group)   #0.00196
#day5
op50.ttest%>%t_test(DAY5 ~ Group)   #0.114
#day6
op50.ttest%>%t_test(DAY6 ~ Group)   #.0446
#day7
op50.ttest%>%t_test(DAY7 ~ Group)   #0.0396
#day8
op50.ttest%>%t_test(DAY8 ~ Group)   #0.0603
#day9
op50.ttest%>%t_test(DAY9 ~ Group)   #0.0603
#day10
op50.ttest%>%t_test(DAY10 ~ Group)   #0.0484

#HQB-5
#day3
H5.ttest%>%t_test(DAY3 ~ Group)   #0.391
#day4
H5.ttest%>%t_test(DAY4 ~ Group)   #0.0171
#day5
H5.ttest%>%t_test(DAY5 ~ Group)   #0.055
#day6
H5.ttest%>%t_test(DAY6 ~ Group)   #0.0932
#day7
H5.ttest%>%t_test(DAY7 ~ Group)   #0.154
#day8
H5.ttest%>%t_test(DAY8 ~ Group)   #0.219
#day9
H5.ttest%>%t_test(DAY9 ~ Group)   #0.219
#day10
H5.ttest%>%t_test(DAY10 ~ Group)   #0.219


############################################################################################
#Fig. 4f
library(tidyverse)
library(ggsurvfit)
library(gghighlight)
library(tidycmprsk)
library(paletteer) 
#F23-HQB5
data <- read_excel("input/Supplementary_N2.xlsx",sheet = 'f23_hqb5_ggsurvfit')

surcical1<-survfit2(Surv(time, status) ~ group, data = data) %>% 
  ggsurvfit(size = 1)+
  add_censor_mark()+
  #add_confidence_interval() +
  #add_risktable(theme=theme_test()+
  #                theme(axis.title = element_blank(),
  #                      axis.text.x = element_blank(),
  #                      axis.ticks.x=element_blank(),
  #                      axis.text.y = element_text(color="black",size=10)),
  #              combine_groups=F)+
  add_quantile(color ="grey80",size=0.8,linetype =5)+theme_classic()+
  labs(y = "Percentage Survival",x='Day')+
  #add_pvalue(caption = "Log-rank {p.value}",location  = "annotation", x = 10,prepend_p=T)+
  guides(colour =guide_legend(ncol=1))+
  scale_color_manual(values = c('#244579FF', '#C6242DFF'))+
  scale_fill_manual(values = c(#paletteer_d("ggthemes::hc_darkunica"),
    paletteer_d("basetheme::brutal")))+
  scale_x_continuous(breaks = seq(0, 30, by = 1))
surcical1

ggsave('sur2.pdf',surcical1,width = 5,height = 2)

#log rank p value
library(survival)
survdiff(Surv(time, status) ~ group, data = data)
# p= 3e-06

#Fig. S16
#N2-HQB5
data <- read_excel("input/Supplementary_N2.xlsx",sheet = 'N2-HQB5_sur')

surcical1<-survfit2(Surv(time, status) ~ group, data = data) %>% 
  ggsurvfit(size = 1)+
  add_censor_mark()+
  #add_confidence_interval() +
  #add_risktable(theme=theme_test()+
  #                theme(axis.title = element_blank(),
  #                      axis.text.x = element_blank(),
  #                      axis.ticks.x=element_blank(),
  #                      axis.text.y = element_text(color="black",size=10)),
  #              combine_groups=F)+
  add_quantile(color ="grey80",size=0.8,linetype =5)+theme_classic()+
  labs(y = "Percentage Survival",x='Day')+
  #add_pvalue(caption = "Log-rank {p.value}",location  = "annotation", x = 10,prepend_p=T)+
  guides(colour =guide_legend(ncol=1))+
  scale_color_manual(values = c('#244579FF', '#C6242DFF'))+
  scale_fill_manual(values = c(#paletteer_d("ggthemes::hc_darkunica"),
    paletteer_d("basetheme::brutal")))+
  scale_x_continuous(breaks = seq(0, 23, by = 1))+
  annotate("text", label = expression(paste('Log-rank ', italic('p'),' = 0.4')), x = 2.5, y = 0.08)
surcical1

ggsave('sur_n2-hqb5.pdf',surcical1,width = 7.5,height = 3)

#log rank p value
library(survival)
survdiff(Surv(time, status) ~ group, data = data)
# p= 0.4

#F23-OP50
data <- read_excel("input/Supplementary_N2.xlsx",sheet = 'F23-OP50_sur')

surcical1<-survfit2(Surv(time, status) ~ group, data = data) %>% 
  ggsurvfit(size = 1)+
  add_censor_mark()+
  #add_confidence_interval() +
  #add_risktable(theme=theme_test()+
  #                theme(axis.title = element_blank(),
  #                      axis.text.x = element_blank(),
  #                      axis.ticks.x=element_blank(),
  #                      axis.text.y = element_text(color="black",size=10)),
  #              combine_groups=F)+
  add_quantile(color ="grey80",size=0.8,linetype =5)+theme_classic()+
  labs(y = "Percentage Survival",x='Day')+
  #add_pvalue(caption = "Log-rank {p.value}",location  = "annotation", x = 10,prepend_p=T)+
  guides(colour =guide_legend(ncol=1))+
  scale_color_manual(values = c('#244579FF', '#C6242DFF'))+
  scale_fill_manual(values = c(#paletteer_d("ggthemes::hc_darkunica"),
    paletteer_d("basetheme::brutal")))+
  scale_x_continuous(breaks = seq(0, 21, by = 1))+
  annotate("text", label = expression(paste('Log-rank ', italic('p'),' = 0.005')), x = 2.5, y = 0.08)
surcical1

ggsave('sur_f23-op50.pdf',surcical1,width = 7.5,height = 3)

#log rank p value
library(survival)
survdiff(Surv(time, status) ~ group, data = data)
# p= p= 0.005

#N2-OP50
data <- read_excel("input/Supplementary_N2.xlsx",sheet = 'N2_OP50_sur')

surcical1<-survfit2(Surv(time, status) ~ group, data = data) %>% 
  ggsurvfit(size = 1)+
  add_censor_mark()+
  #add_confidence_interval() +
  #add_risktable(theme=theme_test()+
  #                theme(axis.title = element_blank(),
  #                      axis.text.x = element_blank(),
  #                      axis.ticks.x=element_blank(),
  #                      axis.text.y = element_text(color="black",size=10)),
  #              combine_groups=F)+
  add_quantile(color ="grey80",size=0.8,linetype =5)+theme_classic()+
  labs(y = "Percentage Survival",x='Day')+
  #add_pvalue(caption = "Log-rank {p.value}",location  = "annotation", x = 10,prepend_p=T)+
  guides(colour =guide_legend(ncol=1))+
  scale_color_manual(values = c('#244579FF', '#C6242DFF'))+
  scale_fill_manual(values = c(#paletteer_d("ggthemes::hc_darkunica"),
    paletteer_d("basetheme::brutal")))+
  scale_x_continuous(breaks = seq(0, 31, by = 1))+
  annotate("text", label = expression(paste('Log-rank ', italic('p'),' = 0.005')), x = 2.5, y = 0.08)
surcical1

ggsave('sur_n2-op50.pdf',surcical1,width = 7.5,height = 3)

#log rank p value
library(survival)
survdiff(Surv(time, status) ~ group, data = data)
# p= 0.005


###########################################################################################
#Fig. 5c-e
library(readxl)
#N2
#HQB-263
data.263<-read_excel("input/Supplementary_N2.xlsx",sheet = 'H263_eeglaying')
#HQB-118
data.118 <- read_excel("input/Supplementary_N2.xlsx",sheet = 'H118_eeglaying')
#HQB-325
data.325 <- read_excel("input/Supplementary_N2.xlsx",sheet = 'H325_eeglaying')
#HQB-123
data.123 <- read_excel("input/Supplementary_N2.xlsx",sheet = 'H123_egglaying')
###########################################################
#Dunnett test
library(DescTools)
#HQB-263
DunnettTest(x=data.263$`Egg laying time`, g=data.263$strain, control = 'HQB-263')
#$`HQB-263`
#diff      lwr.ci      upr.ci   pval    
#CoA-HQB-263      -0.96666667 -1.02251489 -0.91081844 <2e-16 ***
#  CoQ-HQB-263       0.03333333 -0.02251489  0.08918156 0.3358    
#DMSO-HQB-263      0.03333333 -0.02251489  0.08918156 0.3357    
#Heme-HQB-263      0.03333333 -0.02251489  0.08918156 0.3354    
#Pyruvate-HQB-263  0.03333333 -0.02251489  0.08918156 0.3354  

#HQB-118
DunnettTest(x=data.118$`Egg laying time`, g=data.118$strain, control = 'HQB-118')
#$`HQB-118`
#diff     lwr.ci     upr.ci    pval    
#CoA-HQB-118      -1.333333 -1.8918156 -0.7748511 7.8e-05 ***
#  CoQ-HQB-118       0.000000 -0.5584823  0.5584823  1.0000    
#DMSO-HQB-118      0.000000 -0.5584823  0.5584823  1.0000    
#Heme-HQB-118      0.000000 -0.5584823  0.5584823  1.0000    
#Pyruvate-HQB-118 -2.000000 -2.5584823 -1.4415177 1.7e-07 ***

#HQB-325
DunnettTest(x=data.325$`Egg laying time`, g=data.325$strain, control = 'HQB-325')
#$`HQB-325`
#diff      lwr.ci     upr.ci   pval    
#CoA-HQB-325       0.03333333 -0.01861673  0.0852834 0.2976    
#CoQ-HQB-325      -1.96666667 -2.01861673 -1.9147166 <2e-16 ***
#  DMSO-HQB-325     -1.96666667 -2.01861673 -1.9147166 <2e-16 ***
#  Heme-HQB-325      0.03333333 -0.01861673  0.0852834 0.2974    
#Oleate-HQB-325   -3.96666667 -4.01861673 -3.9147166 <2e-16 ***
#  Pyruvate-HQB-325  0.03333333 -0.01861673  0.0852834 0.2974        

#HQB-123
DunnettTest(x=data.123$`Egg laying time`, g=data.123$strain, control = 'HQB-123')

#$`HQB-123`
#diff      lwr.ci      upr.ci   pval    
#CoA-HQB-123      -13.6666667 -14.2251489 -13.1081844 <2e-16 ***
#  CoQ-HQB-123        0.3333333  -0.2251489   0.8918156 0.3357    
#DMSO-HQB-123       0.3333333  -0.2251489   0.8918156 0.3357    
#Heme-HQB-123       0.3333333  -0.2251489   0.8918156 0.3357    
#Pyruvate-HQB-123   0.3333333  -0.2251489   0.8918156 0.3358   

#===============================
library(ggplot2)
library(ggpubr)
library(paletteer) 

#HQB-263
growth.263<-ggboxplot(
  data.263, x = "strain", y = "Egg laying time",
  color = "strain", palette = c('black',paletteer_d("basetheme::minimal")),
  #ylim = c(0,0.8),
  ylab='Egg laying time (h)',xlab = '',size=1,ggtheme=theme_classic()
)+theme(axis.text = element_text(color="black"))+
  #stat_pvalue_manual(p.d5, label = "p.adj.signif")+
  rotate_x_text(45,size=10)#+geom_text(aes(label = p.growth),size=3)
growth.263

#HQB-118
growth.118<-ggboxplot(
  data.118, x = "strain", y = "Egg laying time",
  color = "strain", palette = c('black',paletteer_d("basetheme::minimal")),
  #ylim = c(0,0.8),
  ylab='Egg laying time (h)',xlab = '',size=1,ggtheme=theme_classic()
)+theme(axis.text = element_text(color="black"))+
  #stat_pvalue_manual(p.d5, label = "p.adj.signif")+
  rotate_x_text(45,size=10)#+geom_text(aes(label = p.growth),size=3)
growth.118

#HQB-123
growth.123<-ggboxplot(
  data.123, x = "strain", y = "Egg laying time",
  color = "strain", palette = c('black',paletteer_d("basetheme::minimal")),
  #ylim = c(0,0.8),
  ylab='Egg laying time (h)',xlab = '',size=1,ggtheme=theme_classic()
)+theme(axis.text = element_text(color="black"))+
  #stat_pvalue_manual(p.d5, label = "p.adj.signif")+
  rotate_x_text(45,size=10)#+geom_text(aes(label = p.growth),size=3)
growth.123

merge1<-ggarrange(growth.263, growth.118,growth.123, common.legend = T, legend = 'none',
                  labels = c('a','b', 'c','d'))
merge1
ggsave('N2-adding4.pdf',merge1,height = 6,width = 7)


##########################################################################################
#Fig. S14 and Fig. S38
#OP50 dosage supplementary
library(readxl)
library(ggpubr)
library(tidyverse)
library(ggplot2)
library(ggalt)
library(paletteer)
library(DescTools)
set.seed(666)
#F23
#Fig. S14g
data<-read_excel('input/OP50_dosage.xlsx', sheet= 'F23_D5')

DunnettTest(x=data$D5, g=data$strain, control = 'OP50')
# $OP50
#                             diff     lwr.ci       upr.ci    pval    
# Acetaldehyde-1-OP50  -0.04444444 -0.1919551  0.103066231  0.9910    
# Acetaldehyde-10-OP50  0.01481481 -0.1326959  0.162325491  1.0000    
# Acetaldehyde-5-OP50  -0.05185185 -0.1993625  0.095658824  0.9666    
# acetyl-CoA-0.1-OP50  -0.06666667 -0.2141773  0.080844009  0.8350    
# acetyl-CoA-1-OP50     0.01481481 -0.1326959  0.162325491  1.0000    
# acetyl-CoA-5-OP50    -0.14074074 -0.2882514  0.006769935  0.0694 .  
# CoQ-10-OP50          -0.11851852 -0.2660292  0.028992157  0.1827    
# CoQ-100-OP50         -0.02962963 -0.1771403  0.117881046  0.9999    
# CoQ-50-OP50          -0.01481481 -0.1623255  0.132695861  1.0000    
# Heme-10-OP50         -0.01481481 -0.1623255  0.132695861  1.0000    
# Heme-100-OP50        -0.08888889 -0.2363996  0.058621787  0.5108    
# Heme-500-OP50        -0.20000000 -0.3475107 -0.052489324  0.0030 ** 
# Oleate-1-OP50        -0.40000000 -0.5475107 -0.252489324 1.9e-09 ***
# Oleate-50-OP50       -0.38518519 -0.5326959 -0.237674509 5.3e-09 ***
# Oleate-500-OP50      -0.40740741 -0.5549181 -0.259896732 5.5e-10 ***
# Pyruvate-1-OP50      -0.06666667 -0.2141773  0.080844009  0.8347    
# Pyruvate-2.5-OP50    -0.03703704 -0.1845477  0.110473639  0.9986    
# Pyruvate-5-OP50      -0.40000000 -0.5475107 -0.252489324 1.1e-08 ***

data$group<-factor(data$group,levels=c('OP50','DMSO','CoQ','Heme', 'acetyl-CoA',
                                       'Acetaldehyde','Pyruvate','Oleate'))

plot<-ggboxplot(
  data, x = "strain", y = "D5*100",fill = 'group',
  #color = "Group", 
  palette =paletteer_d("basetheme::minimal"),
  #ylim = c(0,0.8),
  ylab='Proportion of L4 in day 5 (%)',
  xlab = '',
  size=0.5,ggtheme=theme_classic()
)+
  theme(axis.text = element_text(color="black"))+
  #stat_pvalue_manual(p.d5, label = "p.adj.signif")+
  rotate_x_text(45,size=10)
plot
ggsave('op50_various.pdf',plot,width = 8, height = 4)

###################################################################################################
#N2
#Fig, S38
data2<-read_excel('input/OP50_dosage.xlsx', sheet= 'N2')
DunnettTest(x=data2$Time, g=data2$strain, control = 'OP50')
#$OP50
#                          diff     lwr.ci     upr.ci    pval    
#acetyl-CoA-0.1-OP50 -0.6666667 -1.4209322 0.08759885 0.10404    
#acetyl-CoA-1-OP50    0.1666667 -0.5875989 0.92093219 0.99715    
#acetyl-CoA-5-OP50    1.1666667  0.4124011 1.92093219 0.00097 ***
#CoQ-1-OP50           0.1666667 -0.5875989 0.92093219 0.99715     
#CoQ-10-OP50          0.1666667 -0.5875989 0.92093219 0.99715    
#CoQ-50-OP50          0.1666667 -0.5875989 0.92093219 0.99715    
#DMSO-0.25-OP50       0.1666667 -0.5875989 0.92093219 0.99715    
#DMSO-2.5-OP50        0.1666667 -0.5875989 0.92093219 0.99716    
#DMSO-25-OP50         0.1666667 -0.5875989 0.92093219 0.99715    
#Heme-10-OP50         0.1666667 -0.5875989 0.92093219 0.99716    
#Heme-100-OP50        0.1666667 -0.5875989 0.92093219 0.99715    
#Heme-500-OP50        1.3333333  0.5790678 2.08759885 0.00018 ***


plot2<-ggboxplot(
  data2, x = "strain", y = "Time",fill = 'group',
  #color = "Group", 
  palette =paletteer_d("basetheme::minimal"),
  #ylim = c(0,0.8),
  ylab='Proportion of L4 in day 5 (%)',
  xlab = '',
  #size=0.5,
  ggtheme=theme_classic()
)+
  theme(axis.text = element_text(color="black"))+
  #stat_pvalue_manual(p.d5, label = "p.adj.signif")+
  rotate_x_text(45,size=10)
plot2
ggsave('op50_N2_dosage.pdf',plot2,width = 6, height = 4)

####################################################################################
#Fig. S14a-f
#curve plot
#1
library(ggalt)
library(tidyverse)
library(paletteer)
library(readxl)
#coq
curve_coq <-read_excel('input/OP50_dosage.xlsx', sheet='cur_coq')


curve_coq %>% 
  rowwise() %>% 
  mutate(mean_value=mean(c(rep1,rep2,rep3)),
         std_error=plotrix::std.error(c(rep1,rep2,rep3))) %>% 
  ungroup() -> new.curve_coq

new.curve_coq$strain<-factor(new.curve_coq$strain,levels=c('DMSO-0.25','DMSO-2.5','DMSO-25',
                                                           'CoQ-1','CoQ-10','CoQ-50'))
#plot-
cur_coq<- ggplot(data=new.curve_coq,aes(x=Day,y=mean_value,color=strain
))+
  geom_errorbar(aes(ymin=mean_value-std_error,
                    ymax=mean_value+std_error),linetype="solid",
                width=0.1,linewidth=0.2)+
  geom_point(size=0.8)+
  geom_line(linewidth = 0.3)+
  #geom_xspline(spline_shape = -0.5,size =6)+
  scale_linetype_manual(values=c( 'solid',"dotted"))+
  geom_vline(xintercept=5, linewidth= 0.5,color="lightgrey",linetype = "dashed")+
  scale_color_manual(values = c(paletteer_d("awtools::mpalette")),
                     name=NULL)+
  scale_x_continuous(#limits = c(3,10),
    expand = expansion(mult=c(0,0)),
    breaks = seq(10)
  )+
  scale_y_continuous(limits = c(-0.5,90),
                     expand = expansion(mult=c(0,0)),
                     #breaks = seq(),
  )+
  coord_cartesian(clip = "off")+
  theme_classic()+
  theme(legend.text = element_text(hjust=0),
        legend.position = c(0.8,0.35),
        axis.text = element_text(colour = "black")
  )+
  labs(x="Day",
       y=expression(paste("L4 proportion of L. marina (%)")))
cur_coq

#heme
curve_heme <-read_excel('input/OP50_dosage.xlsx', sheet='cur_heme')


curve_heme %>% 
  rowwise() %>% 
  mutate(mean_value=mean(c(rep1,rep2,rep3)),
         std_error=plotrix::std.error(c(rep1,rep2,rep3))) %>% 
  ungroup() -> new.curve_heme

new.curve_heme$strain<-factor(new.curve_heme$strain,levels=c('OP50','Heme-10','Heme-100','Heme-500'))
#plot-
cur_heme<- ggplot(data=new.curve_heme,aes(x=Day,y=mean_value,color=strain
))+
  geom_errorbar(aes(ymin=mean_value-std_error,
                    ymax=mean_value+std_error),linetype="solid",
                width=0.1,linewidth=0.2)+
  geom_point(size=0.8)+
  geom_line(linewidth = 0.3)+
  #geom_xspline(spline_shape = -0.5,size =6)+
  scale_linetype_manual(values=c( 'solid',"dotted"))+
  geom_vline(xintercept=5, linewidth= 0.5,color="lightgrey",linetype = "dashed")+
  scale_color_manual(values = c('black',paletteer_d("awtools::mpalette")),
                     name=NULL)+
  scale_x_continuous(#limits = c(3,10),
    expand = expansion(mult=c(0,0)),
    breaks = seq(10)
  )+
  scale_y_continuous(limits = c(-0.5,90),
                     expand = expansion(mult=c(0,0)),
                     #breaks = seq(),
  )+
  coord_cartesian(clip = "off")+
  theme_classic()+
  theme(legend.text = element_text(hjust=0),
        legend.position = c(0.8,0.3),
        axis.text = element_text(colour = "black")
  )+
  labs(x="Day",
       y=expression(paste("L4 proportion of L. marina (%)")))
cur_heme


#coa
curve_coa <-read_excel('input/OP50_dosage.xlsx', sheet='cur_coa')


curve_coa %>% 
  rowwise() %>% 
  mutate(mean_value=mean(c(rep1,rep2,rep3)),
         std_error=plotrix::std.error(c(rep1,rep2,rep3))) %>% 
  ungroup() -> new.curve_coa

new.curve_coa$strain<-factor(new.curve_coa$strain,levels=c('OP50','Acetyl-CoA-0.1','Acetyl-CoA-1','Acetyl-CoA-5'))
#plot-
cur_coa<- ggplot(data=new.curve_coa,aes(x=Day,y=mean_value,color=strain
))+
  geom_errorbar(aes(ymin=mean_value-std_error,
                    ymax=mean_value+std_error),linetype="solid",
                width=0.1,linewidth=0.2)+
  geom_point(size=0.8)+
  geom_line(linewidth = 0.3)+
  #geom_xspline(spline_shape = -0.5,size =6)+
  scale_linetype_manual(values=c( 'solid',"dotted"))+
  geom_vline(xintercept=5, linewidth= 0.5,color="lightgrey",linetype = "dashed")+
  scale_color_manual(values = c('black',paletteer_d("awtools::mpalette")),
                     name=NULL)+
  scale_x_continuous(#limits = c(3,10),
    expand = expansion(mult=c(0,0)),
    breaks = seq(10)
  )+
  scale_y_continuous(limits = c(-0.5,90),
                     expand = expansion(mult=c(0,0)),
                     #breaks = seq(),
  )+
  coord_cartesian(clip = "off")+
  theme_classic()+
  theme(legend.text = element_text(hjust=0),
        legend.position = c(0.8,0.3),
        axis.text = element_text(colour = "black")
  )+
  labs(x="Day",
       y=expression(paste("L4 proportion of L. marina (%)")))
cur_coa

#Acetaldehyde
curve_Acetaldehyde <-read_excel('input/OP50_dosage.xlsx', sheet='cur_Acetaldehyde')


curve_Acetaldehyde %>% 
  rowwise() %>% 
  mutate(mean_value=mean(c(rep1,rep2,rep3)),
         std_error=plotrix::std.error(c(rep1,rep2,rep3))) %>% 
  ungroup() -> new.curve_Acetaldehyde

new.curve_Acetaldehyde$strain<-factor(new.curve_Acetaldehyde$strain,levels=c('OP50','Acetaldehyde-1','Acetaldehyde-5','Acetaldehyde-10'))

#plot-
cur_Acetaldehyde<- ggplot(data=new.curve_Acetaldehyde,aes(x=Day,y=mean_value,color=strain
))+
  geom_errorbar(aes(ymin=mean_value-std_error,
                    ymax=mean_value+std_error),linetype="solid",
                width=0.1,linewidth=0.2)+
  geom_point(size=0.8)+
  geom_line(linewidth = 0.3)+
  #geom_xspline(spline_shape = -0.5,size =6)+
  scale_linetype_manual(values=c( 'solid',"dotted"))+
  geom_vline(xintercept=5, linewidth= 0.5,color="lightgrey",linetype = "dashed")+
  scale_color_manual(values = c('black',paletteer_d("awtools::mpalette")),
                     name=NULL)+
  scale_x_continuous(#limits = c(3,10),
    expand = expansion(mult=c(0,0)),
    breaks = seq(10)
  )+
  scale_y_continuous(limits = c(-0.5,90),
                     expand = expansion(mult=c(0,0)),
                     #breaks = seq(),
  )+
  coord_cartesian(clip = "off")+
  theme_classic()+
  theme(legend.text = element_text(hjust=0),
        legend.position = c(0.8,0.3),
        axis.text = element_text(colour = "black")
  )+
  labs(x="Day",
       y=expression(paste("L4 proportion of L. marina (%)")))
cur_Acetaldehyde


#Pyruvate
curve_Pyruvate <-read_excel('input/OP50_dosage.xlsx', sheet='cur_Pyruvate')


curve_Pyruvate %>% 
  rowwise() %>% 
  mutate(mean_value=mean(c(rep1,rep2,rep3)),
         std_error=plotrix::std.error(c(rep1,rep2,rep3))) %>% 
  ungroup() -> new.curve_Pyruvate

new.curve_Pyruvate$strain<-factor(new.curve_Pyruvate$strain,levels=c('OP50','Pyruvate-1','Pyruvate-2.5','Pyruvate-5'))

#plot-
cur_Pyruvate<- ggplot(data=new.curve_Pyruvate,aes(x=Day,y=mean_value,color=strain
))+
  geom_errorbar(aes(ymin=mean_value-std_error,
                    ymax=mean_value+std_error),linetype="solid",
                width=0.1,linewidth=0.2)+
  geom_point(size=0.8)+
  geom_line(linewidth = 0.3)+
  #geom_xspline(spline_shape = -0.5,size =6)+
  scale_linetype_manual(values=c( 'solid',"dotted"))+
  geom_vline(xintercept=5, linewidth= 0.5,color="lightgrey",linetype = "dashed")+
  scale_color_manual(values = c('black',paletteer_d("awtools::mpalette")),
                     name=NULL)+
  scale_x_continuous(#limits = c(3,10),
    expand = expansion(mult=c(0,0)),
    breaks = seq(10)
  )+
  scale_y_continuous(limits = c(-0.5,90),
                     expand = expansion(mult=c(0,0)),
                     #breaks = seq(),
  )+
  coord_cartesian(clip = "off")+
  theme_classic()+
  theme(legend.text = element_text(hjust=0),
        legend.position = c(0.8,0.3),
        axis.text = element_text(colour = "black")
  )+
  labs(x="Day",
       y=expression(paste("L4 proportion of L. marina (%)")))
cur_Pyruvate

#oleate
curve_oleate <-read_excel('input/OP50_dosage.xlsx', sheet='cur_oleate')


curve_oleate %>% 
  rowwise() %>% 
  mutate(mean_value=mean(c(rep1,rep2,rep3)),
         std_error=plotrix::std.error(c(rep1,rep2,rep3))) %>% 
  ungroup() -> new.curve_oleate

new.curve_oleate$strain<-factor(new.curve_oleate$strain,levels=c('OP50','Oleate-1','Oleate-50','Oleate-500'))

#plot-
cur_oleate<- ggplot(data=new.curve_oleate,aes(x=Day,y=mean_value,color=strain
))+
  geom_errorbar(aes(ymin=mean_value-std_error,
                    ymax=mean_value+std_error),linetype="solid",
                width=0.1,linewidth=0.2)+
  geom_point(size=0.8)+
  geom_line(linewidth = 0.3)+
  #geom_xspline(spline_shape = -0.5,size =6)+
  scale_linetype_manual(values=c( 'solid',"dotted"))+
  geom_vline(xintercept=5, linewidth= 0.5,color="lightgrey",linetype = "dashed")+
  scale_color_manual(values = c('black',paletteer_d("awtools::mpalette")),
                     name=NULL)+
  scale_x_continuous(#limits = c(3,10),
    expand = expansion(mult=c(0,0)),
    breaks = seq(10)
  )+
  scale_y_continuous(limits = c(-0.5,90),
                     expand = expansion(mult=c(0,0)),
                     #breaks = seq(),
  )+
  coord_cartesian(clip = "off")+
  theme_classic()+
  theme(legend.text = element_text(hjust=0),
        legend.position = c(0.8,0.3),
        axis.text = element_text(colour = "black")
  )+
  labs(x="Day",
       y=expression(paste("L4 proportion of L. marina (%)")))
cur_oleate


#merge
library(ggpubr)
merge1<-ggarrange(cur_coq,cur_heme,cur_coa,
                  ncol = 3,labels = c('a','b','c'))
merge1

merge2<-ggarrange(cur_Acetaldehyde,cur_Pyruvate,cur_oleate,
                  ncol = 3,labels = c('d','e','f'))
merge2

merge12<-ggarrange(merge1,merge2,ncol = 1)
merge12

ggsave('merge_op50_adding.pdf',merge12,height = 5,width = 9)

#########################################################################################################
#Fig. S12 and Fig. S13 and Fig. S19
# The impacts of OP50 mutants on L. marina development with or without dietary supplementation of single bacteria metabolite
library(readxl)
library(rstatix)

#Fig. S13
#OP50

#ISOLATE-1
#Total
OP50_f23<-read_excel("input/mutants.xlsx",sheet = 'OP50_1_plot')
#day3
OP50_D3<-OP50_f23[,c('strain','REP','D3','color','line')]
#only mutants
OP50_D3_mutants<-OP50_D3 %>% filter(strain %in% c('WT','ΔubiG','ΔhemF','ΔadhE'))

#Dunnett test
library(DescTools)
#day3
DunnettTest(x=OP50_D3_mutants$D3, g=OP50_D3_mutants$strain, control = 'WT')
#$WT
# diff     lwr.ci     upr.ci    pval    
# ΔadhE-WT -0.34 -0.4222174 -0.2577826 1.7e-07 ***
#   ΔhemF-WT -0.33 -0.4069073 -0.2530927 2.7e-06 ***
#   ΔubiG-WT -0.40 -0.4822174 -0.3177826 1.4e-09 ***

#

#===============================
library(ggplot2)
library(ggpubr)
library(paletteer) 

OP50_D3$strain<-factor(OP50_D3$strain,levels=c('WT','WT+CoQ','WT+heme','WT+CoA',
                                               'ΔubiG','ΔubiG+CoQ',
                                               'ΔhemF','ΔhemF+heme',
                                               'ΔadhE','ΔadhE+CoA'))
OP50_D3$color<-factor(OP50_D3$color,levels=c('WT','ΔubiG','ΔhemF','ΔadhE'))
OP50_D3$line<-factor(OP50_D3$line,levels=c('solid','dashed'))
#day3
OP50_day3<-ggboxplot(
  OP50_D3, x = "strain", y = "D3",
  color = "color", palette = c('black', '#DA2E20FF','#3870C2FF','#66C84DFF'), linetype = 'line',
  #ylim = c(0,0.8),
  ylab='Proportion of L4 in day 3 (%)',xlab = '',size=0.5,ggtheme=theme_classic()
)+ scale_linetype_manual(values=c( 'solid',"dotted"))+
  #stat_pvalue_manual(p.d5, label = "p.adj.signif")+
  rotate_x_text(45,size=10)+theme(axis.text =  element_text(color="black"))
OP50_day3
ggsave('op50/OP50_day3.pdf',OP50_day3,width = 6,height = 4)

#================================================================
#day4
OP50_D4<-OP50_f23[,c('strain','REP','D4','color','line')]
#only mutants
OP50_D4_mutants<-OP50_D4 %>% filter(strain %in% c('WT','ΔubiG','ΔhemF','ΔadhE'))

#Dunnett test
#day4
DunnettTest(x=OP50_D4_mutants$D4, g=OP50_D4_mutants$strain, control = 'WT')
#$WT
# diff    lwr.ci     upr.ci    pval    
# ΔadhE-WT -27.33333 -48.04302  -6.623649 0.01235 *  
# ΔhemF-WT -20.66667 -40.03880  -1.294531 0.03744 *  
# ΔubiG-WT -45.33333 -66.04302 -24.623649 0.00045 ***

#===============================
#plot
OP50_D4$strain<-factor(OP50_D4$strain,levels=c('WT','WT+CoQ','WT+heme','WT+CoA',
                                               'ΔubiG','ΔubiG+CoQ',
                                               'ΔhemF','ΔhemF+heme',
                                               'ΔadhE','ΔadhE+CoA'))
OP50_D4$color<-factor(OP50_D4$color,levels=c('WT','ΔubiG','ΔhemF','ΔadhE'))
OP50_D4$line<-factor(OP50_D4$line,levels=c('solid','dashed'))
#day4
OP50_day4<-ggboxplot(
  OP50_D4, x = "strain", y = "D4",
  color = "color", palette = c('black', '#DA2E20FF','#3870C2FF','#66C84DFF'), linetype = 'line',
  #ylim = c(0,0.8),
  ylab='Proportion of L4 in day 4 (%)',xlab = '',size=0.5,ggtheme=theme_classic()
)+ scale_linetype_manual(values=c( 'solid',"dotted"))+
  #stat_pvalue_manual(p.d5, label = "p.adj.signif")+
  rotate_x_text(45,size=10)+theme(axis.text =  element_text(color="black"))
OP50_day4

#================================================================
#day5
OP50_D5<-OP50_f23[,c('strain','REP','D5','color','line')]
#only mutants
OP50_D5_mutants<-OP50_D5 %>% filter(strain %in% c('WT','ΔubiG','ΔhemF','ΔadhE'))

#Dunnett test
#day5
DunnettTest(x=OP50_D5_mutants$D5, g=OP50_D5_mutants$strain, control = 'WT')
#$WT
# diff    lwr.ci    upr.ci   pval    
# ΔadhE-WT -18.00000 -43.12169  7.121687 0.1708    
# ΔhemF-WT -12.83333 -36.33252 10.665853 0.3385    
# ΔubiG-WT -32.66667 -57.78835 -7.544980 0.0134 * 

#===============================
#plot
OP50_D5$strain<-factor(OP50_D5$strain,levels=c('WT','WT+CoQ','WT+heme','WT+CoA',
                                               'ΔubiG','ΔubiG+CoQ',
                                               'ΔhemF','ΔhemF+heme',
                                               'ΔadhE','ΔadhE+CoA'))
OP50_D5$color<-factor(OP50_D5$color,levels=c('WT','ΔubiG','ΔhemF','ΔadhE'))
OP50_D5$line<-factor(OP50_D5$line,levels=c('solid','dashed'))
#day5
OP50_day5<-ggboxplot(
  OP50_D5, x = "strain", y = "D5",
  color = "color", palette = c('black', '#DA2E20FF','#3870C2FF','#66C84DFF'), linetype = 'line',
  #ylim = c(0,0.8),
  ylab='Proportion of L4 in day 5 (%)',xlab = '',size=0.5,ggtheme=theme_classic()
)+ scale_linetype_manual(values=c( 'solid',"dotted"))+
  #stat_pvalue_manual(p.d5, label = "p.adj.signif")+
  rotate_x_text(45,size=10)+theme(axis.text =  element_text(color="black"))
OP50_day5

#================================================================
#day6
OP50_D6<-OP50_f23[,c('strain','REP','D6','color','line')]
#only mutants
OP50_D6_mutants<-OP50_D6 %>% filter(strain %in% c('WT','ΔubiG','ΔhemF','ΔadhE'))

#Dunnett test
#day6
DunnettTest(x=OP50_D6_mutants$D6, g=OP50_D6_mutants$strain, control = 'WT')
#$WT
# diff    lwr.ci    upr.ci   pval    
# ΔadhE-WT -13.33333 -34.19720  7.530537 0.2351    
# ΔhemF-WT -10.16667 -29.68303  9.349697 0.3716    
# ΔubiG-WT -25.33333 -46.19720 -4.469463 0.0192 *  
#ubiG t test
OP50_D6_ubiG<-OP50_D6 %>% filter(strain %in% c('ΔubiG','ΔubiG+CoQ'))
OP50_D6_ubiG %>% t_test(D6 ~ line)  # p=0.946
#===============================
#plot
OP50_D6$strain<-factor(OP50_D6$strain,levels=c('WT','WT+CoQ','WT+heme','WT+CoA',
                                               'ΔubiG','ΔubiG+CoQ',
                                               'ΔhemF','ΔhemF+heme',
                                               'ΔadhE','ΔadhE+CoA'))
OP50_D6$color<-factor(OP50_D6$color,levels=c('WT','ΔubiG','ΔhemF','ΔadhE'))
OP50_D6$line<-factor(OP50_D6$line,levels=c('solid','dashed'))

OP50_day6<-ggboxplot(
  OP50_D6, x = "strain", y = "D6",
  color = "color", palette = c('black', '#DA2E20FF','#3870C2FF','#66C84DFF'), linetype = 'line',
  #ylim = c(0,0.8),
  ylab='Proportion of L4 in day 6 (%)',xlab = '',size=0.5,ggtheme=theme_classic()
)+ scale_linetype_manual(values=c( 'solid',"dotted"))+
  #stat_pvalue_manual(p.d5, label = "p.adj.signif")+
  rotate_x_text(45,size=10)+theme(axis.text =  element_text(color="black"))
OP50_day6

#================================================================
#day7
OP50_D7<-OP50_f23[,c('strain','REP','D7','color','line')]
#only mutants
OP50_D7_mutants<-OP50_D7 %>% filter(strain %in% c('WT','ΔubiG','ΔhemF','ΔadhE'))

#Dunnett test
#day7
DunnettTest(x=OP50_D7_mutants$D7, g=OP50_D7_mutants$strain, control = 'WT')
#$WT
# diff    lwr.ci    upr.ci   pval    
# ΔadhE-WT  -8.666667 -29.38104 12.047707 0.5288    
# ΔhemF-WT  -7.166667 -26.54319 12.209855 0.6130    
# ΔubiG-WT -18.000000 -38.71437  2.714373 0.0892 . 
#ubiG t test
OP50_D7_ubiG<-OP50_D7 %>% filter(strain %in% c('ΔubiG','ΔubiG+CoQ'))
OP50_D7_ubiG %>% t_test(D7 ~ line)  # p= 0.71
#===============================
#plot
OP50_D7$strain<-factor(OP50_D7$strain,levels=c('WT','WT+CoQ','WT+heme','WT+CoA',
                                               'ΔubiG','ΔubiG+CoQ',
                                               'ΔhemF','ΔhemF+heme',
                                               'ΔadhE','ΔadhE+CoA'))
OP50_D7$color<-factor(OP50_D7$color,levels=c('WT','ΔubiG','ΔhemF','ΔadhE'))
OP50_D7$line<-factor(OP50_D7$line,levels=c('solid','dashed'))

OP50_day7<-ggboxplot(
  OP50_D7, x = "strain", y = "D7",
  color = "color", palette = c('black', '#DA2E20FF','#3870C2FF','#66C84DFF'), linetype = 'line',
  #ylim = c(0,0.8),
  ylab='Proportion of L4 in day 7 (%)',xlab = '',size=0.5,ggtheme=theme_classic()
)+ scale_linetype_manual(values=c( 'solid',"dotted"))+
  #stat_pvalue_manual(p.d5, label = "p.adj.signif")+
  rotate_x_text(45,size=10)+theme(axis.text =  element_text(color="black"))
OP50_day7

#================================================================
#day8
OP50_D8<-OP50_f23[,c('strain','REP','D8','color','line')]
#only mutants
OP50_D8_mutants<-OP50_D8 %>% filter(strain %in% c('WT','ΔubiG','ΔhemF','ΔadhE'))

#Dunnett test
#day8
DunnettTest(x=OP50_D8_mutants$D8, g=OP50_D8_mutants$strain, control = 'WT')
#$WT
# diff    lwr.ci    upr.ci   pval    
# ΔadhE-WT  -8.666667 -28.69450 11.361167 0.5049    
# ΔhemF-WT  -6.166667 -24.90099 12.567656 0.6855    
# ΔubiG-WT -14.000000 -34.02783  6.027833 0.1842 
#ubiG t test
OP50_D8_ubiG<-OP50_D8 %>% filter(strain %in% c('ΔubiG','ΔubiG+CoQ'))
OP50_D8_ubiG %>% t_test(D8 ~ line)  # p= 0.69

#===============================
#plot
OP50_D8$strain<-factor(OP50_D8$strain,levels=c('WT','WT+CoQ','WT+heme','WT+CoA',
                                               'ΔubiG','ΔubiG+CoQ',
                                               'ΔhemF','ΔhemF+heme',
                                               'ΔadhE','ΔadhE+CoA'))
OP50_D8$color<-factor(OP50_D8$color,levels=c('WT','ΔubiG','ΔhemF','ΔadhE'))
OP50_D8$line<-factor(OP50_D8$line,levels=c('solid','dashed'))

OP50_day8<-ggboxplot(
  OP50_D8, x = "strain", y = "D8",
  color = "color", palette = c('black', '#DA2E20FF','#3870C2FF','#66C84DFF'), linetype = 'line',
  #ylim = c(0,0.8),
  ylab='Proportion of L4 in day 8 (%)',xlab = '',size=0.5,ggtheme=theme_classic()
)+ scale_linetype_manual(values=c( 'solid',"dotted"))+
  #stat_pvalue_manual(p.d5, label = "p.adj.signif")+
  rotate_x_text(45,size=10)+theme(axis.text =  element_text(color="black"))
OP50_day8

#================================================================
#day9
OP50_D9<-OP50_f23[,c('strain','REP','D9','color','line')]
#only mutants
OP50_D9_mutants<-OP50_D9 %>% filter(strain %in% c('WT','ΔubiG','ΔhemF','ΔadhE'))

#Dunnett test
#day9
DunnettTest(x=OP50_D9_mutants$D9, g=OP50_D9_mutants$strain, control = 'WT')
#$WT
# diff    lwr.ci   upr.ci   pval    
# ΔadhE-WT  -9.333333 -27.09978 8.433115 0.3659    
# ΔhemF-WT  -6.833333 -23.45232 9.785658 0.5410    
# ΔubiG-WT -13.333333 -31.09978 4.433115 0.1482    
#ubiG t test
OP50_D9_ubiG<-OP50_D9 %>% filter(strain %in% c('ΔubiG','ΔubiG+CoQ'))
OP50_D9_ubiG %>% t_test(D9 ~ line)  # p= 0.599

#===============================
#plot
OP50_D9$strain<-factor(OP50_D9$strain,levels=c('WT','WT+CoQ','WT+heme','WT+CoA',
                                               'ΔubiG','ΔubiG+CoQ',
                                               'ΔhemF','ΔhemF+heme',
                                               'ΔadhE','ΔadhE+CoA'))
OP50_D9$color<-factor(OP50_D9$color,levels=c('WT','ΔubiG','ΔhemF','ΔadhE'))
OP50_D9$line<-factor(OP50_D9$line,levels=c('solid','dashed'))

OP50_day9<-ggboxplot(
  OP50_D9, x = "strain", y = "D9",
  color = "color", palette = c('black', '#DA2E20FF','#3870C2FF','#66C84DFF'), linetype = 'line',
  #ylim = c(0,0.8),
  ylab='Proportion of L4 in day 9 (%)',xlab = '',size=0.5,ggtheme=theme_classic()
)+ scale_linetype_manual(values=c( 'solid',"dotted"))+
  #stat_pvalue_manual(p.d5, label = "p.adj.signif")+
  rotate_x_text(45,size=10)+theme(axis.text =  element_text(color="black"))
OP50_day9

#================================================================
#day10
OP50_D10<-OP50_f23[,c('strain','REP','D10','color','line')]
#only mutants
OP50_D10_mutants<-OP50_D10 %>% filter(strain %in% c('WT','ΔubiG','ΔhemF','ΔadhE'))

#Dunnett test
#day10
DunnettTest(x=OP50_D10_mutants$D10, g=OP50_D10_mutants$strain, control = 'WT')
#$WT
# diff    lwr.ci   upr.ci   pval    
# ΔadhE-WT  -9.333333 -27.09978 8.433115 0.3657    
# ΔhemF-WT  -6.833333 -23.45232 9.785658 0.5410    
# ΔubiG-WT -13.333333 -31.09978 4.433115 0.1484   
#ubiG t test
OP50_D10_ubiG<-OP50_D10 %>% filter(strain %in% c('ΔubiG','ΔubiG+CoQ'))
OP50_D10_ubiG %>% t_test(D10 ~ line)  # p= 0.599

#===============================
#plot
OP50_D10$strain<-factor(OP50_D10$strain,levels=c('WT','WT+CoQ','WT+heme','WT+CoA',
                                                 'ΔubiG','ΔubiG+CoQ',
                                                 'ΔhemF','ΔhemF+heme',
                                                 'ΔadhE','ΔadhE+CoA'))
OP50_D10$color<-factor(OP50_D10$color,levels=c('WT','ΔubiG','ΔhemF','ΔadhE'))
OP50_D10$line<-factor(OP50_D10$line,levels=c('solid','dashed'))

OP50_day10<-ggboxplot(
  OP50_D10, x = "strain", y = "D10",
  color = "color", palette = c('black', '#DA2E20FF','#3870C2FF','#66C84DFF'), linetype = 'line',
  #ylim = c(0,0.8),
  ylab='Proportion of L4 in day 10 (%)',xlab = '',size=0.5,ggtheme=theme_classic()
)+ scale_linetype_manual(values=c( 'solid',"dotted"))+
  #stat_pvalue_manual(p.d5, label = "p.adj.signif")+
  rotate_x_text(45,size=10)+theme(axis.text =  element_text(color="black"))
OP50_day10

#---------------------------------------------------------------------------------------------
#=================================ISOLATE-2================================
#Total
OP50_f23.2<-read_excel("input/mutants.xlsx",sheet = 'OP50_2_plot')
#day3
OP50_D3.2<-OP50_f23.2[,c('strain','REP','D3','color','line')]
#only mutants
OP50_D3.2_mutants<-OP50_D3.2 %>% filter(strain %in% c('WT','ΔubiG','ΔhemF','ΔadhE'))

#Dunnett test
library(DescTools)
#day3
DunnettTest(x=OP50_D3.2_mutants$D3, g=OP50_D3.2_mutants$strain, control = 'WT')
#$WT
# diff    lwr.ci    upr.ci    pval    
# ΔadhE-WT -29.33333 -35.82571 -22.84096 2.6e-09 ***
# ΔhemF-WT -35.83333 -41.90640 -29.76027 9.4e-11 ***
# ΔubiG-WT -38.66667 -45.15904 -32.17429 1.4e-12 ***

#

#===============================
OP50_D3.2$strain<-factor(OP50_D3.2$strain,levels=c('WT','WT+CoQ','WT+heme','WT+CoA',
                                                   'ΔubiG','ΔubiG+CoQ',
                                                   'ΔhemF','ΔhemF+heme',
                                                   'ΔadhE','ΔadhE+CoA'))
OP50_D3.2$color<-factor(OP50_D3.2$color,levels=c('WT','ΔubiG','ΔhemF','ΔadhE'))
OP50_D3.2$line<-factor(OP50_D3.2$line,levels=c('solid','dashed'))
#day3
OP50_day3.2<-ggboxplot(
  OP50_D3.2, x = "strain", y = "D3",
  color = "color", palette = c('black', '#DA2E20FF','#3870C2FF','#66C84DFF'), linetype = 'line',
  #ylim = c(0,0.8),
  ylab='Proportion of L4 in day 3 (%)',xlab = '',size=0.5,ggtheme=theme_classic()
)+ scale_linetype_manual(values=c( 'solid',"dotted"))+
  #stat_pvalue_manual(p.d5, label = "p.adj.signif")+
  rotate_x_text(45,size=10)+theme(axis.text =  element_text(color="black"))
OP50_day3.2

#=========================================================================================
#day4
OP50_D4.2<-OP50_f23.2[,c('strain','REP','D4','color','line')]
#only mutants
OP50_D4.2_mutants<-OP50_D4.2 %>% filter(strain %in% c('WT','ΔubiG','ΔhemF','ΔadhE'))

#Dunnett test
library(DescTools)
#day4
DunnettTest(x=OP50_D4.2_mutants$D4, g=OP50_D4.2_mutants$strain, control = 'WT')
#$WT
# diff    lwr.ci     upr.ci    pval    
# ΔadhE-WT -14.66667 -32.30143   2.968096 0.10443    
# ΔhemF-WT -32.50000 -48.99581 -16.004190 0.00095 ***
# ΔubiG-WT -52.00000 -69.63476 -34.365238 1.6e-05 ***

#

#===============================
OP50_D4.2$strain<-factor(OP50_D4.2$strain,levels=c('WT','WT+CoQ','WT+heme','WT+CoA',
                                                   'ΔubiG','ΔubiG+CoQ',
                                                   'ΔhemF','ΔhemF+heme',
                                                   'ΔadhE','ΔadhE+CoA'))
OP50_D4.2$color<-factor(OP50_D4.2$color,levels=c('WT','ΔubiG','ΔhemF','ΔadhE'))
OP50_D4.2$line<-factor(OP50_D4.2$line,levels=c('solid','dashed'))
#day4
OP50_day4.2<-ggboxplot(
  OP50_D4.2, x = "strain", y = "D4",
  color = "color", palette = c('black', '#DA2E20FF','#3870C2FF','#66C84DFF'), linetype = 'line',
  #ylim = c(0,0.8),
  ylab='Proportion of L4 in day 4 (%)',xlab = '',size=0.5,ggtheme=theme_classic()
)+ scale_linetype_manual(values=c( 'solid',"dotted"))+
  #stat_pvalue_manual(p.d5, label = "p.adj.signif")+
  rotate_x_text(45,size=10)+theme(axis.text =  element_text(color="black"))
OP50_day4.2


#=========================================================================================
#day5
OP50_D5.2<-OP50_f23.2[,c('strain','REP','D5','color','line')]
#only mutants
OP50_D5.2_mutants<-OP50_D5.2 %>% filter(strain %in% c('WT','ΔubiG','ΔhemF','ΔadhE'))

#Dunnett test
library(DescTools)
#day5
DunnettTest(x=OP50_D5.2_mutants$D5, g=OP50_D5.2_mutants$strain, control = 'WT')
#$WT
# diff    lwr.ci     upr.ci   pval    
# ΔadhE-WT -12.00000 -26.91340   2.913401 0.1176    
# ΔhemF-WT -19.50000 -33.45021  -5.549791 0.0088 ** 
# ΔubiG-WT -36.66667 -51.58007 -21.753266 0.0002 ***

#

#===============================
OP50_D5.2$strain<-factor(OP50_D5.2$strain,levels=c('WT','WT+CoQ','WT+heme','WT+CoA',
                                                   'ΔubiG','ΔubiG+CoQ',
                                                   'ΔhemF','ΔhemF+heme',
                                                   'ΔadhE','ΔadhE+CoA'))
OP50_D5.2$color<-factor(OP50_D5.2$color,levels=c('WT','ΔubiG','ΔhemF','ΔadhE'))
OP50_D5.2$line<-factor(OP50_D5.2$line,levels=c('solid','dashed'))
#day5
OP50_day5.2<-ggboxplot(
  OP50_D5.2, x = "strain", y = "D5",
  color = "color", palette = c('black', '#DA2E20FF','#3870C2FF','#66C84DFF'), linetype = 'line',
  #ylim = c(0,0.8),
  ylab='Proportion of L4 in day 5 (%)',xlab = '',size=0.5,ggtheme=theme_classic()
)+ scale_linetype_manual(values=c( 'solid',"dotted"))+
  #stat_pvalue_manual(p.d5, label = "p.adj.signif")+
  rotate_x_text(45,size=10)+theme(axis.text =  element_text(color="black"))
OP50_day5.2


#=========================================================================================
#day6
OP50_D6.2<-OP50_f23.2[,c('strain','REP','D6','color','line')]
#only mutants
OP50_D6.2_mutants<-OP50_D6.2 %>% filter(strain %in% c('WT','ΔubiG','ΔhemF','ΔadhE'))

#Dunnett test
library(DescTools)
#day6
DunnettTest(x=OP50_D6.2_mutants$D6, g=OP50_D6.2_mutants$strain, control = 'WT')
#$WT
# diff    lwr.ci      upr.ci    pval    
# ΔadhE-WT  -9.333333 -26.26568   7.5990121 0.33212    
# ΔhemF-WT -15.666667 -31.50543   0.1720922 0.05230 .  
# ΔubiG-WT -33.333333 -50.26568 -16.4009879 0.00095 ***

#

#===============================
OP50_D6.2$strain<-factor(OP50_D6.2$strain,levels=c('WT','WT+CoQ','WT+heme','WT+CoA',
                                                   'ΔubiG','ΔubiG+CoQ',
                                                   'ΔhemF','ΔhemF+heme',
                                                   'ΔadhE','ΔadhE+CoA'))
OP50_D6.2$color<-factor(OP50_D6.2$color,levels=c('WT','ΔubiG','ΔhemF','ΔadhE'))
OP50_D6.2$line<-factor(OP50_D6.2$line,levels=c('solid','dashed'))
#day6
OP50_day6.2<-ggboxplot(
  OP50_D6.2, x = "strain", y = "D6",
  color = "color", palette = c('black', '#DA2E20FF','#3870C2FF','#66C84DFF'), linetype = 'line',
  #ylim = c(0,0.8),
  ylab='Proportion of L4 in day 6 (%)',xlab = '',size=0.5,ggtheme=theme_classic()
)+ scale_linetype_manual(values=c( 'solid',"dotted"))+
  #stat_pvalue_manual(p.D6, label = "p.adj.signif")+
  rotate_x_text(45,size=10)+theme(axis.text =  element_text(color="black"))
OP50_day6.2


#=========================================================================================
#day7
OP50_D7.2<-OP50_f23.2[,c('strain','REP','D7','color','line')]
#only mutants
OP50_D7.2_mutants<-OP50_D7.2 %>% filter(strain %in% c('WT','ΔubiG','ΔhemF','ΔadhE'))

#Dunnett test
library(DescTools)
#day7
DunnettTest(x=OP50_D7.2_mutants$D7, g=OP50_D7.2_mutants$strain, control = 'WT')
# $WT
# diff    lwr.ci    upr.ci   pval    
# ΔadhE-WT  -8.666667 -29.19264 11.859305 0.5224    
# ΔhemF-WT -13.833333 -33.03362  5.366955 0.1682    
# ΔubiG-WT -28.666667 -49.19264 -8.140695 0.0091 ** 

#

#===============================
OP50_D7.2$strain<-factor(OP50_D7.2$strain,levels=c('WT','WT+CoQ','WT+heme','WT+CoA',
                                                   'ΔubiG','ΔubiG+CoQ',
                                                   'ΔhemF','ΔhemF+heme',
                                                   'ΔadhE','ΔadhE+CoA'))
OP50_D7.2$color<-factor(OP50_D7.2$color,levels=c('WT','ΔubiG','ΔhemF','ΔadhE'))
OP50_D7.2$line<-factor(OP50_D7.2$line,levels=c('solid','dashed'))
#day7
OP50_day7.2<-ggboxplot(
  OP50_D7.2, x = "strain", y = "D7",
  color = "color", palette = c('black', '#DA2E20FF','#3870C2FF','#66C84DFF'), linetype = 'line',
  #ylim = c(0,0.8),
  ylab='Proportion of L4 in day 7 (%)',xlab = '',size=0.5,ggtheme=theme_classic()
)+ scale_linetype_manual(values=c( 'solid',"dotted"))+
  #stat_pvalue_manual(p.D7, label = "p.adj.signif")+
  rotate_x_text(45,size=10)+theme(axis.text =  element_text(color="black"))
OP50_day7.2


#=========================================================================================
#day8
OP50_D8.2<-OP50_f23.2[,c('strain','REP','D8','color','line')]
#only mutants
OP50_D8.2_mutants<-OP50_D8.2 %>% filter(strain %in% c('WT','ΔubiG','ΔhemF','ΔadhE'))

#Dunnett test
library(DescTools)
#day8
DunnettTest(x=OP50_D8.2_mutants$D8, g=OP50_D8.2_mutants$strain, control = 'WT')
# $WT
# diff    lwr.ci    upr.ci   pval    
# ΔadhE-WT  -8.666667 -27.65907 10.325735 0.4668    
# ΔhemF-WT -13.833333 -31.59910  3.932432 0.1315    
# ΔubiG-WT -26.666667 -45.65907 -7.674265 0.0085 ** 

#

#===============================
OP50_D8.2$strain<-factor(OP50_D8.2$strain,levels=c('WT','WT+CoQ','WT+heme','WT+CoA',
                                                   'ΔubiG','ΔubiG+CoQ',
                                                   'ΔhemF','ΔhemF+heme',
                                                   'ΔadhE','ΔadhE+CoA'))
OP50_D8.2$color<-factor(OP50_D8.2$color,levels=c('WT','ΔubiG','ΔhemF','ΔadhE'))
OP50_D8.2$line<-factor(OP50_D8.2$line,levels=c('solid','dashed'))
#day8
OP50_day8.2<-ggboxplot(
  OP50_D8.2, x = "strain", y = "D8",
  color = "color", palette = c('black', '#DA2E20FF','#3870C2FF','#66C84DFF'), linetype = 'line',
  #ylim = c(0,0.8),
  ylab='Proportion of L4 in day 8 (%)',xlab = '',size=0.5,ggtheme=theme_classic()
)+ scale_linetype_manual(values=c( 'solid',"dotted"))+
  #stat_pvalue_manual(p.D8, label = "p.adj.signif")+
  rotate_x_text(45,size=10)+theme(axis.text =  element_text(color="black"))
OP50_day8.2



#=========================================================================================
#day9
OP50_D9.2<-OP50_f23.2[,c('strain','REP','D9','color','line')]
#only mutants
OP50_D9.2_mutants<-OP50_D9.2 %>% filter(strain %in% c('WT','ΔubiG','ΔhemF','ΔadhE'))

#Dunnett test
library(DescTools)
#day9
DunnettTest(x=OP50_D9.2_mutants$D9, g=OP50_D9.2_mutants$strain, control = 'WT')

# $WT
# diff    lwr.ci    upr.ci   pval    
# ΔadhE-WT  -8.666667 -28.05885 10.725513 0.4819    
# ΔhemF-WT -13.333333 -31.47306  4.806389 0.1583    
# ΔubiG-WT -25.333333 -44.72551 -5.941154 0.0132 * 

#

#===============================
OP50_D9.2$strain<-factor(OP50_D9.2$strain,levels=c('WT','WT+CoQ','WT+heme','WT+CoA',
                                                   'ΔubiG','ΔubiG+CoQ',
                                                   'ΔhemF','ΔhemF+heme',
                                                   'ΔadhE','ΔadhE+CoA'))
OP50_D9.2$color<-factor(OP50_D9.2$color,levels=c('WT','ΔubiG','ΔhemF','ΔadhE'))
OP50_D9.2$line<-factor(OP50_D9.2$line,levels=c('solid','dashed'))
#day9
OP50_day9.2<-ggboxplot(
  OP50_D9.2, x = "strain", y = "D9",
  color = "color", palette = c('black', '#DA2E20FF','#3870C2FF','#66C84DFF'), linetype = 'line',
  #ylim = c(0,0.8),
  ylab='Proportion of L4 in day 9 (%)',xlab = '',size=0.5,ggtheme=theme_classic()
)+ scale_linetype_manual(values=c( 'solid',"dotted"))+
  #stat_pvalue_manual(p.D9, label = "p.adj.signif")+
  rotate_x_text(45,size=10)+theme(axis.text =  element_text(color="black"))
OP50_day9.2


#=========================================================================================
#day10
OP50_D10.2<-OP50_f23.2[,c('strain','REP','D10','color','line')]
#only mutants
OP50_D10.2_mutants<-OP50_D10.2 %>% filter(strain %in% c('WT','ΔubiG','ΔhemF','ΔadhE'))

#Dunnett test
library(DescTools)
#day10
DunnettTest(x=OP50_D10.2_mutants$D10, g=OP50_D10.2_mutants$strain, control = 'WT')

# $WT
# diff    lwr.ci    upr.ci   pval    
# ΔadhE-WT  -8.666667 -28.05885 10.725513 0.4818    
# ΔhemF-WT -13.333333 -31.47306  4.806389 0.1584    
# ΔubiG-WT -25.333333 -44.72551 -5.941154 0.0131 * 

#

#===============================
OP50_D10.2$strain<-factor(OP50_D10.2$strain,levels=c('WT','WT+CoQ','WT+heme','WT+CoA',
                                                     'ΔubiG','ΔubiG+CoQ',
                                                     'ΔhemF','ΔhemF+heme',
                                                     'ΔadhE','ΔadhE+CoA'))
OP50_D10.2$color<-factor(OP50_D10.2$color,levels=c('WT','ΔubiG','ΔhemF','ΔadhE'))
OP50_D10.2$line<-factor(OP50_D10.2$line,levels=c('solid','dashed'))
#day10
OP50_day10.2<-ggboxplot(
  OP50_D10.2, x = "strain", y = "D10",
  color = "color", palette = c('black', '#DA2E20FF','#3870C2FF','#66C84DFF'), linetype = 'line',
  #ylim = c(0,0.8),
  ylab='Proportion of L4 in day 10 (%)',xlab = '',size=0.5,ggtheme=theme_classic()
)+ scale_linetype_manual(values=c( 'solid',"dotted"))+
  #stat_pvalue_manual(p.D10, label = "p.adj.signif")+
  rotate_x_text(45,size=10)+theme(axis.text =  element_text(color="black"))
OP50_day10.2


#merge all
library(ggpubr)
merge1<-ggarrange(OP50_day3,OP50_day4,OP50_day5,OP50_day6,
                  ncol = 4,legend='none')
merge1
merge2<-ggarrange(OP50_day7,OP50_day8,OP50_day9,OP50_day10,
                  ncol = 4,legend='none')
merge2
merge12<-ggarrange(merge1,merge2,ncol = 1,legend='none')
merge12

merge3<-ggarrange(OP50_day3.2,OP50_day4.2,OP50_day5.2,OP50_day6.2,
                  ncol = 4,legend='none')
merge3
merge4<-ggarrange(OP50_day7.2,OP50_day8.2,OP50_day9.2,OP50_day10.2,
                  ncol = 4,legend='none')
merge4
merge34<-ggarrange(merge3,merge4,ncol = 1,legend='none')
merge34
ggsave('box_merge1.pdf',merge12,height = 8,width = 15)
ggsave('box_merge2.pdf',merge34,height = 8,width = 15)

#================================================================================================
#Fig. S19
#N2
#isolate-1
data.n2_1<-read_excel("input/mutants.xlsx",sheet = 'N2_OP50_1')

###########################################################
#Dunnett test
library(DescTools)

mutants_n2.1<-data.n2_1 %>% filter(strain %in% c('WT','ΔubiG','ΔhemF','ΔadhE'))
DunnettTest(x=mutants_n2.1$`Egg-laying time`, g=mutants_n2.1$strain, control = 'WT')


#diff     lwr.ci    upr.ci   pval    
# diff     lwr.ci    upr.ci    pval    
# ΔubiG-WT  1.830000000  1.4903728 2.1696272 2.8e-08 ***
#   ΔhemF-WT  1.996666667  1.6570394 2.3362939 1.0e-11 ***
#   ΔadhE-WT -0.003333333 -0.3429606 0.3362939  1.0000       

#===============================
library(ggplot2)
library(ggpubr)
library(paletteer) 
data.n2_1$strain<-factor(data.n2_1$strain,levels=c('WT','WT+CoQ','WT+heme','WT+CoA',
                                                   'ΔubiG','ΔubiG+CoQ',
                                                   'ΔhemF','ΔhemF+heme',
                                                   'ΔadhE','ΔadhE+CoA'))
data.n2_1$color<-factor(data.n2_1$color,levels = c('WT','ΔubiG','ΔhemF','ΔadhE'))
data.n2_1$line<-factor(data.n2_1$line,levels = c('solid','dashed'))
n2.op50.1<-ggboxplot(
  data.n2_1, x = "strain", y = "`Egg-laying time`",
  color = "color", palette = c('black', '#DA2E20FF','#3870C2FF','#66C84DFF'),
  linetype = 'line',
  #ylim = c(0,0.8),
  ylab='Egg laying time (h)',xlab = '',size=0.5,ggtheme=theme_classic()
)+theme(axis.text = element_text(color="black"))+
  #stat_pvalue_manual(p.d5, label = "p.adj.signif")+
  rotate_x_text(45,size=10)#+geom_text(aes(label = p.growth),size=3)
n2.op50.1
ggsave('OP50/n2.op50.1.pdf',n2.op50.1,width = 6,height = 4)

#---------------------------------------
#isolate-2
data.n2_2<-read_excel("input/mutants.xlsx",sheet = 'N2_OP50_2')

###########################################################
#Dunnett test
library(DescTools)

mutants_n2.2<-data.n2_2 %>% filter(strain %in% c('WT','ΔubiG','ΔhemF','ΔadhE'))
DunnettTest(x=mutants_n2.2$`Egg-laying time`, g=mutants_n2.2$strain, control = 'WT')


#diff     lwr.ci    upr.ci   pval    
# $WT
# diff     lwr.ci    upr.ci    pval    
# ΔadhE-WT -0.003333333 -0.4811839 0.4745172   1e+00    
# ΔhemF-WT  1.830000000  1.3521494 2.3078506   1e-04 ***
#   ΔubiG-WT  1.830000000  1.3521494 2.3078506 3.9e-07 ***


#===============================
library(ggplot2)
library(ggpubr)
library(paletteer) 
data.n2_2$strain<-factor(data.n2_2$strain,levels=c('WT','WT+CoQ','WT+heme','WT+CoA',
                                                   'ΔubiG','ΔubiG+CoQ',
                                                   'ΔhemF','ΔhemF+heme',
                                                   'ΔadhE','ΔadhE+CoA'))
data.n2_2$color<-factor(data.n2_2$color,levels = c('WT','ΔubiG','ΔhemF','ΔadhE'))
data.n2_2$line<-factor(data.n2_2$line,levels = c('solid','dashed'))
n2.op50.2<-ggboxplot(
  data.n2_2, x = "strain", y = "`Egg-laying time`",
  color = "color", palette = c('black', '#DA2E20FF','#3870C2FF','#66C84DFF'),
  linetype = 'line',
  #ylim = c(0,0.8),
  ylab='Egg laying time (h)',xlab = '',size=0.5,ggtheme=theme_classic()
)+theme(axis.text = element_text(color="black"))+
  #stat_pvalue_manual(p.d5, label = "p.adj.signif")+
  rotate_x_text(45,size=10)#+geom_text(aes(label = p.growth),size=3)
n2.op50.2
ggsave('OP50/n2.op50.2.pdf',n2.op50.2,width = 6,height = 4)


n2_merge<-ggarrange(n2.op50.1,n2.op50.2,legend = 'none')
n2_merge
ggsave('OP50/n2.op50.merge.pdf',n2_merge,width = 6,height = 3)

################################################################################################
##########---------------------------------------------------------------------
#Fig. S12
#curve plot
library(ggalt)
library(tidyverse)
library(paletteer)
library(readxl)
#curve_1
curve_1 <-read_excel('input/mutants.xlsx', sheet='OP50_1_curve')


curve_1 %>% 
  rowwise() %>% 
  mutate(mean_value=mean(c(rep1,rep2,rep3)),
         std_error=plotrix::std.error(c(rep1,rep2,rep3))) %>% 
  ungroup() -> new.curve_1

new.curve_1$strain<-factor(new.curve_1$strain,levels=c('WT','ΔubiG','ΔhemF','ΔadhE',
                                                       'ΔubiG+CoQ','ΔhemF+heme','ΔadhE+CoA'))
new.curve_1$color<-factor(new.curve_1$color,levels=c('WT','ΔubiG','ΔhemF','ΔadhE'))
new.curve_1$line<-factor(new.curve_1$line,levels=c('solid','dashed'))

#plot--cur1
op50_1<- ggplot(data=new.curve_1,aes(x=Day,y=mean_value,group=strain, color=color,linetype=line
))+
  geom_errorbar(aes(ymin=mean_value-std_error,
                    ymax=mean_value+std_error),linetype="solid",
                width=0.1,linewidth=0.2)+
  geom_point(size=0.8)+
  geom_line(linewidth =0.3)+
  #geom_xspline(spline_shape = -0.6,size =6)+
  scale_linetype_manual(values=c( 'solid',"dotted"))+
  geom_vline(xintercept=5, linewidth= 0.5,color="lightgrey",linetype = "dashed")+
  scale_color_manual(values = c('black', '#DA2E20FF','#3870C2FF','#66C84DFF'),
                     name=NULL)+
  scale_x_continuous(#limits = c(3,10),
    expand = expansion(mult=c(0,0)),
    breaks = seq(10)
  )+
  scale_y_continuous(limits = c(-0.5,80),
                     expand = expansion(mult=c(0,0)),
                     #breaks = seq(),
  )+
  coord_cartesian(clip = "off")+
  theme_classic()+
  theme(legend.text = element_text(hjust=0),
        # legend.position = c(0.1,0.4)
        axis.text = element_text(colour = "black")
  )+
  labs(x="Day",
       y=expression(paste("L4 proportion of L. marina (%)")))
op50_1

ggsave('OP50/op50_1.pdf', op50_1, width = 5,height = 2.2)

#curve_2
curve_2 <-read_excel('input/mutants.xlsx', sheet='OP50_2_curve')


curve_2 %>% 
  rowwise() %>% 
  mutate(mean_value=mean(c(rep1,rep2,rep3)),
         std_error=plotrix::std.error(c(rep1,rep2,rep3))) %>% 
  ungroup() -> new.curve_2

new.curve_2$strain<-factor(new.curve_2$strain,levels=c('WT','ΔubiG','ΔhemF','ΔadhE',
                                                       'ΔubiG+CoQ','ΔhemF+heme','ΔadhE+CoA'))
new.curve_2$color<-factor(new.curve_2$color,levels=c('WT','ΔubiG','ΔhemF','ΔadhE'))
new.curve_2$line<-factor(new.curve_2$line,levels=c('solid','dashed'))

#plot--cur2
op50_2<- ggplot(data=new.curve_2,aes(x=Day,y=mean_value,group=strain, color=color,linetype=line
))+
  geom_errorbar(aes(ymin=mean_value-std_error,
                    ymax=mean_value+std_error),linetype="solid",
                width=0.1,linewidth=0.2)+
  geom_point(size=0.8)+
  geom_line(linewidth = 0.3)+
  #geom_xspline(spline_shape = -0.5,size =6)+
  scale_linetype_manual(values=c( 'solid',"dotted"))+
  geom_vline(xintercept=5, linewidth= 0.5,color="lightgrey",linetype = "dashed")+
  scale_color_manual(values = c('black', '#DA2E20FF','#3870C2FF','#66C84DFF'),
                     name=NULL)+
  scale_x_continuous(#limits = c(3,10),
    expand = expansion(mult=c(0,0)),
    breaks = seq(10)
  )+
  scale_y_continuous(limits = c(-0.5,80),
                     expand = expansion(mult=c(0,0)),
                     #breaks = seq(),
  )+
  coord_cartesian(clip = "off")+
  theme_classic()+
  theme(legend.text = element_text(hjust=0),
        # legend.position = c(0.1,0.4)
        axis.text = element_text(colour = "black")
  )+
  labs(x="Day",
       y=expression(paste("L4 proportion of L. marina (%)")))
op50_2

ggsave('OP50/op50_1.pdf', op50_1, width = 5,height = 2.2)

#merge
library(ggpubr)
merge<-ggarrange(op50_1,op50_2,ncol = 1)
merge
ggsave('OP50/op50_merge.pdf', merge, width = 5,height = 5)



