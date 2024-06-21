#Fig. 3a was created in iTols using ./input/Fig. 3a/*
#Supplementary Fig. 9 was created in iTols using ./input/Fig. S9/*

#Fig. 3b data from gapseq and manually organized
library(pheatmap)
library(ggpubr)
#plot
library(readxl)
dt <- read_excel("input/metabolic.xlsx", sheet = 1)
#remove first two cols
dt1<-dt[,-c(1,2)]
rownames(dt1)=dt$path_sum
#Replace >0 with 1
library(dplyr)
dt2<-dt1 %>% mutate_if(is.numeric, ~1 * (. > 0))
rownames(dt2)=dt$path_sum
anno_row<-dt[,c(2)]
rownames(anno_row)=dt$path_sum
anno_row<-as.data.frame(anno_row)

anno_col<-read_excel("input/metabolic.xlsx", sheet = 2)
anno_col<-mutate(anno_col,log2.N2=log2(N2))
anno_col1 <- anno_col[,c(20,25,23,5,10)]
rownames(anno_col1)=anno_col$ID
anno_col1<-as.data.frame(anno_col1)
colnames(anno_col1)<-c('Family','Phylum / Class','Gram strain','Log2(N2)','F23')
#add color
anno_colors = list('group'= c('Cell wall'='#8FA87AFF',  'Cofactor'='#244579FF','DNA/RNA'='#9F248FFF',
                              'Lipids'='#C6242DFF', 'Protein'='#F9791EFF'),
                   'Phylum / Class'=c('Actinobacteria'='#EEE685','Bacteroidetes'='#4F94CD',
                                      'Firmicutes'='#EE7600','Alphaproteobacteria'='#00CD00',
                                      'Gammaproteobacteria'='#90EE90'),
                   'Gram strain' = c('Negative'='#7E8BD6FF','Positive'='#E84A5FFF'),
                   'Log2(N2)'=c("#004D40FF","#E0F2F1FF"),'F23'=c('#E3F2FDFF','#0D47A1FF'),
                   'Family'=c('Dermabacteraceae'='#9E9D24FF','Microbacteriaceae'='#AFB42BFF','Micrococcaceae'='#C0CA33FF',
                              'Nocardiaceae'='#DCE775FF', 'Rhodobacteraceae'='#00CD00','Sphingomonadaceae'='#00EE00',
                              'Flavobacteriaceae'='#4F94CD', 'Bacillaceae'='#D95204FF','Planococcaceae'='#F28705FF',
                              'Staphylococcaceae'='#F29F05FF','Streptococcaceae'='#F2B705FF',
                              'Alishewanella'='#C5E1A5FF','Alteromonadaceae'='#AED581FF','Cellvibrionaceae'='#9CCC65FF',
                              'Colwelliaceae'='#8BC34AFF', 'Enterobacteriaceae'='#7CB342FF','Marinobacteraceae'='#689F38FF',
                              'Moraxellaceae'='#558B2FFF','Oceanospirillaceae'='#33691EFF',
                              'Pseudoalteromonadaceae'='#8FBD77FF','Pseudomonadaceae'='#84AD83FF',
                              'Saccharospirillaceae'='#779D8DFF','Shewanellaceae'='#577F9FFF',
                              'Vibrionaceae'='#3E71A8FF','Woeseiaceae'='#698E96FF'
                   )
)

#plot
pheatmap(dt2, cluster_cols = T, cluster_rows = F, 
         clustering_distance_cols = "euclidean",clustering_method='ward.D',
         color= c("#F9FBE7FF","#604A76FF"),legend = F,
         fontsize = 5,fontsize_row = 3, fontsize_col = 3.2,angle_col = 45,
         cutree_cols = 3,
         scale="none",
         border_color= NA,
         border=F,
         annotation_row = anno_row,
         annotation_col = anno_col1,
         annotation_colors = anno_colors,
         gaps_row=c(4,30,35,60),
         annotation_names_row=F,annotation_legend=T,treeheight_col=10,
         cellwidth = 4, cellheight = 3,
         # filename='bsy.pdf',width=9,height=8
)
#darkgreen

#Fig. S30 transports
library(tidyverse)
library(readxl)

transport<- list.files('transport150',pattern = 'Transporter.tbl', full.names = T) %>%
  lapply(read.delim,colClasses=c("NULL", "NULL", "character",
                                 "NULL", "NULL","NULL","NULL",
                                 "NULL", "NULL","NULL","NULL", "NULL","NULL",
                                 "NULL", "NULL"))  %>% 
  rapply( function(x) if(is.character(x)) tolower(x) else x, how="replace")%>% 
  lapply(unique)%>%  
  lapply(cbind,trans='T') %>% 
  reduce(full_join,by=c('sub')) %>% 
  arrange(sub) 

library(stringr)
colnames(transport)<- c('sub',
                        str_replace(list.files('transport150',pattern = '*Transporter.tbl', full.names = F),'-Transporter.tbl', ''))  
library("xlsx")
write.xlsx(transport,'transport.xlsx', sheetName = "transport", append = TRUE)

#Replace T to 1, NA to 0
##---------------------plot-----------------------------------------------------------------
library(pheatmap)
library(ggpubr)
library(readxl)
dt <- read_excel("input/transport.xlsx", sheet = "plot")
#remove first two cols
dt1<-dt[,-c(1,2)]
rownames(dt1)=dt$sub
anno_row<-dt[,c(2)]
rownames(anno_row)=dt$sub
anno_row<-as.data.frame(anno_row)
anno_colors = list('group'= c('Amino acid'='#C0FF3E','Monosaccharides'='#BCD2EE','Fatty acid'='#BA55D3',
                              'Poligo-polysaccharides'='#BDB76B','Vitamin'='#AB82FF','Other'='#7D7D7D'),
                   'Phylum / Class'=c('Actinobacteria'='#EEE685','Bacteroidetes'='#4F94CD',
                                      'Firmicutes'='#EE7600','Alphaproteobacteria'='#00CD00',
                                      'Gammaproteobacteria'='#90EE90'),
                   'Gram strain' = c('Negative'='#7E8BD6FF','Positive'='#E84A5FFF'),
                   'C. elegans growth'=c("#004D40FF","#E0F2F1FF"),'L. marina growth'=c('#E3F2FDFF','#0D47A1FF'),
                   'Family'=c('Dermabacteraceae'='#9E9D24FF','Microbacteriaceae'='#AFB42BFF','Micrococcaceae'='#C0CA33FF','Nocardiaceae'='#DCE775FF',
                              'Rhodobacteraceae'='#00CD00','Sphingomonadaceae'='#00EE00',
                              'Flavobacteriaceae'='#4F94CD',
                              'Bacillaceae'='#D95204FF','Planococcaceae'='#F28705FF','Staphylococcaceae'='#F29F05FF','Streptococcaceae'='#F2B705FF',
                              'Alishewanella'='#C5E1A5FF','Alteromonadaceae'='#AED581FF','Cellvibrionaceae'='#9CCC65FF','Colwelliaceae'='#8BC34AFF',
                              'Enterobacteriaceae'='#7CB342FF','Marinobacteraceae'='#689F38FF','Moraxellaceae'='#558B2FFF','Oceanospirillaceae'='#33691EFF',
                              'Pseudoalteromonadaceae'='#8FBD77FF','Pseudomonadaceae'='#84AD83FF','Saccharospirillaceae'='#779D8DFF','Shewanellaceae'='#577F9FFF',
                              'Vibrionaceae'='#3E71A8FF','Woeseiaceae'='#698E96FF'))
anno_col<-read_excel("transport.xlsx", sheet = "anno")
anno_col<-mutate(anno_col,log2.N2=log2(N2))
anno_col1 <- anno_col[,c(20,25,23,5,10)]
rownames(anno_col1)=anno_col$ID
anno_col1<-as.data.frame(anno_col1)
colnames(anno_col1)<-c('Family','Phylum / Class','Gram strain','C. elegans growth','L. marina growth')
#anno_col<-mutate(anno_col,log2.N2=log2(N2))
#anno_col1 <- anno_col[,c(20,25,23)]
#rownames(anno_col1)=anno_col$ID
#anno_col1<-as.data.frame(anno_col1)
#colnames(anno_col1)<-c('Family','Phylum / Class','Gram strain')

pheatmap(dt1, cluster_cols = T, cluster_rows = F, 
         clustering_distance_cols = "euclidean",clustering_method='ward.D2',
         color= c("#F9FBE7FF","#604A76FF"),legend = F,
         fontsize = 5,fontsize_row = 3, fontsize_col = 3.2,angle_col = 45,
         cutree_cols = 3,scale="none",
         #border_color= 'NA',
         border=F,
         annotation_row = anno_row,annotation_colors = anno_colors,
         annotation_col = anno_col1,
         annotation_names_row=F,annotation_legend=T,
         gaps_row=c(34,41,61,72,78),
         treeheight_col=8,
         cellwidth = 4, cellheight = 3,
         filename="transport.pdf",
)

#Fig. S33  Prediction of carbon sources and metabolic byproducts 
library(tidyverse)
library(readxl)
###------------clear up the outputs from gapseq----------------------------------
#list.files('input/cs_fer_150/',pattern = 'ferm.tbl', full.names = T) %>%
#  lapply(read.delim) %>%
#  reduce(full_join,by=c('id','name')) %>%
#  arrange(id) -> ferm
#library(stringr)
#colnames(ferm)<- c('id','name',
#                   str_replace(list.files('cs_fer_150',pattern = '*ferm.tbl', full.names = F),'-ferm.tbl', ''))  
#
#list.files('cs_fer_150',pattern = 'cs.tbl', full.names = T) %>%
#  lapply(read.delim) %>%
#  reduce(full_join,by=c('id','name')) %>%
#  arrange(id) -> cs
#colnames(cs)<- c('id','name',
#                 str_replace(list.files('cs_fer_150',pattern = '*cs.tbl', full.names = F),'-cs.tbl', ''))  
#library("xlsx")
#write.xlsx(ferm,'sybil.commen.xlsx', sheetName = "ferm", append = TRUE)
#write.xlsx(cs,'sybil.commen.xlsx', sheetName = "cs", append = TRUE)
###-----------------------------------------------------------------------------------------
#plot
library(pheatmap)
library(ggpubr)

#plot-all
library(readxl)
dt <- read_excel("input/sybil.xlsx", sheet = "plot")
# convert NA to 0
dt[is.na(dt)] <- 0
#remove first three cols
dt1<-dt[,-c(1,2,3)]
rownames(dt1)=dt$name
anno_row<-dt[,c(3)]
rownames(anno_row)=dt$name
anno_row<-as.data.frame(anno_row)

anno_col<-read_excel("input/sybil.xlsx", sheet = "anno")
anno_col<-mutate(anno_col,log2.N2=log2(N2))
anno_col1 <- anno_col[,c(20,25,23,5,10)]
rownames(anno_col1)=anno_col$ID
anno_col1<-as.data.frame(anno_col1)
colnames(anno_col1)<-c('Family','Phylum / Class','Gram strain','Log2(N2)','F23')
anno_colors = list('Class'= c('Carbon sources'='#C0FF3E','Metabolic byproducts'='#AEEEEE'),
                   'Phylum / Class'=c('Actinobacteria'='#EEE685','Bacteroidetes'='#4F94CD',
                                      'Firmicutes'='#EE7600','Alphaproteobacteria'='#00CD00',
                                      'Gammaproteobacteria'='#90EE90'),
                   'Gram strain' = c('Negative'='#7E8BD6FF','Positive'='#E84A5FFF'),
                   'Log2(N2)'=c("#004D40FF","#E0F2F1FF"),'F23'=c('#E3F2FDFF','#0D47A1FF'),
                   'Family'=c('Dermabacteraceae'='#9E9D24FF','Microbacteriaceae'='#AFB42BFF','Micrococcaceae'='#C0CA33FF','Nocardiaceae'='#DCE775FF',
                              'Rhodobacteraceae'='#00CD00','Sphingomonadaceae'='#00EE00',
                              'Flavobacteriaceae'='#4F94CD',
                              'Bacillaceae'='#D95204FF','Planococcaceae'='#F28705FF','Staphylococcaceae'='#F29F05FF','Streptococcaceae'='#F2B705FF',
                              'Alishewanella'='#C5E1A5FF','Alteromonadaceae'='#AED581FF','Cellvibrionaceae'='#9CCC65FF','Colwelliaceae'='#8BC34AFF',
                              'Enterobacteriaceae'='#7CB342FF','Marinobacteraceae'='#689F38FF','Moraxellaceae'='#558B2FFF','Oceanospirillaceae'='#33691EFF',
                              'Pseudoalteromonadaceae'='#8FBD77FF','Pseudomonadaceae'='#84AD83FF','Saccharospirillaceae'='#779D8DFF','Shewanellaceae'='#577F9FFF',
                              'Vibrionaceae'='#3E71A8FF','Woeseiaceae'='#698E96FF')
)



pheatmap(dt1, cluster_cols = T, cluster_rows = F, 
         clustering_distance_cols = "euclidean",clustering_method='ward.D',
         color= c("#F9FBE7FF","#604A76FF"),legend = F,
         fontsize = 5,fontsize_row = 3, fontsize_col = 3.2,angle_col = 45,
         cutree_cols = 4,scale="none",
         #border_color= 'NA',
         border=F,
         annotation_row = anno_row,annotation_colors = anno_colors,
         annotation_col = anno_col1,
         annotation_names_row=F,annotation_legend=T,
         gaps_row=55,treeheight_col=8,
         cellwidth = 4, cellheight = 3,
         #filename="ferm_cs.pdf"
)

#Fig. S36 was created using BUSCO

#Fig. S37. Comparison of BUSCO and MEMOTE parameters between complete genomes
##busco score
df.busco<-read.xlsx('input/df.busco.xlsx',1)
##MEmote score
df.memote<-read.xlsx('input/Memote.xlsx',1)
##seq type
df.seq<-read.xlsx('input/BUSCO_MEMOTE.xlsx',1)
merge1<-merge(df.seq,df.busco,by.x='ID',by.y='my_species',sort=F)
merge2<-merge(merge1,df.memote,by='ID',sort=F)
df.test<-merge2[,c(1:7)]

library(rstatix)
library(ggplot2)
library(ggpubr)
#gather
df.test2<-gather(merge2,category, value,S,D,F,M,Total.Score,stoichiometric.consistency,mass.balanced,charge.balanced,metabolite.connectivity)
#df.test%>%wilcox_test(Total.Score ~ Seq)  #0.149
df.test2$category = factor(df.test2$category,levels = c('S','D','F','M','Total.Score',
                                                        'stoichiometric.consistency','mass.balanced',
                                                        'charge.balanced','metabolite.connectivity'))
busco_memote_plot<-ggboxplot(df.test2, x = "Seq", y = "value",
                             color = "Seq", #palette = "jco",
                             add = "jitter",
                             facet.by = "category", short.panel.labs = FALSE,
                             scales = "free_y",ncol=3)+
  stat_compare_means(method = 'wilcox',label.x.npc = 'right',label.y.npc = 'bottom')
busco_memote_plot

ggsave('busco_memote_plot.pdf',busco_memote_plot,width = 10,height = 8)


#Fig. S10
#BIOLOG test in silico
library(data.table)
library(sybil)
library(ggplot2)
library(reshape2)
library(stringr)
library(readxl)
library(xlsx)

gapseq<-readRDS('input/HQbiome.rds')
mysub<-read_excel('input/mysub.xlsx')
set_diet <- function(model, medium, uptake=-100,verbose=TRUE){
  ex <- findExchReact(model) # find exchange reactions
  ex <- ex[grep("^EX_",ex@react_id),] # reduce to real exchange reactions only
  if(is.null(medium)){
    medium <- ex@react_id[which(ex@lowbnd < 0)]
    uptake <- ex@lowbnd[which(ex@lowbnd < 0)]
  } 
  
  ub <- model@uppbnd
  lb <- model@lowbnd
  if(length(uptake) == 1)
    uptake <- rep(uptake, length(medium))
  
  lb[ex@react_pos] <- 0
  idx <- match(medium, react_id(model))
  if(sum(is.na(idx))>0 & verbose) cat("Not found in model:", medium[is.na(idx)], "\n")
  lb[na.omit(idx)] <- uptake[!is.na(idx)]
  model.constrainted <- changeBounds(model, react=1:length(lb), ub=ub, lb=lb)
  
  return(model.constrainted)
}

csource2 <- function(mod, cs, medium=NULL, verbose=FALSE, esource=F, GrowthExNa=0, db="seed"){
  if( sybil::SYBIL_SETTINGS("SOLVER") == "cplexAPI" ) solver_ok=1 else if( sybil::SYBIL_SETTINGS("SOLVER") == "glpkAPI" ) solver_ok=5
  
  if( length(medium) == 0 ){
    if ( verbose ) print("No medium specified, using glucose minimal medium!")
    minmedia <- fread("input/dat/MM_glu.csv")
    medium <- paste0("EX_",minmedia[,compounds],"_e0")
  }
  
  if( esource ){
    if( db=="seed"){
      mql <- "cpd15499[c0]"; mqn  <- "cpd15500[c0]"
      uql <- "cpd15561[c0]"; uqn  <- "cpd15560[c0]"
      nadh<- "cpd00004[c0]"; nad  <- "cpd00003[c0]"
      h   <- "cpd00067[c0]"
    }else{
      mql <- "mql8[c]"; mqn   <- "mqn8[c]"
      uql <- "u8h2[c]"; uqn   <- "u8[c]"
      nadh<- "nadh[c]"; nad   <- "nad[c]"
      h   <- "h[c]"
    }
    mod <- addReact(mod, "ESP1", met=c(mql,h,mqn), Scoef=c(-1,2,1), lb=0, ub=1000) # check if ubiquinone pool can be recycled 
    mod <- addReact(mod, "ESP2", met=c(uql,h,uqn), Scoef=c(-1,2,1), lb=0, ub=1000) # check if menaquinone pool can be recycled 
    mod <- addReact(mod, "ESP3", met=c(nadh,h,nad), Scoef=c(-1,1,1),lb=0, ub=1000) # check if nadh can be recycled
    mod <- changeObjFunc(mod, react=c("ESP1", "ESP2", "ESP3"), obj_coef=c(1,1,1))
  }
  med_withoutCS <- setdiff(medium, "EX_cpd00027_e0") # remove glucose
  if( db=="vmh"){
    dic <- fread("input/dat/SEED2VMH_translation_edited.csv", header=F)
    idx <- match(gsub("_e0","\\(e\\)",med_withoutCS), dic$V1)
    med_withoutCS <- dic$V2[idx]
  }
  dfcs <- data.frame()
  for( carbon in cs){
    if( !carbon %in% react_id(mod) ){
      dfcs <- rbind(dfcs, data.frame(csource=carbon, usage=FALSE, growth=GrowthExNa)) # if no exchange is avalable set set to default value
      next
    }
    med <- c(med_withoutCS, carbon)
    model <- set_diet(mod, med, verbose=verbose)
    sol <- sybil::optimizeProb(model, retOptSol=F)
    usage = sol$stat==solver_ok & round(sol$obj,3)>0
    growth= ifelse(sol$stat==solver_ok, sol$obj, 0)
    dfcs <- rbind(dfcs, data.frame(csource=carbon, usage, growth))
    
    if( verbose & db == "seed" ){
      idx.active <- which(sol$fluxes != 0)
      mod@react_name[idx.active]  
      idx.active.met <- idx.active[which(colSums(abs(mod@S[grep(str_extract(carbon, "cpd[0-9]+"), mod@met_id), idx.active]))!=0)]
      mod@react_name[idx.active.met]
      mod@react_id[idx.active.met]      
    }
  }
  return(dfcs)
}

library(foreach)
library(doParallel)
library(cplexAPI)
SYBIL_SETTINGS("SOLVER", "cplexAPI")
registerDoParallel(3)# 3 cores; detectCores()

op50.nad<-csource2(gapseq$OP50.RDS,cs=mysub$seed,esource=TRUE)  

op50.nad2<-merge(op50.nad,mysub,by.x = 'csource',by.y = 'seed',sort=F)
write.xlsx(op50.nad2,'op50.nad2.xlsx', sheetName = "op50.nad2", append = TRUE)
#Vibrio
HQB20.nad<-csource2(gapseq$`HQB-20.RDS`,cs=mysub$seed,esource=TRUE)
HQB21.nad<-csource2(gapseq$`HQB-21.RDS`,cs=mysub$seed,esource=TRUE)
HQB111.nad<-csource2(gapseq$`HQB-111.RDS`,cs=mysub$seed,esource=TRUE)
HQB114.nad<-csource2(gapseq$`HQB-114.RDS`,cs=mysub$seed,esource=TRUE)
HQB171.nad<-csource2(gapseq$`HQB-171.RDS`,cs=mysub$seed,esource=TRUE)

HQB20.nad2<-merge(HQB20.nad,mysub,by.x = 'csource',by.y = 'seed',sort=F)
HQB21.nad2<-merge(HQB21.nad,mysub,by.x = 'csource',by.y = 'seed',sort=F)
HQB111.nad2<-merge(HQB111.nad,mysub,by.x = 'csource',by.y = 'seed',sort=F)
HQB114.nad2<-merge(HQB114.nad,mysub,by.x = 'csource',by.y = 'seed',sort=F)
HQB171.nad2<-merge(HQB171.nad,mysub,by.x = 'csource',by.y = 'seed',sort=F)

write.xlsx(HQB20.nad2,'vibrio.biolog.xlsx', sheetName = "HQB20.nad2", append = TRUE)
write.xlsx(HQB21.nad2,'vibrio.biolog.xlsx', sheetName = "HQB21.nad2", append = TRUE)
write.xlsx(HQB111.nad2,'vibrio.biolog.xlsx', sheetName = "HQB111.nad2", append = TRUE)
write.xlsx(HQB114.nad2,'vibrio.biolog.xlsx', sheetName = "HQB114.nad2", append = TRUE)
write.xlsx(HQB171.nad2,'vibrio.biolog.xlsx', sheetName = "HQB171.nad2", append = TRUE)

#plot
library(pheatmap)
library(data.table)
library(readxl)
biolog <- read_excel("input/carbon_test_sum.xlsx", sheet = "BIOLOG-SUM")
predict <- read_excel("input/carbon_test_sum.xlsx", sheet = "predict_sum_1.2_150")
biolog1<-biolog[,c(1:6)]
biolog1 = biolog1[,-c(1)]

rownames(biolog1)=biolog$seed_name[]
predict1 = predict[,-c(1)]
rownames(predict1)=predict$seed_name[]
#y = read.csv("Figure4_Data_Metadata.csv",check.names=FALSE)

#y1 = y[,c(4,3,1)]
#rownames(y1)=y$`GENRE Filename`
colfunc1 <- colorRampPalette(c("white","#FFFFCC"))(5)
colfunc11 <- colorRampPalette(c("#FF9966","#CC0000"))(18)
colfunc2 <- colorRampPalette(c("white","black"))
bk <- c(seq(-1,1.5,by=0.5),seq(1.6,15,by=1))
cols<-c(colfunc1,'#FFFFCC',colfunc11)
biolog <- 
  pheatmap(biolog1,
           cellwidth = 15,
           cellheight = 10,
           #annotation_colors = L,
           #annotation_col = biolog1,
           scale="none",
           color = cols,
           border_color = 'gray',
           filename="figure4.pdf",
           show_colnames = T,
           show_rownames = F,
           #labels_col = ColumnLabels,
           #labels_row = rownames(biolog1),
           fontsize_col = 10,
           fontsize_row = 10,
           #clustering_method = "average",
           #clustering_distance_rows = "euclidean",
           #clustering_distance_cols = "euclidean",
           width=10,
           height=15,
           cluster_cols=F,
           cluster_rows=F,
           angle_col= 45,
           breaks=bk,
           legend_breaks=seq(0,16,2))

predict <-
  pheatmap(predict1,
           cellwidth = 15,
           cellheight = 10,
           #annotation_colors = L,
           #annotation_col = biolog1,
           scale="none",
           color = colfunc2(100),
           border_color = 'gray',
           filename="figure3.pdf",
           show_colnames = TRUE,
           #labels_col = ColumnLabels,
           labels_row = rownames(predict1),
           fontsize_col = 10,
           fontsize_row = 10,
           #clustering_method = "average",
           #clustering_distance_rows = "euclidean",
           #clustering_distance_cols = "euclidean",
           width=10,
           height=15,
           cluster_cols=F,
           cluster_rows=F,
           angle_col= 45,
           legend=F,
           fontfamily= "serif")

##Fig. S11
#in silico prediction
library(data.table)
library(sybil)
library(ggplot2)
library(reshape2)
library(stringr)

cs.sub <- fread("input/dat/sub2pwy.csv") # gapseq file containing carbon sources + exchange ids

# file from: http://protraits.irb.hr/data.html
protraits.db <- fread("input/protraits/HQ_bac_db.txt") #Manually select 13 bacteria in the database in the protraits database 
cs.src <- colnames(protraits.db)[3:111] #The first 111 in the database are metabolites
protraits.cs <- protraits.db[protraits.db[, Reduce(`|`, lapply(.SD, `!=`, "?")),.SDcols = cs.src],] # organisms with at least one carbon source prediction

# load models
#getwd()
#setwd('C:/Users/Gemini/Desktop/Code/Fig. 3/input/protraits/fill_RDS_150')
#files <- list.files(path = 'C:/Users/Gemini/Desktop/Code/Fig. 3/input/protraits/fill_RDS_150', pattern = ".RDS$")
#HQBiome_gap1.2_150 <- lapply(files, readRDS)
#names(HQBiome_gap1.2_150)<-c(files)
#
#gapseq <- HQBiome_gap1.2_150
gapseq <- readRDS("input/HQbiome.RDS")

gapseq.id <- gsub(".RDS","",names(gapseq))

#
# matching of carbon source names
#
dt.cs  <- data.table(name=gsub("_"," ",cs.src))
idx <- match(dt.cs$name, tolower(cs.sub$name))
idx[is.na(idx)] <- match( dt.cs$name[is.na(idx)], tolower(cs.sub$altname))
idx[is.na(idx)] <- match( dt.cs$name[is.na(idx)], tolower(str_remove(cs.sub$name, "^.-")))
idx[is.na(idx)] <- match( str_remove(dt.cs$name[is.na(idx)],"^.-"), tolower(cs.sub$name))
idx[is.na(idx)] <- match( str_remove(dt.cs$name[is.na(idx)],"^.-"), tolower(cs.sub$altname))
dt.cs[,seed.name:=cs.sub$name[idx]]
dt.cs[,seed.ex:=cs.sub$exid_seed[idx]]
dt.cs[,vmh.ex:=cs.sub$exid[idx]]
dt.cs <- dt.cs[!is.na(seed.ex) & !seed.ex=="" & !is.na(vmh.ex) & !vmh.ex==""]
# manual modifications
dt.cs[name=="dextrin", seed.ex:="EX_cpd11976_e0"] # dextrin not really distingushable from maltodextrin and dextrin=maltodextrin for seed db


# refseq genomes for protraits organism
#setwd('./../..')
tax2genome <- fread("input/protraits/protraits_genomes.txt", fill = T, header=F, sep="\t")
tax2genome<- tax2genome[1:75,]

set_diet <- function(model, medium, uptake=-100,verbose=TRUE){
  ex <- findExchReact(model)
  ex <- ex[grep("^EX_",ex@react_id),] # reduce to real exchange reactions only
  if(is.null(medium)){
    medium <- ex@react_id[which(ex@lowbnd < 0)]
    uptake <- ex@lowbnd[which(ex@lowbnd < 0)]
  } 
  
  ub <- model@uppbnd
  lb <- model@lowbnd
  if(length(uptake) == 1)
    uptake <- rep(uptake, length(medium))
  
  lb[ex@react_pos] <- 0
  idx <- match(medium, react_id(model))
  if(sum(is.na(idx))>0 & verbose) cat("Not found in model:", medium[is.na(idx)], "\n")
  lb[na.omit(idx)] <- uptake[!is.na(idx)]
  model.constrainted <- changeBounds(model, react=1:length(lb), ub=ub, lb=lb)
  
  return(model.constrainted)
}

csource2 <- function(mod, cs, medium=NULL, verbose=FALSE, esource=F, GrowthExNa=0, db="seed"){
  if( sybil::SYBIL_SETTINGS("SOLVER") == "cplexAPI" ) solver_ok=1 else if( sybil::SYBIL_SETTINGS("SOLVER") == "glpkAPI" ) solver_ok=5
  
  if( length(medium) == 0 ){
    if ( verbose ) print("No medium specified, using glucose minimal medium!")
    minmedia <- fread("input/dat/MM_glu.csv")
    medium <- paste0("EX_",minmedia[,compounds],"_e0")
  }
  
  if( esource ){
    if( db=="seed"){
      mql <- "cpd15499[c0]"; mqn  <- "cpd15500[c0]"
      uql <- "cpd15561[c0]"; uqn  <- "cpd15560[c0]"
      nadh<- "cpd00004[c0]"; nad  <- "cpd00003[c0]"
      h   <- "cpd00067[c0]"
    }else{
      mql <- "mql8[c]"; mqn   <- "mqn8[c]"
      uql <- "u8h2[c]"; uqn   <- "u8[c]"
      nadh<- "nadh[c]"; nad   <- "nad[c]"
      h   <- "h[c]"
    }
    mod <- addReact(mod, "ESP1", met=c(mql,h,mqn), Scoef=c(-1,2,1), lb=0, ub=1000) # check if ubiquinone pool can be recycled 
    mod <- addReact(mod, "ESP2", met=c(uql,h,uqn), Scoef=c(-1,2,1), lb=0, ub=1000) # check if menaquinone pool can be recycled 
    mod <- addReact(mod, "ESP3", met=c(nadh,h,nad), Scoef=c(-1,1,1),lb=0, ub=1000) # check if nadh can be recycled
    mod <- changeObjFunc(mod, react=c("ESP1", "ESP2", "ESP3"), obj_coef=c(1,1,1))
  }
  med_withoutCS <- setdiff(medium, "EX_cpd00027_e0") # remove glucose
  if( db=="vmh"){
    dic <- fread("input/dat/SEED2VMH_translation_edited.csv", header=F)
    idx <- match(gsub("_e0","\\(e\\)",med_withoutCS), dic$V1)
    med_withoutCS <- dic$V2[idx]
  }
  dfcs <- data.frame()
  for( carbon in cs){
    if( !carbon %in% react_id(mod) ){
      dfcs <- rbind(dfcs, data.frame(csource=carbon, usage=FALSE, growth=GrowthExNa)) # if no exchange is avalable set set to default value
      next
    }
    med <- c(med_withoutCS, carbon)
    model <- set_diet(mod, med, verbose=verbose)
    sol <- sybil::optimizeProb(model, retOptSol=F)
    usage = sol$stat==solver_ok & round(sol$obj,3)>0
    growth= ifelse(sol$stat==solver_ok, sol$obj, 0)
    dfcs <- rbind(dfcs, data.frame(csource=carbon, usage, growth))
    
    if( verbose & db == "seed" ){
      idx.active <- which(sol$fluxes != 0)
      mod@react_name[idx.active]  
      idx.active.met <- idx.active[which(colSums(abs(mod@S[grep(str_extract(carbon, "cpd[0-9]+"), mod@met_id), idx.active]))!=0)]
      mod@react_name[idx.active.met]
      mod@react_id[idx.active.met]      
    }
  }
  return(dfcs)
}

library(foreach)
library(doParallel)
library(cplexAPI)
SYBIL_SETTINGS("SOLVER", "cplexAPI")
registerDoParallel(3) # 3 cores; detectCores()

dt.cs.predict <- foreach(i=1:nrow(protraits.cs), .combine=rbind) %dopar%{
  cat("\r",i,"/",nrow(protraits.cs))
  tax.id <- protraits.cs$Tax_ID[i]
  genome.id <- tax2genome[V1==tax.id,V4]
  if( !is.na(genome.id) ){
    idx.gapseq     <- match(genome.id, gsub("GCF","GCA", gapseq.id)) # for old genome files
    #idx.modelseed <- match(genome.id, gsub("GCF","GCA", modelseed.id))
    #idx.carveme   <- match(genome.id, gsub("GCF","GCA", carveme.id)) 
    if( !all(is.na(idx.gapseq))){
      #if( !all(is.na(idx.gapseq)) & !all(is.na(idx.modelseed)) & !all(is.na(idx.carveme)) ){
      mod.gapseq    <- gapseq[[   na.omit(idx.gapseq)[1] ]]
      #mod.modelseed <- modelseed[[na.omit(idx.modelseed)[1] ]]
      #mod.carveme   <- carveme[[  na.omit(idx.carveme)[1] ]]
      
      cs.pot <- protraits.cs[Tax_ID==tax.id, c(1:111)]
      cs <- data.table::melt(cs.pot, id.vars = c("Organism_name", "Tax_ID"))[value!="?" & !is.na(value)]
      cs.sel <- dt.cs[name %in% cs$variable]  
      if( nrow(cs.sel)>0 ){
        growth.gapseq    <- csource2(mod.gapseq   , cs=cs.sel$seed.ex,esource=TRUE)
        #growth.modelseed <- csource2(mod.modelseed, cs=cs.sel$seed.ex,esource=TRUE)
        #growth.carveme   <- csource2(mod.carveme  , cs=cs.sel$vmh.ex,esource=TRUE, db="vmh")
        
        #cs.sel[,org:= genome.id]; 
        cs.sel[,org:= tax.id]; cs.sel[,org.name:= protraits.cs$Organism_name[i]]
        cs.sel[,org.gapseq:= gapseq.id[na.omit(idx.gapseq)[1] ]]
        cs.sel[,protraits:=cs[match(cs.sel$name, cs$variable), value==1]]
        cs.sel[,gapseq:=growth.gapseq$usage]
        #cs.sel[,modelseed:=growth.modelseed$usage]
        #cs.sel[,carveme:=growth.carveme$usage]
        cs.sel
      }
    } else data.table()
  } else data.table()
}

dt.cs.predict.melt <- data.table::melt(dt.cs.predict[,.(org,protraits,gapseq,name,seed.ex)], id.vars = c("org","protraits", "name", "seed.ex"), variable.name = "method")
dt.cs.predict.melt[value==T & protraits==T, validation:="TP"]
dt.cs.predict.melt[value==F & protraits==F, validation:="TN"]
dt.cs.predict.melt[value==T & protraits==F, validation:="FN"]
dt.cs.predict.melt[value==F & protraits==T, validation:="FP"]

dt.cs.predict[gapseq==T & protraits==T, validation:="TP"]
dt.cs.predict[gapseq==F & protraits==F, validation:="TN"]
dt.cs.predict[gapseq==T & protraits==F, validation:="FN"]
dt.cs.predict[gapseq==F & protraits==T, validation:="FP"]

write.csv(dt.cs.predict,'dt.cs.predict.csv')

table(dt.cs.predict.melt[method=="gapseq",validation]); round(.Last.value / dt.cs.predict.melt[method=="gapseq",.N], 2) 

dt.cs.predict.plot <- dt.cs.predict.melt[,.N, by=.(method, validation)]
dt.cs.predict.plot$method <- factor(dt.cs.predict.plot$method, levels = c("carveme", "gapseq", "modelseed"), labels= c("CarveMe", "gapseq", "ModelSEED"), ordered = T)
dt.cs.predict.plot$validation <- factor(dt.cs.predict.plot$validation, levels = c("FN", "FP", "TN", "TP"),labels = c("False negative", "False positive", "True negative", "True positive"), ordered = T)
dt.cs.predict.plot$validation <- factor(dt.cs.predict.plot$validation, levels = c("True positive", "True negative", "False positive", "False negative"))
ggplot(data = dt.cs.predict.plot) + geom_bar(stat="identity",aes(x=validation, y=N,fill=method),position="dodge2", color="black") + 
  #geom_text(aes(x=validation, y=N, label=N), position = position_dodge2(width = 1)) + # wrong order?????
  xlab("") + ylab("Carbon source validation") + labs(fill="") + theme_bw(base_size = 14) + 
  scale_y_continuous(limits=c(0, 30), expand = c(0, 0)) + #scale_x_continuous(expand = c(0, 0)) + 
  scale_fill_discrete() +
  #scale_fill_manual(values=c("#black", "#e41a1c", "#FFB200")) +
  theme(strip.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank())
#ggsave("protraits.pdf", width=6, height=5)


# good/bad predicted carbon sources
dt.quali.all <- dt.cs.predict.melt[,list(FN=sum(validation=="FN"), FP=sum(validation=="FP"),
                                         TN=sum(validation=="TN"), TP=sum(validation=="TP")), by=name]
dt.quali.all[,accuracy:=(TN+TP)/(FN+FP+TN+TP)]
dt.quali.all[,sensitivity:=(TP)/(TP+FN)]
dt.quali.all[,specificity:=(TN)/(TN+FP)]
dt.quali.all[,N:=TN+TP+FP+FN]
dt.quali.all[order(accuracy)][FN+FP+TN+TP>100]

# plot in graphpad using input/carbon_test_sum.xlsx, sheet=1

#Fig. S28a Functional correlation between natural microbiomes and HQbiome isolates.
library(readxl)
library(ggplot2)
library(ggpubr)
library(paletteer) 
library(ggbreak)
db<-read_xlsx('input/path_abun_unstrat.xlsx',sheet = 3)
main_theme <- theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    axis.line=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=10),
                    #legend.position="top",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    #text=element_text(family="sans")
)
db$Group<-factor(db$Group,levels = c('Amino acid metabolism','Aromatic compounds metabolism',
                                     'Carbohydrates metabolism','Cell structure metabolism',
                                     'Cofactor metabolism','Energy metabolism','Lipid metabolism',
                                     'Nucleotide metabolism','Secondary metabolite metabolism',
                                     'Others'))
cor<-ggplot(db,aes(x = gapseq, y = `mean proportion`)) +
  geom_point(aes(color = Group), alpha=1,size = 2) +
  main_theme+
  geom_smooth(method="lm", col="darkgrey", se =F)+ 
  xlab("Pathway proportion of Gapseq (%)")+
  ylab("Pathway proportion of Picrust2 (%)")+
  #scale_y_continuous(breaks = c(1,3,5,7)) +
  #theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5)) +
  #scale_color_identity()+
  #scale_y_break(c(4, 2),expand = F,scales='0.1',space=0,
  #              #expand=expansion(add = c(0, 1))
  #              )+
  scale_colour_manual(name = "Group", values = paletteer_d("ggprism::floral2"))+ 
  #guides(colour = guide_legend(override.aes = list(shape = 15,size=4)))+ 
  stat_cor(aes(x = gapseq, y = `mean proportion`,group = "",label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           method = "spearman", cor.coef.name = "R[s]", 
           label.x = 30, label.y = 1.5, size = 3.5)#+ 
cor
ggsave('cor.pdf',cor,width = 7,height = 5)
#geom_rect(xmin = 75, xmax = 100, ymin = 5, ymax = 7, size = 0.5, fill = "#00000000", color = "darkgrey", linetype=2 ) #画矩形

#Fig. S28b was created by ANIclustermap using ./input/fastANI.txt

#Fig. S32a PCoA plot depicting functional distances (n = 74) between sequenced genomes based on metabolic networks
library(readxl)
library(ggplot2)
library(ggalt)
library(dplyr)
#import
path <- read_excel("input/pcoa.xlsx")
path1<-path[,-1]
rownames(path1)<-path$ID
tax<-read_excel("input/pcoa.xlsx",sheet=2)

#refer to https://github.com/garridoo/atsphere_wgs/blob/master/functional_MDS.R
distance<-dist(path1, method = 'euclidean',p = 2)
distance.dt<-as.matrix(distance)

pcoa <- cmdscale(distance.dt, k=2, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")
points$`Phylum/Class` <- tax$`Phylum/Class`[match(rownames(points), tax$ID)] 
points$`Gram strain` <- tax$`Gram strain`[match(rownames(points), tax$ID)] 

labs=c('Actinobacteria','Bacteroidetes','Firmicutes','Alphaproteobacteria','Gammaproteobacteria')
points$`Phylum/Class` <- factor(points$`Phylum/Class`, levels = labs)
cols <- c(Actinobacteria= '#EEE685',Bacteroidetes = "#4F94CD", Firmicutes= '#EE7600', 
          Alphaproteobacteria = "#00CD00",Gammaproteobacteria = "#90EE90")
#plot
main_theme <- theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    axis.line=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=10),
                    #legend.position="top",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    #text=element_text(family="sans")
)
pcoa <- ggplot(points, aes(x=x, y=y, color=`Phylum/Class`, shape=`Gram strain`))+
  geom_point(alpha=1,size=3)+
  scale_colour_manual(values=cols)+
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+
  main_theme+
  geom_encircle(data = points[points$`Phylum/Class`=='Actinobacteria',],fill='#EEE685',alpha = 0.1, show.legend = F,spread=0.002)+
  geom_encircle(data = points[points$`Phylum/Class`=='Bacteroidetes',],fill='#4F94CD',alpha = 0.1, show.legend = F,spread=0.002)+
  geom_encircle(data = points[points$`Phylum/Class`=='Firmicutes',],fill='#EE7600',alpha = 0.1, show.legend = F,spread=0.002)+
  geom_encircle(data = points[points$`Phylum/Class`=='Alphaproteobacteria',],fill='#00CD00',alpha = 0.1, show.legend = F,spread=0.002)+
  geom_encircle(data = points[points$`Phylum/Class`=='Gammaproteobacteria',],fill='#90EE90',alpha = 0.1, show.legend = F,spread=0.002)+
  theme(legend.position = c(0.16, 0.8))+ theme(legend.key.size = unit(0.03, "inches"))
pcoa

#correlation between pathways and pcoa axis
corr<-as.data.frame(cor(path1,points[,c(1,2)],method = 'spearman'))
library(data.table)
meta <- fread("input/meta_pwy.tbl")
corr$id<-rownames(corr)
corr<-as.data.table(corr)
corr[,name:=meta$name[match(id, meta$id)]]
write.csv(corr,'cor-spearman.csv')


#Fig. S32b
phylum <- c('Actinobacteria',#'Bacteroidetes',
            'Firmicutes',
            'Alphaproteobacteria','Gammaproteobacteria')
levels(phylum)<-c('Actinobacteria',#'Bacteroidetes',
                  'Firmicutes',
                  'Alphaproteobacteria','Gammaproteobacteria')
df <- data.frame(phylum=NULL, distances=NULL)
pb <- txtProgressBar(min=1, max=length(phylum), style=3)

library(reshape2)
rm(df)
df<-setNames(melt(distance.dt), c('x', 'y', 'distance'))
tax.id<-tax[,c(1,33)]
df<-merge(df,tax.id,by.x='x',by.y = 'ID',all.x = T, sort=F)
df<-merge(df,tax.id,by.x='y',by.y = 'ID',all.x = T, sort=F)
colnames(df)<-c('x', 'y', 'distance','tax1','tax2')
rm(df1)
library(dplyr)
df1 = df %>% filter_(~tax1==tax2) %>% filter(distance != 0)
#write.csv(df1,'df1.csv')
library(stringr)
rm(df2)
df2<-mutate(df1,num.x=gsub('HQB-','',df1$x))%>%mutate(num.y=gsub('HQB-','',df1$y))
df2$num.x<-gsub('OP','',df2$num.x)
df2$num.y<-gsub('OP','',df2$num.y)
df2$filt<-as.numeric(df2$num.x)+as.numeric(df2$num.y)
df2$filt2<-as.numeric(df2$num.x)*as.numeric(df2$num.y)
df2$filt3<-as.numeric(df2$filt)*as.numeric(df2$filt2)
df3<-df2[!duplicated(df2$filt3),] #1073
df3<-df3[-which(df3$tax1=="Bacteroidetes"),]  

medians <- aggregate(df3$distance, by=list(df3$tax1), FUN=median)
order <- medians[sort(medians[, 2], index.return=T, decreasing=F)$ix, 1]
df3$tax1 <- factor(df3$tax1, levels=order)
dodge = position_dodge(width=0.8)
jdodge = position_jitterdodge(jitter.width = 0.1, jitter.height = 0, dodge.width = 0.8)
p1<-ggplot(df3, aes(x=tax1, y=distance, color=tax1)) +
  geom_jitter( size=1.7, alpha=0.25,width = 0.1) +
  geom_boxplot(alpha=0.8, #outlier.size=.2, 
               #size=1,#width = 0.25,
               width=0.4,cex=0.8,
               #position=dodge,
               color=c(Alphaproteobacteria = "#00CD00",Actinobacteria= '#EEE685',  Firmicutes= '#EE7600',Gammaproteobacteria = "#90EE90"), 
               fill=NA
  ) +
  geom_violin(alpha = 1,
              fill = NA,
              #colour = NA,
              #position=dodge,
              width=0.9, cex=0.8)+
  labs(x="", y="functional distance") +
  scale_colour_manual(values=cols, guide='none') +
  coord_flip() +
  #scale_y_continuous(limits=c(.25, .75)) +
  main_theme +
  theme(axis.text.y = element_text(size=8))+
  ylab('Pair-wise functional distance')

p1

#sig test
library(rstatix)
#································································································································
#shapiro.test and bartlett.test
shapiro.test(df3$distance)  #W = 0.98811, p-value = 1.23e-07
bartlett.test(distance ~ tax1, data = df3) #Bartlett's K-squared = 7.2628, df = 3, p-value = 0.06398 
#·······························································································································

stat_test <- df3 %>% 
  wilcox_test(distance ~ tax1, p.adjust.method = "fdr",alternative = "two.sided") %>% 
  add_xy_position(x = "tax1")
library(ggpubr)
p1<-p1+stat_pvalue_manual(stat_test, label = "p.adj.signif", step.increase = 0, hide.ns = TRUE, 
                          tip.length = 0, label.size = 5,coord.flip = TRUE,
)+ scale_y_continuous(limits=c(6,25))+
  stat_compare_means(label.x=0.5, label.y = 21,size=3)
p1 

p2<-ggplot(df3, aes(x=distance)) +
  geom_histogram(size=.5, alpha=1, color="grey", fill="grey", binwidth=.01) +
  #labs(title="Functional diversity within Phylum/Class", x="") +
  main_theme +
  theme(legend.position="none", #axis.text.y=element_text(size=8),
        #title=element_text(size=6)
        plot.title = element_text(size = 6),
        axis.text.x= element_blank(),axis.ticks.x=element_blank()
  )+
  ylab('Counts')+xlab(NULL)+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(limits=c(6,25))

p2

library(ggpubr)
b<-ggarrange(p2,p1,ncol=1,align = 'hv',
             heights=c(2, 6)
)
b

#Fig. S32c
#Correlation between pairwise similarities in 16S rRNA sequences and metabolic networks 
merge <- read_excel("input/dist_phy.xlsx",sheet='merge')
scatter <- read_excel("input/dist_phy.xlsx",sheet='scatter')
cor16sMeat<-ggplot(merge,aes(x = metabolic,y = phylogene)) +  
  geom_point(alpha = 0.2,color='#404B55')+ 
  scale_x_continuous(name='Pairwise metabolic identity', limits = c(58,100)) + 
  scale_y_continuous(name='Pairwise 16S rRNA identity',limits = c(58,100)) + 
  theme(axis.text=element_text(size=10, color = 'black', face="bold"), 
        axis.title=element_text(size=20,face="bold"), 
        panel.grid.minor = element_line(colour = "#FFCC00", size=0.25, linetype = "dashed"), 
        panel.grid.major = element_line(colour = "#FFCC00", size=0.4))+ 
  geom_abline(color = "black", lty=2, size = 0.4) +  
  geom_point(data=scatter,  mapping=aes(metabolic,phylogene), color = '#D94D26',alpha = 0.5)+
  theme_own()+
  stat_cor(aes(x = metabolic, y = phylogene,label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           method = "spearman", cor.coef.name = "R[s]", size = 4.5,label.x = 71, label.y = 69)+
  theme(panel.border = element_blank())+ 
  theme(axis.ticks=element_blank())  +
  annotate(geom = "rect", xmin = 71, xmax = 89,
           ymin = 96,ymax = 100, alpha=0,
           color="#D94D26",
           lty="dashed")+ 
  theme(aspect.ratio=1)  

#geom_smooth(method="lm", col="#E6E6FA", se =F) 
cor16sMeat

cor16sMeat2<-ggMarginal(cor16sMeat,type='density',xparams=list(fill='grey',size=0.5,alpha = 0.5),
                        yparams=list(fill='grey',size=0.5,alpha = 0.5))
cor16sMeat2


#Fig. S32d
#pvclust
library(pvclust)
library(dendextend)
pv.db <- read_excel("input/dist_phy.xlsx",sheet='merge150')
pv.db<-pv.db[,-1]
meta_clust <- pvclust(pv.db, method.dist="euclidean", method.hclust="average", nboot=1000, parallel=TRUE)
meta_clust1<-meta_clust

#origin plot
plot(meta_clust,print.num = F,print.pv=c('au','bp'),
     cex.pv=0.7, font.pv=2, main='',xlab='',sub='',cex=0.8, font=2,
     font.lab =2, cex.lab = 1,tip.color=1)
#pvrect(meta_clust, alpha=0.97)

#main plot
plot(meta_clust,print.num = F,print.pv=c('au'),
     col.pv = c(au=7),cex.pv=0.8,font.pv=2,main='',xlab='',sub='',cex=0.8,font=2,
     font.lab =2, cex.lab = 1)
#colours were change in AI

#Fig. S27
#Fig. S27a
library(tidyverse)
library(readxl)
COG<-list.files('input/cogfinal',pattern = 'final.csv', full.names = T) %>%
  lapply(read.csv,header = FALSE) %>%
  reduce(full_join,by=c('V1'))

#change row names
library(stringr)
colnames(COG)<- c('COG',
                  str_replace(list.files('input/cogfinal/',pattern = '*csv', full.names = F),'.emapper.annotations-COG-final.csv', ''))  
library("xlsx")
#write.xlsx(COG,'cog.xlsx')
cog2=as.data.frame(COG)

cog2<-cog2[,-1]
rownames(cog2)<-COG$COG
cog2<-cog2[-24,]
cog3<-as.data.frame(100 * proportions(as.matrix(cog2), 2))
#write.xlsx(cog3,'cog_sum.xlsx')
library(tidyverse)
library(matrixStats)

cog4<-cog3%>%mutate(mean=rowMeans(.))
row_sd <- apply(cog3, 1, sd)
cog4<- cbind(cog4, sd = row_sd)
#write.xlsx(cog4,'cog_sum.xlsx')
cog5<-cog4[,c('mean','sd')]

p.cog<-ggplot(cog5,aes(row.names(cog5),mean))+
  geom_col(aes(fill='black'))+
  geom_errorbar(aes(row.names(cog5),ymin=mean-sd,ymax=mean+sd,width=.25))+
  ylab('Average percentage (%)')+xlab(NULL)+theme_classic()+
  scale_y_continuous(expand = c(0,0),limits = c(0,25))+
  theme(legend.position = "none")+
  theme(panel.grid.major.y=element_line(colour='grey'),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values='#8B1A1A')
p.cog
ggsave('cog.pdf',p.cog,width = 5, height = 2.5)

#Fig. S27b
#-----merge the number of reaction, gene and metabolites from 74 metabolic model-------------------------------
library(data.table)
i=1
sum<-data.frame()
gapseq <- readRDS("input/HQbiome.RDS")
for (i in 1:74) {
  
  mod<-gapseq[[i]]
  sum<-rbind(sum,data.frame(name=mod@mod_name,reactions=mod@react_num,
                            metabolites=mod@met_num,genes=length(mod@allGenes)))
  
}

write.csv(sum,'modelSum.csv')
##-------------reload excel from modelSum.csv------------------------------------------------------------------------
library(readxl)
library(ggplot2)
library(ggpubr)   
library(rstatix)  
model.sum<-read_xlsx('input/modelSum.xlsx')
cols <- c(Actinobacteria= '#EEE685',Bacteroidetes = "#4F94CD", Firmicutes= '#EE7600', 
          Alphaproteobacteria = "#00CD00",Gammaproteobacteria = "#90EE90")
labs=c('Actinobacteria','Bacteroidetes','Firmicutes','Alphaproteobacteria','Gammaproteobacteria')
model.sum$Tax <- factor(model.sum$Tax, levels = labs)
main_theme <- theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    axis.line=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=10),
                    #legend.position="top",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    #text=element_text(family="sans")
)
#plot
#genomes
p1<-ggplot(model.sum, aes(x=Tax, y=genome, color=Tax)) +
  geom_jitter( size=1.7, alpha=0.25,width = 0.1) +
  geom_boxplot(alpha=0.3, #outlier.size=1, 
               #size=1,#width = 0.25,
               width=0.4,cex=0.5,
               #position=dodge,
               #color=c(Alphaproteobacteria = "#00CD00",Actinobacteria= '#EEE685',Gammaproteobacteria = "#90EE90", Firmicutes= '#EE7600'), 
               fill=cols, 
  ) +
  #geom_violin(alpha = 1,
  #            fill = NA,
  #            #colour = NA,
  #            #position=dodge,
  #            width=0.9, cex=0.8)+
  #labs(x="", y="Reactions") +
  scale_colour_manual(values=cols, guide='none') +
  #coord_flip() +
  #scale_y_continuous(limits=c(.25, .75)) +
  main_theme +
  theme(axis.text.y = element_text(size=8))+
  rremove("x.text")+rremove("xlab")+
  geom_hline(yintercept = mean(model.sum$genome), color="grey",linetype = "dashed")+
  ylab('Genome size')
p1

#sig test
#shapiro.test and bartlett.test
shapiro.test(model.sum$growth)   #W = 0.78244, p-value = 4.237e-09
bartlett.test(genome ~ Tax, data = model.sum)  #Bartlett's K-squared = 15.937, df = 4, p-value = 0.003104

model.sum %>% kruskal_test(genome ~ Tax)  # p=0.00704
stat_test.genome<-model.sum %>% 
  dunn_test(genome ~ Tax, p.adjust.method = "fdr")%>% 
  remove_ns() %>% add_xy_position(x = "Tax",step.increase = 0.02,dodge=0.01)
# Actinobacteria      Gammaproteobacteria     7    44        53 0.004 0.021 * 
# Alphaproteobacteria Gammaproteobacteria     8    44        58 0.002 0.018 * 
p1<-p1+stat_compare_means(label.x.npc= "left", label.y.npc = "bottom",size=2)+
  stat_pvalue_manual(stat_test.genome, label = "p.adj.signif", step.increase = 0, hide.ns = TRUE,
                     tip.length = 0, label.size = 3)
p1

#genes
p2<-ggplot(model.sum, aes(x=Tax, y=genes, color=Tax)) +
  geom_jitter( size=1.7, alpha=0.25,width = 0.1) +
  geom_boxplot(alpha=0.3, #outlier.size=1, 
               #size=1,#width = 0.25,
               width=0.4,cex=0.5,
               #position=dodge,
               #color=c(Alphaproteobacteria = "#00CD00",Actinobacteria= '#EEE685',Gammaproteobacteria = "#90EE90", Firmicutes= '#EE7600'), 
               fill=cols, 
  ) +
  #geom_violin(alpha = 1,
  #            fill = NA,
  #            #colour = NA,
  #            #position=dodge,
  #            width=0.9, cex=0.8)+
  #labs(x="", y="Reactions") +
  scale_colour_manual(values=cols, guide='none') +
  #coord_flip() +
  #scale_y_continuous(limits=c(.25, .75)) +
  main_theme +
  theme(axis.text.y = element_text(size=8))+
  rremove("x.text")+rremove("xlab")+
  geom_hline(yintercept = mean(model.sum$genes), color="grey",linetype = "dashed")+
  ylab('Number of genes')
p2

#sig test
#shapiro.test and bartlett.test
shapiro.test(model.sum$genes)  # p=0.049  
bartlett.test(genes ~ Tax, data = model.sum) ##p=0.02592,  
#genes
model.sum %>% kruskal_test(genes ~ Tax)  # p=0.00471
stat_test.genes<-model.sum %>%
  dunn_test(genes ~ Tax, p.adjust.method = "fdr")%>% 
  remove_ns() %>% add_xy_position(x = "Tax",step.increase = 0.02,dodge=0.01)

p2<-p2+stat_compare_means(label.x.npc= "left", label.y.npc = "bottom",size=2)+
  stat_pvalue_manual(stat_test.genes, label = "p.adj.signif", step.increase = 0, hide.ns = TRUE,
                     tip.length = 0, label.size = 3)
p2

#reactions
p3<-ggplot(model.sum, aes(x=Tax, y=reactions, color=Tax)) +
  geom_jitter( size=1.7, alpha=0.25,width = 0.1) +
  geom_boxplot(alpha=0.3, #outlier.size=1, 
               #size=1,#width = 0.25,
               width=0.4,cex=0.5,
               #position=dodge,
               #color=c(Alphaproteobacteria = "#00CD00",Actinobacteria= '#EEE685',Gammaproteobacteria = "#90EE90", Firmicutes= '#EE7600'), 
               fill=cols, 
  ) +
  #geom_violin(alpha = 1,
  #            fill = NA,
  #            #colour = NA,
  #            #position=dodge,
  #            width=0.9, cex=0.8)+
  #labs(x="", y="Reactions") +
  scale_colour_manual(values=cols, guide='none') +
  #coord_flip() +
  #scale_y_continuous(limits=c(.25, .75)) +
  main_theme +
  theme(axis.text.y = element_text(size=8))+
  rremove("x.text")+rremove("xlab")+
  geom_hline(yintercept = mean(model.sum$reactions), color="grey",linetype = "dashed")+
  ylab('Number of reactions')
p3
#sig test 
#shapiro.test and bartlett.test
shapiro.test(model.sum$reactions)  #  p-value = 0.3603  
bartlett.test(reactions ~ Tax, data = model.sum) #p=0.0117
#genes
model.sum %>% kruskal_test(reactions ~ Tax)  # p=0.000506 
stat_test.reactions<-model.sum %>%
  dunn_test(reactions ~ Tax, p.adjust.method = "fdr")%>% 
  remove_ns() %>% add_xy_position(x = "Tax",step.increase = 0.02)

p3<-p3+stat_compare_means(label.x.npc= "left", label.y.npc = "bottom",size=2)+
  stat_pvalue_manual(stat_test.reactions, label = "p.adj.signif", step.increase = 0, hide.ns = TRUE,
                     tip.length = 0, label.size = 3)
p3

#metabolites
p4<-ggplot(model.sum, aes(x=Tax, y=metabolites, color=Tax)) +
  geom_jitter( size=1.7, alpha=0.25,width = 0.1) +
  geom_boxplot(alpha=0.3, #outlier.size=1, 
               #size=1,#width = 0.25,
               width=0.4,cex=0.5,
               #position=dodge,
               #color=c(Alphaproteobacteria = "#00CD00",Actinobacteria= '#EEE685',Gammaproteobacteria = "#90EE90", Firmicutes= '#EE7600'), 
               fill=cols, 
  ) +
  #geom_violin(alpha = 1,
  #            fill = NA,
  #            #colour = NA,
  #            #position=dodge,
  #            width=0.9, cex=0.8)+
  #labs(x="", y="Reactions") +
  scale_colour_manual(values=cols, guide='none') +
  #coord_flip() +
  #scale_y_continuous(limits=c(.25, .75)) +
  main_theme +
  theme(axis.text.y = element_text(size=8))+
  rremove("x.text")+
  rremove("xlab")+
  geom_hline(yintercept = mean(model.sum$metabolites), color="grey",linetype = "dashed")+
  ylab('Number of metabolites')
p4
#sig test
#shapiro.test and bartlett.test
shapiro.test(model.sum$metabolites)  #  p-value = 0.06107
bartlett.test(metabolites ~ Tax, data = model.sum) ##p=0.004917
#metabolites
model.sum %>% kruskal_test(metabolites ~ Tax)  # p=0.00108  
stat_test.metabolites<-model.sum %>%
  dunn_test(metabolites ~ Tax, p.adjust.method = "fdr")%>% 
  remove_ns() %>% add_xy_position(x = "Tax",step.increase = 0.02)

p4<-p4+stat_compare_means(label.x.npc= "left", label.y.npc = "bottom",size=2)+
  stat_pvalue_manual(stat_test.metabolites, label = "p.adj.signif", step.increase = 0, hide.ns = TRUE,
                     tip.length = 0, label.size = 3)
p4

#bacterial growth
p5<-ggplot(model.sum, aes(x=Tax, y=growth, color=Tax)) +
  geom_jitter( size=1.7, alpha=0.25,width = 0.1) +
  geom_boxplot(alpha=0.3, #outlier.size=1, 
               #size=1,#width = 0.25,
               width=0.4,cex=0.5,
               #position=dodge,
               #color=c(Alphaproteobacteria = "#00CD00",Actinobacteria= '#EEE685',Gammaproteobacteria = "#90EE90", Firmicutes= '#EE7600'), 
               fill=cols, 
  ) +
  #geom_violin(alpha = 1,
  #            fill = NA,
  #            #colour = NA,
  #            #position=dodge,
  #            width=0.9, cex=0.8)+
  #labs(x="", y="Reactions") +
  scale_colour_manual(values=cols, guide='none') +
  #coord_flip() +
  #scale_y_continuous(limits=c(.25, .75)) +
  main_theme +
  theme(axis.text.y = element_text(size=8))+
  rremove("x.text")+
  rremove("xlab")+
  geom_hline(yintercept = mean(model.sum$growth), color="grey",linetype = "dashed")+
  ylab('Growth rate')+
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      #vjust = 0.5
    ))
p5
#sig test
#shapiro.test and bartlett.test
shapiro.test(model.sum$growth)  #  p-value = 4.237e-09
bartlett.test(growth ~ Tax, data = model.sum) ## p=8.712e-10
#metabolites
model.sum %>% kruskal_test(growth ~ Tax)  # p=0.000536   
stat_test.growth<-model.sum %>%
  dunn_test(growth ~ Tax, p.adjust.method = "fdr")%>% 
  remove_ns() %>% add_xy_position(x = "Tax",step.increase = 0.02)

p5<-p5+stat_compare_means(label.x.npc= "left", label.y.npc = "bottom",size=2)+
  stat_pvalue_manual(stat_test.growth, label = "p.adj.signif", step.increase = 0, hide.ns = TRUE,
                     tip.length = 0, label.size = 3)
p5

#arrange
abcde<-ggarrange(p1,p2,p3,p4,p5,ncol = 1,heights =c(1,1,1,1,1.5),
                 align ='v')
abcde

ggsave('abcde.pdf',abcde, width = 7,height = 11)

#Fig. S27c
library(readxl)
library(ggtern)

uast <- read_excel("input/UAST.xlsx", sheet = "UAST b=150!")
uast$Strategy<-factor(uast$Strategy,levels = c('Stress tolerator','Competitor','Ruderal','Mix'))
uast$`Phylum/Class`<-factor(uast$`Phylum/Class`,levels = c('Actinobacteria','Bacteroidetes','Firmicutes',
                                                           'Alphaproteobacteria','Gammaproteobacteria'))

p3<-ggtern(data=uast, aes(x=R, y=S, z=C))+
  geom_mask()+
  geom_point(aes(shape=`Phylum`,
                 color=Strategy,
                 size=mean_5
  ),alpha=0.5,
  #position=position_jitter(h=0.01,w=0.01)
  )+
  scale_size(range = c(3, 10))+
  scale_color_manual(values  = c('#E31A1C','#228B22','#1F78B4','#8B658B'))+
  scale_shape_manual(values = c(15,18,17,16,12))+
  guides(size="none") +#theme_bw() +
  #theme_custom(
  #  base_size = 12,
  #  base_family = "",
  #  tern.plot.background = 'white',
  #  tern.panel.background = NULL,
  #  col.T = "#228B22",
  #  col.L = "#E31A1C",
  #  col.R = "#1F78B4",
  #  col.grid.minor = "white"
  #)+
  theme_rgbw()+
  theme(axis.text=element_blank(), axis.ticks=element_blank())
p3
ggsave('p3.pdf',p1,width = 7,height = 6)


#Fig. S27d
library(tidyverse)
library(tidyr)
library(corrplot)
library(RColorBrewer)
library(paletteer)
cor_df<-model.sum[,c(1,3:15,19:21,28:30,33,38,44)]%>%column_to_rownames(var ='name')
colnames(cor_df)<-c('Genome','tRNA','rRNA','Growth','Auxotrophy','Eps','Antibiotic','Catabolism',
                    'Siderophore','Degradation','SS','CS','RS',
                    'Reactions','Metabolites','Genes','ATP','Biomass','Essential Reactions','N2 growth','F23 growth',
                    'F23 survival')
M<-cor(cor_df, method = 'spearman')
res1 <- cor.mtest(cor_df, conf.level = .95)
corrplot.mixed(M)
corrplot(M,
         type = "upper",
         method = 'pie',
         order = "FPC",#hclust.method = "average",addrect=2
         col = paletteer_d("RColorBrewer::BrBG"),
         tl.cex = 0.7, 
         tl.col = "black", 
         tl.pos='td',
         tl.srt = 45,
         cl.align.text = "l",
         cl.pos='r',
         cl.length=9,
         p.mat = res1$p,insig='blank',
         diag=F,
         outline='lightgrey',
         mar = c(2,0,0,2),
         addgrid.col = "grey",
         
)

corrplot(M, type = "lower",
         method = 'square',
         order = "FPC",
         col = paletteer_d("RColorBrewer::BrBG"),
         tl.cex = 0.7,tl.col = "black",
         cl.align.text = "l",
         cl.pos='n',
         #cl.ratio = 0.1,
         outline='white',
         p.mat = res1$p,insig='blank',
         add=T,
         tl.pos='n',
         diag=F,
         addgrid.col = "grey"
)
#################################################################################3
#Fig. S8 wea created using BRIG with fasta files in ./input/BRIG/*
###########
#=======================================================================
#Fig. S23 Predicted metabolism of the HQbiome-1 and HQbiome-2 communities using BacArena.
library(BacArena)
library(data.table)
SYBIL_SETTINGS("SOLVER", "cplexAPI")
library(Biostrings)
library(stringr)
library(ggplot2)


models.hqbiome1    <- readRDS("input/HQbiome1.RDS")
models.hqbiome2    <- readRDS("input/HQbiome2.RDS")

medium.2216 <- fread("input/2216e.csv")
medium.2216[, ex.rxn := paste0("EX_", compounds, "_e0")]

get_sim <- function(models, tsteps=7){
  arena <- Arena(n=40, m=40)
  for(i in seq_along(models)){
    mod <- models[[i]]
    if( "EX_cpd00221_e0" %in% mod@react_id ) mod <- rmReact(mod, react="EX_cpd00221_e0") # d-lactate
    bac   <- Bac(mod, setAllExInf=TRUE) 
    arena <- addOrg(arena, bac, amount = 5)
  }
  
  arena <- addSubs(arena,  smax = medium.2216$maxFlux/100, mediac = medium.2216$ex.rxn,unit = "mM")
  sim <- simEnv(arena, time=tsteps, sec_obj='mtf')
  
  return(sim)
}

sim.hqbiome1 <- get_sim(models.hqbiome1,  tsteps=7)
sim.hqbiome2 <- get_sim(models.hqbiome2,  tsteps=7)


saveRDS(sim.hqbiome1, "input/sim.HQbiome1-40_40.RDS", compress = "xz")
saveRDS(sim.hqbiome2, "input/sim.HQbiome2-40_40.RDS", compress = "xz")


######HQbiome-1  
sim.hqbiome1 <- readRDS("input/sim.HQbiome1-40_40.RDS")
sim.hqbiome2 <- readRDS("input/sim.HQbiome2-40_40.RDS")

#top30 subs
dt.sim.hqbiome1<-sim.hqbiome1@exchangeslist[[4]]
plotGrowthCurve(sim.hqbiome1)[[2]]
VarSubs.hqbiome1<-getVarSubs(sim.hqbiome1)
subs <- names(head(getVarSubs(sim.hqbiome1),30))
#Fig. S24
plotSubCurve(sim.hqbiome1, mediac=subs,useNames=T)[[1]]
plotSubCurve(sim.hqbiome2, mediac=subs,useNames=T)[[1]]
plotSpecActivity(sim.hqbiome1, subs=subs)[[4]]


seed <- fread("input/seed_metabolites_edited.tsv")

subs <- names(head(getVarSubs(sim.hqbiome1),30))
subs
subs1<-str_replace(subs,'EX_', '')
subs2<-str_replace(subs1,'_e0', '')
subs3<-str_replace(subs2,'_c0', '')
sub.ferm.seed <-subs3
sub.ferm.seed

###plot
library(data.table)
library(ggplot2)
library(dplyr)
my_breaks <- c(0.01,0.1,1,10,100,1000)

# HQbiome-1
a.gs <- ex.hebiomq1
a.gs$recon <- "HQbiome-1"
a.gs <- a.gs[time == 4]
a.gs[, spec.name := gsub("^ ","", spec.name)]
a.gs[, sub := gsub("-e0$","",sub)]
a.gs <- a.gs[abs(mflux)>0]
#a[stat == "uptake", mflux := mflux*100]
a.gs[, mflux2 := abs(mflux)]


# plotting ( High: )
sub.incl <- "CO2|O2|Acetate|H2O|Succinate|L-Glutamate|L-Proline|L-Threonine|Formate|L-Alanine|Citrate|L-Isoleucine|L-Cysteine|H2S|L-Aspartate|Propionate|L-Serine|Glycine|L-Arginine|L-Valine|L-Leucine|Urea
Phosphate|H2|L-Lysine"
p1 <- ggplot(a.gs[grepl(sub.incl, sub)], aes(sub, spec)) + 
  #geom_tile(data = tmp[!grepl(sub.excl, V1)], aes(V1, V2, fill = N)) +
  geom_tile(aes(fill=mflux2), colour = "white", size = 1.5) + 
  scale_fill_gradient(low = "white", high = "#FF5500FF", name = "Production\n [unit]", trans = "log",
                      breaks = my_breaks, labels = my_breaks) +
  theme_bw() +
  facet_grid(stat~recon) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1
ggsave(filename = "plots/HQbiome-1_red.pdf", plot = p1, width = 6.25, height = 8)

p2 <- ggplot(a.gs[grepl(sub.incl, sub)], aes(sub, spec)) + 
  #geom_tile(data = tmp[!grepl(sub.excl, V1)], aes(V1, V2, fill = N)) +
  geom_tile(aes(fill=mflux2), colour = "white", size = 1.5) + 
  scale_fill_gradient(low = "white", high = "#008B00FF", name = "Production\n [unit]", trans = "log",
                      breaks = my_breaks, labels = my_breaks) +
  theme_bw() +
  facet_grid(stat~recon) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p2
ggsave(filename = "plots/HQbiome-1_green.pdf", plot = p2, width = 6.25, height = 8)


#####################################################################################
########HQbiome-2
sim.hqbiome2 <- readRDS("dat/sim.HQbiome2-40_40.RDS")

#top30 subs
dt.sim.hqbiome2<-sim.hqbiome2@exchangeslist[[4]]
plotGrowthCurve(sim.hqbiome2)[[2]]
VarSubs.hqbiome2<-getVarSubs(sim.hqbiome2)
subs <- names(head(getVarSubs(sim.hqbiome2),30))

plotSubCurve(sim.hqbiome2, mediac=subs,useNames=T)[[1]]

plotSpecActivity(sim.hqbiome1, subs=subs)[[4]]

seed <- fread("dat/seed_metabolites_edited.tsv")

seedEX2bigg <- function(l){
  is.ex <- str_starts(l[[1]], "EX_")
  dict <- fread("dat/SEED2VMH_translation_edited.csv", header = F)
  if( is.ex ){
    idx  <- match(gsub("_e0$","",l), gsub("\\(e\\)$","",dict$V1))  
    bigg <- gsub("_((L|D|R)\\(e\\))","__\\1", dict$V2[idx])
  }else{
    idx  <- match(l, str_extract(dict$V1, "cpd[0-9]+"))
    bigg <- str_extract(dict$V2[idx], "(?<=EX_).*?(?=\\(e\\))")
    bigg <- gsub("_((L|D|R)$)","__\\1", bigg)
  }
  return(bigg)
}


subs <- names(head(getVarSubs(sim.hqbiome2),30))
subs
subs1<-str_replace(subs,'EX_', '')
subs2<-str_replace(subs1,'_e0', '')
subs3<-str_replace(subs2,'_c0', '')
sub.ferm.seed <-subs3
sub.ferm.seed

sub.ferm.bigg <- seedEX2bigg(sub.ferm.seed)

sim.activity <- function(sim, namespace.seed, sub.ferm#, agora=F
)                         {
  if( namespace.seed ){
    dat.sub <- unique(c(setdiff(sub.ferm, "cpd00011"), "cpd00239"))
    dat <- data.table(BacArena::plotSpecActivity(sim, subs=paste0("EX_",dat.sub,"_e0"),useNames = T, rm_unused = T, ret_data = T))
    dat[,id:=seed$id[match(tolower(gsub("-e0$","",sub)), tolower(seed$name))]]
    dat[sub=="raffinose", id:="cpd00382"]
    #   dat[,spec.name:=trimws(genome.desc$name[match(gsub(".fna.sbml","", spec), gsub(".fna.gz","",genome.desc$file))])]
  }else{
    dat.sub <- unique(c(setdiff(sub.ferm, seedEX2bigg("cpd00011")), seedEX2bigg("cpd00239")))
    dat <- data.table(plotSpecActivity(sim, subs=paste0("EX_",dat.sub,"(e)"),useNames = F, rm_unused = T, ret_data = T))
    dat[,id:=str_extract(sub, "(?<=EX_).*?(?=\\(e\\))")]
    
  }
  dat[,stat:=ifelse(mflux>0, "production", ifelse(mflux<0, "uptake", NA))]
  dat[spec=="barkeri_iAF692", spec.name:="Methanosarcina barkeri"]
}

# plot growth curve
growth.dt <- data.table()
growth.dt <- rbind(growth.dt, data.table(BacArena::plotGrowthCurve(sim.hqbiome1, use_biomass = F, ret_data = T))[, method:="gapseq"])
p <- ggplot(growth.dt, aes(x=time, y=value)) + 
  stat_summary(fun = "sum", geom="line") + 
  facet_wrap(~method) + 
  ylab("Total number of organisms") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p

ex.hebiomq2    <- sim.activity(sim.hqbiome2, namespace.seed = T, sub.ferm=sub.ferm.seed)

###plot
library(data.table)
library(ggplot2)
library(dplyr)
my_breaks <- c(0.01,0.1,1,10,100,1000)

# HQbiome-2
a.gs <- ex.hebiomq2
a.gs$recon <- "HQbiome-2"
a.gs <- a.gs[time == 4]

a.gs[, spec.name := gsub("^ ","", spec.name)]
a.gs[, sub := gsub("-e0$","",sub)]

a.gs <- a.gs[abs(mflux)>0]
#a[stat == "uptake", mflux := mflux*100]
a.gs[, mflux2 := abs(mflux)]


# plotting ( High: )
sub.incl <- "CO2|O2|Acetate|Succinate|Citrate|Formate|L-Glutamate
|L-Threonine|L-Alanine|L-Cysteine|L-Serine|H2S|L-Proline|Nitrate|Nitrite|L-Arginine
|H2O|L-Isoleucine|Urea|L-Aspartate|Propionate|L-Histidine|L-Valine
|Phosphate|Glycine|L-Leucine"
p1 <- ggplot(a.gs[grepl(sub.incl, sub)], aes(sub, spec)) + 
  #geom_tile(data = tmp[!grepl(sub.excl, V1)], aes(V1, V2, fill = N)) +
  geom_tile(aes(fill=mflux2), colour = "white", size = 1.5) + 
  scale_fill_gradient(low = "white", high = "#FF5500FF", name = "Production\n [unit]", trans = "log",
                      breaks = my_breaks, labels = my_breaks) +
  theme_bw() +
  facet_grid(stat~recon) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1
ggsave(filename = "plots/HQbiome-2_red.pdf", plot = p1, width = 6.25, height = 8)

p2 <- ggplot(a.gs[grepl(sub.incl, sub)], aes(sub, spec)) + 
  #geom_tile(data = tmp[!grepl(sub.excl, V1)], aes(V1, V2, fill = N)) +
  geom_tile(aes(fill=mflux2), colour = "white", size = 1.5) + 
  scale_fill_gradient(low = "white", high = "#008B00FF", name = "Production\n [unit]", trans = "log",
                      breaks = my_breaks, labels = my_breaks) +
  theme_bw() +
  facet_grid(stat~recon) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p2
ggsave(filename = "plots/HQbiome-2_green.pdf", plot = p2, width = 6.25, height = 8)

###################################################################################################
#Fig. S31. The effects of bacteria with different adaptive strategies on nematodes development
library(readxl)
library(ggtern)

uast <- read_excel("input/UAST_GROUP.xlsx", 
                   sheet = "UAST b=150!")
uast$Strategy<-factor(uast$Strategy,levels = c('Stress tolerator','Competitor','Ruderal','Mix'))
uast$`Phylum/Class`<-factor(uast$`Phylum/Class`,levels = c('Actinobacteria','Bacteroidetes','Firmicutes',
                                                           'Alphaproteobacteria','Gammaproteobacteria'))

#Fig. S31 bc

library(rstatix)
library(ggpubr)
main_theme <- theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    axis.line=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=10),
                    #legend.position="top",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    #text=element_text(family="sans")
)
cols <- c(Actinobacteria= '#EEE685',Bacteroidetes = "#4F94CD", Firmicutes= '#EE7600', 
          Alphaproteobacteria = "#00CD00",Gammaproteobacteria = "#90EE90")

####L. marina
b<-ggplot(uast, aes(x=Strategy, y=mean_5)) +
  geom_jitter( size=1.7, alpha=0.9,width = 0.1) +
  geom_boxplot(alpha=0.8, #outlier.size=.2, 
               #size=1,#width = 0.25,
               width=0.4,cex=0.8,
               #position=dodge,
               color=c('#E31A1C','#228B22','#1F78B4','#8B658B'), 
               fill=NA
  ) +
  #geom_violin(aes(group = Strategy),alpha = 0.1,
  #            fill = NA,
  #            #colour = NA,
  #            #position=dodge,
  #            width=0.5, cex=0.8)+
  labs(x="", y="Perportion of L4 in day 5 (%)") +
  scale_colour_manual(values=cols, guide='none') +
  coord_flip() +
  #scale_y_continuous(limits=c(-10, 70)) +
  main_theme +
  theme(axis.text.y = element_text(size=8))
b
#ggsave('p2.pdf',P2,width = 5,height = 4)
#test
kruskal_test(uast, strategy ~ mean_5)  #0.552
uast %>% wilcox_test (mean_5 ~ Strategy)


###C. elegans
c<-ggplot(uast, aes(x=Strategy, y=N2)) +
  geom_jitter( size=1.7, alpha=0.9,width = 0.1) +
  geom_boxplot(alpha=0.8, #outlier.size=.2, 
               #size=1,#width = 0.25,
               width=0.4,cex=0.8,
               #position=dodge,
               color=c('#E31A1C','#228B22','#1F78B4','#8B658B'), 
               fill=NA
  ) +
  #geom_violin(aes(group = Strategy),alpha = 0.1,
  #            fill = NA,
  #            #colour = NA,
  #            #position=dodge,
  #            width=0.5, cex=0.8)+
  labs(x="", y="Egg laying time of C. elegans (%)") +
  scale_colour_manual(values=cols, guide='none') +
  coord_flip() +
  #scale_y_continuous(limits=c(-10, 70)) +
  main_theme +
  theme(axis.text.y = element_text(size=8))
c
#ggsave('p2.pdf',P2,width = 5,height = 4)
#test
kruskal_test(uast, strategy ~ N2)   #0.7
uast %>% wilcox_test (N2 ~ Strategy)

uast2<-uast %>% filter( Strategy == c("Competitor",'Ruderal'))

wilcox.test(N2 ~ Strategy, data = uast2)

bc<- ggarrange(b,c,ncol=2,widths = c(1,1),common.legend =F,
               labels = c('a','b'))
bc
ggsave('bc.pdf',bc,width = 8,height = 4)


#cor test
#f34
cor_test(uast,vars = mean_5,vars2 =c(S,C,R),method='spearman') #NS
#var1   var2     cor statistic     p method  
#<chr>  <chr>  <dbl>     <dbl> <dbl> <chr>   
#  1 mean_5 S     -0.054    71189. 0.646 Spearman
#2 mean_5 C      0.086    61697. 0.465 Spearman
#3 mean_5 R     -0.15     77883. 0.192 Spearman

#n2
cor_test(uast,vars = N2,vars2 =c(S,C,R),method='spearman')
#var1  var2     cor statistic     p method  
#<chr> <chr>  <dbl>     <dbl> <dbl> <chr>   
#  1 N2    S     -0.048    70752. 0.686 Spearman
#2 N2    C     -0.082    73044. 0.489 Spearman
#3 N2    R      0.13     58977. 0.282 Spearman

#sig test
library(dplyr)
library(rstatix)

#································································································································
shapiro.test(uast$mean_5)  #W = 0.82593, p-value = 6.738e-08    
bartlett.test(mean_5 ~ Strategy, data = uast) #Bartlett's K-squared = 0.85633, df = 3, p-value = 0.836 
# p<0.05为方差不齐，大于为齐#·············································· ·················································································
uast %>% kruskal_test(mean_5 ~ Strategy)  #0.97

#F23
uast %>% 
  wilcox_test(mean_5 ~ Strategy, p.adjust.method = "fdr") %>% 
  add_xy_position(x = "Strategy")  
#N2
uast %>% 
  wilcox_test(mean_5 ~ Strategy, p.adjust.method = "fdr") %>% 
  add_xy_position(x = "Strategy")


###############################################################################################################
###Fig. S31a group--strategies 
group<-read_excel('input/UAST_GROUP.xlsx',sheet='Group')
uast2<- uast %>% merge(group, by.x='org',by.y = 'ID', sort = F) %>% select(#'org','strategy',
  'Strategy','Group')#%>% 
#filter(Strategy != c("Mix"))

# chisq.test
uast_1_2 <- uast2 %>% filter(Group != c("Intermediate"))
chisq_test(table(uast_1_2))  #0.737     

uast_1_3 <- uast2 %>% filter(Group != c("HQbiome-2"))
chisq_test(table(uast_1_3))  #0.91     

uast_2_3 <- uast2 %>% filter(Group != c("HQbiome-1"))
chisq_test(table(uast_2_3))  #0.99      

#plot
library(reshape2)
color<-c("#E31A1C", "#228B22", "#1F78B4", "#8B658B")
plot.df<-read_xlsx('input/uast_group2.xlsx',sheet=2) %>% melt
plot.df$value<-plot.df$value*100
plot.df$Group<-factor(plot.df$Group,levels=c('Stress tolerator','Competitor','Ruderal','Mix'))
plot.df$variable<-factor(plot.df$variable,levels=c('HQbiome-1','Intermediate','HQbiome-2'))
a<-ggplot(plot.df, aes(x=variable, y=value*100, fill=Group))+
  labs(y="Proportion (%)",x='')+
  main_theme+
  geom_bar(position="fill",
           stat="identity",
  )+
  scale_fill_manual(values=color)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
a
ggsave('stack plot.pdf',P6,width = 4,height = 3)





