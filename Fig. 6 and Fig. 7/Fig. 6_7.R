library(tidyverse) # Data processing
library(ggplot2) # Plot figures
library(qiime2R) # QIIME2 artifacts to phyloseq object 
library(vegan) # Ecology analysis
library(phyloseq) # Base microbiome data structure
library(microbiome) # Microbiome data analysis and visualization
library(phylosmith) # Microbiome data analysis and visualization
library(microbiomeutilities) # Microbiome data analysis and visualization
library(ggordiplots) # Vegan ordination plots for ggplot2
library(lubridate) # Date formatting
library(readxl) # Read in excel .xlsx
library(ggrepel) # Prevent overlapping text in figures
library(eulerr) # Venn diagram visualization for core microbiome
library(treemap) # Tree map visualization
library(ggpubr) # Data visualization - wrapper
library(rstatix) # Tidy statistical tests
library(RColorBrewer) # Color selection for graphics
library(tidytext) # Text repel on graphics 
library(dendextend) # Dendrogram plotting
library(ggsci) # color
library(paletteer)
#HQbiome analisis
physeq_gut <- qza_to_phyloseq(
  features = "input/HQbiome_qza/table_deblur.qza", # ASV/OTU table
  #tree = "tree_rooted.qza", # Phylogenetic tree
  taxonomy = "input/HQbiome_qza/tax.qza", # Taxonomy file
  metadata = "input/HQbiome_qza/meta.txt" # Sample metadata
)
physeq_gut
#remove back
physeq_gut1<-subset_samples(physeq_gut, !(group %in% c('BACK', 'Inoculum','lawn')))
physeq_gut2<-subset_samples(physeq_gut1, !(name %in% c('ND1-3','FD1-1'
                                                       ,'ND7-1','ND7-3'
)))
#add ND7-1
physeq_N7.1 <- qza_to_phyloseq(
  features = "input/N7.1_qza/table_deblur.qza", # ASV/OTU table
  #tree = "tree_rooted.qza", # Phylogenetic tree
  taxonomy = "input/N7.1_qza/tax.qza", # Taxonomy file
  metadata = "input/N7.1_qza/meta.txt" # Sample metadata
)
sample_names(physeq_N7.1)<-'ND7-1'    #correct name to ND7-1 
sample_data(physeq_N7.1)$total_reads <- sample_sums(physeq_N7.1)
#MERGE
physeq_gut_all<-merge_phyloseq(physeq_gut2, physeq_N7.1)
physeq_gut_all
#lawn
physeq_lawn <- qza_to_phyloseq(
  features = "input/qza-lawn/table_deblur.qza", # ASV/OTU table
  taxonomy = "input/qza-lawn/tax.qza", # Taxonomy file
  metadata = "input/qza-lawn/meta.txt" # Sample metadata
)
physeq_lawn
physeq_lawn1<-subset_samples(physeq_lawn, !(name %in% c('F23-lawn1.3','F23-lawn3.3','F23-lawn5.2','F23-lawn7.3',
                                                        'N2-lawn1.3','N2-lawn3.3','N2-lawn5.3','N2-lawn7.3',
                                                        'F23-lawn7.1','F23-lawn7.2')))
#merge phyloseq
physeq_all<-merge_phyloseq(physeq_gut_all, physeq_lawn1)
physeq_all

#merge genus
physeq_all <- physeq_all %>%
  tax_glom(taxrank = "Genus", NArm = TRUE) # agglomerate on Genus level
tax_table(physeq_all) <- tax_table(physeq_all)[,1:6]

#order
physeq_all@sam_data$group<-factor(physeq_all@sam_data$group,levels = c('Gut-F23','Gut-N2','lawn-F23','lawn-N2','Inoculum'))

physeq_all@sam_data$name <- factor(physeq_all@sam_data$name, 
                                   levels = c("FD1-2","FD1-3","FD3-2","FD3-1","FD5-1","FD5-2",
                                              "ND1-1","ND1-2","ND3-1","ND3-2","ND5-1",'ND5-2','ND7-2','ND7-1',
                                              "F23-lawn1.1","F23-lawn1.2","F23-lawn3.1","F23-lawn3.2",
                                              "F23-lawn5.1","F23-lawn5.3","N2-lawn1.1","N2-lawn1.2",
                                              "N2-lawn3.1","N2-lawn3.2","N2-lawn5.1","N2-lawn5.2",
                                              "N2-lawn7.1","N2-lawn7.2","Inoculum1","Inoculum2","Inoculum3"
                                   ))

physeq_all@sam_data$name

#merge samples
phyloseq_merge<-merge_samples(physeq_all,'merge')

#Draw graphs of F23 or N2 separately 
#Fig. 6a
#F23
F23_phyloseq<-subset_samples(physeq_all,physeq_all@sam_data$group%in%c('Gut-F23','lawn-F23'))
#merge samples
F23_phyloseq_merge<-merge_samples(F23_phyloseq,'merge')
#F23  plot
F23_physeq_other_Genus <- F23_phyloseq_merge %>%
  tax_glom(taxrank = "Genus") %>%    # agglomerate on Class level
  transform_sample_counts(function(x) {100 * x/sum(x)}) %>%  # transform to relative abundance
  psmelt() %>%  # convert from phyloseq object to a long data frame
  mutate(Taxa_Order = Genus)# %>%
#mutate(Taxa_Order = replace(Taxa_Order, Abundance < 1, "<1%"))

F23_Genus_BP <- 
  ggplot(data=F23_physeq_other_Genus, 
         aes(x=Sample, # order by date for each library season
             y=Abundance, 
             fill=fct_reorder(Taxa_Order, Abundance, .desc = TRUE))) + # order taxa within bars by abundance
  geom_bar(stat="identity", 
           position = position_stack(reverse = TRUE), # order taxa in reverse
           width = 0.9) + 
  #geom_col(width = 1)+
  scale_fill_manual(values = c('Shewanella'='#3A82E4FF','Vibrio'='#60BD68FF','Bacillus'='#FF9800FF',
                               'Lactococcus'='#E9A820FF','Pseudoalteromonas'='#D25D38FF','Pseudomonas'='#C6242DFF',
                               'Aliivibrio'='#244579FF','Galactobacter'='#3F51B5FF','Marinomonas'='#8A6842FF',
                               'Neptunomonas'='#66E3D9FF','<1%'='#999999FF'
                               
  )
  #paletteer_d("ggsci::default_igv",50,direction=-1)
  , # assign legend colors
  guide = guide_legend(reverse = TRUE)) + # make legend match order of taxa in barplot
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), # Y axis by 10% increments
                     expand = c(0,0)) + # remove whitespace in plot
  facet_grid(~group, # draw facets by date
             scales = "free_x", 
             space = "free_x",
             drop = F
  ) + #switch="x" 
  theme(legend.direction="vertical", 
        axis.text.x = element_text(color = 'black',angle = 45, hjust = 1),
        #element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y = element_text(#face = "bold", 
          size = 10,color = 'black'),
        strip.text = element_text(face = "bold",#, size = 15
        ),
        strip.background = element_rect(color = "white", fill = "#DDEC7EFF"
        ),
        legend.title = element_text( face = "bold"),
        plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
        panel.spacing = unit(0, "lines"),
        legend.text = element_text(face="italic"),
        panel.background = element_rect(fill = 'white', color = 'white'),  
        panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
  labs(  face = "bold",
         y = "Relative Abundance %", 
         fill = "Genus")+
  guides(fill = guide_legend(ncol = 1))
F23_Genus_BP

#Fig. 7a
#N2
N2_phyloseq<-subset_samples(physeq_all,physeq_all@sam_data$group%in%c('Gut-N2','lawn-N2'))
#merge samples
N2_phyloseq_merge<-merge_samples(N2_phyloseq,'merge')
#N2  plot
N2_physeq_other_Genus <- N2_phyloseq_merge %>%
  tax_glom(taxrank = "Genus") %>%    # agglomerate on Class level
  transform_sample_counts(function(x) {100 * x/sum(x)}) %>%  # transform to relative abundance
  psmelt() %>%  # convert from phyloseq object to a long data frame
  mutate(Taxa_Order = Genus)# %>%
#mutate(Taxa_Order = replace(Taxa_Order, Abundance < 1, "<1%"))

N2_Genus_BP <- 
  ggplot(data=N2_physeq_other_Genus, 
         aes(x=Sample, # order by date for each library season
             y=Abundance, 
             fill=fct_reorder(Taxa_Order, Abundance, .desc = TRUE))) + # order taxa within bars by abundance
  geom_bar(stat="identity", 
           position = position_stack(reverse = TRUE), # order taxa in reverse
           width = 0.9) + 
  #geom_col(width = 1)+
  scale_fill_manual(values = c('Shewanella'='#3A82E4FF','Pseudomonas'='#C6242DFF','Lactococcus'='#E9A820FF',
                               'Psychrobacter'='#8D4B08FF','Neptunomonas'='#66E3D9FF',
                               'Parasphingorhabdus'='#F17CB0FF','Paenisporosarcina'='#00FF7FFF',
                               'Paracoccus'='#006064FF','Vibrio'='#60BD68FF','Postechiella'='#F15854FF',
                               'Pseudoclavibacter'='#56B4E9FF','Alkalihalobacillus'='#F0E442FF',
                               'Sulfitobacter'='#009E73FF','Metabacillus'='#CC79A7FF',
                               'Bacillus'='#FF9800FF',
                               '<1%'='#999999FF')
                    , # assign legend colors
                    guide = guide_legend(reverse = TRUE)) + # make legend match order of taxa in barplot
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), # Y axis by 10% increments
                     expand = c(0,0)) + # remove whitespace in plot
  facet_grid(~group, # draw facets by date
             scales = "free_x", 
             space = "free_x",
             drop = F
  ) + #switch="x" 
  theme(legend.direction="vertical", 
        axis.text.x = element_text(color = 'black',angle = 45, hjust = 1),
        #element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y = element_text(#face = "bold", 
          size = 10,color = 'black'),
        strip.text = element_text(face = "bold",#, size = 15
        ),
        strip.background = element_rect(color = "white", fill = "#DDEC7EFF"
        ),
        legend.title = element_text( face = "bold"),
        plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
        panel.spacing = unit(0, "lines"),
        legend.text = element_text(face="italic"),
        panel.background = element_rect(fill = 'white', color = 'white'), 
        panel.border = element_rect(colour = "black", fill=NA, size=1.5)
  ) +
  labs(  face = "bold",
         y = "Relative Abundance %", 
         fill = "Genus")+
  guides(fill = guide_legend(ncol = 1))+
  theme()
N2_Genus_BP

#Fig. S20a
#Inoculum
ino_phyloseq<-subset_samples(physeq_all,physeq_all@sam_data$group%in%c('Inoculum'))
#merge samples
N2_phyloseq_merge<-merge_samples(N2_phyloseq,'merge')

ino_physeq_other_Genus <- ino_phyloseq %>%
  tax_glom(taxrank = "Genus") %>%    # agglomerate on Class level
  transform_sample_counts(function(x) {100 * x/sum(x)}) %>%  # transform to relative abundance
  psmelt() %>%  # convert from phyloseq object to a long data frame
  mutate(Taxa_Order = Genus)# %>%
#mutate(Taxa_Order = replace(Taxa_Order, Abundance < 1, "<1%"))

ino_Genus_BP <- 
  ggplot(data=ino_physeq_other_Genus, 
         aes(x=Sample, # order by date for each library season
             y=Abundance, 
             fill=fct_reorder(Taxa_Order, Abundance, .desc = TRUE))) + # order taxa within bars by abundance
  geom_bar(stat="identity", 
           position = position_stack(reverse = TRUE), # order taxa in reverse
           width = 0.9) + 
  #geom_col(width = 1)+
  scale_fill_manual(values = c('Vibrio'='#60BD68FF', 'Shewanella'='#3A82E4FF','Psychrobacter'='#8D4B08FF',
                               'Rheinheimera'='#004949FF','Parasphingorhabdus'='#F17CB0FF','Aliivibrio'='#244579FF',
                               'Pseudomonas'='#C6242DFF', 'Enterovibrio'='#DECF3FFF','Lacinutrix'='#F17CB0FF',
                               'Sulfitobacter'='#009E73FF','Pseudoalteromonas'='#D25D38FF','Marinobacter'='#66E3D9FF',
                               'Thalassotalea'='#00FF7FFF','HQ336491_g'='#FCA3B7FF','Neptunomonas'='#66E3D9FF',
                               'Litoreibacter'='#98FB98FF','Saccharospirillum'='#D070B9FF','Yoonia'='#66CDAAFF',
                               'Bacillus'='#FF9800FF','Alteromonas'='#F49538FF','Marinagarivorans'='#7D5329FF',
                               'Postechiella'='#F15854FF','Marinomonas'='#8A6842FF','Salinimonas'='#E7D202FF',
                               'Staphylococcus'='#BF616AFF','Alkalihalobacillus'='#F0E442FF','Lactococcus'='#E9A820FF',
                               'Priestia'='#8FA87AFF','Thaumasiovibrio'='#5E81ACFF','Cytobacillus'='#FFFF6DFF',
                               'Galactobacter'='#3F51B5FF','Planococcus'='#24FF24FF','Paracoccus'='#006064FF',
                               'Halobacillus'='#24FF24FF','Metabacillus'='#CC79A7FF','Niallia'='#DB6D00FF',
                               'Rhodococcus'='#924900FF','Paenisporosarcina'='#00FF7FFF','Colwellia'='#920000FF',
                               'Microbacterium'='#B6DBFFFF','Loktanella'='#6DB6FFFF','Arenivirga'='#B66DFFFF'
  )
  #paletteer_d("ggsci::default_igv",50,direction=-1)
  , # assign legend colors
  guide = guide_legend(reverse = TRUE)) + # make legend match order of taxa in barplot
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), # Y axis by 10% increments
                     expand = c(0,0)) + # remove whitespace in plot
  facet_grid(~group, # draw facets by date
             scales = "free_x", 
             space = "free_x",
             drop = F
  ) + #switch="x" 
  theme(legend.direction="vertical", 
        axis.text.x = element_text(color = 'black',angle = 45, hjust = 1),
        #element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y = element_text(#face = "bold", 
          size = 10,color = 'black'),
        strip.text = element_text(face = "bold",#, size = 15
        ),
        strip.background = element_rect(color = "white", fill = "#DDEC7EFF"
        ),
        legend.title = element_text( face = "bold"),
        plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
        panel.spacing = unit(0, "lines"),
        legend.text = element_text(face="italic"),
        panel.background = element_rect(fill = 'white', color = 'white'), 
        panel.border = element_rect(colour = "black", fill=NA, size=1.5)
  ) +
  labs(  face = "bold",
         y = "Relative Abundance %", 
         fill = "Genus")+
  guides(fill = guide_legend(ncol = 3))
ino_Genus_BP

#Fig. S22a
#PcoA
library(MicrobiotaProcess)
library(patchwork)

#genus level
genus_db<-read_excel('input/tax-HQbiome.xlsx',sheet = 'pcoa1')
genus_db1<-genus_db[,-1]
row.names(genus_db1)<-genus_db$otu
distMatrix <- vegdist(genus_db1,method = "bray")
pCoa <- cmdscale(distMatrix, eig = T,k = 2 )
varExp <- (eigenvals(pCoa)/sum(eigenvals(pCoa)))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)

pCoaVecs <- as.data.frame(pCoa$points)
colnames(pCoaVecs) <- paste0("PCo",c(1:2))
pCoaVecs$ID <- row.names(pCoaVecs)

meta<-read_excel('input/tax-HQbiome.xlsx',sheet = 'meta')
pCoaVecs2<-merge(pCoaVecs,meta,by='ID')

#PERMANOVA adonis2
adonis_result <- adonis2(genus_db1 ~ group, data=meta, permutations=999)

library(paletteer)
pcoa.p<-ggplot(pCoaVecs2,aes(x=PCo1,y=PCo2,color=group)) + 
  geom_point(size=2) + 
  theme_classic() + 
  #scale_color_gradientn(colours = rainbow(5)) + 
  xlab(paste0('PCoA1 (',round(xVar,2),' %)')) + 
  ylab(paste0('PCoA2 (',round(yVar,2),' %)')) + 
  theme(axis.text = element_text(color = 'black'))+
  scale_colour_manual(name = "Group", values = paletteer_d("basetheme::void"))+
  stat_ellipse(data = pCoaVecs2, aes(x = PCo1, y = PCo2),
               level = 0.68, type='norm',size=0.6)+
  annotate("text",x=-0.2,y=-0.5,label=paste("p= ", adonis_result$`Pr(>F)`[1]),size=3)
#ggforce::geom_mark_ellipse(aes(group=group),expand = unit(1, "mm"),tol =0.05)
#save size =
pcoa.p

##############################################################################################################33
#####################################################################################################
#differential abundance analysis Fig. S21a-d
library(DESeq2)
library(xlsx)
physeq_all
#F23
F23_phyloseq<-subset_samples(physeq_all,physeq_all@sam_data$group%in%c('Gut-F23','lawn-F23'))
#merge samples
#F23_phyloseq_merge<-merge_samples(F23_phyloseq,'merge')

##======deseq2   gut and lawn
#day 1
F23_phyloseq_day1<-subset_samples(F23_phyloseq,F23_phyloseq@sam_data$merge%in%c('F23-Day1','F23-lawn-Day1'))
head(sample_data(F23_phyloseq_day1)$group, n=10)
diagdds = phyloseq_to_deseq2(F23_phyloseq_day1, ~ group)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #lawn.F23 vs Gut.F23
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(F23_phyloseq_day1)[rownames(sigtab), ], "matrix"))
head(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}


# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
d1.lawn.vs.gut<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
d1.lawn.vs.gut
ggsave('differential abundance/lawn.vs.gut_F23_day1.pdf',d1.lawn.vs.gut,width = 5,height = 3)
write.xlsx(sigtab2,'differential abundance/sigtab_lawn.vs.gut_F23_day1.xlsx')

#day 3
F23_phyloseq_day3<-subset_samples(F23_phyloseq,F23_phyloseq@sam_data$merge%in%c('F23-Day3','F23-lawn-Day3'))
head(sample_data(F23_phyloseq_day3)$group, n=10)
diagdds = phyloseq_to_deseq2(F23_phyloseq_day3, ~ group)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #lawn.F23 vs Gut.F23
sigtab
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(F23_phyloseq_day3)[rownames(sigtab), ], "matrix"))
head(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}


# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
d3.lawn.vs.gut<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
d3.lawn.vs.gut
ggsave('differential abundance/lawn.vs.gut_F23_day3.pdf',d3.lawn.vs.gut,width = 5,height = 3)
write.xlsx(sigtab2,'differential abundance/sigtab_lawn.vs.gut_F23_day3.xlsx')

#day 5
F23_phyloseq_day5<-subset_samples(F23_phyloseq,F23_phyloseq@sam_data$merge%in%c('F23-Day5','F23-lawn-Day5'))
head(sample_data(F23_phyloseq_day5)$group, n=10)
diagdds = phyloseq_to_deseq2(F23_phyloseq_day5, ~ group)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #lawn.F23 vs Gut.F23
sigtab
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(F23_phyloseq_day5)[rownames(sigtab), ], "matrix"))
head(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}


# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
d5.lawn.vs.gut<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
d5.lawn.vs.gut
ggsave('differential abundance/lawn.vs.gut_F23_day5.pdf',d5.lawn.vs.gut,width = 5,height = 3)
write.xlsx(sigtab2,'differential abundance/sigtab_lawn.vs.gut_F23_day5.xlsx')


###gut  day1 vs day3
F23_gut_phyloseq<-subset_samples(F23_phyloseq,F23_phyloseq@sam_data$group%in%c('Gut-F23'))
F23_gut_phyloseq_13<-subset_samples(F23_gut_phyloseq,F23_gut_phyloseq@sam_data$merge %in%c('F23-Day1','F23-Day3'))

head(sample_data(F23_gut_phyloseq_13)$merge, n=15)

diagdds = phyloseq_to_deseq2(F23_gut_phyloseq_13, ~ merge)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #F23.Day3 vs F23.Day1 
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(F23_gut_phyloseq_13)[rownames(sigtab), ], "matrix"))
head(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}


# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
d3.vs.d1<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
d3.vs.d1
ggsave('differential abundance/d3.vs.d1t.pdf',d3.vs.d1,width = 3,height = 2.5)
write.xlsx(sigtab2,'differential abundance/sigtab_d3.vs.d1.xlsx')

###gut  day1 vs day5
F23_gut_phyloseq<-subset_samples(F23_phyloseq,F23_phyloseq@sam_data$group%in%c('Gut-F23'))
F23_gut_phyloseq_15<-subset_samples(F23_gut_phyloseq,F23_gut_phyloseq@sam_data$merge %in%c('F23-Day1','F23-Day5'))

head(sample_data(F23_gut_phyloseq_15)$merge, n=15)

diagdds = phyloseq_to_deseq2(F23_gut_phyloseq_15, ~ merge)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]   #F23.Day5 vs F23.Day1 
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(F23_gut_phyloseq_15)[rownames(sigtab), ], "matrix"))
head(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}


# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
d5.vs.d1<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
d5.vs.d1
ggsave('differential abundance/d5.vs.d1.pdf',d5.vs.d1,width = 3,height = 2.5)
write.xlsx(sigtab2,'differential abundance/sigtab_d5.vs.d1.xlsx')


###gut  day3 vs day5  No difference
F23_gut_phyloseq<-subset_samples(F23_phyloseq,F23_phyloseq@sam_data$group%in%c('Gut-F23'))
F23_gut_phyloseq_35<-subset_samples(F23_gut_phyloseq,F23_gut_phyloseq@sam_data$merge %in%c('F23-Day3','F23-Day5'))

head(sample_data(F23_gut_phyloseq_35)$merge, n=15)

diagdds = phyloseq_to_deseq2(F23_gut_phyloseq_35, ~ merge)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]   #F23.Day5 vs F23.Day3  
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(F23_gut_phyloseq_35)[rownames(sigtab), ], "matrix"))
head(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}


# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
d5.vs.d3<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
d5.vs.d3
ggsave('differential abundance/d5.vs.d3.pdf',d5.vs.d3,width = 3,height = 2.5)
write.xlsx(sigtab2,'differential abundance/sigtab_d5.vs.d3.xlsx')


####################################################################################
#deseq2   N2

N2_phyloseq<-subset_samples(physeq_all,physeq_all@sam_data$group%in%c('Gut-N2','lawn-N2'))
#merge samples
#F23_phyloseq_merge<-merge_samples(F23_phyloseq,'merge')

##======deseq2   gut and lawn
#day 1
N2_phyloseq_day1<-subset_samples(N2_phyloseq,N2_phyloseq@sam_data$merge%in%c('N2-Day1','N2-lawn-Day1'))
head(sample_data(N2_phyloseq_day1)$group, n=10)
diagdds = phyloseq_to_deseq2(N2_phyloseq_day1, ~ group)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #lawn.N2 vs Gut.N2
sigtab
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(N2_phyloseq_day1)[rownames(sigtab), ], "matrix"))
head(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}


# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
d1.lawn.vs.gut<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
d1.lawn.vs.gut
ggsave('differential abundance/lawn.vs.gut_N2_day1.pdf',d1.lawn.vs.gut,width = 5,height = 3)
write.xlsx(sigtab2,'differential abundance/sigtab_lawn.vs.gut_N2_day1.xlsx')

#day 3
N2_phyloseq_day3<-subset_samples(N2_phyloseq,N2_phyloseq@sam_data$merge%in%c('N2-Day3','N2-lawn-Day3'))
head(sample_data(N2_phyloseq_day3)$group, n=10)
diagdds = phyloseq_to_deseq2(N2_phyloseq_day3, ~ group)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #lawn.N2 vs Gut.N2
sigtab
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(N2_phyloseq_day3)[rownames(sigtab), ], "matrix"))
head(sigtab)

# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
d3.lawn.vs.gut<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
d3.lawn.vs.gut
ggsave('differential abundance/lawn.vs.gut_N2_day3.pdf',d3.lawn.vs.gut,width = 5,height = 3)
write.xlsx(sigtab2,'differential abundance/sigtab_lawn.vs.gut_N2_day3.xlsx')

#day 5
N2_phyloseq_day5<-subset_samples(N2_phyloseq,N2_phyloseq@sam_data$merge%in%c('N2-Day5','N2-lawn-Day5'))
head(sample_data(N2_phyloseq_day5)$group, n=10)
diagdds = phyloseq_to_deseq2(N2_phyloseq_day5, ~ group)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #lawn.N2 vs Gut.N2
sigtab
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(N2_phyloseq_day5)[rownames(sigtab), ], "matrix"))
head(sigtab)

# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
d5.lawn.vs.gut<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
d5.lawn.vs.gut
ggsave('differential abundance/lawn.vs.gut_N2_day5.pdf',d5.lawn.vs.gut,width = 5,height = 3)
write.xlsx(sigtab2,'differential abundance/sigtab_lawn.vs.gut_N2_day5.xlsx')

#day 7
N2_phyloseq_day7<-subset_samples(N2_phyloseq,N2_phyloseq@sam_data$merge%in%c('N2-Day7','N2-lawn-Day7'))
head(sample_data(N2_phyloseq_day7)$group, n=10)
diagdds = phyloseq_to_deseq2(N2_phyloseq_day7, ~ group)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #lawn.N2 vs Gut.N2
sigtab
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(N2_phyloseq_day7)[rownames(sigtab), ], "matrix"))
head(sigtab)

# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
d7.lawn.vs.gut<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
d7.lawn.vs.gut
ggsave('differential abundance/lawn.vs.gut_N2_day7.pdf',d7.lawn.vs.gut,width = 5,height = 3)
write.xlsx(sigtab2,'differential abundance/sigtab_lawn.vs.gut_N2_day7.xlsx')



###gut  day1 vs day3   No difference
N2_gut_phyloseq<-subset_samples(N2_phyloseq,N2_phyloseq@sam_data$group%in%c('Gut-N2'))
N2_gut_phyloseq_13<-subset_samples(N2_gut_phyloseq,N2_gut_phyloseq@sam_data$merge %in%c('N2-Day1','N2-Day3'))

head(sample_data(N2_gut_phyloseq_13)$merge, n=15)

diagdds = phyloseq_to_deseq2(N2_gut_phyloseq_13, ~ merge)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #N2.Day3 vs N2.Day1 


###gut  day1 vs day5  N2.Day5 vs N2.Day1    
N2_gut_phyloseq<-subset_samples(N2_phyloseq,N2_phyloseq@sam_data$group%in%c('Gut-N2'))
N2_gut_phyloseq_15<-subset_samples(N2_gut_phyloseq,N2_gut_phyloseq@sam_data$merge %in%c('N2-Day1','N2-Day5'))

head(sample_data(N2_gut_phyloseq_15)$merge, n=15)

diagdds = phyloseq_to_deseq2(N2_gut_phyloseq_15, ~ merge)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #N2.Day5 vs N2.Day1 
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(N2_phyloseq_2)[rownames(sigtab), ], "matrix"))
head(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}


# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
day5.vs.day1<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
day5.vs.day1
ggsave('differential abundance/N2.day5.vs.day1.pdf',day5.vs.day1,width = 5,height = 3)
write.xlsx(sigtab2,'differential abundance/N2.sigtab_day5.vs.day1.xlsx')


###gut  day1 vs day7  
N2_gut_phyloseq<-subset_samples(N2_phyloseq,N2_phyloseq@sam_data$group%in%c('Gut-N2'))
N2_gut_phyloseq_17<-subset_samples(N2_gut_phyloseq,N2_gut_phyloseq@sam_data$merge %in%c('N2-Day1','N2-Day7'))

head(sample_data(N2_gut_phyloseq_17)$merge, n=15)

diagdds = phyloseq_to_deseq2(N2_gut_phyloseq_17, ~ merge)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #N2.Day7 vs N2.Day1 
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(N2_gut_phyloseq_17)[rownames(sigtab), ], "matrix"))
head(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}


# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
N2.d7.vs.d1<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
N2.d7.vs.d1
ggsave('differential abundance/N2.d7.vs.d1.pdf',N2.d7.vs.d1,width = 3,height = 2.5)
write.xlsx(sigtab2,'differential abundance/N2.sigtab_d7.vs.d1.xlsx')

###gut  day5 vs day3  
N2_gut_phyloseq<-subset_samples(N2_phyloseq,N2_phyloseq@sam_data$group%in%c('Gut-N2'))
N2_gut_phyloseq_35<-subset_samples(N2_gut_phyloseq,N2_gut_phyloseq@sam_data$merge %in%c('N2-Day3','N2-Day5'))

head(sample_data(N2_gut_phyloseq_35)$merge, n=15)

diagdds = phyloseq_to_deseq2(N2_gut_phyloseq_35, ~ merge)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #N2.Day5 vs N2.Day3 
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(N2_gut_phyloseq_35)[rownames(sigtab), ], "matrix"))
head(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}


# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
N2.d5.vs.d3<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
N2.d5.vs.d3
ggsave('differential abundance/N2.d5.vs.d3.pdf',N2.d5.vs.d3,width = 3,height = 2.5)
write.xlsx(sigtab2,'differential abundance/N2.sigtab_d5.vs.d3.xlsx')


###gut  day7 vs day3  
N2_gut_phyloseq<-subset_samples(N2_phyloseq,N2_phyloseq@sam_data$group%in%c('Gut-N2'))
N2_gut_phyloseq_37<-subset_samples(N2_gut_phyloseq,N2_gut_phyloseq@sam_data$merge %in%c('N2-Day3','N2-Day7'))

head(sample_data(N2_gut_phyloseq_37)$merge, n=15)

diagdds = phyloseq_to_deseq2(N2_gut_phyloseq_37, ~ merge)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #N2.Day7 vs N2.Day3 
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(N2_gut_phyloseq_37)[rownames(sigtab), ], "matrix"))
head(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}


# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
N2.d7.vs.d3<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
N2.d7.vs.d3
ggsave('differential abundance/N2.d7.vs.d3.pdf',N2.d7.vs.d3,width = 3,height = 2.5)
write.xlsx(sigtab2,'differential abundance/N2.sigtab_d7.vs.d3.xlsx')


###gut  day5 vs day7  no differnece
N2_gut_phyloseq<-subset_samples(N2_phyloseq,N2_phyloseq@sam_data$group%in%c('Gut-N2'))
N2_gut_phyloseq_57<-subset_samples(N2_gut_phyloseq,N2_gut_phyloseq@sam_data$merge %in%c('N2-Day5','N2-Day7'))

head(sample_data(N2_gut_phyloseq_57)$merge, n=15)

diagdds = phyloseq_to_deseq2(N2_gut_phyloseq_57, ~ merge)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #N2.Day7 vs N2.Day5  



###===========================
#heatmap plot
library(pheatmap)
library(ggpubr)
library(readxl)

#####F23  DDA Fig. S21a-b
#gut
dt <- read_excel("./input/F23_deseq2.xlsx", sheet = 'gut')
#remove first col
dt1<-dt[,-c(1)]
rownames(dt1)=dt$`Genus`
library(pheatmap)
#plot  
bk <- c(seq(-10,-0.1,by=0.1),seq(0,10,by=0.1))
pheatmap(dt1, cluster_cols = F, cluster_rows = T, 
         clustering_distance_cols = "euclidean",clustering_method='average',
         color= c(colorRampPalette(colors = c("#510051FF","white"))(length(bk)/2),
                  colorRampPalette(colors = c("white","#005100FF"))(length(bk)/2)),
         legend_breaks=seq(-10,10,2),
         breaks=bk,
         fontsize = 5,fontsize_row = 3, fontsize_col = 3.2,angle_col = 45,
         scale="none",
         border_color= 'grey',
         border=F,
         annotation_names_row=F,annotation_legend=T,treeheight_col=10,
         cellwidth = 8, cellheight = 6,
        # filename='differential abundance/F23_gut_DDA.pdf',width=3,height=3
)

#gut vs lawn
dt <- read_excel("./input/F23_deseq2.xlsx", sheet = 'gut_vs_lawn')
#remove first col
dt1<-dt[,-c(1)]
rownames(dt1)=dt$`Genus`

#plot  
bk <- c(seq(-12,-0.1,by=0.1),seq(0,12,by=0.1))
pheatmap(dt1, cluster_cols = F, cluster_rows = T, 
         clustering_distance_cols = "euclidean",clustering_method='average',
         color= c(colorRampPalette(colors = c("#510051FF","white"))(length(bk)/2),
                  colorRampPalette(colors = c("white","#005100FF"))(length(bk)/2)),
         legend_breaks=seq(-10,10,2),
         breaks=bk,
         fontsize = 5,fontsize_row = 3, fontsize_col = 3.2,angle_col = 45,
         scale="none",
         border_color= 'grey',
         border=F,
         annotation_names_row=F,annotation_legend=T,treeheight_col=10,
         cellwidth = 8, cellheight = 6,
         #filename='differential abundance/F23_gut_vs_lawn_DDA.pdf',width=3,height=3
)



#Fig. S21cd
#####N2  DDA
#gut
dt <- read_excel("./input/N2_deseq2.xlsx", sheet = 'gut')
#remove first col
dt1<-dt[,-c(1)]
rownames(dt1)=dt$`Genus`

#plot  
bk <- c(seq(-10,-0.1,by=0.1),seq(0,10,by=0.1))
pheatmap(dt1, cluster_cols = F, cluster_rows = T, 
         clustering_distance_cols = "euclidean",clustering_method='average',
         color= c(colorRampPalette(colors = c("#510051FF","white"))(length(bk)/2),
                  colorRampPalette(colors = c("white","#005100FF"))(length(bk)/2)),
         legend_breaks=seq(-10,10,2),
         breaks=bk,
         fontsize = 5,fontsize_row = 3, fontsize_col = 3.2,angle_col = 45,
         scale="none",
         border_color= 'grey',
         border=F,
         annotation_names_row=F,annotation_legend=T,treeheight_col=10,
         cellwidth = 8, cellheight = 6,
         #filename='differential abundance/N2_gut_DDA.pdf',width=3,height=3
)

#gut vs lawn
dt <- read_excel("./input/N2_deseq2.xlsx", sheet = 'gut_vs_lawn')
#remove first col
dt1<-dt[,-c(1)]
rownames(dt1)=dt$`Genus`

#plot  
bk <- c(seq(-12,-0.1,by=0.1),seq(0,12,by=0.1))
pheatmap(dt1, cluster_cols = F, cluster_rows = T, 
         clustering_distance_cols = "euclidean",clustering_method='average',
         color= c(colorRampPalette(colors = c("#510051FF","white"))(length(bk)/2),
                  colorRampPalette(colors = c("white","#005100FF"))(length(bk)/2)),
         legend_breaks=seq(-10,10,2),
         breaks=bk,
         fontsize = 5,fontsize_row = 3, fontsize_col = 3.2,angle_col = 45,
         scale="none",
         border_color= 'grey',
         border=F,
         annotation_names_row=F,annotation_legend=T,treeheight_col=10,
         cellwidth = 8, cellheight = 6,
         #filename='differential abundance/N2_gut_vs_lawn_DDA.pdf',width=3,height=3
)


#====================================================================================================
#HQbiome-1 analysis
physeq <- qza_to_phyloseq(
  features = "input/HQbiome-1_qza/table_deblur.qza", # ASV/OTU table
  taxonomy = "input/HQbiome-1_qza/tax.qza", # Taxonomy file
  metadata = "input/HQbiome-1_qza/meta.txt" # Sample metadata
)
physeq

#DAY1 AND DAY5 MICROBIOTA
physeq_d1d5 <- qza_to_phyloseq(
  features = "./input/d1d5_qza/table_deblur.qza", # ASV/OTU table
  taxonomy = "./input/d1d5_qza/tax.qza", # Taxonomy file
  metadata = "./input/d1d5_qza/meta.txt" # Sample metadata
)
physeq_d1d5

#merge all and d1d5
physeq_all2<-merge_phyloseq(physeq, physeq_d1d5)
#select samples
physeq_all<-subset_samples(physeq_all2, !(name %in% c('FL3.1','FL3.2','NL3.2','NL3.3','NL3.4',
                                                      'NL5.1','NL5.2','NL5.3',
                                                      'NL7.1','NL7.2','NL7.3','NL7.4',
                                                      'F1.1-all','F1.2-all', 'F5.4',
                                                      'N3.4','N5.1','N5.5','N7.2','N7.4','N7.5')))

#merge genus
physeq_all <- physeq_all %>%
  tax_glom(taxrank = "Genus", NArm = TRUE) # agglomerate on Genus level
tax_table(physeq_all) <- tax_table(physeq_all)[,1:6]
#order
physeq_all@sam_data$group<-factor(physeq_all@sam_data$group,levels = c('Microbiota-F23','Microbiota-N2','Lawn-F23','Lawn-N2','Inoculum'))

physeq_all@sam_data$name <- factor(physeq_all@sam_data$name, 
                                   levels = c( "F1.4","F1.5","F3.5","F3.6", "F3.7","F5.5","F5.6","F5.7",
                                               "N1.2","N1.3","N3.1","N3.5","N3.6","N5.2","N5.4","N5.6","N7.1","N7.3",
                                               "FL1.1","FL1.2","FL1.3","FL3.3","FL3.5","FL3.4","FL5.1","FL5.2","FL5.3",
                                               "NL1.1", "NL1.2","NL1.3","NL3.1","NL3.5","NL3.6", "NL5.4","NL5.5","NL5.6","NL7.5","NL7.6","NL7.7",
                                               "Inoculum-1","Inoculum-2","Inoculum-3"))
#Fig. 6b
#F23
F23_phyloseq<-subset_samples(physeq_all,physeq_all@sam_data$group%in%c('Microbiota-F23','Lawn-F23'))
#merge samples
F23_phyloseq_merge<-merge_samples(F23_phyloseq,'merge')
#F23  plot
F23_physeq_other_Genus <- F23_phyloseq_merge %>%
  tax_glom(taxrank = "Genus") %>%    # agglomerate on Class level
  transform_sample_counts(function(x) {100 * x/sum(x)}) %>%  # transform to relative abundance
  psmelt() %>%  # convert from phyloseq object to a long data frame
  mutate(Taxa_Order = Genus)# %>%
#mutate(Taxa_Order = replace(Taxa_Order, Abundance < 1, "<1%"))


#colourCount = length(unique(physeq_other_Genus$Taxa_Order))  # obtain number of colors needed
#getPalette = colorRampPalette(brewer.pal(9, "Accent")) # create palette with Set1 ,3

F23_Genus_BP <- 
  ggplot(data=F23_physeq_other_Genus, 
         aes(x=Sample, # order by date for each library season
             y=Abundance, 
             fill=fct_reorder(Taxa_Order, Abundance, .desc = TRUE))) + # order taxa within bars by abundance
  geom_bar(stat="identity", 
           position = position_stack(reverse = TRUE), # order taxa in reverse
           width = 0.9) + 
  #geom_col(width = 1)+
  scale_fill_manual(values = c('Vibrio'='#60BD68FF','Paracoccus'='#006064FF','Parasphingorhabdus'='#F17CB0FF',
                               'Pseudomonas'='#C6242DFF','Rheinheimera'='#004949FF','Enterovibrio'='#DECF3FFF',
                               'Shewanella'='#3A82E4FF','Neptunomonas'='#66E3D9FF','<1%'='#999999FF')
                    #paletteer_d("ggsci::default_igv",50,direction=-1)
                    , # assign legend colors
                    guide = guide_legend(reverse = TRUE)) + # make legend match order of taxa in barplot
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), # Y axis by 10% increments
                     expand = c(0,0)) + # remove whitespace in plot
  facet_grid(~group, # draw facets by date
             scales = "free_x", 
             space = "free_x",
             drop = F
  ) + 
  theme(legend.direction="vertical", 
        axis.text.x = element_text(color = 'black',angle = 45, hjust = 1),
        #element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y = element_text(#face = "bold", 
          size = 10,color = 'black'),
        strip.text = element_text(face = "bold",#, size = 15
        ),
        strip.background = element_rect(color = "white", fill = "#DDEC7EFF"
        ),
        legend.title = element_text( face = "bold"),
        plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
        panel.spacing = unit(0, "lines"),
        legend.text = element_text(face="italic"),
        panel.background = element_rect(fill = 'white', color = 'white'), 
        panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
  labs(  face = "bold",
         y = "Relative Abundance %", 
         fill = "Genus")+
  guides(fill = guide_legend(ncol = 1))
F23_Genus_BP

#Fig. 7b
#N2
N2_phyloseq<-subset_samples(physeq_all,physeq_all@sam_data$group%in%c('Microbiota-N2','Lawn-N2'))
#merge samples
N2_phyloseq_merge<-merge_samples(N2_phyloseq,'merge')
#F23  plot
N2_physeq_other_Genus <- N2_phyloseq_merge %>%
  tax_glom(taxrank = "Genus") %>%    # agglomerate on Class level
  transform_sample_counts(function(x) {100 * x/sum(x)}) %>%  # transform to relative abundance
  psmelt() %>%  # convert from phyloseq object to a long data frame
  mutate(Taxa_Order = Genus)# %>%
#mutate(Taxa_Order = replace(Taxa_Order, Abundance < 1, "<1%"))


#colourCount = length(unique(physeq_other_Genus$Taxa_Order))  # obtain number of colors needed
#getPalette = colorRampPalette(brewer.pal(9, "Accent")) # create palette with Set1 ,3

N2_Genus_BP <- 
  ggplot(data=N2_physeq_other_Genus, 
         aes(x=Sample, # order by date for each library season
             y=Abundance, 
             fill=fct_reorder(Taxa_Order, Abundance, .desc = TRUE))) + # order taxa within bars by abundance
  geom_bar(stat="identity", 
           position = position_stack(reverse = TRUE), # order taxa in reverse
           width = 0.9) + 
  #geom_col(width = 1)+
  scale_fill_manual(values = c('Pseudomonas'='#C6242DFF', 'Enterovibrio'='#DECF3FFF',
                               'Vibrio'='#60BD68FF','Paracoccus'='#006064FF','Alteromonas'='#F49538FF',
                               'Rheinheimera'='#004949FF','Shewanella'='#3A82E4FF','Neptunomonas'='#66E3D9FF',
                               'Parasphingorhabdus'='#F17CB0FF','Loktanella'='#006DDBFF','<1%'='#999999FF')
                    , # assign legend colors
                    guide = guide_legend(reverse = TRUE)) + # make legend match order of taxa in barplot
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), # Y axis by 10% increments
                     expand = c(0,0)) + # remove whitespace in plot
  facet_grid(~group, # draw facets by date
             scales = "free_x", 
             space = "free_x",
             drop = F
  ) + #switch="x" 
  theme(legend.direction="vertical", 
        axis.text.x = element_text(color = 'black',angle = 45, hjust = 1),
        #element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y = element_text(#face = "bold", 
          size = 10,color = 'black'),
        strip.text = element_text(face = "bold",#, size = 15
        ),
        strip.background = element_rect(color = "white", fill = "#DDEC7EFF"
        ),
        legend.title = element_text( face = "bold"),
        plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
        panel.spacing = unit(0, "lines"),
        legend.text = element_text(face="italic"),
        panel.background = element_rect(fill = 'white', color = 'white'), 
        panel.border = element_rect(colour = "black", fill=NA, size=1.5)
  ) +
  labs(  face = "bold",
         y = "Relative Abundance %", 
         fill = "Genus")+
  guides(fill = guide_legend(ncol = 1))+
  theme()
N2_Genus_BP

#Fig. S20b
#Inoculum  plot
ino_phyloseq<-subset_samples(physeq_all,physeq_all@sam_data$group%in%c('Inoculum'))
ino_physeq_other_Genus <- ino_phyloseq %>%
  tax_glom(taxrank = "Genus") %>%    # agglomerate on Class level
  transform_sample_counts(function(x) {100 * x/sum(x)}) %>%  # transform to relative abundance
  psmelt() %>%  # convert from phyloseq object to a long data frame
  mutate(Taxa_Order = Genus)# %>%
#mutate(Taxa_Order = replace(Taxa_Order, Abundance < 1, "<1%"))


ino_Genus_BP <- 
  ggplot(data=ino_physeq_other_Genus, 
         aes(x=Sample, # order by date for each library season
             y=Abundance, 
             fill=fct_reorder(Taxa_Order, Abundance, .desc = TRUE))) + # order taxa within bars by abundance
  geom_bar(stat="identity", 
           position = position_stack(reverse = TRUE), # order taxa in reverse
           width = 0.9) + 
  #geom_col(width = 1)+
  scale_fill_manual(values = c('Vibrio'='#60BD68FF', 'Shewanella'='#3A82E4FF','Psychrobacter'='#8D4B08FF',
                               'Rheinheimera'='#004949FF','Parasphingorhabdus'='#F17CB0FF','Aliivibrio'='#244579FF',
                               'Pseudomonas'='#C6242DFF', 'Enterovibrio'='#DECF3FFF','Lacinutrix'='#F17CB0FF',
                               'Sulfitobacter'='#009E73FF','Pseudoalteromonas'='#D25D38FF','Marinobacter'='#66E3D9FF',
                               'Thalassotalea'='#00FF7FFF','HQ336491_g'='#FCA3B7FF','Neptunomonas'='#66E3D9FF',
                               'Litoreibacter'='#98FB98FF','Saccharospirillum'='#D070B9FF','Yoonia'='#66CDAAFF',
                               'Bacillus'='#FF9800FF','Alteromonas'='#F49538FF','Marinagarivorans'='#7D5329FF',
                               'Postechiella'='#F15854FF','Marinomonas'='#8A6842FF','Salinimonas'='#E7D202FF',
                               'Staphylococcus'='#BF616AFF','Alkalihalobacillus'='#F0E442FF','Lactococcus'='#E9A820FF',
                               'Priestia'='#8FA87AFF','Thaumasiovibrio'='#5E81ACFF','Cytobacillus'='#FFFF6DFF',
                               'Galactobacter'='#3F51B5FF','Planococcus'='#24FF24FF','Paracoccus'='#006064FF',
                               'Halobacillus'='#24FF24FF','Metabacillus'='#CC79A7FF','Niallia'='#DB6D00FF',
                               'Rhodococcus'='#924900FF','Paenisporosarcina'='#00FF7FFF','Colwellia'='#920000FF',
                               'Microbacterium'='#B6DBFFFF','Loktanella'='#6DB6FFFF','Arenivirga'='#B66DFFFF'
  )
  #paletteer_d("ggsci::default_igv",50,direction=-1)
  , # assign legend colors
  guide = guide_legend(reverse = TRUE)) + # make legend match order of taxa in barplot
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), # Y axis by 10% increments
                     expand = c(0,0)) + # remove whitespace in plot
  facet_grid(~group, # draw facets by date
             scales = "free_x", 
             space = "free_x",
             drop = F
  ) + #switch="x" 
  theme(legend.direction="vertical", 
        axis.text.x = element_text(color = 'black',angle = 45, hjust = 1),
        #element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y = element_text(#face = "bold", 
          size = 10,color = 'black'),
        strip.text = element_text(face = "bold",#, size = 15
        ),
        strip.background = element_rect(color = "white", fill = "#DDEC7EFF"
        ),
        legend.title = element_text( face = "bold"),
        plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
        panel.spacing = unit(0, "lines"),
        legend.text = element_text(face="italic"),
        panel.background = element_rect(fill = 'white', color = 'white'),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5)
  ) +
  labs(  face = "bold",
         y = "Relative Abundance %", 
         fill = "Genus")+
  guides(fill = guide_legend(ncol = 3))
ino_Genus_BP

#Fig. S22b
a_phyloseq<-subset_samples(physeq_all, physeq_all@sam_data$merge%in%c('F23-gut-day3','F23-lawn-day3', 
                                                                      'F23-gut-day5','F23-lawn-day5',
                                                                      'N2-gut-day5','N2-lawn-day5',
                                                                      "N2-gut-day7","N2-lawn-day7",
                                                                      'Inoculum'))
a_phyloseq@sam_data$group<-factor(a_phyloseq@sam_data$group,levels = c('Microbiota-F23','Lawn-F23','Microbiota-N2','Lawn-N2','Inoculum'))
genus_db<-read_excel('input/HQbiome-1-tax.xlsx',sheet = 'pcoa2')
genus_db1<-genus_db[,-1]
row.names(genus_db1)<-genus_db$ID
distMatrix <- vegdist(genus_db1,method = "bray")
pCoa <- cmdscale(distMatrix, eig = T,k = 2 )
varExp <- (eigenvals(pCoa)/sum(eigenvals(pCoa)))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)

pCoaVecs <- as.data.frame(pCoa$points)
colnames(pCoaVecs) <- paste0("PCo",c(1:2))
pCoaVecs$ID <- row.names(pCoaVecs)

meta<-read_excel('input/HQbiome-1-tax.xlsx',sheet = 'meta')
pCoaVecs2<-merge(pCoaVecs,meta,by='ID')

#PERMANOVA adonis2
adonis_result <- adonis2(genus_db1 ~ group, data=meta, permutations=999)


library(paletteer) 
pcoa.p2<-ggplot(pCoaVecs2,aes(x=PCo1,y=PCo2,color=group)) + 
  geom_point(size=2) + 
  theme_classic() + 
  #scale_color_gradientn(colours = rainbow(5)) + 
  xlab(paste0('PCoA1 (',round(xVar,2),' %)')) + 
  ylab(paste0('PCoA2 (',round(yVar,2),' %)')) + 
  theme(axis.text = element_text(color = 'black'))+
  scale_colour_manual(name = "Group", values = paletteer_d("basetheme::void"))+
  stat_ellipse(data = pCoaVecs2, aes(x = PCo1, y = PCo2),
               level = 0.68, type='norm',size=0.6)+
  annotate("text",x=-0.2,y=-0.5,label=paste("p= ", adonis_result$`Pr(>F)`[1]),size=3)
#ggforce::geom_mark_ellipse(aes(group=group),expand = unit(1, "mm"),tol =0.05)
#save size =
pcoa.p2

################################################################################################3
#Fig. S21e-g
#differential abundance analysis
library(DESeq2)
library(xlsx)

#F23
F23_phyloseq<-subset_samples(physeq_all,physeq_all@sam_data$group%in%c('Microbiota-F23','Lawn-F23'))
#merge samples
#F23_phyloseq_merge<-merge_samples(F23_phyloseq,'merge')

##======deseq2   gut and lawn
#day 1
F23_phyloseq_day1<-subset_samples(F23_phyloseq,F23_phyloseq@sam_data$merge%in%c('F23-gut-day1','F23-lawn-day1'))
head(sample_data(F23_phyloseq_day1)$group, n=10)
diagdds = phyloseq_to_deseq2(F23_phyloseq_day1, ~ group)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #lawn.F23 vs Gut.F23
sigtab
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(F23_phyloseq_day1)[rownames(sigtab), ], "matrix"))
head(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}


# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
d1.lawn.vs.gut<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
d1.lawn.vs.gut
ggsave('differential abundance/lawn.vs.gut_F23_day1.pdf',d1.lawn.vs.gut,width = 5,height = 3)
write.xlsx(sigtab2,'differential abundance/sigtab_lawn.vs.gut_F23_day1.xlsx')

#day 3   #NS
F23_phyloseq_day3<-subset_samples(F23_phyloseq,F23_phyloseq@sam_data$merge%in%c('F23-gut-day3','F23-lawn-day3'))
head(sample_data(F23_phyloseq_day3)$group, n=10)
diagdds = phyloseq_to_deseq2(F23_phyloseq_day3, ~ group)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #lawn.F23 vs Gut.F23
sigtab

##NS

#day 5
F23_phyloseq_day5<-subset_samples(F23_phyloseq,F23_phyloseq@sam_data$merge%in%c('F23-gut-day5','F23-lawn-day5'))
head(sample_data(F23_phyloseq_day5)$group, n=10)
diagdds = phyloseq_to_deseq2(F23_phyloseq_day5, ~ group)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #lawn.F23 vs Gut.F23
sigtab
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(F23_phyloseq_day5)[rownames(sigtab), ], "matrix"))
head(sigtab)

# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
d5.lawn.vs.gut<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
d5.lawn.vs.gut
ggsave('differential abundance/lawn.vs.gut_F23_day5.pdf',d5.lawn.vs.gut,width = 5,height = 3)
write.xlsx(sigtab2,'differential abundance/sigtab_lawn.vs.gut_F23_day5.xlsx')


###gut  day1 vs day3
F23_gut_phyloseq<-subset_samples(F23_phyloseq,F23_phyloseq@sam_data$group%in%c('Microbiota-F23'))
F23_gut_phyloseq_13<-subset_samples(F23_gut_phyloseq,F23_gut_phyloseq@sam_data$merge %in%c('F23-gut-day1','F23-gut-day3'))

head(sample_data(F23_gut_phyloseq_13)$merge, n=15)

diagdds = phyloseq_to_deseq2(F23_gut_phyloseq_13, ~ merge)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #F23.Day3 vs F23.Day1 
sigtab
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(F23_gut_phyloseq_13)[rownames(sigtab), ], "matrix"))
head(sigtab)

# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
d3.vs.d1<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
d3.vs.d1
ggsave('differential abundance/F23.d3.vs.d1.pdf',d3.vs.d1,width = 3,height = 2.5)
write.xlsx(sigtab2,'differential abundance/F23_sigtab_d3.vs.d1.xlsx')

###gut  day1 vs day5
F23_gut_phyloseq<-subset_samples(F23_phyloseq,F23_phyloseq@sam_data$group%in%c('Microbiota-F23'))
F23_gut_phyloseq_15<-subset_samples(F23_gut_phyloseq,F23_gut_phyloseq@sam_data$merge %in%c('F23-gut-day1','F23-gut-day5'))

head(sample_data(F23_gut_phyloseq_15)$merge, n=15)

diagdds = phyloseq_to_deseq2(F23_gut_phyloseq_15, ~ merge)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]   #F23.Day5 vs F23.Day1 
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(F23_gut_phyloseq_15)[rownames(sigtab), ], "matrix"))
head(sigtab)

# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
d5.vs.d1<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
d5.vs.d1
ggsave('differential abundance/f23.d5.vs.d1.pdf',d5.vs.d1,width = 3,height = 2.5)
write.xlsx(sigtab2,'differential abundance/F23_sigtab_d5.vs.d1.xlsx')


###gut  day3 vs day5 
F23_gut_phyloseq<-subset_samples(F23_phyloseq,F23_phyloseq@sam_data$group%in%c('Microbiota-F23'))
F23_gut_phyloseq_35<-subset_samples(F23_gut_phyloseq,F23_gut_phyloseq@sam_data$merge %in%c('F23-gut-day3','F23-gut-day5'))

head(sample_data(F23_gut_phyloseq_35)$merge, n=15)

diagdds = phyloseq_to_deseq2(F23_gut_phyloseq_35, ~ merge)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]   #F23.Day5 vs F23.Day3  
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(F23_gut_phyloseq_35)[rownames(sigtab), ], "matrix"))
head(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}


# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
d5.vs.d3<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
d5.vs.d3
ggsave('differential abundance/F23.d5.vs.d3.pdf',d5.vs.d3,width = 3,height = 2.5)
write.xlsx(sigtab2,'differential abundance/F23_sigtab_d5.vs.d3.xlsx')


####################################################################################
#deseq2   N2

N2_phyloseq<-subset_samples(physeq_all,physeq_all@sam_data$group%in%c('Microbiota-N2','Lawn-N2'))
#merge samples
#F23_phyloseq_merge<-merge_samples(F23_phyloseq,'merge')

##======deseq2   gut and lawn
#day 1
N2_phyloseq_day1<-subset_samples(N2_phyloseq,N2_phyloseq@sam_data$merge%in%c('N2-gut-day1','N2-lawn-day1'))
head(sample_data(N2_phyloseq_day1)$group, n=10)
diagdds = phyloseq_to_deseq2(N2_phyloseq_day1, ~ group)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #lawn.N2 vs Gut.N2
sigtab
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(N2_phyloseq_day1)[rownames(sigtab), ], "matrix"))
head(sigtab)

# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
d1.lawn.vs.gut<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
d1.lawn.vs.gut
ggsave('differential abundance/lawn.vs.gut_N2_day1.pdf',d1.lawn.vs.gut,width = 5,height = 3)
write.xlsx(sigtab2,'differential abundance/sigtab_lawn.vs.gut_N2_day1.xlsx')

#day 3
N2_phyloseq_day3<-subset_samples(N2_phyloseq,N2_phyloseq@sam_data$merge%in%c('N2-gut-day3','N2-lawn-day3'))
head(sample_data(N2_phyloseq_day3)$group, n=10)
diagdds = phyloseq_to_deseq2(N2_phyloseq_day3, ~ group)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #lawn.N2 vs Gut.N2
sigtab
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(N2_phyloseq_day3)[rownames(sigtab), ], "matrix"))
head(sigtab)

# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
d3.lawn.vs.gut<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
d3.lawn.vs.gut
ggsave('differential abundance/lawn.vs.gut_N2_day3.pdf',d3.lawn.vs.gut,width = 5,height = 3)
write.xlsx(sigtab2,'differential abundance/sigtab_lawn.vs.gut_N2_day3.xlsx')

#day 5
N2_phyloseq_day5<-subset_samples(N2_phyloseq,N2_phyloseq@sam_data$merge%in%c('N2-gut-day5','N2-lawn-day5'))
head(sample_data(N2_phyloseq_day5)$group, n=10)
diagdds = phyloseq_to_deseq2(N2_phyloseq_day5, ~ group)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #lawn.N2 vs Gut.N2
sigtab
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(N2_phyloseq_day5)[rownames(sigtab), ], "matrix"))
head(sigtab)

# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
d5.lawn.vs.gut<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
d5.lawn.vs.gut
ggsave('differential abundance/lawn.vs.gut_N2_day5.pdf',d5.lawn.vs.gut,width = 5,height = 3)
write.xlsx(sigtab2,'differential abundance/sigtab_lawn.vs.gut_N2_day5.xlsx')

#day 7
N2_phyloseq_day7<-subset_samples(N2_phyloseq,N2_phyloseq@sam_data$merge%in%c('N2-gut-day7','N2-lawn-day7'))
head(sample_data(N2_phyloseq_day7)$group, n=10)
diagdds = phyloseq_to_deseq2(N2_phyloseq_day7, ~ group)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #lawn.N2 vs Gut.N2
sigtab
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(N2_phyloseq_day7)[rownames(sigtab), ], "matrix"))
head(sigtab)

# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
d7.lawn.vs.gut<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
d7.lawn.vs.gut
ggsave('differential abundance/lawn.vs.gut_N2_day7.pdf',d7.lawn.vs.gut,width = 5,height = 3)
write.xlsx(sigtab2,'differential abundance/sigtab_lawn.vs.gut_N2_day7.xlsx')

###gut  day1 vs day3     NS
N2_gut_phyloseq<-subset_samples(N2_phyloseq, N2_phyloseq@sam_data$group%in%c('Microbiota-N2'))
N2_gut_phyloseq_13<-subset_samples(N2_gut_phyloseq, N2_gut_phyloseq@sam_data$merge %in%c('N2-gut-day1','N2-gut-day3'))

head(sample_data(N2_gut_phyloseq_13)$merge, n=15)

diagdds = phyloseq_to_deseq2(N2_gut_phyloseq_13, ~ merge)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #N2.Day3 vs N2.Day1 
sigtab
#NS

###gut  day1 vs day5  N2.Day5 vs N2.Day1   no difference
N2_gut_phyloseq<-subset_samples(N2_phyloseq,N2_phyloseq@sam_data$group%in%c('Microbiota-N2'))
N2_gut_phyloseq_15<-subset_samples(N2_gut_phyloseq,N2_gut_phyloseq@sam_data$merge %in%c('N2-gut-day1','N2-gut-day5'))

head(sample_data(N2_gut_phyloseq_15)$merge, n=15)

diagdds = phyloseq_to_deseq2(N2_gut_phyloseq_15, ~ merge)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #N2.Day5 vs N2.Day1 
N2-gut-day  #NS


###gut  day1 vs day7  
N2_gut_phyloseq<-subset_samples(N2_phyloseq,N2_phyloseq@sam_data$group%in%c('Microbiota-N2'))
N2_gut_phyloseq_17<-subset_samples(N2_gut_phyloseq, N2_gut_phyloseq@sam_data$merge %in%c('N2-gut-day1','N2-gut-day7'))

head(sample_data(N2_gut_phyloseq_17)$merge, n=15)

diagdds = phyloseq_to_deseq2(N2_gut_phyloseq_17, ~ merge)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #N2.Day7 vs N2.Day1 

#NS
###gut  day5 vs day3  


N2_gut_phyloseq_35<-subset_samples(N2_gut_phyloseq,N2_gut_phyloseq@sam_data$merge %in%c('N2-gut-day3','N2-gut-day5'))

head(sample_data(N2_gut_phyloseq_35)$merge, n=15)

diagdds = phyloseq_to_deseq2(N2_gut_phyloseq_35, ~ merge)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #N2.Day5 vs N2.Day3 
sigtab  #NS


###gut  day7 vs day3  

N2_gut_phyloseq_37<-subset_samples(N2_gut_phyloseq,N2_gut_phyloseq@sam_data$merge %in%c('N2-gut-day3','N2-gut-day7'))

head(sample_data(N2_gut_phyloseq_37)$merge, n=15)

diagdds = phyloseq_to_deseq2(N2_gut_phyloseq_37, ~ merge)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #N2.Day7 vs N2.Day3 
sigtab  #NS

###gut  day5 vs day7  no differnece

N2_gut_phyloseq_57<-subset_samples(N2_gut_phyloseq,N2_gut_phyloseq@sam_data$merge %in%c('N2-gut-day5','N2-gut-day7'))

head(sample_data(N2_gut_phyloseq_57)$merge, n=15)

diagdds = phyloseq_to_deseq2(N2_gut_phyloseq_57, ~ merge)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #N2.Day7 vs N2.Day5  
sigtab
#NS

###===========================
#heatmap plot Fig. S21e-g
library(pheatmap)
library(ggpubr)
library(readxl)

#####F23  DDA
#gut
dt <- read_excel("./input/F23_HQbiome-1_deseq2.xlsx", sheet = 'gut')
#remove first col
dt1<-dt[,-c(1)]
rownames(dt1)=dt$`Genus`

#plot  
bk <- c(seq(-10,-0.1,by=0.1),seq(0,10,by=0.1))
pheatmap(dt1, cluster_cols = F, cluster_rows = T, 
         clustering_distance_cols = "euclidean",clustering_method='average',
         color= c(colorRampPalette(colors = c("#510051FF","white"))(length(bk)/2),
                  colorRampPalette(colors = c("white","#005100FF"))(length(bk)/2)),
         legend_breaks=seq(-10,10,2),
         breaks=bk,
         fontsize = 5,fontsize_row = 3, fontsize_col = 3.2,angle_col = 45,
         scale="none",
         border_color= 'grey',
         border=F,
         annotation_names_row=F,annotation_legend=T,treeheight_col=10,
         cellwidth = 8, cellheight = 6,
         #filename='differential abundance/F23_gut_DDA.pdf',width=3,height=3
)

#gut vs lawn
dt <- read_excel("./input/F23_HQbiome-1_deseq2.xlsx", sheet = 'gut_vs_lawn')
#remove first col
dt1<-dt[,-c(1)]
rownames(dt1)=dt$`Genus`

#plot  
bk <- c(seq(-12,-0.1,by=0.1),seq(0,12,by=0.1))
pheatmap(dt1, cluster_cols = F, cluster_rows = T, 
         clustering_distance_cols = "euclidean",clustering_method='average',
         color= c(colorRampPalette(colors = c("#510051FF","white"))(length(bk)/2),
                  colorRampPalette(colors = c("white","#005100FF"))(length(bk)/2)),
         legend_breaks=seq(-10,10,2),
         breaks=bk,
         fontsize = 5,fontsize_row = 3, fontsize_col = 3.2,angle_col = 45,
         scale="none",
         border_color= 'grey',
         border=F,
         annotation_names_row=F,annotation_legend=T,treeheight_col=10,
         cellwidth = 8, cellheight = 6,
        # filename='differential abundance/F23_gut_vs_lawn_DDA.pdf',width=3,height=3
)




#####N2  DDA
#gut
#none significance

#gut vs lawn
dt <- read_excel("./input/N2_HQbiome-1_deseq2.xlsx", sheet = 'gut_vs_lawn')
#remove first col
dt1<-dt[,-c(1)]
rownames(dt1)=dt$`Genus`

#plot  
bk <- c(seq(-12,-0.1,by=0.1),seq(0,12,by=0.1))
pheatmap(dt1, cluster_cols = F, cluster_rows = T, 
         clustering_distance_cols = "euclidean",clustering_method='average',
         color= c(colorRampPalette(colors = c("#510051FF","white"))(length(bk)/2),
                  colorRampPalette(colors = c("white","#005100FF"))(length(bk)/2)),
         legend_breaks=seq(-10,10,2),
         breaks=bk,
         fontsize = 5,fontsize_row = 3, fontsize_col = 3.2,angle_col = 45,
         scale="none",
         border_color= 'grey',
         border=F,
         annotation_names_row=F,annotation_legend=T,treeheight_col=10,
         cellwidth = 8, cellheight = 6,
         #filename='differential abundance/N2_gut_vs_lawn_DDA.pdf',width=3,height=3
)


#Fig. 6d
library(readxl)
library(rstatix)
data <- read_excel("input/HQBiome_lifespan_F23.xlsx",sheet = '!3GROUP-f23')
#lifespan
shapiro.test(data$lifespan) 
bartlett.test(lifespan ~ GROUP, data = data)   
#growth-f23 
shapiro.test(data$d5) 
bartlett.test(d5 ~ GROUP, data = data)   
#ANOVA 
data %>% anova_test(lifespan ~ GROUP)  # p=0.007
data %>% anova_test(d5 ~ GROUP)  #  p=4.03e-05
library(ggplot2)
library(ggpubr)
p.d5<-data %>% 
  tukey_hsd(d5 ~ GROUP, paired = T) %>% add_y_position()

p.lifespan<- data %>%  
  tukey_hsd(lifespan ~ GROUP, paired = T)  %>% add_y_position(step.increase = 0.48)
#Fig. 5c
F23<-ggboxplot(
  data, x = "GROUP", y = "d5",
  color = "GROUP", palette = c('#004586FF','#579D1CFF','#FF420EFF'), ylim = c(0,0.8),
  ylab='Growth rate (%)',xlab = '',size=1
)+
  stat_pvalue_manual(p.d5, label = "p.adj.signif")+rotate_x_text(45,size=10)
F23

lifespan<-ggboxplot(
  data, x = "GROUP", y = "lifespan",
  color = "GROUP", palette = c('#004586FF','#579D1CFF','#FF420EFF'), ylim = c(5,12),
  ylab='Mean lifespan (day)',xlab = '',size=1
)+
  stat_pvalue_manual(p.lifespan, label = "p.adj.signif",hide.ns=T)+rotate_x_text(45,size=10)
lifespan

#life curve
library(tidyverse)
library(ggsurvfit)
library(gghighlight)
library(tidycmprsk)
library(paletteer) 
data2 <- read_excel("input/HQBiome_lifespan_F23.xlsx",sheet = '!ggsurvfit-f23')
data2$group<-factor(data2$group,levels=c('HQBiome','Group1','Group2'))
#group1 and group2 and HQbiome
#Fig. 6f
surcical1<-
  survfit2(Surv(time, status) ~ group, data = data2) %>% 
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
  guides(colour =guide_legend(ncol=1))+
  scale_color_manual(values = c('#004586FF','#579D1CFF','#FF420EFF'))+
  scale_fill_manual(values = c('#004586FF','#579D1CFF','#FF420EFF'))+
  scale_x_continuous(breaks = seq(0, 18, by = 1))+
  theme(axis.text = element_text(color="black"))
#theme(legend.position="none")
surcical1

#log rank p value
library(survival)
#HQBiome-G1  
G1_ALL<- subset(data2,group!='Group2')
survdiff(Surv(time, status) ~ group, data = G1_ALL)  
#Chisq= 4.2  on 1 degrees of freedom, p= 0.04 

#HQBiome-G2
G2_ALL<- subset(data2,group!='Group1')
survdiff(Surv(time, status) ~ group, data = G2_ALL)
#Chisq= 11.6  on 1 degrees of freedom, p= 7e-04

#G1-G2
G1_G2<- subset(data2,group!='HQBiome')
survdiff(Surv(time, status) ~ group, data = G1_G2)
# Chisq= 1.8  on 1 degrees of freedom, p= 0.2 

#==========================================================
###  N2  Fig. 7de
data <- read_excel("input/HQbiome_lifespan_N2.xlsx",sheet = 'plot_mean')


#lifespan
shapiro.test(data$lifespan)  
bartlett.test(lifespan ~ GROUP, data = data)   

#growth-f23 
shapiro.test(data$egg_laying_time)  #    
bartlett.test(egg_laying_time ~ GROUP, data = data)   


##non para --growth
library(rstatix)
data %>% kruskal_test(egg_laying_time ~ GROUP)  # P=0.0226

#wilcox 
p.growth_n2 <- data %>% 
  wilcox_test(egg_laying_time ~ GROUP, p.adjust.method = "fdr")

#dunn_test  
p.growth_n2 <- data %>% 
  dunn_test(egg_laying_time ~ GROUP, p.adjust.method = "fdr")

write.csv(p.lifespan_n2,'./N2/p.lifespan_n2_dunn_test.csv')

#ANOVA 
data %>% anova_test(lifespan ~ GROUP)  #p=0.004

library(ggplot2)
library(ggpubr)
p.growth.n2<-data %>% 
  tukey_hsd(egg_laying_time ~ GROUP, paired = T) %>% add_y_position()

p.lifespan.n2 <- data %>%  
  tukey_hsd(lifespan ~ GROUP, paired = T)  %>% add_y_position()
#Fig. 7de
n2_growth<-ggboxplot(
  data, x = "GROUP", y = "egg_laying_time",
  color = "GROUP", palette = c('#004586FF','#579D1CFF','#FF420EFF'), #ylim = c(55,70),
  ylab='Egg laying time (h)',xlab = '',size=1
)+
  stat_pvalue_manual(p.growth.n2, label = "p.adj.signif")+rotate_x_text(45,size=10)
n2_growth

n2_lifespan<-ggboxplot(
  data, x = "GROUP", y = "lifespan",
  color = "GROUP", palette = c('#004586FF','#579D1CFF','#FF420EFF'),# ylim = c(5,12),
  ylab='Mean lifespan (day)',xlab = '',size=1
)+
  stat_pvalue_manual(p.lifespan.n2, label = "p.adj.signif",hide.ns=T)+rotate_x_text(45,size=10)
n2_lifespan

merge<-ggarrange(n2_growth,n2_lifespan,common.legend = T,legend = 'top')
merge
ggsave('./N2/merge.pdf',merge,height = 3,width = 5)

# survival curve
#Fig. 7f
library(tidyverse)
library(ggsurvfit)
library(gghighlight)
library(tidycmprsk)
library(paletteer) 
data2 <- read_excel("input/HQbiome_lifespan_N2.xlsx",sheet = 'survival')
data2$group<-factor(data2$group,levels=c('HQbiome','HQbiome-1','HQbiome-2'))

surcical1<-
  survfit2(Surv(time, status) ~ group, data = data2) %>% 
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
  guides(colour =guide_legend(ncol=1))+
  scale_color_manual(values = c('#004586FF','#579D1CFF','#FF420EFF'))+
  scale_fill_manual(values = c('#004586FF','#579D1CFF','#FF420EFF'))+
  scale_x_continuous(breaks = seq(0, 29, by = 1))+
  theme(axis.text = element_text(color="black"))
#theme(legend.position="none")
surcical1
ggsave('./N2/sur_all_n2.pdf',surcical1,width = 7,height = 3)


#log rank p value
#pairwise p values
library(survminer)
pair<-pairwise_survdiff(Surv(time, status) ~ group,
                        data = data2)
pair.frame<-as.data.frame(pair$p.value)
pair.frame.star<-symnum(pair$p.value, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
                        symbols = c("****", "***", "**", "*", "+", " "),
                        abbr.colnames = FALSE, na = "")

#            HQbiome HQbiome-1
#  HQbiome-1 0.47    -        
#  HQbiome-2 6.4e-05 1.2e-05

#====================================================================================================
#Fig. 6c and Fig. 7c  gut microbiota and lawn composition grown on HQbiome-3
#HQBiome-3 F23
physeq <- qza_to_phyloseq(
  features = "input/F23-HQbiome-3_qza/table_deblur.qza", # ASV/OTU table
  #tree = "tree_rooted.qza", # Phylogenetic tree
  taxonomy = "input/F23-HQbiome-3_qza/tax.qza", # Taxonomy file
  metadata = "input/F23-HQbiome-3_qza/meta.txt" # Sample metadata
)
physeq


physeq_all<-subset_samples(physeq_all, !(name %in% c('FL3.1','FL3.2','NL3.2','NL3.3','NL3.4',
                                                     'NL5.1','NL5.2','NL5.3',
                                                     'NL7.1','NL7.2','NL7.3','NL7.4',
                                                     'F1.1-all','F1.2-all', 'F5.4',
                                                     'N3.4','N5.1','N5.5','N7.2','N7.4','N7.5')))

#merge genus
physeq <- physeq %>%
  tax_glom(taxrank = "Genus", NArm = TRUE) # agglomerate on Genus level
tax_table(physeq) <- tax_table(physeq)[,1:6]
sample_data(physeq)
#ÅÅÐò
physeq@sam_data$group<-factor(physeq@sam_data$group,levels = c('Microbiota','Lawn','Inoculum'))

physeq@sam_data$name <- factor(physeq@sam_data$name, 
                               levels = c( "gut1"))

physeq_all@sam_data$name
#Add total read counts to each library¡¯s sample data
#sample_data(physeq)$total_reads <- sample_sums(physeq)

#otu<-as.data.frame(tax_table(physeq))
#Genus
library(paletteer)
physeq_other_Genus <- physeq %>%
  tax_glom(taxrank = "Genus") %>%    # agglomerate on Class level
  transform_sample_counts(function(x) {100 * x/sum(x)}) %>%  # transform to relative abundance
  psmelt() %>%  # convert from phyloseq object to a long data frame
  mutate(Taxa_Order = Genus)# %>%
#mutate(Taxa_Order = replace(Taxa_Order, Abundance < 1, "<1%"))


#colourCount = length(unique(physeq_other_Genus$Taxa_Order))  # obtain number of colors needed
#getPalette = colorRampPalette(brewer.pal(9, "Accent")) # create palette with Set1 ,3

#-----------------------------------------------------------------------------------------
#microbiota and lawn
#Fig. 6c
F23_phyloseq<-subset_samples(physeq,physeq@sam_data$group%in%c('Microbiota','Lawn'))
#merge samples
F23_phyloseq_merge<-merge_samples(F23_phyloseq,'merge')
#F23  plot
F23_physeq_other_Genus <- F23_phyloseq_merge %>%
  tax_glom(taxrank = "Genus") %>%    # agglomerate on Class level
  transform_sample_counts(function(x) {100 * x/sum(x)}) %>%  # transform to relative abundance
  psmelt() %>%  # convert from phyloseq object to a long data frame
  mutate(Taxa_Order = Genus) %>%
  mutate(Taxa_Order = replace(Taxa_Order, Abundance < 1, "<1%"))


#colourCount = length(unique(physeq_other_Genus$Taxa_Order))  # obtain number of colors needed
#getPalette = colorRampPalette(brewer.pal(9, "Accent")) # create palette with Set1 ,3

F23_Genus_BP <- 
  ggplot(data=F23_physeq_other_Genus, 
         aes(x=Sample, # order by date for each library season
             y=Abundance, 
             fill=fct_reorder(Taxa_Order, Abundance, .desc = TRUE))) + # order taxa within bars by abundance
  geom_bar(stat="identity", 
           position = position_stack(reverse = TRUE), # order taxa in reverse
           width = 0.9) + 
  #geom_col(width = 1)+
  scale_fill_manual(values = c('Vibrio'='#60BD68FF','Paracoccus'='#006064FF','Parasphingorhabdus'='#F17CB0FF',
                               'Pseudomonas'='#C6242DFF','Rheinheimera'='#004949FF','Enterovibrio'='#DECF3FFF',
                               'Shewanella'='#3A82E4FF','Neptunomonas'='#66E3D9FF','<1%'='#999999FF')
                    #paletteer_d("ggsci::default_igv",50,direction=-1)
                    , # assign legend colors
                    guide = guide_legend(reverse = TRUE)) + # make legend match order of taxa in barplot
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), # Y axis by 10% increments
                     expand = c(0,0)) + # remove whitespace in plot
  facet_grid(~group, # draw facets by date
             scales = "free_x", 
             space = "free_x",
             drop = F
  ) +
  theme(legend.direction="vertical", 
        axis.text.x = element_text(color = 'black',angle = 45, hjust = 1),
        #element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y = element_text(#face = "bold", 
          size = 10,color = 'black'),
        strip.text = element_text(face = "bold",#, size = 15
        ),
        strip.background = element_rect(color = "white", fill = "#DDEC7EFF"
        ),
        legend.title = element_text( face = "bold"),
        plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
        panel.spacing = unit(0, "lines"),
        legend.text = element_text(face="italic"),
        panel.background = element_rect(fill = 'white', color = 'white'),  
        panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
  labs(  face = "bold",
         y = "Relative Abundance %", 
         fill = "Genus")+
  guides(fill = guide_legend(ncol = 1))
F23_Genus_BP


ggsave('HQbiome-3.pdf',F23_Genus_BP,width = 3,height = 4)

#Fig. 20c
#Inoculum
ino_phyloseq<-subset_samples(physeq,physeq@sam_data$group%in%c('Inoculum'))

#Inoculum  plot
ino_physeq_other_Genus <- ino_phyloseq %>%
  tax_glom(taxrank = "Genus") %>%    # agglomerate on Class level
  transform_sample_counts(function(x) {100 * x/sum(x)}) %>%  # transform to relative abundance
  psmelt() %>%  # convert from phyloseq object to a long data frame
  mutate(Taxa_Order = Genus)# %>%
#mutate(Taxa_Order = replace(Taxa_Order, Abundance < 1, "<1%"))


#colourCount = length(unique(physeq_other_Genus$Taxa_Order))  # obtain number of colors needed
#getPalette = colorRampPalette(brewer.pal(9, "Accent")) # create palette with Set1 ,3

ino_Genus_BP <- 
  ggplot(data=ino_physeq_other_Genus, 
         aes(x=Sample, # order by date for each library season
             y=Abundance, 
             fill=fct_reorder(Taxa_Order, Abundance, .desc = TRUE))) + # order taxa within bars by abundance
  geom_bar(stat="identity", 
           position = position_stack(reverse = TRUE), # order taxa in reverse
           width = 0.9) + 
  #geom_col(width = 1)+
  scale_fill_manual(values = c('Vibrio'='#60BD68FF', 'Shewanella'='#3A82E4FF','Psychrobacter'='#8D4B08FF',
                               'Rheinheimera'='#004949FF','Parasphingorhabdus'='#F17CB0FF','Aliivibrio'='#244579FF',
                               'Pseudomonas'='#C6242DFF', 'Enterovibrio'='#DECF3FFF','Lacinutrix'='#F17CB0FF',
                               'Sulfitobacter'='#009E73FF','Pseudoalteromonas'='#D25D38FF','Marinobacter'='#66E3D9FF',
                               'Thalassotalea'='#00FF7FFF','HQ336491_g'='#FCA3B7FF','Neptunomonas'='#66E3D9FF',
                               'Litoreibacter'='#98FB98FF','Saccharospirillum'='#D070B9FF','Yoonia'='#66CDAAFF',
                               'Bacillus'='#FF9800FF','Alteromonas'='#F49538FF','Marinagarivorans'='#7D5329FF',
                               'Postechiella'='#F15854FF','Marinomonas'='#8A6842FF','Salinimonas'='#E7D202FF',
                               'Staphylococcus'='#BF616AFF','Alkalihalobacillus'='#F0E442FF','Lactococcus'='#E9A820FF',
                               'Priestia'='#8FA87AFF','Thaumasiovibrio'='#5E81ACFF','Cytobacillus'='#FFFF6DFF',
                               'Galactobacter'='#3F51B5FF','Planococcus'='#24FF24FF','Paracoccus'='#006064FF',
                               'Halobacillus'='#24FF24FF','Metabacillus'='#CC79A7FF','Niallia'='#DB6D00FF',
                               'Rhodococcus'='#924900FF','Paenisporosarcina'='#00FF7FFF','Colwellia'='#920000FF',
                               'Microbacterium'='#B6DBFFFF','Loktanella'='#6DB6FFFF','Arenivirga'='#B66DFFFF'
  )
  #paletteer_d("ggsci::default_igv",50,direction=-1)
  , # assign legend colors
  guide = guide_legend(reverse = TRUE)) + # make legend match order of taxa in barplot
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), # Y axis by 10% increments
                     expand = c(0,0)) + # remove whitespace in plot
  facet_grid(~group, # draw facets by date
             scales = "free_x", 
             space = "free_x",
             drop = F
  ) + 
  theme(legend.direction="vertical", 
        axis.text.x = element_text(color = 'black',angle = 45, hjust = 1),
        #element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y = element_text(#face = "bold", 
          size = 10,color = 'black'),
        strip.text = element_text(face = "bold",#, size = 15
        ),
        strip.background = element_rect(color = "white", fill = "#DDEC7EFF"
        ),
        legend.title = element_text( face = "bold"),
        plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
        panel.spacing = unit(0, "lines"),
        legend.text = element_text(face="italic"),
        panel.background = element_rect(fill = 'white', color = 'white'),  
        panel.border = element_rect(colour = "black", fill=NA, size=1.5)
  ) +
  labs(  face = "bold",
         y = "Relative Abundance %", 
         fill = "Genus")+
  guides(fill = guide_legend(ncol = 1))
ino_Genus_BP
ggsave('inoculum.pdf',ino_Genus_BP,width = 3,height = 4)

#========================================================================================================
#Fig. 7c
#HQBiome-3 N2
physeq <- qza_to_phyloseq(
  features = "input/N2-HQbiome-3_qza/table_deblur.qza", # ASV/OTU table
  #tree = "tree_rooted.qza", # Phylogenetic tree
  taxonomy = "input/N2-HQbiome-3_qza/tax.qza", # Taxonomy file
  metadata = "input/N2-HQbiome-3_qza/meta.txt" # Sample metadata
)
physeq

#merge genus
physeq <- physeq %>%
  tax_glom(taxrank = "Genus", NArm = TRUE) # agglomerate on Genus level
tax_table(physeq) <- tax_table(physeq)[,1:6]
sample_data(physeq)

physeq@sam_data$group<-factor(physeq@sam_data$group,levels = c('Microbiota','Lawn'))

#-----------------------------------------------------------------------------------------
#microbiota  lawn
N2_phyloseq<-subset_samples(physeq,physeq@sam_data$group%in%c('Microbiota','Lawn'))
#merge samples
N2_phyloseq_merge<-merge_samples(N2_phyloseq,'merge')
#F23  plot
N2_physeq_other_Genus <- N2_phyloseq_merge %>%
  tax_glom(taxrank = "Genus") %>%    # agglomerate on Class level
  transform_sample_counts(function(x) {100 * x/sum(x)}) %>%  # transform to relative abundance
  psmelt() %>%  # convert from phyloseq object to a long data frame
  mutate(Taxa_Order = Genus) %>%
  mutate(Taxa_Order = replace(Taxa_Order, Abundance < 1, "<1%"))


#colourCount = length(unique(physeq_other_Genus$Taxa_Order))  # obtain number of colors needed
#getPalette = colorRampPalette(brewer.pal(9, "Accent")) # create palette with Set1 ,3

N2_Genus_BP <- 
  ggplot(data=N2_physeq_other_Genus, 
         aes(x=Sample, # order by date for each library season
             y=Abundance, 
             fill=fct_reorder(Taxa_Order, Abundance, .desc = TRUE))) + # order taxa within bars by abundance
  geom_bar(stat="identity", 
           position = position_stack(reverse = TRUE), # order taxa in reverse
           width = 0.9) + 
  #geom_col(width = 1)+
  scale_fill_manual(values = c('Vibrio'='#60BD68FF','Paracoccus'='#006064FF','Parasphingorhabdus'='#F17CB0FF',
                               'Pseudomonas'='#C6242DFF','Rheinheimera'='#004949FF','Enterovibrio'='#DECF3FFF',
                               'Shewanella'='#3A82E4FF','Neptunomonas'='#66E3D9FF','<1%'='#999999FF')
                    #paletteer_d("ggsci::default_igv",50,direction=-1)
                    , # assign legend colors
                    guide = guide_legend(reverse = TRUE)) + # make legend match order of taxa in barplot
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4), # Y axis by 10% increments
                     expand = c(0,0)) + # remove whitespace in plot
  facet_grid(~group, # draw facets by date
             scales = "free_x", 
             space = "free_x",
             drop = F
  ) + 
  theme(legend.direction="vertical", 
        axis.text.x = element_text(color = 'black',angle = 45, hjust = 1),
        #element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y = element_text(#face = "bold", 
          size = 10,color = 'black'),
        strip.text = element_text(face = "bold",#, size = 15
        ),
        strip.background = element_rect(color = "white", fill = "#DDEC7EFF"
        ),
        legend.title = element_text( face = "bold"),
        plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
        panel.spacing = unit(0, "lines"),
        legend.text = element_text(face="italic"),
        panel.background = element_rect(fill = 'white', color = 'white'),  
        panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
  labs(  face = "bold",
         y = "Relative Abundance %", 
         fill = "Genus")+
  guides(fill = guide_legend(ncol = 1))
N2_Genus_BP


#Fig. S21h
##########################################################################################################
#differential abundance analysis
library(DESeq2)
library(xlsx)
physeq_all
#F23
F23_phyloseq<-subset_samples(physeq,physeq@sam_data$group%in%c('Microbiota','Lawn'))
#merge samples
#F23_phyloseq_merge<-merge_samples(F23_phyloseq,'merge')

##======deseq2   gut and lawn
#day 3
F23_phyloseq_day3<-subset_samples(F23_phyloseq,F23_phyloseq@sam_data$merge%in%c('F23-gut-day3','F23-lawn-day3'))
head(sample_data(F23_phyloseq_day3)$group, n=10)
diagdds = phyloseq_to_deseq2(F23_phyloseq_day3, ~ group)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #lawn.F23 vs Gut.F23
sigtab
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(F23_phyloseq_day3)[rownames(sigtab), ], "matrix"))
head(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}


# remove records with NA at Genus
sigtab2 = subset(sigtab, !is.na(Genus))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
d3.lawn.vs.gut<-ggplot(sigtab2, aes(x=Genus, y=log2FoldChange#, color=Phylum
)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
d3.lawn.vs.gut
ggsave('differential abundance/lawn.vs.gut_F23_day3.pdf',d3.lawn.vs.gut,width = 5,height = 3)
write.xlsx(sigtab2,'differential abundance/sigtab_lawn.vs.gut_F23_day3.xlsx')

#heatmap
library(pheatmap)
library(ggpubr)
library(readxl)

#gut vs lawn
dt <- read_excel("./input/F23_HQbiome-3_deseq2.xlsx", sheet = 1)
#remove first col
dt1<-dt[,-c(1)]
rownames(dt1)=dt$`Genus`

#plot  
bk <- c(seq(-12,-0.1,by=0.1),seq(0,12,by=0.1))
pheatmap(dt1, cluster_cols = F, cluster_rows = F, 
         #clustering_distance_cols = "euclidean",clustering_method='average',
         color= c(colorRampPalette(colors = c("#510051FF","white"))(length(bk)/2),
                  colorRampPalette(colors = c("white","#005100FF"))(length(bk)/2)),
         legend_breaks=seq(-10,10,2),
         breaks=bk,
         fontsize = 5,fontsize_row = 3, fontsize_col = 3.2,angle_col = 45,
         scale="none",
         border_color= 'grey',
         border=F,
         annotation_names_row=F,annotation_legend=T,treeheight_col=10,
         cellwidth = 8, cellheight = 6,
         #filename='differential abundance/N2_gut_vs_lawn_DDA.pdf',width=3,height=3
)

#N2- HQbiome-3  no difference
N2_phyloseq<-subset_samples(physeq,physeq@sam_data$group%in%c('Microbiota','Lawn'))
#merge samples
#F23_phyloseq_merge<-merge_samples(F23_phyloseq,'merge')

##======deseq2   gut and lawn
#day 3
N2_phyloseq_day3<-subset_samples(N2_phyloseq,N2_phyloseq@sam_data$merge%in%c('N2-gut-day5','N2-lawn-day5'))
head(sample_data(N2_phyloseq_day3)$group, n=10)
diagdds = phyloseq_to_deseq2(N2_phyloseq_day3, ~ group)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
diagdds_geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans=diagdds_geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]  #lawn.F23 vs Gut.F23
sigtab

#NS




