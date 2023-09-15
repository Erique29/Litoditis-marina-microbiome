#Fig. 1a
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

#load data
physeq <- qza_to_phyloseq(
  features = "./input/table_filter.qza", # ASV/OTU table
  tree = "./input/tree_rooted.qza", # Phylogenetic tree
  taxonomy = "./input/ez_vsearch_tax.qza", # Taxonomy file
  metadata = "./input/meta.txt" # Sample metadata
)

#Add total read counts to each library¡¯s sample data
sample_data(physeq)$total_reads <- sample_sums(physeq)
#Convert ¡°date¡± from character to an R factor
physeq@sam_data$Date <- factor(physeq@sam_data$Date, 
                               levels = c("Jun-20", "Aug-20", "Oct-20","Dec-20","Jan-21","Feb-21","Apr-21","Jun-21"))

#remove 'Unassigned' in phylum
physeq<- subset_taxa(physeq, Kingdom != "Unassigned")    #5783 OTUs

#Species labels are unreliable and many taxa are labeled as ¡°uncultured¡± at species level.
physeq <- physeq %>%
  tax_glom(taxrank = "Genus", NArm = TRUE) # agglomerate on Genus level
tax_table(physeq) <- tax_table(physeq)[,1:6]
sample_sums(physeq)

#Bar Charts at phylum and class levels
physeq2 = physeq
plot.tax = as.data.frame(tax_table(physeq2))
plot.tax = data.frame(lapply(plot.tax, as.character), stringsAsFactors = F)
plot.tax$Phylum[plot.tax$Phylum=="Proteobacteria"] = plot.tax$Class[plot.tax$Phylum=="Proteobacteria"] 
plot.tax$Phylum<- gsub(" ", "", plot.tax$Phylum, fixed = TRUE)
plot.tax[] = lapply(plot.tax, factor)
plot.tax.table = tax_table(plot.tax)
rownames(plot.tax.table) = rownames(tax_table(physeq))
colnames(plot.tax.table) = colnames(tax_table(physeq))
tax_table(physeq2) = plot.tax.table
#plot
physeq_other_phylum2 <- physeq2 %>%
  tax_glom(taxrank = "Phylum") %>%    # agglomerate on Class level
  transform_sample_counts(function(x) {100 * x/sum(x)}) %>%  # transform to relative abundance
  psmelt() %>%  # convert from phyloseq object to a long data frame
  mutate(Taxa_Order = Phylum) %>%
  mutate(Taxa_Order = replace(Taxa_Order, Abundance < 2, "<2%")) %>%
  mutate(Taxa_Order = factor(Taxa_Order)) %>%
  unite("year-month", year:month, sep = "-")

colourCount = length(unique(physeq_other_phylum2$Taxa_Order))  # obtain number of colors needed
getPalette = colorRampPalette(brewer.pal(9, "Set1")) # create palette with Set1
cols<-c('Gammaproteobacteria'='#90EE90',"Bacteroidetes" = "#4F94CD","Firmicutes"="#EE7600",'Deltaproteobacteria'='#71C671',
        "Actinobacteria" = "#EEE685","Cyanobacteria" = "#99991EFF",'Alphaproteobacteria'='#00CD00',
        'Betaproteobacteria'='#43CD80',"Acidobacteria" = "#6699FFFF",'Verrucomicrobia'='#8B3E2F',
        'Peregrinibacteria'='#8AC1D4FF','Latescibacteria_WS3'='#DA2E20FF','Planctomycetes'='#99CCFFFF',
        "Gemmatimonadetes"="#FFFF00FF",'Rhodothermaeota'='#CCCC99FF','Saccharibacteria_TM7'='#666600FF',
        '<2%'='#919191')


phylum_BP2 <- 
  ggplot(data=physeq_other_phylum2, 
         aes(x=Sample, # order by date for each library season
             y=Abundance, 
             fill=fct_reorder(Taxa_Order, Abundance, .desc = TRUE))) + # order taxa within bars by abundance
  geom_bar(stat="identity", 
           position = position_stack(reverse = TRUE), # order taxa in reverse
           width = 1) + 
  scale_fill_manual(values = cols, # assign legend colors
                    guide = guide_legend(reverse = TRUE)
  ) + # make legend match order of taxa in barplot
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), # Y axis by 10% increments
                     expand = c(0,0)) + # remove whitespace in plot
  facet_grid(~Date, # draw facets by date
             scales = "free_x", 
             space = "free_x",switch="x") + 
  theme(legend.direction="vertical", 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        #axis.text.y = element_text(face = "bold", size = 10#, family = "serif"
        #                           ),
        strip.text = element_text( size = 10,color='black'#, family = "serif"
        ),
        #legend.text = element_text(family = "serif"), 
        legend.title = element_text(face = "bold", #family = "serif"
                                    size=10),
        plot.margin = unit(c(0.1,0,0.1,0.1), "cm"),
        panel.spacing = unit(0, "lines")) +
  theme(legend.key.size = unit(0.1, "inches"))+
  labs(y = "Relative Frequency (%)", 
       fill = "Phylum/Class")+
  theme(axis.text = element_text(color = 'black'))+
  theme(
    strip.background = element_rect(
      color = "white", fill = "#EEE8CD"),
    panel.grid = element_blank(),
    panel.border = element_blank())
phylum_BP2
ggsave('a.pdf',phylum_BP2,width = 10, height = 4)

#Fig. 1b was created in Graphpad using data 'top 15 genus' in ./input/Fig. 1b.xlsx

#Fig. 1c was created using iTol (https://itol.embl.de/)  based on data in ./input/iTols_CFU/*

#Fig. 1d
library(reshape2)
cp<- read_excel("input/cfu_vs_otu.xlsx")
dt.cp<-melt(cp,id.vars='Phylum/Class')
dt.cp$`Phylum/Class`<-factor(dt.cp$`Phylum/Class`,
                             levels =c('Acidobacteria','Actinobacteria','Bacteroidetes',
                                       'Chlorobi','Cyanobacteria','Deinococcus-Thermus',
                                       'Firmicutes','Gemmatimonadetes','Planctomycetes',
                                       'Alphaproteobacteria','Betaproteobacteria','Deltaproteobacteria',
                                       'Gammaproteobacteria','Epsilonproteobacteria','Rhodothermaeota',
                                       'Verrucomicrobia'))
cols<-c("Acidobacteria" = "#6699FFFF", "Actinobacteria" = "#EEE685", "Bacteroidetes" = "#4F94CD", 
        "Chlorobi" = "#CC33FFFF", "Cyanobacteria" = "#99991EFF",
        "Deinococcus-Thermus"="#FFCCCCFF","Firmicutes"="#EE7600","Gemmatimonadetes"="#FFFF00FF",
        'Planctomycetes'='#99CCFFFF','Alphaproteobacteria'='#00CD00','Betaproteobacteria'='#43CD80',
        'Deltaproteobacteria'='#71C671','Gammaproteobacteria'='#90EE90','Epsilonproteobacteria'='#4EEE94',
        'Rhodothermaeota'='#CCCC99FF','Verrucomicrobia'='#8B3E2F')

p4<-ggplot(dt.cp,aes(variable,value,fill=`Phylum/Class`)) + 
  geom_col(width = 0.8, size = 0.6) +
  #geom_hline(yintercept = 0, color = "red1", size = 0.7)+
  theme_classic()+xlab(NULL)+ylab('Relative abundance (%)')+
  theme(axis.text = element_text(color = 'black'))+
  scale_fill_manual(values = cols)+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.key.size = unit(0.1, "inches"),
        legend.title = element_text(face = "bold", size=9),
        axis.text.x = element_text(angle = 45, hjust = 1))
p4


#Fig. S1
library(maps) 
library(mapdata) 
library(ggplot2)
library(plyr)
library(rgdal)
library(magrittr)
library(dplyr)
library(ggmap)
library(maptools)
library(reshape2)
library(readxl)
library(ggplot2)
library(dplyr)
library(maps)
library(sf)
library(ggspatial)
library(ggthemes)

china <- map_data("china")
ggplot(data = china, aes(x = long, y = lat, group = group))   +
  geom_polygon(colour="grey")+
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.05, "in"), pad_y = unit(0.05, "in"),
                         style = north_arrow_nautical)+
  theme_bw()+ylim(15,55)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                               panel.background = element_rect(fill = 'white', colour = 'white'), 
                               axis.line = element_line(colour = "white"), legend.position="none",
                               axis.ticks=element_blank(), axis.text.x=element_blank(),
                               axis.text.y=element_blank())

china<-map("china", col = "black", ylim = c(18, 54), panel.first = grid())

meta <- read_xlsx("input/sampledata.xlsx")
meta$Date <- factor(meta$Date,
                    levels = c('Aug-19','Oct-19','Dec-19',
                               "Jun-20", "Aug-20", "Oct-20","Dec-20","Jan-21","Feb-21",
                               "Apr-21","Jun-21",'Aug-21','Oct-21','Dec-21','Feb-22','Apr-22'))
qd<-get_stamenmap(bbox = c(left = 120.325, bottom = 36.044, right = 120.342, top = 36.060),
                  zoom = 16, maptype = "watercolor", crop = T, color = 'bw')	
ggmap(qd)

#show_point_shapes()
map<-ggmap(qd)+ geom_point(data = meta, aes(x = Lon, y = Lat), size = 2, color = 'black',shape=8)+
  #scale_shape_manual(values = c(2,3,4,5,6,7,8,10))+
  theme(legend.text = element_text(face = "bold",size=12), legend.title = element_text( face = "bold",size=14),
        legend.background = element_rect(color = "black", linetype = "solid", size = 0.4))+
  labs(x = 'Longitude', y = 'Latitude')+theme(axis.title = element_text( face = "bold",size = 20),
                                              axis.text = element_text( face = "bold"))+
  ggtitle('Sampling sites')+theme(plot.title = element_text(hjust = 0.5,  face = "bold",size=20))+
  annotate('text', x = 120.332, y = 36.051, label='Huiquan Bay', size = 12)+
  annotation_north_arrow(location = "tl", which_north = "true", 
                         pad_x = unit(0.05, "in"), pad_y = unit(0.05, "in"),
                         style = north_arrow_nautical)
map  
ggsave('map.pdf',map,width=8,height = 8)
#merge two plots in AI.


#Fig. S2 was created using iTol (https://itol.embl.de/)  based on data in ./input/iTols_OTU/* 


#Fig. S3
#Alpha Diversity
set.seed(1980)
physeq_rarefy <- rarefy_even_depth(physeq, sample.size = min(sample_sums(physeq)),
                                   rngseed = FALSE , replace = FALSE, trimOTUs = TRUE, verbose = FALSE)

# Calculate alpha diversity metrics
physeq_alpha_div <- microbiome::alpha(physeq_rarefy, index = "all")
# separate out the metadata from physeq
physeq_meta <- meta(physeq_rarefy)
#add the sample names for merging
physeq_meta$sample_name <- rownames(physeq_meta)
# add the sample names to diversity results
physeq_alpha_div$sample_name <- rownames(physeq_alpha_div)
# merge
physeq_alpha_div <- merge(physeq_alpha_div,physeq_meta, by = "sample_name")

#write.csv(physeq_alpha_div,'physeq_alpha_div.csv')

#plot_richness(physeq_7_rarefy, measures = c("Observed", "Chao1", "Shannon", "InvSimpson"), 
#              nrow = 4, color = "Date") +
#  theme(axis.text.x = element_text(angle = 60, face = "plain", hjust = .5, vjust = .5, size = 10,)) +
#  scale_x_discrete(labels = ) + 
#  labs(color = "Date") +
#  xlab("")

#significant analysis
#Observed Taxa
# perform Wilcoxon test for boxplot
#¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤
#shapiro.test and bartlett.test
shapiro.test(physeq_alpha_div$observed)  #W = 0.94708, p-value = 0.002353
bartlett.test(observed ~ Date, data = physeq_alpha_div) #Bartlett's K-squared = 11.347, df = 7, p-value = 0.1242  
#¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤
#stat_test_obs <- physeq_alpha_div %>% 
#  wilcox_test(observed ~ Date, p.adjust.method = "fdr") %>% 
#  add_xy_position(x = "Date")
##add letters in plots
library(multcompView)
wilcox_obs <-  pairwise.wilcox.test(physeq_alpha_div$observed, physeq_alpha_div$Date, 
                                    p.adjust.method = "BH") 

library(rcompanion)
wilcox_obs1<-fullPTable(wilcox_obs$p.value)  #convert data to full matrix
obs_letter <- multcompLetters(wilcox_obs1,
                              compare="<",
                              threshold=0.05,
                              Letters=letters,
                              reversed = FALSE)
print(list(wilcox_pairwise = wilcox_obs, obs_letter))

#multcompLetters2(observed ~ Date, stat_test_obs$p.adj, physeq_alpha_div)


text<-physeq_alpha_div %>% 
  group_by(Date) %>% 
  slice_max(observed, n = 1, with_ties=FALSE)   
text2<-text[,c(24,2)]
#text2<-text2[order(-text2$observed),]
text2$labels<-obs_letter[['Letters']]


##          Jun-20 Aug-20 Oct-20 Dec-20 Jan-21 Feb-21 Apr-21 Jun-21 
##          "a"   "bc"   "ad"   "be"    "f"    "e"  "bce"   "cd" 


#cld(object = stat_test_obs$group1,
#    Letters = letters)
# create boxplot for graphing later
alpha_obs <- ggboxplot(physeq_alpha_div, 
                       x = "Date", 
                       y = "observed",
                       color = "Date", 
                       palette = pal_aaas("default")(10),
                       legend = "none",
                       add = "jitter") + 
  stat_compare_means(label.x.npc= "left", label.y.npc = "bottom",size=2) + # add Kruskal-Wallis test to boxplot
  #stat_pvalue_manual(stat_test_obs, label = "p.adj.signif", step.increase = 0, hide.ns = TRUE, 
  #                   tip.length = 0, label.size = 5) + # add Wilcoxon test to boxplot,hide.ns
  ylab("Observed Taxa") +
  xlab("") +
  labs(color = "Date")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=6),
        axis.title.y=element_text(size=8))+
  geom_text(data = text2, aes(x = Date, y = observed, label = labels,color=Date,vjust =-1))
alpha_obs
#pdf output in 10*8   y.position = c(158,162,166,170,174,178, 182),



#Chao 1
#¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤
#shapiro.test and bartlett.test
shapiro.test(physeq_alpha_div$chao1)  #W = 0.9681, p-value = 0.04295
bartlett.test(chao1 ~ Date, data = physeq_alpha_div) #Bartlett's K-squared = 10.164, df = 7, p-value = 0.1795  
#¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤

stat_test_chao1 <- physeq_alpha_div %>% 
  wilcox_test(chao1 ~ Date, p.adjust.method = "fdr") %>% 
  add_xy_position(x = "Date")

alpha_Chao1 <- ggboxplot(physeq_alpha_div, 
                         x = "Date", 
                         y = "chao1",
                         color = "Date", 
                         palette = pal_aaas("default")(10),
                         legend = "none",
                         add = "jitter") + 
  stat_compare_means(label.x.npc= "left", label.y.npc = "bottom",size=2) +
  # stat_pvalue_manual(stat_test_chao1, label = "p.adj.signif", step.increase = 0, hide.ns = TRUE,
  #                    tip.length = 0, y.position = c(250,258,266,274,282,290,298),label.size = 5) +
  ylab("Chao1 Index") +
  xlab("") +
  labs(color = "Date")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=6),
        axis.title.y=element_text(size=8))
alpha_Chao1
#Shannon Diversity
#¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤
#shapiro.test and bartlett.test
shapiro.test(physeq_alpha_div$diversity_shannon)  #W = 0.91472, p-value = 5.402e-05
bartlett.test(diversity_shannon ~ Date, data = physeq_alpha_div) #Bartlett's K-squared = 20.395, df = 7, p-value = 0.004776  

#¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤

#stat_test_shannon <- physeq_alpha_div %>% 
#  wilcox_test(diversity_shannon ~ Date, p.adjust.method = "fdr") %>%
#  add_xy_position(x = "Date")
##add letters in plots
#library(multcompView)
wilcox_shannon <-  pairwise.wilcox.test(physeq_alpha_div$diversity_shannon, physeq_alpha_div$Date, 
                                        p.adjust.method = "BH") 

#library(rcompanion)
wilcox_shannon1<-fullPTable(wilcox_shannon$p.value)  #concert data to full matrix
shannon_letter <- multcompLetters(wilcox_shannon1,
                                  compare="<",
                                  threshold=0.05,
                                  Letters=letters,
                                  reversed = FALSE)
print(list(wilcox_pairwise = wilcox_shannon, shannon_letter))

#multcompLetters2(observed ~ Date, stat_test_obs$p.adj, physeq_alpha_div)


text_shannon<-physeq_alpha_div %>% 
  group_by(Date) %>% 
  slice_max(diversity_shannon, n = 1, with_ties=FALSE)  
text_shannon2<-text_shannon[,c(24,6)]  
#text2<-text2[order(-text2$observed),]
text_shannon2$labels<-shannon_letter[['Letters']]

alpha_shannon <- ggboxplot(physeq_alpha_div, 
                           x = "Date", 
                           y = "diversity_shannon",
                           color = "Date", 
                           palette = pal_aaas("default")(10),
                           legend = "none",
                           add = "jitter") + 
  stat_compare_means(label.x.npc= "left", label.y.npc = "bottom",size=2) +
  #stat_pvalue_manual(stat_test_shannon, label = "p.adj.signif", step.increase = 0, hide.ns = TRUE,
  #                   tip.length = 0, y.position = c(4.5,4.6,4.7,4.8,4.9),label.size = 5) +
  ylab("Shannon Diversity") +
  xlab("") +
  labs(color = "Date")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=6),
        axis.title.y=element_text(size=8))+
  geom_text(data = text_shannon2, aes(x = Date, y = diversity_shannon, label = labels,color=Date,vjust =-1))

alpha_shannon
#simpson
#¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤
#shapiro.test and bartlett.test
shapiro.test(physeq_alpha_div$diversity_inverse_simpson)  #W = 0.97487, p-value = 0.1161
bartlett.test(diversity_inverse_simpson ~ Date, data = physeq_alpha_div) #Bartlett's K-squared = 36.353, df = 7, p-value = 6.217e-06 
#¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤¡¤

stat_test_simpson <- physeq_alpha_div %>% 
  wilcox_test(diversity_inverse_simpson ~ Date, p.adjust.method = "fdr") %>% 
  add_xy_position(x = "Date")

alpha_simpson <- ggboxplot(physeq_alpha_div, 
                           x = "Date", 
                           y = "diversity_inverse_simpson",
                           color = "Date", 
                           palette = pal_aaas("default")(10),
                           legend = "none",
                           add = "jitter") + 
  stat_compare_means(label.x.npc= "left", label.y.npc = "bottom",size=2) +
  #stat_pvalue_manual(stat_test_simpson, label = "p.adj.signif", step.increase = 0, hide.ns = TRUE,
  #                   tip.length = 0,y.position = c(45,47,49,51),label.size = 5) +
  ylab("Inverse Simpson Diversity") +
  xlab("") +
  labs(color = "Date")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=6),
        axis.title.y=element_text(size=8))

alpha_simpson
#output significant results

library(writexl)
write_xlsx(as.data.frame(stat_test_obs),'obs.xlsx')
write_xlsx(as.data.frame(stat_test_chao1),'chao1.xlsx')
write_xlsx(as.data.frame(stat_test_shannon),'shannon.xlsx')
write_xlsx(as.data.frame(stat_test_simpson),'simpson.xlsx')

#merge plots
#bcde<-ggarrange(alpha_obs, alpha_Chao1, alpha_shannon, alpha_simpson,ncol=4,nrow=1,
#                common.legend = TRUE, legend = "none", align = "hv",labels=c('b','c','d','e'))
#bcde

bc<-ggarrange(alpha_obs, alpha_shannon, ncol=2,nrow=1,
              common.legend = TRUE, legend = "none", align = "hv",labels=c('b','c'))
bc
ggsave('bc.pdf',bc,width =8 ,height = 4)


# #Supplementary Fig. 4 has been removed
# #Occupancy¨Cabundance curves
# library(UpSetR)
# library(tidyverse)
# library(reshape2)
# library(vegan)
# 
# theme_set(theme_light())
# 
# tax<-readxl::read_excel("input/curves.xlsx", sheet = "OTUs")
# #tax_sub<-subset(tax,Phylum %in% c('Proteobacteria','Bacteroidetes','Cyanobacteria','Actinobacteria',
# #                                          'Verrucomicrobia','Acidobacteria','Planctomycetes',
# #                                          'Gemmatimonadetes','Firmicutes','Rhodothermaeota'))
# cols<-c('Proteobacteria'='#00CD00',"Bacteroidetes" = "#4F94CD","Cyanobacteria" = "#99991EFF",
#         "Actinobacteria" = "#EEE685",'Verrucomicrobia'='#8B3E2F', "Acidobacteria" = "#6699FFFF",
#         'Planctomycetes'='#99CCFFFF',"Gemmatimonadetes"="#FFFF00FF","Firmicutes"="#EE7600",
#         'Rhodothermaeota'='#CCCC99FF','Others'='#838B8B')
# otu<-tax[,-c(82:88)]%>%column_to_rownames('OTU')
# tax1<-tax[,c(1,82:88)]%>%column_to_rownames('OTU')
# 
# otu_PA <- 1*((otu>0)==1)     # presence-absence data
# otu_occ <- rowSums(otu_PA)/ncol(otu_PA)    # occupancy calculation
# otu_rel <- apply(decostand(otu, method="total", MARGIN=2),1, mean)     # relative abundance 
# occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') # combining occupancy and abundance 
# merge_occ_abun<-merge(occ_abun,tax1,by.x = 'otu',by.y='row.names')
# fill.marker<-merge(otu_PA, merge_occ_abun, by.x = 'row.names', by.y = 'otu')
# #write.csv(fill.marker,'fill.marker.csv')
# 
# plot.dt<-read.csv('./input/fill.marker.csv')
# plot.dt2<-plot.dt[,c(88,84:86)]
# #colnames(plot.dt2)<-c('Phylum','fill','otu_occ','otu_rel')
# plot.dt2$fill<-factor(plot.dt2$fill,levels=c('Shared among ¡Ý70 samples',
#                                              'Shared among 40-70 samples',
#                                              'Shared among <40 samples',
#                                              'Non-shared OTUs'))
# curve<-ggplot(data=plot.dt2, aes(x=log10(otu_rel), y=otu_occ, color=Phylum,shape=fill))+
#   geom_point(size = 2)+scale_color_manual(values=cols)+ 
#   scale_shape_manual(values = c(19,10,3,15))+
#   labs(x=paste('log(mean relative abundace per OTU)\n (n=',nrow(occ_abun),' OTUs)',sep=''), y=paste('Occupancy (n=',ncol(otu),')', sep=''), fill=NULL)+
#   geom_smooth(aes(x=log10(otu_rel), y=otu_occ,group = "" ),
#               #method = "loess",
#               #formula = "y ~ s(x) ",
#               se=T,
#               #stat = "identity",
#               #group='',
#               color="grey",
#               span = 0.9)+   
#   ylim(0, 1.)+
#   #mytheme = 
#   theme_classic() + theme(axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 8))+
#   theme(axis.title.y= element_text(size=12))+theme(axis.title.x = element_text(size = 12))+
#   theme(legend.title=element_blank(),legend.text=element_text(size=7,color = 'black'),
#         legend.position = c(0.2, 0.6))+
#   guides(color=guide_legend(ncol=2),
#          fill =guide_legend(ncol=2) ,
#          shape = guide_legend(ncol=1))
# curve
# 
# ggsave('curve.pdf',curve,width = 6,height = 5)
# 
# #-----------------draw plots at phylum
# library(microbiome)
# #top 10 phylums
# physeq_top10Phylum<- subset_taxa(physeq, Phylum== 'Proteobacteria' | Phylum=='Bacteroidetes' | 
#                                    Phylum=='Cyanobacteria' | Phylum=='Actinobacteria' | Phylum=='Verrucomicrobia' | Phylum==
#                                    'Acidobacteria' | Phylum=='Planctomycetes' | Phylum== 'Gemmatimonadetes' | Phylum==
#                                    'Firmicutes' | Phylum=='Rhodothermaeota')
# cols2<-c('Proteobacteria'='#00CD00',"Bacteroidetes" = "#4F94CD","Cyanobacteria" = "#99991EFF",
#          "Actinobacteria" = "#EEE685",'Verrucomicrobia'='#8B3E2F', "Acidobacteria" = "#6699FFFF",
#          'Planctomycetes'='#99CCFFFF',"Gemmatimonadetes"="#FFFF00FF","Firmicutes"="#EE7600",
#          'Rhodothermaeota'='#CCCC99FF')
# phylum_curve<-plot_taxa_prevalence(physeq_top10Phylum, "Phylum")+
#   scale_color_manual(values=cols2)+
#   labs(x='log(mean relative abundace per OTU)', y='Occupancy')+
#   geom_smooth(aes(x=abundance, y=prevalence,group = "" ),
#               #method = "loess",
#               #formula = 'y ~ s(x, bs = "cs")',
#               se=T,
#               #stat = "identity",
#               #group='',
#               color="grey",
#               span = 0.9)+
#   theme(axis.text.x = element_text(angle = 0))+
#   theme(legend.position = "none")
# phylum_curve
# ggsave('phylum_curve.pdf',phylum_curve,width = 7.5,height = 4)


