#packages
library(readxl)
library(ggpubr)
library(tidyverse)
library(ggplot2)
library(ggalt)
library(paletteer)

#k-means cluster
library(data.table)
library(dplyr)
library(RColorBrewer)
library(indicspecies)
library(superheat)
library(factoextra)
library(cluster)

#-----------------all collections---------------------------------------------------------------------
tax<-read_excel('input/cluster.xlsx',sheet = 1)
phyno_f23n2<-tax[,c(16,20)]
phyno_f23<-tax[,20] 
phyno_n2<-tax[,16]
rownames(phyno_f23n2)<-tax$Strain
rownames(phyno_f23)<-tax$Strain
rownames(phyno_n2)<-tax$Strain

#f23-------k-means
#k value
a<-fviz_nbclust(phyno_f23, kmeans, method = "wss", linecolor  = "black",nboot=1000)#+
#geom_vline(xintercept = 4, linetype = 2,colour = 'lightgrey')  #k=3
a
#fviz_nbclust(phyno, kmeans, method = "gap_stat") #k=5
b<-fviz_nbclust(phyno_f23, kmeans, method = "silhouette",
                linecolor  = "black",nboot=1000)#+
#geom_vline(xintercept = 4, linetype = 2,colour = 'lightgrey')# k=2
b

#gap_stat <- clusGap(phyno,
#                    FUN = kmeans,
#                   nstart = 100,
#                    K.max = 10,
#                    B = 100)
#fviz_gap_stat(gap_stat)   #k=6
# kmeans-clust
set.seed(666)
km <- kmeans(phyno_f23, centers = 3, nstart = 1000)
aggregate(phyno_f23, by=list(cluster=km$cluster), mean)
final_f23 <- cbind(phyno_f23, cluster = km$cluster)

#cluster        d5
#1       1  2.395787
#2       2 58.118657
#3       3 32.471555

write.csv(final_f23,'km_final_f23.csv')
#plot-cluster
plot(silhouette(km$cluster, daisy(phyno_f23)),  
     border=NA,
     col=c("#EE2617FF", "#F2A241FF", "#558934FF"),
     main = 'Silhouette width of each cluster',
     xlab='Silhouette width',
     #asp = 0.01
)

#N2  k-means
#k value
a<-fviz_nbclust(phyno_n2, kmeans, method = "wss", linecolor  = "black",nboot=1000)#+
#geom_vline(xintercept = 4, linetype = 2,colour = 'lightgrey')  #k=2
a
#fviz_nbclust(phyno, kmeans, method = "gap_stat") #k=5
b<-fviz_nbclust(phyno_f23n2, kmeans, method = "silhouette",
                linecolor  = "black",nboot=1000)+
  geom_vline(xintercept = 4, linetype = 2,colour = 'lightgrey')# k=2
b

#gap_stat <- clusGap(phyno,
#                    FUN = kmeans,
#                   nstart = 100,
#                    K.max = 10,
#                    B = 100)
#fviz_gap_stat(gap_stat)   #k=6
# kmeans-clust
set.seed(666)
km <- kmeans(phyno_n2, centers = 3, nstart = 1000)

aggregate(phyno_n2, by=list(cluster=km$cluster), mean)
#cluster        N2
#1       1 166.55128
#2       2  57.63669
#3       3  92.26013

final_n2 <- cbind(phyno_n2, cluster = km$cluster)

write.csv(final_n2,'km_final_n2.csv')
#plot-cluster
plot(silhouette(km$cluster, daisy(phyno_n2)),  
     border=NA,
     col=c("#EE2617FF", "#F2A241FF", "#558934FF"),
     main = 'Silhouette width of each cluster',
     xlab='Silhouette width',
     #asp = 0.01
)


#-----------------------HQbiome------------------------------------------------------------
tax<-read_excel('input/cluster.xlsx',sheet = 2)
phyno_f23n2<-tax[,c(4,9)]
phyno_f23<-tax[,9] 
phyno_n2<-tax[,5]
rownames(phyno_f23n2)<-tax$ID
rownames(phyno_f23)<-tax$ID
rownames(phyno_n2)<-tax$ID

#f23
#k value
fviz_nbclust(phyno_f23, kmeans, method = "wss", linecolor  = "black",nboot=1000)#+
#geom_vline(xintercept = 4, linetype = 2,colour = 'lightgrey')  #k=2

#fviz_nbclust(phyno, kmeans, method = "gap_stat") #k=5
fviz_nbclust(phyno_f23n2, kmeans, method = "silhouette",
             linecolor  = "black",nboot=1000)+
  geom_vline(xintercept = 4, linetype = 2,colour = 'lightgrey')# k=2

set.seed(666)
km <- kmeans(phyno_f23, centers = 3, nstart = 1000)

aggregate(phyno_f23, by=list(cluster=km$cluster), mean)
#cluster    mean_5
#1       1  1.075838
#2       2 60.555556
#3       3 35.476190
final_f23 <- cbind(phyno_f23, cluster = km$cluster)

write.csv(final_f23,'km_final_f23.csv')
#plot-cluster
plot(silhouette(km$cluster, daisy(phyno_n2)),  
     border=NA,
     col=c("#EE2617FF", "#F2A241FF", "#558934FF"),
     main = 'Silhouette width of each cluster',
     xlab='Silhouette width',
     #asp = 0.01
)

#N2
set.seed(666)
km_n2 <- kmeans(phyno_n2, centers = 4, nstart = 1000)

aggregate(phyno_n2, by=list(cluster=km_n2$cluster), mean)
final_n2 <- cbind(tax[,c(4,5)], cluster = km_n2$cluster)

write.csv(final_n2,'km_final_n2.csv')

# Manually organize the above data into the source data for Fig.2.
#Fig. 2a,b were created in graphpad using ./input/Fig. 2ab.xlsx

#================================================================================================
#Fig. 2c
#revised version
library(readxl)
library(ggpubr)
library(tidyverse)
library(ggplot2)
library(ggalt)
library(paletteer)

HQbiome <-read_excel('input/Fig. 2 subset_plot.xlsx', sheet='error_bar')
plotrix::std.error(c(1,2,3))


HQbiome %>% 
  rowwise() %>% 
  mutate(mean_value=mean(c(rep1_100,rep2_100,rep3_100)),
         std_error=plotrix::std.error(c(rep1_100,rep2_100,rep3_100))) %>% 
  ungroup() -> new.HQbiome

new.HQbiome$group<-factor(new.HQbiome$group,levels=c('Beneficial','Intermediate','Detrimental'))
new.HQbiome$strain<-factor(new.HQbiome$strain,levels=c('OP50','HQB-281','HQB-154','HQB-476',
                                                       'HQB-55', 'HQB-372', 'HQB-361',
                                                       'HQB-498','HQB-99','HQB-181'))

#plot
fig2c<- ggplot(data=new.HQbiome,aes(x=Day,y=mean_value, color=strain,linetype=group))+
  geom_errorbar(aes(ymin=mean_value-std_error,
                    ymax=mean_value+std_error),linetype="solid",
                width=0.1,linewidth=0.5)+
  geom_point(size=1)+
  #geom_line(linewidth =0.8)+
  geom_xspline(spline_shape = -0.5)+
  scale_linetype_manual(values=c("solid", "dashed",'dotted'))+
  geom_vline(xintercept=5, linewidth= 0.5,color="lightgrey",linetype = "dashed")+
  scale_color_manual(values = c('OP50'='black','HQB-281'='#E27069FF','HQB-154'='#B33941FF','HQB-476'='#7852A9FF',
                                'HQB-55'='#1175BBFF', 'HQB-372'='#2E8B57FF', 'HQB-361'='#DBA520FF',
                                'HQB-498'='#EF7215FF','HQB-99'='#4AC6AEFF','HQB-181'='#997950FF'
  ),
  #labels=c("Vector","mRHIM2","mRHIM1",expression(paste("mZ",alpha,"2")),"WT"),
  #breaks = c("A","F","E","D","B"),
  name=NULL)+
  scale_x_continuous(#limits = c(3,10),
    expand = expansion(mult=c(0,0)),
    breaks = seq(10)
  )+
  scale_y_continuous(limits = c(0,80),
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
fig2c
ggsave('fig2c.pdf', fig2c, width = 5,height = 2.2)


#old version
# HQbiome <-read_excel('input/Fig. 2 subset_plot.xlsx', sheet='HQbiome')%>%
#   #subset(Family %in% 'Vibrionaceae')%>%
#   select(ID, Group, Day3,Day4,Day5,Day6,Day7,Day8,Day9,Day10)%>%
#   gather(Day, Proportion, -c(ID, Group))%>%
#   subset(ID%in%c('OP50','HQB-281','HQB-154','HQB-476',  #beneficial
#                  'HQB-55', 'HQB-372', 'HQB-361',   #Intermediate
#                  'HQB-498','HQB-99','HQB-181'       #detrimental
#   ))
# HQbiome$Day<-factor(HQbiome$Day,levels = c('Day3','Day4','Day5','Day6','Day7','Day8','Day9','Day10'))  
# HQbiome$Group<-factor(HQbiome$Group,levels=c('Beneficial','Intermediate','Detrimental'))
# 
# color_V<-c('OP50'='black','HQB-281'='#E27069FF','HQB-154'='#B33941FF','HQB-476'='#7852A9FF',
#            'HQB-55'='#1175BBFF', 'HQB-372'='#2E8B57FF', 'HQB-361'='#DBA520FF',
#            'HQB-498'='#EF7215FF','HQB-99'='#4AC6AEFF','HQB-181'='#997950FF')
# #paletteer_d("basetheme::royal")
# #plot
# HQ<-ggplot(data=HQbiome,
#            aes(x=Day,y=Proportion,
#                group=ID,colour=ID,linetype=Group))+
#   #geom_point(size=2)+
#   scale_linetype_manual(values=c("solid", "dashed",'dotted'))+
#   geom_xspline(spline_shape = -0.5,size=1.5)+
#   labs(#x="Vibrionaceae", 
#     y="L4 proportion of F23 (%)")+
#   scale_colour_manual(values=color_V)+
#   geom_vline(xintercept='Day5', color="lightgrey",linetype = "dashed")+
#   theme_bw()+
#   theme(axis.text = element_text(color = 'black'),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.title.x = element_text(face="italic"),
#         legend.box = "horizontal",
#         panel.grid.major=element_line(colour=NA),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent",colour = NA),
#         panel.grid.minor = element_blank()   
#   )+
#   guides(fill=guide_legend(ncol=2))+
#   scale_y_continuous(expand = c(0.0, 0), limits=c(0,75),breaks=c(0,25,50,75))+  
#   scale_x_discrete(expand = c(0,0)) 
# 
# HQ
# ggsave('HQ.pdf', HQ, width = 7,height = 4)


#Fig. S4
library(readxl)
library(ggpubr)
library(tidyverse)

#F23 stack
growth<-read_excel('input/growth-all.xlsx', sheet=1)
cols2 <- c(Actinobacteria= '#EEE685',Bacteroidetes = "#4F94CD", Firmicutes= '#EE7600', 
           Alphaproteobacteria = "#00CD00",Gammaproteobacteria = "#90EE90")

group1<-subset(growth,group1%in%1)  #detrimental
group2<-subset(growth,group1%in%2)  #benefitial
group3<-subset(growth,group1%in%3)  #intermidiate

group1.2<-group1[,c('ID','Day3','Day4','Day5','Day6','Day7','Day8','Day9','Day10','Remains','Survival')]%>%
  mutate(ID=fct_reorder(ID, Survival)
  ) %>% 
  gather(day, value, -ID)%>%subset(day!='Survival')
group1.2$day<-factor(group1.2$day,levels =c('Remains','Day10','Day9','Day8','Day7','Day6','Day5','Day4','Day3'))

group2.2<-group2[,c('ID','Day3','Day4','Day5','Day6','Day7','Day8','Day9','Day10','Remains','Survival')]%>%
  mutate(ID=fct_reorder(ID, Survival)
  ) %>% 
  gather(day, value, -ID)%>%subset(day!='Survival')
group2.2$day<-factor(group2.2$day,levels =c('Remains','Day10','Day9','Day8','Day7','Day6','Day5','Day4','Day3'))

group3.2<-group3[,c('ID','Day3','Day4','Day5','Day6','Day7','Day8','Day9','Day10','Remains','Survival')]%>%
  mutate(ID=fct_reorder(ID, Survival)
  ) %>% 
  gather(day, value, -ID)%>%subset(day!='Survival')
group3.2$day<-factor(group3.2$day,levels =c('Remains','Day10','Day9','Day8','Day7','Day6','Day5','Day4','Day3'))

#color
color1<-group1[,c('ID','color')]
color2<-group2[,c('ID','color')]
color3<-group3[,c('ID','color')]

color1$ID<-factor(color1$ID,levels=(levels(group1.2$ID)))
color1<-color1[order(color1$ID),]
color2$ID<-factor(color2$ID,levels=(levels(group2.2$ID)))
color2<-color2[order(color2$ID),]
color3$ID<-factor(color3$ID,levels=(levels(group3.2$ID)))
color3<-color3[order(color3$ID),]
#plot
#group1
group1.p<-ggplot(group1.2,aes(ID, value))+
  geom_col(aes(fill = day), width = 1, color = "white", size = 0.6
  )+
  coord_flip(clip = "off")+
  scale_y_continuous(expand = c(0, 0), limits=c(0,100),"F23 developmental rate (%)")+
  scale_fill_manual(values = c("Day3" = "#008ECEFF", "Day4" = "#00A9E0FF", "Day5" = "#59C7EBFF", 
                               "Day6" = "#FC8D59FF", "Day7" = "#EF6548FF",
                               "Day8"="#D7301FFF","Day9"="#B30000FF","Day10"="#7F0000FF",
                               'Remains'='#696969'))+
  theme_classic()+
  xlab(NULL)+theme(#axis.text.y = element_blank(),
    axis.text.y = element_text(size = 5, colour=color1$color
    ),
    axis.text = element_text(color = 'black'))

group1.p
#ggsave('1.pdf',group1.p,width = 10,height = 12)
#group2
group2.p<-ggplot(group2.2,aes(ID, value))+
  geom_col(aes(fill = day), width = 1, color = "white", size = 0.6
  )+
  coord_flip(clip = "off")+
  scale_y_continuous(expand = c(0, 0), limits=c(0,100),"F23 developmental rate (%)")+
  scale_fill_manual(values = c("Day3" = "#008ECEFF", "Day4" = "#00A9E0FF", "Day5" = "#59C7EBFF", 
                               "Day6" = "#FC8D59FF", "Day7" = "#EF6548FF",
                               "Day8"="#D7301FFF","Day9"="#B30000FF","Day10"="#7F0000FF",
                               'Remains'='#696969'))+
  theme_classic()+
  xlab(NULL)+theme(#axis.text.y = element_blank(),
    axis.text.y = element_text(size = 5, colour=color2$color
    ),
    axis.text = element_text(color = 'black'))
group2.p

#group3
group3.p<-ggplot(group3.2,aes(ID, value))+
  geom_col(aes(fill = day), width = 1, color = "white", size = 0.6
  )+
  coord_flip(clip = "off")+
  scale_y_continuous(expand = c(0, 0), limits=c(0,100),"F23 developmental rate (%)")+
  scale_fill_manual(values = c("Day3" = "#008ECEFF", "Day4" = "#00A9E0FF", "Day5" = "#59C7EBFF", 
                               "Day6" = "#FC8D59FF", "Day7" = "#EF6548FF",
                               "Day8"="#D7301FFF","Day9"="#B30000FF","Day10"="#7F0000FF",
                               'Remains'='#696969'))+
  theme_classic()+
  xlab(NULL)+theme(#axis.text.y = element_blank(),
    axis.text.y = element_text(size = 5, colour=color3$color),
    axis.text = element_text(color = 'black'))
group3.p

group2_3<-ggarrange(group2.p,group3.p, common.legend = T, nrow=2,legend='none')
group2_3
all<-ggarrange(group2_3,group1.p,common.legend = T, nrow=1,legend='right')
all
ggsave('all.pdf',all,width = 15,height = 14)


#Fig. S5 was created in graphpad using data in ./input/Supplementary Fig. 5.xlsx

#Fig. S6
#Flavobacteriaceae
Flavobacteriaceae <-read_excel('./input/Fig. 2 subset_plot.xlsx', sheet='Flavobacteriaceae')%>%
  select(Strain,`Group-F23`,Day3,Day4,Day5,Day6,Day7,Day8,Day9,Day10)%>%
  gather(Day, Proportion, -c(Strain,`Group-F23`))%>%
  subset(Strain%in%c('OP50','HQB-134','HQB-156','HQB-178',  #beneficial
                     'HQB-109', 'HQB-113', 'HQB-310',   #Intermediate
                     'HQB-220','HQB-250','HQB-416'       #detrimental
  ))
Flavobacteriaceae$Day<-factor(Flavobacteriaceae$Day,levels = c('Day3','Day4','Day5','Day6','Day7','Day8','Day9','Day10'))  
Flavobacteriaceae$`Group-F23`<-factor(Flavobacteriaceae$`Group-F23`,levels=c('Beneficial','Intermediate','Detrimental'))
color_F<-c('OP50'='black','HQB-134'='#5E81ACFF','HQB-156'='#8FA87AFF','HQB-178'='#BF616AFF',
           'HQB-109'='#E7D202FF', 'HQB-113'='#7D5329FF', 'HQB-310'='#F49538FF',
           'HQB-220'='#66CDAAFF','HQB-250'='#D070B9FF','HQB-416'='#98FB98FF')
#paletteer_d("basetheme::brutal")
#plot
p.F<-ggplot(data=Flavobacteriaceae,
            aes(x=Day,y=Proportion,
                group=Strain,colour=Strain,linetype=`Group-F23`))+
  #geom_point(size=2)+
  scale_linetype_manual(values=c("solid", "dashed",'dotted'))+
  geom_xspline(spline_shape = -0.5,size=1.5)+
  labs(x="Flavobacteriaceae", y="L4 proportion of F23 (%)")+
  scale_colour_manual(values=color_F)+
  geom_vline(xintercept='Day5', color="lightgrey",linetype = "dashed")+
  theme_bw()+
  theme(axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(face="italic"),
        legend.box = "horizontal",
        panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank()   
  )+
  scale_y_continuous(expand = c(0.0, 0), limits=c(0,75),breaks=c(0,25,50,75))+  
  scale_x_discrete(expand = c(0,0)) 

p.F

#=============================================================================================
#Rhodobacteraceae
Rhodobacteraceae <-read_excel('./input/Fig. 2 subset_plot.xlsx', sheet='Rhodobacteraceae')%>%
  select(Strain,`Group-F23`,Day3,Day4,Day5,Day6,Day7,Day8,Day9,Day10)%>%
  gather(Day, Proportion, -c(Strain,`Group-F23`))%>%
  subset(Strain%in%c('OP50','HQB-154','HQB-177','HQB-267',  #beneficial
                     'HQB-161', 'HQB-446', 'HQB-462',   #Intermediate
                     'HQB-74','HQB-139','HQB-163'       #detrimental
  ))
Rhodobacteraceae$Day<-factor(Rhodobacteraceae$Day,levels = c('Day3','Day4','Day5','Day6','Day7','Day8','Day9','Day10'))  
Rhodobacteraceae$`Group-F23`<-factor(Rhodobacteraceae$`Group-F23`,levels=c('Beneficial','Intermediate','Detrimental'))

color_R<-c('OP50'='black','HQB-154'='#42BA90FF','HQB-177'='#BEE948FF','HQB-267'='#3A82E4FF',
           'HQB-161'='#D25D38FF', 'HQB-446'='#E9A820FF', 'HQB-462'='#9041BAFF',
           'HQB-74'='#795C32FF','HQB-139'='#EC5578FF','HQB-163'='#00B7EBFF')
#paletteer_d("basetheme::deepblue")
#plot
p.R<-ggplot(data=Rhodobacteraceae,
            aes(x=Day,y=Proportion,
                group=Strain,colour=Strain,linetype=`Group-F23`))+
  #geom_point(size=2)+
  scale_linetype_manual(values=c("solid", "dashed",'dotted'))+
  geom_xspline(spline_shape = -0.5,size=1.5)+
  labs(x="Rhodobacteraceae", y="L4 proportion of F23 (%)")+
  scale_colour_manual(values=color_R)+
  geom_vline(xintercept='Day5', color="lightgrey",linetype = "dashed")+
  theme_bw()+
  theme(axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(face="italic"),
        legend.box = "horizontal",
        panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank()   
  )+
  guides(fill=guide_legend(ncol=2))+
  scale_y_continuous(expand = c(0.0, 0), limits=c(0,75),breaks=c(0,25,50,75))+  
  scale_x_discrete(expand = c(0,0)) 
p.R

#=============================================================================================
#Shewanellaceae
Shewanellaceae <-read_excel('input/Fig. 2 subset_plot.xlsx', sheet='Shewanellaceae')%>%
  select(Strain,`Group-F23`,Day3,Day4,Day5,Day6,Day7,Day8,Day9,Day10)%>%
  gather(Day, Proportion, -c(Strain,`Group-F23`))%>%
  subset(Strain%in%c('OP50','HQB-76','HQB-86','HQB-162',  #beneficial
                     'HQB-5', 'HQB-62', 'HQB-453',   #Intermediate
                     'HQB-370','HQB-433','HQB-469'       #detrimental
  ))
Shewanellaceae$Day<-factor(Shewanellaceae$Day,levels = c('Day3','Day4','Day5','Day6','Day7','Day8','Day9','Day10'))  
Shewanellaceae$`Group-F23`<-factor(Shewanellaceae$`Group-F23`,levels=c('Beneficial','Intermediate','Detrimental'))

color_S<-c('OP50'='black','HQB-76'='#5DA5DAFF','HQB-86'='#FAA43AFF','HQB-162'='#60BD68FF',
           'HQB-5'='#F15854FF', 'HQB-62'='#B276B2FF', 'HQB-453'='#8D4B08FF',
           'HQB-370'='#DECF3FFF','HQB-433'='#F17CB0FF','HQB-469'='#66E3D9FF')
#paletteer_d("basetheme::void")
#plot
p.S<-ggplot(data=Shewanellaceae,
            aes(x=Day,y=Proportion,
                group=Strain,colour=Strain,linetype=`Group-F23`))+
  #geom_point(size=2)+
  scale_linetype_manual(values=c("solid", "dashed",'dotted'))+
  geom_xspline(spline_shape = -0.5,size=1.5)+
  labs(x="Shewanellaceae", y="L4 proportion of F23 (%)")+
  scale_colour_manual(values=color_S)+
  geom_vline(xintercept='Day5', color="lightgrey",linetype = "dashed")+
  theme_bw()+
  theme(axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(face="italic"),
        legend.box = "horizontal",
        panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank()   
  )+
  guides(fill=guide_legend(ncol=2))+
  scale_y_continuous(expand = c(0.0, 0), limits=c(0,75),breaks=c(0,25,50,75))+  
  scale_x_discrete(expand = c(0,0)) 
p.S

#================================================================================================
#Vibrionaceae
Vibrionaceae <-read_excel('input/Fig. 2 subset_plot.xlsx', sheet='Vibrionaceae')%>%
  select(Strain,`Group-F23`,Day3,Day4,Day5,Day6,Day7,Day8,Day9,Day10)%>%
  gather(Day, Proportion, -c(Strain,`Group-F23`))%>%
  subset(Strain%in%c('OP50','HQB-114','HQB-410','HQB-22',  #beneficial
                     'HQB-19', 'HQB-20', 'HQB-91',   #Intermediate
                     'HQB-21','HQB-39','HQB-150'       #detrimental
  ))
Vibrionaceae$Day<-factor(Vibrionaceae$Day,levels = c('Day3','Day4','Day5','Day6','Day7','Day8','Day9','Day10'))  
Vibrionaceae$`Group-F23`<-factor(Vibrionaceae$`Group-F23`,levels=c('Beneficial','Intermediate','Detrimental'))

color_V<-c('OP50'='black','HQB-114'='#E27069FF','HQB-410'='#B33941FF','HQB-22'='#7852A9FF',
           'HQB-19'='#1175BBFF', 'HQB-20'='#2E8B57FF', 'HQB-91'='#DBA520FF',
           'HQB-21'='#EF7215FF','HQB-39'='#4AC6AEFF','HQB-150'='#997950FF')
#paletteer_d("basetheme::royal")
#plot
p.V<-ggplot(data=Vibrionaceae,
            aes(x=Day,y=Proportion,
                group=Strain,colour=Strain,linetype=`Group-F23`))+
  #geom_point(size=2)+
  scale_linetype_manual(values=c("solid", "dashed",'dotted'))+
  geom_xspline(spline_shape = -0.5,size=1.5)+
  labs(x="Vibrionaceae", y="L4 proportion of F23 (%)")+
  scale_colour_manual(values=color_V)+
  geom_vline(xintercept='Day5', color="lightgrey",linetype = "dashed")+
  theme_bw()+
  theme(axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(face="italic"),
        legend.box = "horizontal",
        panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank()   
  )+
  guides(fill=guide_legend(ncol=2))+
  scale_y_continuous(expand = c(0.0, 0), limits=c(0,75),breaks=c(0,25,50,75))+  
  scale_x_discrete(expand = c(0,0)) 

p.V

Fig.2<-ggarrange(p.F,p.R,p.S,p.V)
Fig.2
ggsave('S6.pdf',Fig.2,width = 8.5,height = 5)


#Fig. S25 was created in graphpad using ./input/Supplementary Fig. 25.xlsx

#Fig. S7
library(tidyverse)
library(ggpubr)
library(ape)
library(phytools)
library(geiger)
library(data.table)
library(readxl)

# Load Data ---------------------------------------------------------------
all.strain<-read_excel('input/signal/signal.xlsx',sheet=3)
Tree = read.tree('input/signal/all_trim.fasta.raxml.bestTree') 
str(all.strain)
# check whehter the Treepheno data and the ProtPhenoDataframe names for the bacs are ok
name.check(Tree, all.strain, data.names = all.strain$Strain)
PhenoDataVec_HQ <- setNames(all.strain$log2N2, all.strain$Strain)

# only keep the ones in the tree for which phenotypic data and phylogenetic data is available and then midpoint.root the tree
Treepheno <- treedata(Tree, PhenoDataVec_HQ)
Treepheno$phy <- midpoint.root(Treepheno$phy)

#Pagel and Abouheif including all strains

# Subset strains based on a column split and if a split contains more than size entries, randomly select size entries from the split. 
# Returns list of x subsetted based on the selection above.
SubsetStrains <- function(df, x, split, size){
  # df is data frame containing x (Strain_ID) and the split (e.g. Family) as columns
  # size -> number of entries that are ok per split and to which number it will be subset
  SplitLists <- split(df[[x]], df[[split]])
  SplitLengths <- as.data.frame(lengths(SplitLists))
  ToSample <- row.names(SplitLengths)[SplitLengths>=size]
  UseAllofSplit <- df[[x]][!(df[[split]] %in% ToSample)]
  Selected_Sub <- c()
  for(Split in ToSample){
    SplitSample <- df[x][df[split] == Split]
    Selected <- sample(SplitSample, size, replace = FALSE)
    Selected_Sub <- c(Selected_Sub, Selected)
  }
  StrainSubset <- c(UseAllofSplit, Selected_Sub)
  return(StrainSubset)
  
}

Strains <- SubsetStrains(df = all.strain, x = "Strain", split = "family", size = 500)
phenodata_sub_N2 <- setNames(all.strain[["log10N2"]][all.strain[["Strain"]] %in% Strains], all.strain[["Strain"]][all.strain[["Strain"]] %in% Strains])
phenodata_sub_F23 <- setNames(all.strain[["d5"]][all.strain[["Strain"]] %in% Strains], all.strain[["Strain"]][all.strain[["Strain"]] %in% Strains])
# subset the tree and drop all tips that don't have phenotypic data
Treesub_N2 <- treedata(Treepheno$phy, phenodata_sub_N2, warning=F)
Treesub_F23 <- treedata(Treepheno$phy, phenodata_sub_F23, warning=F)

TreesubSorted_N2 <- Treesub_N2$data[match(Treesub_N2$phy$tip.label, row.names(Treesub_N2$data)),]
TreesubSorted_F23 <- Treesub_F23$data[match(Treesub_F23$phy$tip.label, row.names(Treesub_F23$data)),]


#N2-Pagel  lambda : 0.693686   P-value: 1.58612e-15
phylosig(Treepheno$phy, TreesubSorted_N2, method = "lambda", test = TRUE)
#F23-d5-Pagel  lambda : 0.410681  P-value: 3.15433e-13 
phylosig(Treepheno$phy, TreesubSorted_F23, method = "lambda", test = TRUE)

# Then Abouheif's Cmean
library(adephylo)
W1_N2 <-proxTips(Treesub_N2$phy,method="oriAbouheif") 
Abouheif_N2 <- abouheif.moran(TreesubSorted_N2,W1_N2)
Abouheif_N2
# N2 Abouheif's Cmean 0.3044989 p=0.001
W1_F23 <- proxTips(Treesub_F23$phy,method="oriAbouheif") 
Abouheif_F23 <- abouheif.moran(TreesubSorted_F23,W1_F23)
Abouheif_F23
#F23 Abouheif's Cmean 0.1918774 p=0.001

#value=1 would mean that all the observed variation is explained by phylogeny



############

# Subset strains based on a column split and if a split contains more than size entries, randomly select size entries from the split. 
# Returns list of x subsetted based on the selection above.
SubsetStrains <- function(df, x, split, size){
  # df is data frame containing x (Strain_ID) and the split (e.g. Family) as columns
  # size -> number of entries that are ok per split and to which number it will be subset
  SplitLists <- split(df[[x]], df[[split]])
  SplitLengths <- as.data.frame(lengths(SplitLists))
  ToSample <- row.names(SplitLengths)[SplitLengths>=size]
  UseAllofSplit <- df[[x]][!(df[[split]] %in% ToSample)]
  Selected_Sub <- c()
  for(Split in ToSample){
    SplitSample <- df[x][df[split] == Split]
    Selected <- sample(SplitSample, size, replace = FALSE)
    Selected_Sub <- c(Selected_Sub, Selected)
  }
  StrainSubset <- c(UseAllofSplit, Selected_Sub)
  return(StrainSubset)
  
}

# # Doing random subsets of the data and calculate phylogenetic signal using Pagel's lambda
#N2
PagelsLambda_subsetlist_persplit_N2 <- function( phenodata, phylotree, x, split, size){
  # Note phenodata is a dataframe containing a column x corresponding to the names in the tip.labels of phylotree, 
  # a column for the split (e.g. Family) and at least a column "N2.log", which contains the phenotypic data
  # check whehter the Treepheno data and the ProtPhenoDataframe names for the bacs are ok
  if(name.check(phylotree, phenodata, data.names= phenodata[[x]]) != "OK") warning("Not all taxa in phylogenetic tree and phenotypic data are equal or present in both data.")
  x <- x
  split <- split
  size <- size
  Strains <- SubsetStrains(phenodata, x, split, size)
  phenodata_sub <- setNames(phenodata[["N2"]][phenodata[[x]] %in% Strains], phenodata[[x]][phenodata[[x]] %in% Strains])
  # subset the tree and drop all tips that don't have phenotypic data
  Treesub <- treedata(phylotree, phenodata_sub, warning=F)
  
  # using Pagel's lambda. Note that we could also add se  = standard error -> could make sense
  Lambda = phylosig(Treesub$phy, phenodata_sub, method = "lambda", test = TRUE)
  return(Lambda[1:4])
}

# Function to do this for rep (number of replications) of size (number of strains per rep) automatically and return data frame
PagelsLambda_replic_Split <- function(phenodata, phylotree, x, split, size, rep){
  x <- x
  split <- split
  phenodata <- phenodata
  phylotree <- phylotree
  size <- size
  rep <- rep
  # subset
  out1 <- as.data.frame(t(replicate( expr = PagelsLambda_subsetlist_persplit_N2(phenodata, phylotree, x, split, size),  n = rep, simplify = "array")))
  out1_df <- as.data.frame(matrix(unlist(out1), nrow = length(unlist(out1[1]))))
  colnames(out1_df) <- c("lambda", "logL", "logL0", "P")
  out1_df$size <- size
  out1_df$split <- split
  
  return(out1_df)
}

# Combine
set.seed(1111)
Pagel_FamilySplit_9_1000rep_N2 <- PagelsLambda_replic_Split(all.strain, Treepheno$phy, "Strain", "family", 5, 1000)

#plot
Pagel_Fam9_plot_N2 <-  ggplot(Pagel_FamilySplit_9_1000rep_N2, aes(x = lambda, y = P))+
  geom_point(alpha = 0.2) +
  theme_bw()+
  scale_y_log10()+
  xlab(expression(paste("Pagel's")~lambda))+
  ylab(expression(italic("P")))+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank())
Pagel_Fam9_plot_N2

write.table(Pagel_FamilySplit_9_1000rep_N2, "Pagel_FamilySplit_5_1000rep_N2.txt", sep="\t", row.names=F)



#F23
PagelsLambda_subsetlist_persplit_F23 <- function( phenodata, phylotree, x, split, size){
  # Note phenodata is a dataframe containing a column x corresponding to the names in the tip.labels of phylotree, 
  # a column for the split (e.g. Family) and at least a column "d5", which contains the phenotypic data
  # check whehter the Treepheno data and the ProtPhenoDataframe names for the bacs are ok
  if(name.check(phylotree, phenodata, data.names= phenodata[[x]]) != "OK") warning("Not all taxa in phylogenetic tree and phenotypic data are equal or present in both data.")
  x <- x
  split <- split
  size <- size
  Strains <- SubsetStrains(phenodata, x, split, size)
  phenodata_sub <- setNames(phenodata[["d5"]][phenodata[[x]] %in% Strains], phenodata[[x]][phenodata[[x]] %in% Strains])
  # subset the tree and drop all tips that don't have phenotypic data
  Treesub <- treedata(phylotree, phenodata_sub, warning=F)
  
  # using Pagel's lambda. Note that we could also add se  = standard error -> could make sense
  Lambda = phylosig(Treesub$phy, phenodata_sub, method = "lambda", test = TRUE)
  return(Lambda[1:4])
}


# Function to do this for rep (number of replications) of size (number of strains per rep) automatically and return data frame
PagelsLambda_replic_Split_F23 <- function(phenodata, phylotree, x, split, size, rep){
  x <- x
  split <- split
  phenodata <- phenodata
  phylotree <- phylotree
  size <- size
  rep <- rep
  # subset
  out1 <- as.data.frame(t(replicate( expr = PagelsLambda_subsetlist_persplit_F23(phenodata, phylotree, x, split, size),  n = rep, simplify = "array")))
  out1_df <- as.data.frame(matrix(unlist(out1), nrow = length(unlist(out1[1]))))
  colnames(out1_df) <- c("lambda", "logL", "logL0", "P")
  out1_df$size <- size
  out1_df$split <- split
  
  return(out1_df)
}

set.seed(1111)
Pagel_FamilySplit_9_1000rep_F23 <- PagelsLambda_replic_Split_F23(all.strain, Treepheno$phy, "Strain", "family", 5, 1000)

Pagel_Fam9_plot_F23 <-  ggplot(Pagel_FamilySplit_9_1000rep_F23, aes(x = lambda, y = P))+
  geom_point(alpha = 0.2) +
  theme_bw()+
  scale_y_log10()+
  xlab(expression(paste("Pagel's ")~lambda))+
  ylab(expression(italic("P")))+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank())
Pagel_Fam9_plot_F23
write.table(Pagel_FamilySplit_9_1000rep_F23, "Pagel_FamilySplit_5_1000rep_F23.txt", sep="\t", row.names=F)

###############################-----------------------------------------------------------------
# Abouheif's Cmean


AbouheifFam9Split_N2 <- data.frame(Obs = numeric(), Std.Obs = numeric(), Expectation = numeric(), Variance = numeric(), 
                                   Pvalue = numeric(), rep = numeric())
AbouheifFam9Split_F23 <- data.frame(Obs = numeric(), Std.Obs = numeric(), Expectation = numeric(), Variance = numeric(), 
                                    Pvalue = numeric(), rep = numeric())

for(i in c(1:500)){
  
  Strains <- SubsetStrains(df = all.strain, x = "Strain", split = "family", size = 5)
  phenodata_sub_N2 <- setNames(all.strain[["N2"]][all.strain[["Strain"]] %in% Strains], all.strain[["Strain"]][all.strain[["Strain"]] %in% Strains])
  phenodata_sub_F23 <- setNames(all.strain[["d5"]][all.strain[["Strain"]] %in% Strains], all.strain[["Strain"]][all.strain[["Strain"]] %in% Strains])
  # subset the tree and drop all tips that don't have phenotypic data
  Treesub_N2 <- treedata(Treepheno$phy, phenodata_sub_N2, warning=F)
  Treesub_F23 <- treedata(Treepheno$phy, phenodata_sub_F23, warning=F)
  
  TreesubSorted_N2 <- Treesub_N2$data[match(Treesub_N2$phy$tip.label, row.names(Treesub_N2$data)),]
  TreesubSorted_F23 <- Treesub_F23$data[match(Treesub_F23$phy$tip.label, row.names(Treesub_F23$data)),]
  TreesubSorted_N2_df<-as.data.frame(TreesubSorted_N2)
  
  # Then Abouheif's Cmean
  W1_N2 <-proxTips(Treesub_N2$phy,method="oriAbouheif") 
  Abouheif_N2 <- abouheif.moran(TreesubSorted_N2_df,W1_N2)
  W1_F23 <- proxTips(Treesub_F23$phy,method="oriAbouheif") 
  Abouheif_F23 <- abouheif.moran(TreesubSorted_F23,W1_F23)
  
  outputAbou_N2 <- setNames(c(Abouheif_N2$obs,Abouheif_N2$expvar$Std.Obs, Abouheif_N2$expvar$Expectation, Abouheif_N2$expvar$Variance, 
                              Abouheif_N2$pvalue, Abouheif_N2$rep), nm = c("Obs", "Std.Obs" , "Expectation", "Variance" , 
                                                                           "Pvalue", "rep"))
  AbouheifFam9Split_N2[i,] <- outputAbou_N2
  
  outputAbou_F23 <- setNames(c(Abouheif_F23$obs,Abouheif_F23$expvar$Std.Obs, Abouheif_F23$expvar$Expectation, Abouheif_F23$expvar$Variance, 
                               Abouheif_F23$pvalue, Abouheif_F23$rep), nm = c("Obs", "Std.Obs" , "Expectation", "Variance" , 
                                                                              "Pvalue", "rep"))
  AbouheifFam9Split_F23[i,] <- outputAbou_F23
  
  
}

write.table(AbouheifFam9Split_N2, file = "AbouheifFam5Split_N2.txt", sep="\t", row.names=F)
write.table(AbouheifFam9Split_F23, file = "AbouheifFam5Split_F23.txt", sep="\t", row.names=F)

#PLOT
Abouheif_Fam9_density_N2 <-  ggplot(AbouheifFam9Split_N2, aes(x = Obs))+
  geom_density()+ 
  theme_bw()+
  xlab(expression(paste("Abouheif's ", C[mean])))+
  scale_x_continuous()+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank())#+
#theme(axis.line = element_line(colour = "black"))
Abouheif_Fam9_density_N2
Abouheif_Fam9_density_F23 <-  ggplot(AbouheifFam9Split_F23, aes(x = Obs))+
  geom_density()+ 
  theme_bw()+
  xlab(expression(paste("Abouheif's ", C[mean])))+
  scale_x_continuous()+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border = element_blank())
Abouheif_Fam9_density_F23
# no y labels
phylogeneticSignal<- ggarrange(Pagel_Fam9_plot_N2, Abouheif_Fam9_density_N2,
                               Pagel_Fam9_plot_F23, Abouheif_Fam9_density_F23, 
                               nrow=2, ncol=2, labels = "auto")
#add y labels
require(grid)
phylogeneticSignal_N2<-ggarrange(Pagel_Fam9_plot_N2, Abouheif_Fam9_density_N2, 
                                 nrow=1, ncol=2, labels = c('c','d'))
phylogeneticSignal_N2.lab<-annotate_figure(phylogeneticSignal_N2, 
                                           left = textGrob("Effects on N2 growth phenotype", rot = 90#, hjust=-0.3, gp = gpar(cex = 1.3)
                                           ))

phylogeneticSignal_F23<-ggarrange(Pagel_Fam9_plot_F23, Abouheif_Fam9_density_F23, 
                                  nrow=1, ncol=2, labels = c('a','b'))
phylogeneticSignal_F23.lab<-annotate_figure(phylogeneticSignal_F23, 
                                            left = textGrob("Effects on F23 growth phenotype", rot = 90#, hjust=-0.2, gp = gpar(cex = 1.3)
                                            ))
phylogeneticSignal_all<-ggarrange( phylogeneticSignal_F23.lab, phylogeneticSignal_N2.lab,
                                   nrow=2)
phylogeneticSignal_all

ggsave(filename="phylogeneticSignal_all.pdf", plot = phylogeneticSignal_all, width = 7, height = 7, useDingbats=F)

######################################################################################
#Fig. S8 were generated using BRIG with FASTA files located in ./BRIG/*

#####################################################################################
#Fig. S29ab
library(tidyverse)
library(ggthemes)
library(cowplot)
library(patchwork)
library(emoGG)
library(ggplot2)
library(readxl)
library(ggpubr)

growth<-read_excel('input/growth.74_plot.xlsx')
growth_d5<-growth[,c('ID','N2','d5')]
growth_d5.1<-
  mutate(growth_d5,
         N2 = N2, 
         d5 = -d5,
  ) %>% 
  gather(cat, value, -ID)

#add tree
library(ggtree)
library(paletteer)
tax<-read_excel('input/growth.74_plot.xlsx',sheet = 'tax')
#16s
#tree<-read.tree("input/74trim.fasta.raxml.bestTree")
#whole genome
tree<-read.tree("input/eztree-HQBiome.nwk")


cols2 <- c(Actinobacteria= '#EEE685',Bacteroidetes = "#4F94CD", Firmicutes= '#EE7600', 
           Alphaproteobacteria = "#00CD00",Gammaproteobacteria = "#90EE90")

#cols2<-c(Bacilli='#E69F00',Actinomycetia='#674ea7', Alphaproteobacteria='#56B4E9',Flavobacteriaceae='#ffe599',
#         Chromatiales='#0000ff', Alteromonadales='#6d9eeb', Oceanospirillales = '#a2c4c9',
#         Pseudomonas='#CC79A7', Marinagarivorans='#0c343d', Psychrobacter = '#0072B2',
#         Shewanella = '#D55E00', Vibrionaceae ='#9ACD32',Escherichia='#31890b')

p  = ggtree(tree, branch.length='none',layout = 'roundrect') %<+% tax #roundrect
labs=c('Actinobacteria','Bacteroidetes','Firmicutes','Alphaproteobacteria','Gammaproteobacteria')
tax$`Phylum/Class`<-factor(tax$`Phylum/Class`, levels = labs)
p1<-p +
  #geom_tiplab(align = T, size = 3#,hjust=-0.1
  #            ) +
  geom_tippoint(aes(color=`Phylum/Class`), size=4, shape = 15)+
  scale_color_manual(values = cols2#,guide = guide_legend(ncol = 7)
  )+
  theme_tree(legend.position = 'none')
# geom_cladelabel(node=5,label = '11')+
#geom_hilight(node=7, fill="darkgreen", alpha=.6)

#xlim(0,74)
p1

theme_own = function(base_size = 11, base_family = "",base_line_size = base_size/22, 
                     base_rect_size = base_size/22)
{
  theme_grey(base_size = base_size, base_family = base_family, 
             base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
    theme(panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA, colour = "grey20"), 
          panel.grid = element_line(colour = "grey92"), 
          panel.grid.minor = element_blank(), #panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "grey85", colour = "grey20"), 
          legend.key = element_rect(fill = "white", colour = NA,size = 1), complete = TRUE,
          #axis.text.y = element_blank(),
          #legend.position = "top",
          axis.ticks.y = element_blank(), 
          axis.text.y = element_text(hjust = 0, color='black'),
          axis.text = element_text(family=''),
          legend.text = element_text(size = 10,family=''))  
}


#F23 stack
growth.f23<-read_excel('input/growth.74_plot.xlsx', sheet='stack.f23')
growth.f23.2<-growth.f23[,c('ID','d3','d4','d5','d6','d7','d8','d9','d10','Remains')]%>%
  mutate() %>% 
  gather(day, value, -ID)
growth.f23.2$day<-factor(growth.f23.2$day,levels =c('Remains','d10','d9','d8','d7','d6','d5','d4','d3'))
#Sort by tree
growth.f23.2$ID<-
  factor(growth.f23.2$ID,
         levels = p$data %>% na.omit() %>% arrange(y) %>% pull(label))
f23<-ggplot(growth.f23.2,aes(ID, value))+
  geom_col(aes(fill = day), width = 1, color = "white", size = 0.6)+
  coord_flip(clip = "off")+
  scale_y_continuous(expand = c(0, 0),"F23 developmental rate")+
  scale_fill_manual(values = c("d3" = "#008ECEFF", "d4" = "#00A9E0FF", "d5" = "#59C7EBFF", 
                               "d6" = "#FC8D59FF", "d7" = "#EF6548FF",
                               "d8"="#D7301FFF","d9"="#B30000FF","d10"="#7F0000FF",
                               'Remains'='#696969'))+
  theme_classic()+
  xlab(NULL)+theme(axis.text.y = element_text(size = 10,family='',color = 'black',hjust = 0),axis.text.x = element_text(color = 'black'))

f23
#N2 box
growth.n2<-read_excel('input/growth.74_plot.xlsx', sheet='box.n2')
growth.n2.2<-growth.n2[,-c(14)] %>%
  mutate() %>% 
  gather(rep, value, -ID) %>% na.omit
#Sort by tree
growth.n2.2$ID<-
  factor(growth.n2.2$ID,
         levels = p$data %>% na.omit() %>% arrange(y) %>% pull(label))
n2<-ggplot(growth.n2.2,aes(ID, value)) + geom_boxplot(color = "#B34660",size=0.7)+#geom_point()+
  coord_flip(clip = "off")+
  scale_y_log10(#expand = c(0, 1),
    "Egg laying time of N2 (h)",breaks = c(50,60,70,100,130,170))+
  labs(x = NULL, y = NULL#, 
       #title = "",
       # subtitle = "\n"
  ) +
  theme(axis.text = element_text(size = 6),
        axis.ticks.length = unit(1, "pt"))+
  theme_classic() + 
  xlab(NULL)+
  theme(axis.text.y =  element_blank(),
        axis.text.x = element_text(color = 'black')
        #,axis.text.y = element_blank()
  )

n2
abc<-ggarrange(p1, f23,n2,nrow=1,align ='h',widths =c(0.4,0.6,0.4),labels = c("a", "b",'c')) 
abc
ggsave(filename="S19ab.pdf", plot = abc, width = 9, height = 9, useDingbats=F)

#Fig. S29c
library(readxl)
library(ggplot2)
library(ggalt)
library(dplyr)
#import
all <- read_excel("input/growth-pcoa.xlsx")
growth <- all[,c(1,4,7:15)]
growth<-growth[,-1]
rownames(growth)<-all$ID
growth.log<-log2(growth+1)
#distance calculation
#using correlation as distance
growth1<-t(growth)
d <- 1 - cor(growth1,method = 'spearman')

distance<-as.matrix(dist(growth.log, method = 'maximum',p = 2))

pcoa <- cmdscale(distance, k=2, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")
points$`Phylum/Class` <- all$`Phylum/Class`[match(rownames(points), all$ID)] 
points$`Gram strain` <- all$`Gram strain`[match(rownames(points), all$ID)] 

labs=c('Actinobacteria','Bacteroidetes','Firmicutes','Alphaproteobacteria','Gammaproteobacteria')
points$`Phylum/Class` <- factor(points$`Phylum/Class`, levels = labs)
cols <- c(Actinobacteria= '#EEE685',Bacteroidetes = "#4F94CD", Firmicutes= '#EE7600', 
          Alphaproteobacteria = "#00CD00",Gammaproteobacteria = "#90EE90")
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

pcoa.plot<-ggplot(points, aes(x=x, y=y, color=`Phylum/Class`, shape=`Gram strain`))+
  geom_point(alpha=1,size=3)+
  scale_colour_manual(values=cols)+
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+
  main_theme+
  #geom_encircle(data = points[points$`Phylum/Class`=='Actinobacteria',],fill='#EEE685',alpha = 0.1, show.legend = F,spread=0.002)+
  #geom_encircle(data = points[points$`Phylum/Class`=='Bacteroidetes',],fill='#4F94CD',alpha = 0.1, show.legend = F,spread=0.002)+
  #geom_encircle(data = points[points$`Phylum/Class`=='Firmicutes',],fill='#EE7600',alpha = 0.1, show.legend = F,spread=0.002)+
  #geom_encircle(data = points[points$`Phylum/Class`=='Alphaproteobacteria',],fill='#00CD00',alpha = 0.1, show.legend = F,spread=0.002)+
  #geom_encircle(data = points[points$`Phylum/Class`=='Gammaproteobacteria',],fill='#90EE90',alpha = 0.1, show.legend = F,spread=0.002)+
  theme(legend.position = c(0.8, 0.2))+ 
  theme(legend.key.size = unit(0.03, "inches"))+
  stat_ellipse(data = points, aes(color = `Gram strain`),level=0.68) 

pcoa.plot
corr<-as.data.frame(cor(growth.log,points[,c(1,2)],method = 'spearman')) #pcoa1:mean5, r=0.9546475 p= 2.2e-16; pcoa2:N2 r=0.100027  p= 0.3965
corr
#柴麻p value
cor.test(growth.log$mean_5,points[,c(1)],method = 'spearman',exact=FALSE)
cor.test(growth.log$N2,points[,c(2)],method = 'spearman',exact=FALSE)
cor_get_pval(corr)

#Fig. S29d
library(reshape2)
#rm(df)
df<-setNames(melt(distance), c('x', 'y', 'distance'))
tax.id<-all[,c(1,33)]
df<-merge(df,tax.id,by.x='x',by.y = 'ID',all.x = T, sort=F)
df<-merge(df,tax.id,by.x='y',by.y = 'ID',all.x = T, sort=F)
colnames(df)<-c('x', 'y', 'distance','tax1','tax2')
#rm(df1)
library(dplyr)
df1 = df %>% filter_(~tax1==tax2) %>% filter(distance != 0)
#write.csv(df1,'df1.csv')
library(stringr)
#rm(df2)
df2<-mutate(df1,num.x=gsub('HQB-','',df1$x))%>%mutate(num.y=gsub('HQB-','',df1$y))
df2$num.x<-gsub('OP','',df2$num.x)
df2$num.y<-gsub('OP','',df2$num.y)
df2$filt<-as.numeric(df2$num.x)+as.numeric(df2$num.y)
df2$filt2<-as.numeric(df2$num.x)*as.numeric(df2$num.y)
df2$filt3<-as.numeric(df2$filt)*(as.numeric(df2$filt2))^2
df3<-df2[!duplicated(df2$filt3),]  #1072
df3<-df3[-which(df3$tax1=="Bacteroidetes"),]  #remove Bacteroidetes, only two species

medians <- aggregate(df3$distance, by=list(df3$tax1), FUN=median)
order <- medians[sort(medians[, 2], index.return=T, decreasing=F)$ix, 1]
df3$tax1 <- factor(df3$tax1, levels=order)
#dodge = position_dodge(width=0.8)
#jdodge = position_jitterdodge(jitter.width = 0.1, jitter.height = 0, dodge.width = 0.8)
p1<-ggplot(df3, aes(x=tax1, y=distance, color=tax1)) +
  geom_jitter( size=1.7, alpha=0.25,width = 0.1) +
  geom_boxplot(alpha=0.8, #outlier.size=.2, 
               #size=1,#width = 0.25,
               width=0.4,cex=0.8,
               #position=dodge,
               color=c(Alphaproteobacteria = "#00CD00",Gammaproteobacteria ='#EEE685',Actinobacteria=  "#90EE90", Firmicutes= '#EE7600'), 
               fill=NA
  ) +
  #geom_violin(alpha = 1,
  #            fill = NA,
  #            #colour = NA,
  #            #position=dodge,
  #            width=0.9, cex=0.8)+
  labs(x="", y="functional distance") +
  scale_colour_manual(values=cols, guide='none') +
  coord_flip() +
  #scale_y_continuous(limits=c(.25, .75)) +
  main_theme +
  theme(axis.text.y = element_text(size=8))+
  ylab('Pair-wise phenotypic distance')

p1

library(rstatix)
#，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，
shapiro.test(df3$distance)  #W = 0.88302, p-value < 2.2e-16 
bartlett.test(distance ~ tax1, data = df3) #Bartlett's K-squared = 33.009, df = 3, p-value = 3.206e-07 
#，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，，

stat_test <- df3 %>% 
  wilcox_test(distance ~ tax1, p.adjust.method = "fdr",alternative = "two.sided") %>% 
  remove_ns()%>%add_xy_position(x = "tax1")
library(ggpubr)
p1<-p1+stat_pvalue_manual(stat_test, label = "p.adj.signif", step.increase = 0, hide.ns = TRUE, 
                          tip.length = 0, label.size = 5,coord.flip = TRUE,
)+scale_y_continuous(limits=c(0,10))+
  stat_compare_means(label.x=0.5, label.y = 7,size=3)
p1 

p2<-ggplot(df3, aes(x=distance)) +
  geom_histogram(size=.5, alpha=1, color="grey", fill="grey", binwidth=.01) +
  #labs(title="Functional diversity within Phylum/Class", x="") +
  main_theme +
  theme(legend.position="none", #axis.text.y=element_text(size=8),
        #title=element_text(size=6)
        plot.title = element_text(size = 6),
        axis.text.x= element_blank(),axis.ticks.x = element_blank()
  )+
  ylab('Counts')+xlab(NULL)+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(limits=c(0,10))

p2
library(ggpubr)

b<-ggarrange(p2,p1,ncol=1,align = 'hv',
             heights=c(2, 6)
)
b
ab<-ggarrange(pcoa.plot,b,ncol=2,widths = c(1.1,1),common.legend =F,
              labels = c('a','b'))
ab
ggsave('pcoa3.pdf',ab,height=5,width=10)


#Fig. S17
library(tidyverse)
library(ggpubr)
library(readxl)
library(ggforce)
library(RColorBrewer)
library(gridExtra)
theme_own = function(base_size = 11, base_family = "",base_line_size = base_size/22, 
                     base_rect_size = base_size/22)
{
  theme_grey(base_size = base_size, base_family = base_family, 
             base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
    theme(panel.background = element_rect(fill = "white", 
                                          colour = NA), panel.border = element_rect(fill = NA, 
                                                                                    colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "grey85", 
                                          colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                        colour = NA), complete = TRUE,
          axis.text.x = element_text(angle = 90, hjust = 1))
}
#all strains
all.growth_Info <- read_excel("input/growth_cor.xlsx",sheet = 1)
#remove Epsilonproteobacteria , only one species
#day 5 growth rate of L. marina
all.growth_Info <-all.growth_Info[-241,]

#cols <- c(Alphaproteobacteria = "#7CCD7C",Betaproteobacteria = "#308014",Gammaproteobacteria = "#9ACD32",
#          Proteobacteria='#d9d9d9',Bacteroidetes = "#1C86EE", Actinobac='#980000',
#          Actinobacteria= '#EE6363', Firmicutes= '#FFD700' )
cols <- c(Actinobacteria= '#EEE685',Bacteroidetes = "#4F94CD", Firmicutes= '#EE7600', 
          Alphaproteobacteria = "#00CD00",Betaproteobacteria='#43CD80',
          Gammaproteobacteria = "#90EE90"#,Epsilonproteobacteria='#4EEE94'
)
figall =ggplot(all.growth_Info,aes(x = d5, y = N2.log)) +
  geom_point(aes(color = `Phylum/Class`,shape= Gram), alpha=1,size = 2) +
  theme_own()+
  geom_smooth(method="lm", col="darkgrey", se =F)+ 
  xlab("Proportion of L4 of L. marina in day 5 (%)")+
  #scale_y_continuous(breaks = c(1,3,5,7)) +
  theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5)) +
  #scale_color_identity()+
  scale_colour_manual(name = "Phylum / Class", values = cols)+ 
  guides(colour = guide_legend(override.aes = list(shape = 15,size=4)))+ 
  ylab("Egg laying time of C. elegans (-log(h))")+
  stat_cor(aes(x = N2.log, y = d5,label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           method = "spearman", cor.coef.name = "R[s]", 
           label.x = 10, label.y = -2, size = 4.5)#+ 
#geom_rect(xmin = 75, xmax = 100, ymin = 5, ymax = 7, size = 0.5, fill = "#00000000", color = "darkgrey", linetype=2 ) #鮫裳侘
figall
ggsave(filename="growth_correlation_all.pdf", plot = figall, width = 6.75, height = 4.25, useDingbats=F)

#day 10 growth of L. marina
figall_10 =ggplot(all.growth_Info,aes(y = N2.log, x = d10)) +
  geom_point(aes(color = `Phylum/Class`,shape= Gram), alpha=1,size = 2) +
  theme_own()+
  geom_smooth(method="lm", col="darkgrey", se =F)+ 
  xlab("Proportion of L4 of L. marina in day 10 (%)")+
  #scale_y_continuous(breaks = c(1,3,5,7)) +
  theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5)) +
  #scale_color_identity()+
  scale_colour_manual(name = "Phylum / Class", values = cols)+ 弼
  guides(colour = guide_legend(override.aes = list(shape = 15,size=4)))+ 
  ylab("Egg laying time of C. elegans (-log(h))")+
  scale_x_continuous(breaks = c(0,20,40,60,80))+
  stat_cor(aes(x = N2.log, y = d10,label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           method = "spearman", cor.coef.name = "R[s]", 
           label.x = 10, label.y = -2, size = 4.5)#+ 
#geom_rect(xmin = 75, xmax = 100, ymin = 5, ymax = 7, size = 0.5, fill = "#00000000", color = "darkgrey", linetype=2 ) #鮫裳侘
figall_10

#day5 and day10
ab<-ggarrange(figall,figall_10,ncol=1,heights =c(1,1),widths = c(1,1), 
              labels = c("a", "b"),common.legend = T,legend='right')
ab

#col.class-d5   Fig. S17c
all.growth_Info$`Phylum/Class` <- factor(all.growth_Info$`Phylum/Class`, levels = unique(all.growth_Info$`Phylum/Class`))
all.growth_Info$`Phylum/Class` <- factor(all.growth_Info$`Phylum/Class`, levels = 
                                           c('Actinobacteria','Bacteroidetes','Firmicutes','Alphaproteobacteria',
                                             'Betaproteobacteria','Gammaproteobacteria','Epsilonproteobacteria'))

col.class = 
  ggplot(all.growth_Info, aes(x = d5, y = N2.log)) + 
  geom_point(aes(color = color), alpha=0.8) +
  scale_color_identity()+
  theme_own()+
  facet_wrap(~`Phylum/Class`)+
  xlab("Proportion of L4 of L. marina in day 5 (%)")+
  scale_y_continuous("Egg laying time of C. elegans (-log(h))") +
  theme(legend.position = "bottom")+
  #theme(strip.text = element_text(face = "italic"))+
  geom_smooth(method="lm", col="darkgrey", se =F,size=0.4)+
  stat_cor(aes(x = N2.log, y = d10,label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           method = "spearman", cor.coef.name = "R[s]", 
           label.x = 10, label.y = -2, size = 4.5)#+ 
col.class
Figabc <- ggarrange(ab, col.class, nrow=2,labels = c("", "c"),heights =c(1,0.8),widths = c(2.5,1))
Figabc
ggsave('Fig S22.pdf',Figabc,width = 8, height = 12, useDingbats=F)

##############################################################################
#Fig. S26 correlation of the effects of HQbiome isolates on L. marina and C. elegans
growth_Info <- read_excel("input/growth_cor.xlsx", sheet = '74')
theme_own = function(base_size = 11, base_family = "",base_line_size = base_size/22, 
                     base_rect_size = base_size/22)
{
  theme_grey(base_size = base_size, base_family = base_family, 
             base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
    theme(panel.background = element_rect(fill = "white", 
                                          colour = NA), panel.border = element_rect(fill = NA, 
                                                                                    colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "grey85", 
                                          colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                        colour = NA), complete = TRUE,
          axis.text.x = element_text(angle = 90, hjust = 1))
}
cols <- c(Alphaproteobacteria = "#7CCD7C",Betaproteobacteria = "#308014",Gammaproteobacteria = "#9ACD32",
          Proteobacteria='#d9d9d9',Bacteroidetes = "#1C86EE", Actinobac='#980000',
          Actinobacteria= '#EE6363', Firmicutes= '#FFD700' )

fig74=ggplot(growth_Info,aes(x = mean_5, y = log_N2)) +
  geom_point(aes(color = `Phylum/Class`,shape= `Gram strain`), alpha=1,size = 2) +
  theme_own()+
  geom_smooth(method="lm", col="darkgrey", se =F)+ #指拷��
  xlab("Adult rate of F23 in day 5 (%)")+
  #scale_y_continuous(breaks = c(1,3,5,7)) +
  theme(axis.text.x = element_text(angle=0, hjust=0.5, vjust=0.5)) +
  #scale_color_identity()+
  scale_colour_manual(name = "Taxanomy", values = cols)+
  ylab("Egg laying time of N2 (1/log(h))")+
  stat_cor(aes(x = log_N2, y = mean_5,label = paste(..r.label.., ..p.label.., sep = "~`,`~")), method = "spearman", cor.coef.name = "rho", 
           label.x = 30, label.y = 0.5, size = 4.5)#+ #指拷狼方
  #geom_rect(xmin = 75, xmax = 100, ymin = 5, ymax = 7, size = 0.5, fill = "#00000000", color = "darkgrey", linetype=2 ) #鮫裳侘
fig74


