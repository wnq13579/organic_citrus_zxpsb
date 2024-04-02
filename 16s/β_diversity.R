#### 1. loading required libraries ####
library(ggplot2)#作图 plot
library(ggpubr)#添加显著性标记, Add the significance marker
library(ggsignif)#添加显著性标记, Add the significance marker
library(dplyr)#数据清洗，Data cleaning
library(plyr)#数据清洗，Data cleaning
library(ggthemes)#ggplot所用主题，Themes for ggplot2
library(readxl)#读入 excel, read excel
library(ggsci)#配色，color scheme
library(showtext)#字体设置, font setting
library(extrafont)#使用系统字体，Using the system fonts
library(sysfonts)#加载系统字体，loading the system fonts
library(Cairo)#抗锯齿,anti-aliasing
library(ape)
library(phyloseq)#microbiome analysis
library(vegan)#Adonis analysis
#### 2. Setting themes and working dictionary path ####
loadfonts()
Sys.setenv(R_GSCMD = "C:/Program Files (x86)/gs/gs9.50/bin/gswin32c.exe")

mytheme1 <- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#4d4d4d",size=0.2),
                              panel.background = element_rect(colour = "#4d4d4d",size=0.2),
                              text = element_text(family = "Arial"),
                              strip.text = element_text(size = 6,hjust = 0.5),
                              plot.title = element_text(size = 6,hjust = 0.5),
                              axis.text=element_text(size=6,color = "#4D4D4D"),
                              axis.title=element_text(size = 6),
                              legend.text = element_text(size = 6),
                              legend.title = element_text(size = 6),
                              legend.background = element_blank(),
                              panel.border = element_rect(colour = NA),
                              axis.line = element_line(color = "#4D4D4D",size=0.2),
                              axis.ticks.length = unit(0.8, "mm"))#移除整体的边???

mytheme_bigfonts <- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#4d4d4d",size=0.2),
                                      panel.background = element_rect(colour = "#4d4d4d",size=0.2),
                                      text = element_text(family = "Arial"),
                                      strip.text = element_text(size = 9,hjust = 0.5),
                                      plot.title = element_text(size = 9,hjust = 0.5),
                                      axis.text=element_text(size=9,color = "#4D4D4D"),
                                      axis.title=element_text(size = 9),
                                      legend.text = element_text(size = 9),
                                      legend.title = element_text(size = 9),
                                      legend.background = element_blank(),
                                      panel.border = element_rect(colour = NA),
                                      axis.line = element_line(color = "#4D4D4D",size=0.2),
                                      axis.ticks.length = unit(0.8, "mm"))#移除整体的边???

wdImport <- c("E:/whx_24samples_filtered1/01.import")
wdOutput <- c("E:/whx_24samples_filtered1/02.diversity/β_diversity")
data.set.name = '_Filtered1' 
#### 3. All samples ####
### 3.1 Import and process data ###
setwd(wdImport)
ASVCount <- read.table("feature_table_24samples_filtered1.txt",header=T,row.names = 1,na.strings = c("NA"))
SampleData <- read.table("24SampleData1.txt",header=T,row.names = 1,na.strings = c("NA"))
colnames(ASVCount)<- row.names(SampleData)
RootedTree_FourStages_rhizosphere_Filtered1<- read.tree("tree.nwk")
Df <- phyloseq(otu_table(ASVCount, taxa_are_rows = T),sample_data(SampleData),phy_tree(RootedTree_FourStages_rhizosphere_Filtered1))
Df
Dfr <-transform_sample_counts(Df, function(x) x / sum(x) )
setwd(wdOutput)
sample_data(Dfr)$group<- factor(sample_data(Dfr)$group,level=c("CK_B","CZB_B","NF_B","CK_R","CZB_R","NF_R"))
sample_data(Dfr)$treatment<- factor(sample_data(Dfr)$treatment,level=c("CK","CZB","NF"))
Dfr_rhizo<-subset_samples(Dfr,compartment=="rhizosphere")
Dfr_Bulk<-subset_samples(Dfr,compartment=="bulk")

#### 4.bray ####
### 4.1 Rhizo ###
Dfr_rhizosphere <- subset_samples(Dfr_rhizo)
Dfr_rhizosphere
SampleData_rhizosphere<-filter(SampleData,compartment=="rhizosphere")
bray_dis_rhizosphere <- distance(Dfr_rhizosphere, method = 'bray') 
adonis2(formula=bray_dis_rhizosphere ~ group, 
        data=SampleData_rhizosphere,
        permutations=999)

#plots#
Dfr_rhizosphere
sample_data(Dfr_rhizosphere)$group<-factor(sample_data(Dfr_rhizosphere)$group)
bray_PCoA_rhizosphere <- ordinate(Dfr_rhizosphere, method="PCoA", distance="bray")
sample_data(Dfr_rhizosphere)$treatment<-factor(sample_data(Dfr_rhizosphere)$treatment)
bray_PCoAPoints_rhizosphere <- plot_ordination(Dfr_rhizosphere, bray_PCoA_rhizosphere,type = "samples", color="treatment")+
  guides(color=guide_legend(title=NULL),shape=guide_legend(title=NULL))+geom_point(size=1)+
  xlab("PCoA1(32.9%)")+ylab("PCoA2(22.1%)")+
  scale_x_continuous(breaks=c(-0.2,0.0,0.2))+scale_y_continuous(breaks=c(-0.2,0.0,0.2))+
  stat_ellipse(geom = "polygon",level=0.682,linetype =2,size=0,aes(fill=treatment,color=treatment),alpha=0.2)+
  scale_color_manual(values=c("#5C4A46","#DA9464","#888e4a"))+
  mytheme1
bray_PCoAPoints_rhizosphere
#setwd(wdOutput)
#getwd()
#ggsave(paste("bray_PCoAPoints_rhizosphere_filter",data.set.name,".pdf",sep=""),bray_PCoAPoints_rhizosphere,
#       device=cairo_pdf,width=90,height=80,dpi = 600,units = "mm")

### 4.2 bulk ###
Dfr_bulk <- subset_samples(Dfr_Bulk)
Dfr_bulk
SampleData_bulk<-filter(SampleData,compartment=="bulk")
#adonis
bray_dis_bulk <- distance(Dfr_bulk, method = 'bray') 
adonis2(formula=bray_dis_bulk ~ group, 
        data=SampleData_bulk,
        permutations=999)

#plots#
Dfr_bulk
sample_data(Dfr_bulk)$group<-factor(sample_data(Dfr_bulk)$group)
bray_PCoA_bulk <- ordinate(Dfr_bulk, method="PCoA", distance="bray")
sample_data(Dfr_bulk)$treatment<-factor(sample_data(Dfr_bulk)$treatment)
bray_PCoAPoints_bulk <- plot_ordination(Dfr_bulk, bray_PCoA_bulk,type = "samples", color="treatment")+
  guides(color=guide_legend(title=NULL),shape=guide_legend(title=NULL))+geom_point(size=1)+
  xlab("PCoA1(23.8%)")+ylab("PCoA2(20.9%)")+
  scale_x_continuous(breaks=c(-0.2,0.0,0.2))+scale_y_continuous(breaks=c(-0.2,0.0,0.2))+
  scale_color_manual(values=c("#5C4A46","#DA9464","#888e4a"))+
  stat_ellipse(geom = "polygon",level=0.682,linetype =2,size=0,aes(fill=treatment),alpha=0.2)+
  mytheme1
bray_PCoAPoints_bulk
#setwd(wdOutput)
#getwd()
#ggsave(paste("bray_PCoAPoints_bulk_filter",data.set.name,".pdf",sep=""),bray_PCoAPoints_bulk,
#       device=cairo_pdf,width=90,height=80,dpi = 600,units = "mm")
### 4.3 ALL ###
Dfr_all <- subset_samples(Dfr)
Dfr_all
SampleData_all<-filter(SampleData)
#adonis
bray_dis_all <- distance(Dfr_all, method = 'bray') 
adonis2(formula=bray_dis_all ~ compartment, 
        data=SampleData_all,
        permutations=999)

#plots#
Dfr_all
sample_data(Dfr_all)$compartment<-factor(sample_data(Dfr_all)$compartment)
bray_PCoA_all<- ordinate(Dfr_all, method="PCoA", distance="bray")
sample_data(Dfr_all)$treatment<-factor(sample_data(Dfr_all)$treatment)
bray_PCoAPoints_all <- plot_ordination(Dfr_all, bray_PCoA_all,type = "samples", color="compartment")+
  guides(color=guide_legend(title=NULL),shape=guide_legend(title=NULL))+geom_point(size=1)+
  xlab("PCoA1(23.6%)")+ylab("PCoA2(18.9%)")+
  scale_x_continuous(breaks=c(-0.2,0.0,0.2))+scale_y_continuous(breaks=c(-0.2,0.0,0.2))+
  scale_color_manual(values=c("#5C4A46","#809A54"))+
  stat_ellipse(geom = "polygon",level=0.682,linetype =2,size=0,aes(fill=compartment),alpha=0.2)+
  mytheme1
bray_PCoAPoints_all
#setwd(wdOutput)
#getwd()
#ggsave(paste("bray_PCoAPoints_filter",data.set.name,".pdf",sep=""),bray_PCoAPoints_all,
#       device=cairo_pdf,width=90,height=80,dpi = 600,units = "mm")


####ALL####
bray_PCoAPoints_all_plot <- ggarrange(bray_PCoAPoints_all,bray_PCoAPoints_bulk,bray_PCoAPoints_rhizosphere,ncol = 3, labels = c("A", "B", "C"),legend = "top")
bray_PCoAPoints_all_plot
setwd(wdOutput)
getwd()
ggsave(paste("bray_PCoAPoints_all_plot",".pdf",sep=""),bray_PCoAPoints_all_plot,
       device=cairo_pdf,width=270,height=120,dpi = 600,units = "mm")
