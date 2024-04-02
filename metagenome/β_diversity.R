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

wdImport <- c("E:/microbiome/01.import")
wdOutput <- c("E:/microbiome/02.diversity/β_diversity")
data.set.name = '_Filtered1' 
#### 3. All samples ####
### 3.1 Import and process data ###
setwd(wdImport)
ASVCount <- read.table("24Sample_SpeciesPercent25_filtered1.txt",header=T,row.names = 1,na.strings = c("NA"))
SampleData <- read.table("24SampleData1.txt",header=T,row.names = 1,na.strings = c("NA"))
colnames(ASVCount)<- row.names(SampleData)
Df <- phyloseq(otu_table(ASVCount, taxa_are_rows = T),sample_data(SampleData))
Df
Dfr <-transform_sample_counts(Df, function(x) x / sum(x) )
setwd(wdOutput)
sample_data(Dfr)$group<- factor(sample_data(Dfr)$group,level=c("CK_B","CZB_B","NF_B","CK_R","CZB_R","NF_R"))
sample_data(Dfr)$treatment<- factor(sample_data(Dfr)$treatment,level=c("CK","CZB","NF"))
Dfr_rhizo<-subset_samples(Dfr,compartment=="rhizosphere")
Dfr_Bulk<-subset_samples(Dfr,compartment=="bulk")


#### 4.bray czb####
### 4.1 Rhizo ###
Dfr_rhizosphere_CZB <- subset_samples(Dfr_rhizo,treatment=="CK"|treatment=="CZB")
Dfr_rhizosphere_CZB
SampleData_rhizosphere_CZB<-filter(SampleData,compartment=="rhizosphere",treatment=="CK"|treatment=="CZB")
bray_dis_rhizosphere <- distance(Dfr_rhizosphere_CZB, method = 'bray') 
adonis2(formula=bray_dis_rhizosphere ~ group, 
        data=SampleData_rhizosphere_CZB,
        permutations=999)

#plots#
Dfr_rhizosphere_CZB
sample_data(Dfr_rhizosphere_CZB)$group<-factor(sample_data(Dfr_rhizosphere_CZB)$group)
bray_PCoA_rhizosphere <- ordinate(Dfr_rhizosphere_CZB, method="PCoA", distance="bray")
sample_data(Dfr_rhizosphere_CZB)$treatment<-factor(sample_data(Dfr_rhizosphere_CZB)$treatment)
bray_PCoAPoints_rhizosphere_CZB <- plot_ordination(Dfr_rhizosphere_CZB, bray_PCoA_rhizosphere,type = "samples", color="treatment")+
  guides(color=guide_legend(title=NULL),shape=guide_legend(title=NULL))+geom_point(size=1)+
  
  scale_x_continuous(breaks=c(-0.2,0.0,0.2))+scale_y_continuous(breaks=c(-0.2,0.0,0.2))+
  stat_ellipse(geom = "polygon",level=0.682,linetype =2,size=0,aes(fill=treatment,color=treatment),alpha=0.2)+
  scale_color_manual(values=c("#EB9752","#46C580"))+
  mytheme1
bray_PCoAPoints_rhizosphere_CZB
#setwd(wdOutput)
#getwd()
#ggsave(paste("bray_PCoAPoints_rhizosphere_CZB_filter",data.set.name,".pdf",sep=""),bray_PCoAPoints_rhizosphere_CZB,
#       device=cairo_pdf,width=90,height=80,dpi = 600,units = "mm")


#### 4.bray NF####
### 4.1 Rhizo ###
Dfr_rhizosphere_NF <- subset_samples(Dfr_rhizo,treatment=="CK"|treatment=="NF")
Dfr_rhizosphere_NF
SampleData_rhizosphere_NF<-filter(SampleData,compartment=="rhizosphere",treatment=="CK"|treatment=="NF")
bray_dis_rhizosphere <- distance(Dfr_rhizosphere_NF, method = 'bray') 
adonis2(formula=bray_dis_rhizosphere ~ group, 
        data=SampleData_rhizosphere_NF,
        permutations=999)

#plots#
Dfr_rhizosphere_NF
sample_data(Dfr_rhizosphere_NF)$group<-factor(sample_data(Dfr_rhizosphere_NF)$group)
bray_PCoA_rhizosphere <- ordinate(Dfr_rhizosphere_NF, method="PCoA", distance="bray")
sample_data(Dfr_rhizosphere_NF)$treatment<-factor(sample_data(Dfr_rhizosphere_NF)$treatment)
bray_PCoAPoints_rhizosphere_NF <- plot_ordination(Dfr_rhizosphere_NF, bray_PCoA_rhizosphere,type = "samples", color="treatment")+
  guides(color=guide_legend(title=NULL),shape=guide_legend(title=NULL))+geom_point(size=1)+
  
  scale_x_continuous(breaks=c(-0.2,0.0,0.2))+scale_y_continuous(breaks=c(-0.2,0.0,0.2))+
  stat_ellipse(geom = "polygon",level=0.682,linetype =2,size=0,aes(fill=treatment,color=treatment),alpha=0.2)+
  scale_color_manual(values=c("#EB9752","#75DAFF"))+
  mytheme1
bray_PCoAPoints_rhizosphere_NF
#setwd(wdOutput)
#getwd()
#ggsave(paste("bray_PCoAPoints_rhizosphere_NF_filter",data.set.name,".pdf",sep=""),bray_PCoAPoints_rhizosphere_NF,
#       device=cairo_pdf,width=90,height=80,dpi = 600,units = "mm")

#### 4.2 bulk CZB####
Dfr_bulk_CZB <- subset_samples(Dfr_Bulk,treatment=="CK"|treatment=="CZB")
Dfr_bulk_CZB
SampleData_bulk_CZB<-filter(SampleData,compartment=="bulk",treatment=="CK"|treatment=="CZB")
#adonis
bray_dis_bulk <- distance(Dfr_bulk_CZB, method = 'bray') 
adonis2(formula=bray_dis_bulk ~ group, 
        data=SampleData_bulk_CZB,
        permutations=999)

#plots#
Dfr_bulk_CZB
sample_data(Dfr_bulk_CZB)$group<-factor(sample_data(Dfr_bulk_CZB)$group)
bray_PCoA_bulk <- ordinate(Dfr_bulk_CZB, method="PCoA", distance="bray")
sample_data(Dfr_bulk_CZB)$treatment<-factor(sample_data(Dfr_bulk_CZB)$treatment)
bray_PCoAPoints_bulk_CZB <- plot_ordination(Dfr_bulk_CZB, bray_PCoA_bulk,type = "samples", color="treatment")+
  guides(color=guide_legend(title=NULL),shape=guide_legend(title=NULL))+geom_point(size=1)+
  
  scale_x_continuous(breaks=c(-0.2,0.0,0.2))+scale_y_continuous(breaks=c(-0.2,0.0,0.2))+
  scale_color_manual(values=c("#EB9752","#46C580"))+
  stat_ellipse(geom = "polygon",level=0.682,linetype =2,size=0,aes(fill=treatment),alpha=0.2)+
  mytheme1
bray_PCoAPoints_bulk_CZB
#setwd(wdOutput)
#getwd()
#ggsave(paste("bray_PCoAPoints_bulk_CZB_filter",data.set.name,".pdf",sep=""),bray_PCoAPoints_bulk_CZB,
#       device=cairo_pdf,width=90,height=80,dpi = 600,units = "mm")
### 4.2 bulk NF###
Dfr_bulk_NF <- subset_samples(Dfr_Bulk,treatment=="CK"|treatment=="NF")
Dfr_bulk_NF
SampleData_bulk_NF<-filter(SampleData,compartment=="bulk",treatment=="CK"|treatment=="NF")
#adonis
bray_dis_bulk <- distance(Dfr_bulk_NF, method = 'bray') 
adonis2(formula=bray_dis_bulk ~ group, 
        data=SampleData_bulk_NF,
        permutations=999)

#plots#
Dfr_bulk_NF
sample_data(Dfr_bulk_NF)$group<-factor(sample_data(Dfr_bulk_NF)$group)
bray_PCoA_bulk <- ordinate(Dfr_bulk, method="PCoA", distance="bray")
sample_data(Dfr_bulk_NF)$treatment<-factor(sample_data(Dfr_bulk_NF)$treatment)
bray_PCoAPoints_bulk_NF <- plot_ordination(Dfr_bulk_NF, bray_PCoA_bulk,type = "samples", color="treatment")+
  guides(color=guide_legend(title=NULL),shape=guide_legend(title=NULL))+geom_point(size=1)+
  
  scale_x_continuous(breaks=c(-0.2,0.0,0.2))+scale_y_continuous(breaks=c(-0.2,0.0,0.2))+
  scale_color_manual(values=c("#EB9752","#75DAFF"))+
  stat_ellipse(geom = "polygon",level=0.682,linetype =2,size=0,aes(fill=treatment),alpha=0.2)+
  mytheme1
bray_PCoAPoints_bulk_NF
#setwd(wdOutput)
#getwd()
#ggsave(paste("bray_PCoAPoints_bulk_NF_filter",data.set.name,".pdf",sep=""),bray_PCoAPoints_bulk_NF,
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
  
  scale_x_continuous(breaks=c(-0.2,0.0,0.2))+scale_y_continuous(breaks=c(-0.2,0.0,0.2))+
  scale_color_manual(values=c("#EB9752","#75DAFF"))+
  stat_ellipse(geom = "polygon",level=0.682,linetype =2,size=0,aes(fill=compartment),alpha=0.2)+
  mytheme1
bray_PCoAPoints_all
#setwd(wdOutput)
#getwd()
#ggsave(paste("bray_PCoAPoints_filter",data.set.name,".pdf",sep=""),bray_PCoAPoints_all,
#       device=cairo_pdf,width=90,height=80,dpi = 600,units = "mm")


####ALL####
bray_PCoAPoints_all_plot <- ggarrange(bray_PCoAPoints_all,bray_PCoAPoints_bulk_CZB,bray_PCoAPoints_bulk_NF,bray_PCoAPoints_rhizosphere_CZB,bray_PCoAPoints_rhizosphere_NF,ncol = 3,nrow = 2, labels = c("A", "B", "C", "D" ,"E"),legend = "top")
bray_PCoAPoints_all_plot
setwd(wdOutput)
getwd()
ggsave(paste("bray_PCoAPoints_all_plot_111",".pdf",sep=""),bray_PCoAPoints_all_plot,
       device=cairo_pdf,width=270,height=200,dpi = 600,units = "mm")
