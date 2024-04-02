#### 1. loading required libraries ####
library(ggplot2)#作图 plot
library(ggpubr)#添加显著性标记, Add the significance marker
library(ggsignif)#添加显著性标记, Add the significance marker
library(dplyr)#数据清洗，Data cleaning
library(plyr)#数据清洗，Data cleaning
library(ggthemes)#ggplot所用主题，Themes for ggplot2
library(readxl)#读入 excel, read excel
library(readr)#读入tsv
library(ggsci)#配色，color scheme
library(showtext)#字体设置, font setting
library(extrafont)#使用系统字体，Using the system fonts
library(sysfonts)#加载系统字体，loading the system fonts
library(Cairo)#抗锯齿,anti-aliasing
library(ape)
library(phyloseq)#microbiome analysis
library(vegan)#Adonis analysis
library(RColorBrewer)# 修改配色方案
library(amplicon)#扩增子R包

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
wdOutput <- c("E:/whx_24samples_filtered1/03.taxonomy")
data.set.name = '_Filtered1' 

#### 3. All samples ####
### 3.1 Import and process data ###
setwd(wdImport)
SampleData <- read.table("24SampleData.txt",header=T,row.names = 1,na.strings = c("NA"))
tax_phylum<- read.table("24Sample_PhylumPercent_filtered1.txt",header=T,sep=";",row.names = 1,na.strings = c("NA"))
SampleData$treatment<-factor(SampleData$treatment,levels=c("CK","CZB","NF"))
SampleData$group<-factor(SampleData$group,levels=c("CK_B","CZB_B","NF_B","CK_R","CZB_R","NF_R"))


#### 4 ####
### 4.1 all###
SampleData<-filter(SampleData)
(stackplot=tax_stackplot(tax_phylum, SampleData, groupID="group", topN=11, style="group")+
    scale_fill_manual(values=c("#809A54","#888e4a","#a4994f","#8d8243","#907640","#F2D8AE","#DA9464","#bc7f57","#8C5D42","#74564e","#5C4A46"))+
    mytheme1 )
zoom=1.5
stackplot
setwd(wdOutput)
getwd()
ggsave(paste("stackplot",".pdf",sep=""),
       stackplot,device=cairo_pdf,width=120,height=90,dpi = 300,units = "mm")
### 4.2 rhizo###
SampleData_rhizo<-filter(SampleData,compartment=="rhizosphere")
(stackplot_rhizo=tax_stackplot(tax_phylum, SampleData_rhizo, groupID="treatment", topN=11, style="group")
  + scale_fill_brewer(palette="Set3"))
zoom=1.5
stackplot_rhizo
setwd(wdOutput)
getwd()
ggsave(paste("stackplot_rhizo",".pdf",sep=""),
       stackplot_rhizo,device=cairo_pdf,width=180,height=120,dpi = 300,units = "mm")
### 4.3 bulk###
SampleData_bulk<-filter(SampleData,compartment=="bulk")
(stackplot_bulk=tax_stackplot(tax_phylum, SampleData_bulk, groupID="treatment", topN=10, style="group")
  + scale_fill_brewer(palette="Set3"))
zoom=1.5
stackplot_bulk
setwd(wdOutput)
getwd()
ggsave(paste("stackplot_bulk",".pdf",sep=""),
       stackplot_bulk,device=cairo_pdf,width=180,height=120,dpi = 300,units = "mm")

