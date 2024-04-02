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
library(vegan)#Adonis analysis
library(ComplexHeatmap)#heatmapping
library(pheatmap)#heatmapping
library(cluster)#hierarchical clustering
library(psych)#correlationship analysis#
#### 2. Setting themes and working dictionary path ####
loadfonts()
Sys.setenv(R_GSCMD = "C:/Program Files (x86)/gs/gs9.50/bin/gswin32c.exe")

mytheme1 <- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#4d4d4d",linewidth=0.2),
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
                              axis.line = element_line(color = "#4D4D4D",linewidth=0.2),
                              axis.ticks.length = unit(0.8, "mm"))#移除整体的边???

mytheme_bigfonts <- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#4d4d4d",linewidth=0.2),
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
                                      axis.line = element_line(color = "#4D4D4D",linewidth=0.2),
                                      axis.ticks.length = unit(0.8, "mm"))#移除整体的边???
wdImport<-("E:/microbiome/05.path")
wdOutput<- ("E:/microbiome/05.path")
####P####
setwd(wdImport)
getwd()
p_stamp_pathabundance_relab<- read.table("p_stamp_pathabundance_relab_mean.txt",header=T,row.names = 1,na.strings = c("NA"))
p_stamp_pathabundance_relab<-scale(t(p_stamp_pathabundance_relab),center = T)
p_stamp_pathabundance_relab<-t(p_stamp_pathabundance_relab)
#write.table(genus_biomarker, "genus_biomarker.txt", row.names = T, sep = ',', quote = FALSE)


p_stamp_pathabundance_relab_pheatmap<-pheatmap(p_stamp_pathabundance_relab,fontsize=10,border_color = "black",
                              fontsize_row =15,fontsize_col = 15,
                              fontsize_number = 22,
                              #显著性标记的符号大小
                              cluster_rows=T,clustering_distance_rows="correlation",
                              #指明行聚类，聚类依据的距离
                              cluster_cols=F,clustering_distance_cols="euclidean",
                              #指明列聚类，聚类依据的距离
                              clustering_method="centroid")
#聚类方法


p_stamp_pathabundance_relab_pheatmap

setwd(wdOutput)
getwd()
ggsave(paste("p_stamp_pathabundance_relab_pheatmap",".pdf",sep=""),p_stamp_pathabundance_relab_pheatmap,
       device=cairo_pdf,width=300,height=120,dpi = 600,units = "mm")
####ACID####
setwd(wdImport)
getwd()
acid_stamp_pathabundance_relab<- read.table("acid_r.txt",header=T,row.names = 1,na.strings = c("NA"))
acid_stamp_pathabundance_relab<-scale(t(acid_stamp_pathabundance_relab),center = T)
acid_stamp_pathabundance_relab<-t(acid_stamp_pathabundance_relab)
#write.table(genus_biomarker, "genus_biomarker.txt", row.names = T, sep = ',', quote = FALSE)


acid_stamp_pathabundance_relab_pheatmap<-pheatmap(acid_stamp_pathabundance_relab,fontsize=10,border_color = "black",
                                             fontsize_row =15,fontsize_col = 15,
                                             fontsize_number = 22,
                                             #显著性标记的符号大小
                                             cluster_rows=T,clustering_distance_rows="correlation",
                                             #指明行聚类，聚类依据的距离
                                             cluster_cols=F,clustering_distance_cols="euclidean",
                                             #指明列聚类，聚类依据的距离
                                             clustering_method="centroid")
#聚类方法


acid_stamp_pathabundance_relab_pheatmap

setwd(wdOutput)
getwd()
ggsave(paste("acid_stamp_pathabundance_relab_pheatmap",".pdf",sep=""),acid_stamp_pathabundance_relab_pheatmap,
       device=cairo_pdf,width=390,height=120,dpi = 600,units = "mm")
