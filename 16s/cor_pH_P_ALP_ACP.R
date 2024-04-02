#### 1. loading required libraries ####
library(ggplot2)#作图 plot
library(ggpubr)#添加显著性标记, Add the significance marker
library(ggsignif)#添加显著性标记, Add the significance marker
library(dplyr)#数据清洗，Data cleaning
library(plyr)#数据清洗，Data cleaning
library(reshape2)#数据清洗，Data cleaning
library(ggthemes)#ggplot所用主题，Themes for ggplot2
library(ggsci)#TOP期刊配色方案, color plates from top journal
library(readxl)#读入 excel, read excel
library(showtext)#字体设置, font setting
library(extrafont)#使用系统字体，Using the system fonts
library(sysfonts)#加载系统字体，loading the system fonts
library(Cairo)#抗锯齿,anti-aliasing
library(pheatmap)#heatmapping
library(cluster)#hierarchical clustering
library(psych)#correlationship analysis#
library(VennDiagram)#Venn plot

#### 2. Setting themes and working dictionary path ####
loadfonts()
Sys.setenv(R_GSCMD = "C:/Program Files (x86)/gs/gs9.50/bin/gswin32c.exe")

mytheme <- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#4d4d4d"),
                             text = element_text(family = "Times New Roman"),
                             strip.text = element_text(size = 6,hjust = 0.5),
                             plot.title = element_text(size = 6,hjust = 0.5),
                             axis.text=element_text(size=6,color = "#4D4D4D"),
                             axis.title=element_text(size = 6),
                             legend.text = element_text(size = 6),
                             legend.title = element_text(size = 6),
                             legend.background = element_blank(),
                             panel.border = element_rect(colour = NA),
                             axis.line = element_line(color = "#4D4D4D",size=0.2),
                             axis.ticks.length = unit(0.8, "mm"))

mytheme1 <- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#4d4d4d",size=0.2),
                              panel.background = element_rect(colour = "#4d4d4d",size=0.2),
                              text = element_text(family = "Times New Roman"),
                              strip.text = element_text(size = 5,hjust = 0.5),
                              plot.title = element_text(size = 5,hjust = 0.5),
                              axis.text=element_text(size=5,color = "#4D4D4D"),
                              axis.title=element_text(size = 5),
                              legend.text = element_text(size = 5),
                              legend.title = element_text(size = 5),
                              legend.background = element_blank(),
                              axis.line = element_line(color = "#4D4D4D",size=0.2),
                              axis.ticks.length = unit(0.8, "mm"))
wdImport<-("F:/whx_24samples_filtered1/01.import")
wdOutput <- ("F:/whx_24samples_filtered1/06.env/soil_nutrition")

#### 3 Import and process data ####
setwd(wdImport)
getwd()
soil <- read_excel("env_soil1.xlsx",
                             sheet = "Sheet1")
setwd(wdImport)
getwd()
soil1 <- read.table("env_soil1.txt",sep="\t",na.strings="",header = T,comment.char = "",check.names = F,stringsAsFactors = F)
SampleData <- read.table("24SampleData.txt",header=T,na.strings = c("NA"))
soil_all <- merge(SampleData, soil1, by = 'sampleid')
soil_all$compartment<-factor(soil_all$compartment,levels=c("bulk","rhizosphere"))
####pH_P####
pH_P_cor<-corr.test(soil$pH,soil$P_bulk,
                                              method="pearson",adjust="BH",minlength=5)
pH_P_cor
pH_P_cor$r
pH_P_cor$p

##plots
pH_P_cor_line <- ggplot(soil, aes(x= pH, y= P_bulk))+
  geom_point(color="#E64B35")+
  geom_smooth(method = "lm",color="#E64B35")+
  stat_regline_equation(label.x = 3, label.y = 14)+
  mytheme
pH_P_cor_line

formula <- y ~ x
pH_P_cor_line <- ggplot(soil, aes(x= pH, y= P_bulk))+
  geom_point(color="#E64B35")+
  geom_smooth(method = "lm",color="#E64B35")+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    formula = formula,label.x = 8,label.y = 90) +
  stat_cor (method = "pearson",
            label.x = 8,label.y = 95)+
  mytheme
pH_P_cor_line
#setwd(wdOutput)
#getwd()
#ggsave(paste("pH_P_cor_line",".pdf",sep=""),
       #pH_P_cor_line,
       #device=cairo_pdf,width=80,height=60,dpi = 600,units = "mm")

####ACP_P####
#all
ACP_P_cor<-corr.test(soil_all$ACP,soil_all$P,
                    method="pearson",adjust="BH",minlength=5)
ACP_P_cor
ACP_P_cor$r
ACP_P_cor$p


##plots
ACP_P_cor_line <- ggplot(soil_all, aes(x= ACP, y= P))+
  geom_point(aes(color=compartment))+
  geom_smooth(method = "lm",color="#E64B35")+
  stat_regline_equation(label.x = 3, label.y = 14)+
  mytheme
ACP_P_cor_line

#bulk
formula <- y ~ x
ACP_P_cor_line <- ggplot(soil, aes(x= ACP_bulk, y= P_bulk))+
  geom_point(color="#E64B35")+
  geom_smooth(method = "lm",color="#E64B35")+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    formula = formula,label.x = 7,label.y = 85) +
  stat_cor (method = "pearson",
            label.x = 7,label.y = 80)+
  mytheme
ACP_P_cor_line
#rhizo
formula <- y ~ x
ACP_P_cor_line <- ggplot(soil, aes(x= ACP_rhizo, y= P_rhizo))+
  geom_point(color="#E64B35")+
  geom_smooth(method = "lm",color="#E64B35")+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    formula = formula,label.x = 7,label.y = 85) +
  stat_cor (method = "pearson",
            label.x = 7,label.y = 80)+
  mytheme
ACP_P_cor_line
#all
formula <- y ~ x
ACP_P_cor_line <- ggplot(soil_all, aes(x= ACP, y= P))+
  geom_point(aes(color=compartment))+
  geom_smooth(method = "lm",color="#E64B35")+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    formula = formula,label.x = 7,label.y = 85) +
  stat_cor (method = "pearson",
            label.x = 7,label.y = 80)+
  mytheme
ACP_P_cor_line
#setwd(wdOutput)
#getwd()
#ggsave(paste("ACP_P_cor_line",".pdf",sep=""),
      # ACP_P_cor_line,
       #device=cairo_pdf,width=80,height=60,dpi = 600,units = "mm")


####ALP_P####
ALP_P_cor<-corr.test(soil_all$ALP,soil_all$P,
                     method="pearson",adjust="BH",minlength=5)
ALP_P_cor
ALP_P_cor$r
ALP_P_cor$p

##plots
ALP_P_cor_line <- ggplot(soil_all, aes(x= ALP, y= P))+
  geom_point(aes(color=compartment))+
  geom_smooth(method = "lm",color="#E64B35")+
  stat_regline_equation(label.x = 3, label.y = 14)+
  mytheme
ALP_P_cor_line

#bulk
formula <- y ~ x
ALP_P_cor_line <- ggplot(soil, aes(x= ALP_bulk, y= P_bulk))+
  geom_point(color="#E64B35")+
  geom_smooth(method = "lm",color="#E64B35")+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    formula = formula,label.x = 7,label.y = 85) +
  stat_cor (method = "pearson",
            label.x = 7,label.y = 80)+
  mytheme
ALP_P_cor_line
#rhizo
formula <- y ~ x
ALP_P_cor_line <- ggplot(soil, aes(x= ALP_rhizo, y= P_rhizo))+
  geom_point(color="#E64B35")+
  geom_smooth(method = "lm",color="#E64B35")+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    formula = formula,label.x = 7,label.y = 85) +
  stat_cor (method = "pearson",
            label.x = 7,label.y = 80)+
  mytheme
ALP_P_cor_line

#all
formula <- y ~ x
ALP_P_cor_line <- ggplot(soil_all, aes(x= ALP, y= P))+
  geom_point(aes(color=compartment))+
  geom_smooth(method = "lm",color="#E64B35")+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    formula = formula,label.x = 13.5,label.y = 85) +
  stat_cor (method = "pearson",
            label.x = 13.5,label.y = 80)+
  mytheme
ALP_P_cor_line

#setwd(wdOutput)
#getwd()
#ggsave(paste("ALP_P_cor_line",".pdf",sep=""),
       #ALP_P_cor_line,
       #device=cairo_pdf,width=80,height=60,dpi = 600,units = "mm")




####Fe_P_bulk####
Fe_P_bulk_cor<-corr.test(soil$Fe_bulk,soil$P_bulk,
                    method="pearson",adjust="BH",minlength=5)
Fe_P_bulk_cor
Fe_P_bulk_cor$r
Fe_P_bulk_cor$p

##plots
Fe_P_bulk_cor_line <- ggplot(soil, aes(x= Fe_bulk, y= P_bulk))+
  geom_point(color="#E64B35")+
  geom_smooth(method = "lm",color="#E64B35")+
  stat_regline_equation(label.x = 3, label.y = 14)+
  mytheme
Fe_P_bulk_cor_line

formula <- y ~ x
Fe_P_bulk_cor_line <- ggplot(soil, aes(x= Fe_bulk, y= P_bulk))+
  geom_point(color="#E64B35")+
  geom_smooth(method = "lm",color="#E64B35")+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    formula = formula,label.x = 8,label.y = 90) +
  stat_cor (method = "pearson",
            label.x = 8,label.y = 95)+
  mytheme
Fe_P_bulk_cor_line

setwd(wdOutput)
getwd()
ggsave(paste("Fe_P_bulk_cor_line",".pdf",sep=""),
Fe_P_bulk_cor_line,
device=cairo_pdf,width=80,height=60,dpi = 600,units = "mm")


####Fe_P_rhizo####
Fe_P_rhizo_cor<-corr.test(soil$Fe_rhizo,soil$P_rhizo,
                         method="pearson",adjust="BH",minlength=5)
Fe_P_rhizo_cor
Fe_P_rhizo_cor$r
Fe_P_rhizo_cor$p

##plots
Fe_P_rhizo_cor_line <- ggplot(soil, aes(x= Fe_rhizo, y= P_rhizo))+
  geom_point(color="#E64B35")+
  geom_smooth(method = "lm",color="#E64B35")+
  stat_regline_equation(label.x = 3, label.y = 14)+
  mytheme
Fe_P_rhizo_cor_line

formula <- y ~ x
Fe_P_rhizo_cor_line <- ggplot(soil, aes(x= Fe_rhizo, y= P_rhizo))+
  geom_point(color="#E64B35")+
  geom_smooth(method = "lm",color="#E64B35")+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    formula = formula,label.x = 8,label.y = 90) +
  stat_cor (method = "pearson",
            label.x = 8,label.y = 95)+
  mytheme
Fe_P_rhizo_cor_line

setwd(wdOutput)
getwd()
ggsave(paste("Fe_P_rhizo_cor_line",".pdf",sep=""),
       Fe_P_rhizo_cor_line,
       device=cairo_pdf,width=80,height=60,dpi = 600,units = "mm")












#### 4.3 all ####
cor_plot <- ggarrange(pH_P_cor_line,ACP_P_cor_line,ALP_P_cor_line,
                    ncol = 3, labels = c("A", "B","C"),common.legend = T,legend = "right")
cor_plot
setwd(wdOutput)
getwd()
ggsave(paste("cor_plot",".pdf",sep=""),cor_plot,
       device=cairo_pdf,width=240,height=90,dpi = 300,units = "mm")