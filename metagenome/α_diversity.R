#### 1. Loading package####
library(ggplot2)#作图 plot
library(ggpubr)#添加显著性标记, Add the significance marKer
library(ggsignif)#添加显著性标记, Add the significance marKer
library(dplyr)#数据清洗，Data cleaning
library(plyr)#数据清洗，Data cleaning
library(reshape2)#数据清洗，Data cleaning
library(ggthemes)#ggplot所用主题，Themes for ggplot2
library(grid)#分面和嵌合图，facet and Mosaic graK
library(agricolae)#多重比较，Multiple comparisons.
library(readxl)#读入 excel, read excel
library(readr)#读入 tsv
library(tidyr)#数据处理
library(ggsci)#配色，color scheme
library(survminer)#正态分布检验和方差齐性检验
library(showtext)#字体设置, font setting
library(car)#方差齐性检验，homogeneity test of variance, levene test
library(extrafont)#使用系统字体，Using the system fonts
library(sysfonts)#加载系统字体，loading the system fonts
library(Cairo)#抗锯齿,anti-aliasing
library(stringr)#字符串处理.string manipulation
library(graphics)#坐标轴表达式，expression for axis
library(cowplot)#拼图
library(forecast)#BoxCox转换 
library(PMCMRplus)#非参数检验


#### 2. setting theme and filepath ####
loadfonts()
Sys.setenv(R_GSCMD = "C:/Program Files (x86)/gs/gs9.50/bin/gswin32c.exe")

mytheme <- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#000000"),
                             text = element_text(family = "Arial"),
                             strip.text = element_text(size=7,hjust = 0.5),
                             plot.title = element_text(size=7,hjust = 0.5),
                             axis.text=element_text(size=7,color = "#808080"),
                             axis.title=element_text(size=7),
                             legend.text = element_text(size=7),
                             legend.title = element_text(size=7),
                             legend.background = element_blank(),
                             panel.border = element_rect(colour = NA),
                             axis.line = element_line(color = "black",size=0.4))#移除整体的边???

FacetTheme <- theme_few()+theme(strip.background = element_rect(fill="gray72",colour ="#000000"),
                                text = element_text(family = "Arial"),
                                strip.text = element_text(size=8,hjust = 0.5),
                                plot.title = element_text(size=8,hjust = 0.5),
                                axis.text =element_text(size=8,color = "black"),
                                axis.title =element_text(size=8,color = "black"),
                                legend.text = element_text(size=8,color = "black"),
                                legend.title = element_text(size=8,color = "black"),
                                legend.background = element_blank(),
                                axis.line = element_line(color = "black",size=0.4))#移除整体的边???
wdImport<- c("E:/microbiome/01.import")
wdOutput <- c("E:/microbiome/02.diversity/α_diversity")




#### 3. α_diversity ####
### 3.1 Import and process data ###
### 3.1.1 chao1 ###
setwd(wdImport)
SampleData <- read.table("24SampleData1.txt",header=T,na.strings = c("NA"))
chao1 <- read.table("chao1.txt",header=T,na.strings = c("NA"))
chao1 <- merge(SampleData, chao1, by = 'sampleid')
chao1$group<-factor(chao1$group,levels=c("CK_B","CZB_B","NF_B","CK_R","CZB_R","NF_R"))
chao1$treatment<-factor(chao1$treatment,levels=c("CK","CZB","NF"))
chao1$compartment<-factor(chao1$compartment,levels=c("bulk","rhizosphere"))
chao1_bulk<-filter(chao1,compartment=="bulk")
chao1_rhizo<-filter(chao1,compartment=="rhizosphere")
### 3.1.2 shannon ###
setwd(wdImport)
SampleData <- read.table("24SampleData1.txt",header=T,na.strings = c("NA"))
shannon <- read.table("shannon.txt",header=T,na.strings = c("NA"))
shannon <- merge(SampleData,shannon, by = 'sampleid')
shannon$group<-factor(shannon$group,levels=c("CK_B","CZB_B","NF_B","CK_R","CZB_R","NF_R"))
shannon$treatment<-factor(shannon$treatment,levels=c("CK","CZB","NF"))
shannon$compartment<-factor(shannon$compartment,levels=c("bulk","rhizosphere"))
shannon_bulk<-filter(shannon,compartment=="bulk")
shannon_rhizo<-filter(shannon,compartment=="rhizosphere")
#### 4 all####
#### 4.1 chao1###
## 4.1.1 statistical analysis##
chao1$treatment<-factor(chao1$treatment,levels=c("CK","CZB","NF"))
chao1$compartment<-factor(chao1$compartment,levels=c("bulk","rhizosphere"))
chao1_mean <- aggregate(chao1$chao1, by=list(chao1$treatment, chao1$compartment), FUN=mean)
chao1_sd <- aggregate(chao1$chao1, by=list(chao1$treatment, chao1$compartment), FUN=sd)
chao1_len <- aggregate(chao1$chao1, by=list(chao1$treatment, chao1$compartment), FUN=length)
df_res <- data.frame(chao1_mean, sd=chao1_sd$x, len=chao1_len$x)
colnames(df_res) = c("treatment", "compartment", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#bulk
leveneTest(chao1 ~ treatment, data = chao1_bulk)#p>0.05，则满足方差齐性
shapiro.test(chao1_bulk$chao1)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
chao1_bulk<-chao1_bulk%>%mutate(boxcox_chao1 =BoxCox(chao1_bulk$chao1,lambda="auto"))


leveneTest(boxcox_chao1 ~ treatment, data = chao1_bulk)#p>0.05，则满足方差齐性
shapiro.test(chao1_bulk$boxcox_chao1)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_chao1_bulk<-aov(data=chao1_bulk,chao1~treatment)
summary(aov_model_chao1_bulk)
LSD_model_chao1_bulk<-LSD.test(aov_model_chao1_bulk,"treatment",p.adj = "BH")
LSD_model_chao1_bulk

#Non-parameter test
kruskal.test(chao1~treatment, data = chao1_bulk)
aov_model_chao1_bulk<-aov(data=chao1_bulk,boxcox_chao1~treatment)
dunnettT3Test(aov_model_chao1_bulk,p.adjust.method = "BH")

#rhizo
leveneTest(chao1 ~ treatment, data = chao1_rhizo)#p>0.05，则满足方差齐性
shapiro.test(chao1_rhizo$chao1)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
chao1_rhizo<-chao1_rhizo%>%mutate(boxcox_chao1 =BoxCox(chao1_rhizo$chao1,lambda="auto"))


leveneTest(boxcox_chao1 ~ treatment, data = chao1_rhizo)#p>0.05，则满足方差齐性
shapiro.test(chao1_rhizo$boxcox_chao1)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_chao1_rhizo<-aov(data=chao1_rhizo,chao1~treatment)
summary(aov_model_chao1_rhizo)
LSD_model_chao1_rhizo<-LSD.test(aov_model_chao1_rhizo,"treatment",p.adj = "BH")
LSD_model_chao1_rhizo

#Non-parameter test
kruskal.test(chao1~treatment, data = chao1_rhizo)
aov_model_chao1_rhizo<-aov(data=chao1_rhizo,boxcox_chao1~treatment)
dunnettT3Test(aov_model_chao1_rhizo,p.adjust.method = "BH")

###test 
wilcox.test(chao1~compartment, data = chao1, var.equal = TRUE)




# 5.1.2 Plot #
chao1_plot <- ggplot(df_res, aes(x=compartment, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "", y = "chao1", fill = "treatment")+
  scale_fill_manual(values=c("#5C4A46","#DA9464","#888e4a"))+
  FacetTheme
chao1_plot
setwd(wdOutput)
getwd()
ggsave(paste("chao1_plot",".pdf",sep=""),
       chao1_plot,device=cairo_pdf,width=90,height=60,dpi = 300,units = "mm")


### 4.2 shannon ###
## 4.1.1 statistical analysis##
shannon$treatment<-factor(shannon$treatment,levels=c("CK","CZB","NF"))
shannon$compartment<-factor(shannon$compartment,levels=c("bulk","rhizosphere"))
shannon_mean <- aggregate(shannon$shannon, by=list(shannon$treatment, shannon$compartment), FUN=mean)
shannon_sd <- aggregate(shannon$shannon, by=list(shannon$treatment, shannon$compartment), FUN=sd)
shannon_len <- aggregate(shannon$shannon, by=list(shannon$treatment, shannon$compartment), FUN=length)
df_res <- data.frame(shannon_mean, sd=shannon_sd$x, len=shannon_len$x)
colnames(df_res) = c("treatment", "compartment", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#bulk
leveneTest(shannon ~ treatment, data = shannon_bulk)#p>0.05，则满足方差齐性
shapiro.test(shannon_bulk$shannon)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
shannon_bulk<-shannon_bulk%>%mutate(boxcox_shannon =BoxCox(shannon_bulk$shannon,lambda="auto"))


leveneTest(boxcox_shannon ~ treatment, data = shannon_bulk)#p>0.05，则满足方差齐性
shapiro.test(shannon_bulk$boxcox_shannon)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_shannon_bulk<-aov(data=shannon_bulk,shannon~treatment)
summary(aov_model_shannon_bulk)
LSD_model_shannon_bulk<-LSD.test(aov_model_shannon_bulk,"treatment",p.adj = "BH")
LSD_model_shannon_bulk

#Non-parameter test
kruskal.test(shannon~treatment, data = shannon_bulk)
aov_model_shannon_bulk<-aov(data=shannon_bulk,boxcox_shannon~treatment)
dunnettT3Test(aov_model_shannon_bulk,p.adjust.method = "BH")

#rhizo
leveneTest(shannon ~ treatment, data = shannon_rhizo)#p>0.05，则满足方差齐性
shapiro.test(shannon_rhizo$shannon)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
shannon_rhizo<-shannon_rhizo%>%mutate(boxcox_shannon =BoxCox(shannon_rhizo$shannon,lambda="auto"))


leveneTest(boxcox_shannon ~ treatment, data = shannon_rhizo)#p>0.05，则满足方差齐性
shapiro.test(shannon_rhizo$boxcox_shannon)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_shannon_rhizo<-aov(data=shannon_rhizo,shannon~treatment)
summary(aov_model_shannon_rhizo)
LSD_model_shannon_rhizo<-LSD.test(aov_model_shannon_rhizo,"treatment",p.adj = "BH")
LSD_model_shannon_rhizo

#Non-parameter test
kruskal.test(shannon~treatment, data = shannon_rhizo)
aov_model_shannon_rhizo<-aov(data=shannon_rhizo,boxcox_shannon~treatment)
dunnettT3Test(aov_model_shannon_rhizo,p.adjust.method = "BH")

###test 
wilcox.test(shannon~compartment, data = shannon, var.equal = TRUE)

# 5.1.2 Plot #
shannon_plot <- ggplot(df_res, aes(x=compartment, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "", y = "shannon", fill = "treatment")+
  scale_fill_manual(values=c("#5C4A46","#DA9464","#888e4a"))+
  FacetTheme
shannon_plot
setwd(wdOutput)
getwd()
ggsave(paste("shannon_plot",".pdf",sep=""),
       shannon_plot,device=cairo_pdf,width=90,height=60,dpi = 300,units = "mm")

### 4.3 all ###
α_plot <- ggarrange(chao1_plot,shannon_plot,
                    ncol = 2, labels = c("A", "B"),common.legend = T,legend = "right")
α_plot
setwd(wdOutput)
getwd()
ggsave(paste("α_plot",".pdf",sep=""),α_plot,
       device=cairo_pdf,width=240,height=120,dpi = 300,units = "mm")