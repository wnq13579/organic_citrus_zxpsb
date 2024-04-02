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
                                text = element_text(family = "Times New Roman"),
                                strip.text = element_text(size=8,hjust = 0.5),
                                plot.title = element_text(size=8,hjust = 0.5),
                                axis.text =element_text(size=8,color = "black"),
                                axis.title =element_text(size=8,color = "black"),
                                legend.text = element_text(size=8,color = "black"),
                                legend.title = element_text(size=8,color = "black"),
                                legend.background = element_blank(),
                                axis.line = element_line(color = "black",size=0.4))#移除整体的边???



wdImport<- c("E:/whx_24samples_filtered1/01.import")
wdOutput <- c("E:/whx_24samples_filtered1/06.env/soil_enzyme_activity")

#### 3. soil_enzyme_activity ####
### 3.1 Import and process data ###
setwd(wdImport)
SampleData <- read.table("24SampleData.txt",header=T,na.strings = c("NA"))
soil_enzyme_activity<- read.table("soil_enzyme_activity.txt",header=T,na.strings = c("NA"))
soil_enzyme_activity <- merge(SampleData, soil_enzyme_activity, by = 'sampleid')

ACP<-soil_enzyme_activity
ACP$group<-factor(ACP$group,levels=c("CK_B","CZB_B","NF_B","CK_R","CZB_R","NF_R"))
ACP$treatment<-factor(ACP$treatment,levels=c("CK","CZB","NF"))
ACP$compartment<-factor(ACP$compartment,levels=c("bulk","rhizosphere"))
ACP_bulk<-filter(ACP,compartment=="bulk")
ACP_rhizo<-filter(ACP,compartment=="rhizosphere")
### 3.1.2 ALP ###
setwd(wdImport)
ALP<-soil_enzyme_activity
ALP$group<-factor(ALP$group,levels=c("CK_B","CZB_B","NF_B","CK_R","CZB_R","NF_R"))
ALP$treatment<-factor(ALP$treatment,levels=c("CK","CZB","NF"))
ALP$compartment<-factor(ALP$compartment,levels=c("bulk","rhizosphere"))
ALP_bulk<-filter(ALP,compartment=="bulk")
ALP_rhizo<-filter(ALP,compartment=="rhizosphere")
####4 ACP####
### 4.1 ACP ###
## 4.1.1 statistical analysis##
ACP_mean <- aggregate(ACP$ACP, by=list(ACP$treatment, ACP$compartment), FUN=mean)
ACP_sd <- aggregate(ACP$ACP, by=list(ACP$treatment, ACP$compartment), FUN=sd)
ACP_len <- aggregate(ACP$ACP, by=list(ACP$treatment, ACP$compartment), FUN=length)
df_res <- data.frame(ACP_mean, sd=ACP_sd$x, len=ACP_len$x)
colnames(df_res) = c("treatment", "compartment", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#bulk
leveneTest(ACP ~ treatment, data = ACP_bulk)#p>0.05，则满足方差齐性
shapiro.test(ACP_bulk$ACP)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
ACP_bulk<-ACP_bulk%>%mutate(boxcox_ACP =BoxCox(ACP_bulk$ACP,lambda="auto"))


leveneTest(boxcox_ACP ~ treatment, data = ACP_bulk)#p>0.05，则满足方差齐性
shapiro.test(ACP_bulk$boxcox_ACP)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_ACP_bulk<-aov(data=ACP_bulk,ACP~treatment)
summary(aov_model_ACP_bulk)
LSD_model_ACP_bulk<-LSD.test(aov_model_ACP_bulk,"treatment",p.adj = "BH")
LSD_model_ACP_bulk
#DUNCAN
aov_model_ACP_bulk<-aov(data=ACP_bulk,ACP~treatment1)
summary(aov_model_ACP_bulk)
compare_means(data=ACP_bulk,ACP~treatment1,method = "anova")
duncan_result_ACP_bulk<- duncan.test(aov_model_ACP_bulk,"treatment1")
duncan_result_ACP_bulk
#Non-parameter test
kruskal.test(ACP~treatment, data = ACP_bulk)
aov_model_ACP_bulk<-aov(data=ACP_bulk,boxcox_ACP~treatment)
dunnettT3Test(aov_model_ACP_bulk,p.adjust.method = "BH")

#rhizo
leveneTest(ACP ~ treatment, data = ACP_rhizo)#p>0.05，则满足方差齐性
shapiro.test(ACP_rhizo$ACP)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
ACP_rhizo<-ACP_rhizo%>%mutate(boxcox_ACP =BoxCox(ACP_rhizo$ACP,lambda="auto"))


leveneTest(boxcox_ACP ~ treatment, data = ACP_rhizo)#p>0.05，则满足方差齐性
shapiro.test(ACP_rhizo$boxcox_ACP)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_ACP_rhizo<-aov(data=ACP_rhizo,ACP~treatment)
summary(aov_model_ACP_rhizo)
LSD_model_ACP_rhizo<-LSD.test(aov_model_ACP_rhizo,"treatment",p.adj = "BH")
LSD_model_ACP_rhizo
#DUNCAN
aov_model_ACP_rhizo<-aov(data=ACP_rhizo,ACP~treatment1)
summary(aov_model_ACP_rhizo)
compare_means(data=ACP_rhizo,ACP~treatment1,method = "anova")
duncan_result_ACP_rhizo<- duncan.test(aov_model_ACP_rhizo,"treatment1")
duncan_result_ACP_rhizo
#Non-parameter test
kruskal.test(ACP~treatment, data = ACP_rhizo)
aov_model_ACP_rhizo<-aov(data=ACP_rhizo,boxcox_ACP~treatment)
dunnettT3Test(aov_model_ACP_rhizo,p.adjust.method = "BH")
###t.test 
var.test(ACP_bulk$ACP,ACP_rhizo$ACP)#>0.05表示方差齐性
t.test(ACP_bulk$ACP,ACP_rhizo$ACP,var.equal = T)#参数检验
compare_means(ACP~compartment, data = ACP, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验
# 5.1.2 Plot #
ACP_plot <- ggplot(df_res, aes(x=compartment, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "", y = "ACP", fill = "treatment")+
  scale_fill_manual(values=c("#5C4A46","#DA9464","#888e4a"))+
  FacetTheme
ACP_plot
#setwd(wdOutput)
#getwd()
#ggsave(paste("ACP_plot",".pdf",sep=""),
       #ACP_plot,device=cairo_pdf,width=90,height=60,dpi = 300,units = "mm")


#### 5 ALP ####
## 4.1.1 statistical analysis##
ALP$treatment<-factor(ALP$treatment,levels=c("CK","CZB","NF"))
ALP$compartment<-factor(ALP$compartment,levels=c("bulk","rhizosphere"))
ALP_mean <- aggregate(ALP$ALP, by=list(ALP$treatment, ALP$compartment), FUN=mean)
ALP_sd <- aggregate(ALP$ALP, by=list(ALP$treatment, ALP$compartment), FUN=sd)
ALP_len <- aggregate(ALP$ALP, by=list(ALP$treatment, ALP$compartment), FUN=length)
df_res <- data.frame(ALP_mean, sd=ALP_sd$x, len=ALP_len$x)
colnames(df_res) = c("treatment", "compartment", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#bulk
leveneTest(ALP ~ treatment, data = ALP_bulk)#p>0.05，则满足方差齐性
shapiro.test(ALP_bulk$ALP)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
ALP_bulk<-ALP_bulk%>%mutate(boxcox_ALP =BoxCox(ALP_bulk$ALP,lambda="auto"))


leveneTest(boxcox_ALP ~ treatment, data = ALP_bulk)#p>0.05，则满足方差齐性
shapiro.test(ALP_bulk$boxcox_ALP)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_ALP_bulk<-aov(data=ALP_bulk,ALP~treatment)
summary(aov_model_ALP_bulk)
LSD_model_ALP_bulk<-LSD.test(aov_model_ALP_bulk,"treatment",p.adj = "BH")
LSD_model_ALP_bulk
#DUNCAN
aov_model_ALP_bulk<-aov(data=ALP_bulk,ALP~treatment1)
summary(aov_model_ALP_bulk)
compare_means(data=ALP_bulk,ALP~treatment1,method = "anova")
duncan_result_ALP_bulk<- duncan.test(aov_model_ALP_bulk,"treatment1")
duncan_result_ALP_bulk
#Non-parameter test
kruskal.test(ALP~treatment, data = ALP_bulk)
aov_model_ALP_bulk<-aov(data=ALP_bulk,boxcox_ALP~treatment)
dunnettT3Test(aov_model_ALP_bulk,p.adjust.method = "BH")

#rhizo
leveneTest(ALP ~ treatment, data = ALP_rhizo)#p>0.05，则满足方差齐性
shapiro.test(ALP_rhizo$ALP)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
ALP_rhizo<-ALP_rhizo%>%mutate(boxcox_ALP =BoxCox(ALP_rhizo$ALP,lambda="auto"))


leveneTest(boxcox_ALP ~ treatment, data = ALP_rhizo)#p>0.05，则满足方差齐性
shapiro.test(ALP_rhizo$boxcox_ALP)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_ALP_rhizo<-aov(data=ALP_rhizo,ALP~treatment)
summary(aov_model_ALP_rhizo)
LSD_model_ALP_rhizo<-LSD.test(aov_model_ALP_rhizo,"treatment",p.adj = "BH")
LSD_model_ALP_rhizo
#DUNCAN
aov_model_ALP_rhizo<-aov(data=ALP_rhizo,ALP~treatment1)
summary(aov_model_ALP_rhizo)
compare_means(data=ALP_rhizo,ALP~treatment1,method = "anova")
duncan_result_ALP_rhizo<- duncan.test(aov_model_ALP_rhizo,"treatment1")
duncan_result_ALP_rhizo
#Non-parameter test
kruskal.test(ALP~treatment, data = ALP_rhizo)
aov_model_ALP_rhizo<-aov(data=ALP_rhizo,boxcox_ALP~treatment)
dunnettT3Test(aov_model_ALP_rhizo,p.adjust.method = "BH")
###t.test 
var.test(ALP_bulk$ALP,ALP_rhizo$ALP)#>0.05表示方差齐性
t.test(ALP_bulk$ALP,ALP_rhizo$ALP,var.equal = T)#参数检验
compare_means(ALP~compartment, data = ALP, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验


# 5.1.2 Plot #
ALP_plot <- ggplot(df_res, aes(x=compartment, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "", y = "ALP", fill = "treatment")+
  scale_fill_manual(values=c("#5C4A46","#DA9464","#888e4a"))+
  FacetTheme
ALP_plot
#setwd(wdOutput)
#getwd()
#ggsave(paste("ALP_plot",".pdf",sep=""),
      # ALP_plot,device=cairo_pdf,width=90,height=60,dpi = 300,units = "mm")

### 4.3 all ###
α_plot <- ggarrange(ACP_plot,ALP_plot,
                    ncol = 2, labels = c("A", "B"),common.legend = T,legend = "right")
α_plot
setwd(wdOutput)
getwd()
ggsave(paste("α_plot",".pdf",sep=""),α_plot,
       device=cairo_pdf,width=240,height=120,dpi = 300,units = "mm")