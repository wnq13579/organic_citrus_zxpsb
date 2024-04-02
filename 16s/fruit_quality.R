#### 1. Loading package ####
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
wdOutput <- c("E:/whx_24samples_filtered1/06.env/fruit_quality")

#### 3. fruit ####
### 3.1 Import and process data ###
setwd(wdImport)
SampleData <- read.table("24SampleData.txt",header=T,na.strings = c("NA"))
fruit_quality<- read.table("fruit_quality.txt",header=T,na.strings = c("NA"))
fruit_quality <- merge(SampleData, fruit_quality, by = 'sampleid')
fruit_quality$treatment<-factor(fruit_quality$treatment,levels=c("NF","CZB","CK"))
fruit_quality$time<-factor(fruit_quality$time,levels=c("2022","2021"))
fruit_quality$treatment1<-factor(fruit_quality$treatment1,levels=c("CK_21","CZB_21","NF_21","CK_22","CZB_22","NF_22"))
fruit_quality_2021<-filter(fruit_quality,time=="2021")
fruit_quality_2022<-filter(fruit_quality,time=="2022")
#### 4 single_fruit_weight####
### 4.1 single_fruit_weight ###
## 4.1.1 statistical analysis##
single_fruit_weight_mean <- aggregate(fruit_quality$single_fruit_weight, by=list(fruit_quality$treatment, fruit_quality$time), FUN=mean)
single_fruit_weight_sd <- aggregate(fruit_quality$single_fruit_weight, by=list(fruit_quality$treatment, fruit_quality$time), FUN=sd)
single_fruit_weight_len <- aggregate(fruit_quality$single_fruit_weight, by=list(fruit_quality$treatment, fruit_quality$time), FUN=length)
df_res <- data.frame(single_fruit_weight_mean, sd=single_fruit_weight_sd$x, len=single_fruit_weight_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(single_fruit_weight ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$single_fruit_weight)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2021<-fruit_quality_2021%>%mutate(boxcox_single_fruit_weight =BoxCox(fruit_quality_2021$single_fruit_weight,lambda="auto"))


leveneTest(boxcox_single_fruit_weight ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$boxcox_single_fruit_weight)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,single_fruit_weight~treatment1)
summary(aov_model_fruit_quality_2021)
LSD_model_fruit_quality_2021<-LSD.test(aov_model_fruit_quality_2021,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2021
#DUNCAN
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,single_fruit_weight~treatment1)
summary(aov_model_fruit_quality_2021)
compare_means(data=fruit_quality_2021,single_fruit_weight~treatment1,method = "anova")
duncan_result_soil_single_fruit_weight<- duncan.test(aov_model_fruit_quality_2021,"treatment1")
duncan_result_soil_single_fruit_weight
#Non-parameter test
kruskal.test(single_fruit_weight~treatment1, data = fruit_quality_2021)
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,boxcox_single_fruit_weight~treatment1)
dunnettT3Test(aov_model_fruit_quality_2021,p.adjust.method = "BH")

#2022
leveneTest(single_fruit_weight ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$single_fruit_weight)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2022<-fruit_quality_2022%>%mutate(boxcox_single_fruit_weight =BoxCox(fruit_quality_2022$single_fruit_weight,lambda="auto"))


leveneTest(boxcox_single_fruit_weight ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$boxcox_single_fruit_weight)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,single_fruit_weight~treatment1)
summary(aov_model_fruit_quality_2022)
LSD_model_fruit_quality_2022<-LSD.test(aov_model_fruit_quality_2022,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2022
#DUNCAN
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,single_fruit_weight~treatment1)
summary(aov_model_fruit_quality_2022)
compare_means(data=fruit_quality_2022,single_fruit_weight~treatment1,method = "anova")
duncan_result_soil_single_fruit_weight<- duncan.test(aov_model_fruit_quality_2022,"treatment1")
duncan_result_soil_single_fruit_weight
#Non-parameter test
kruskal.test(single_fruit_weight~treatment1, data = fruit_quality_2022)
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,boxcox_single_fruit_weight~treatment1)
dunnettT3Test(aov_model_fruit_quality_2022,p.adjust.method = "BH")

###t.test 
var.test(fruit_quality_2021$single_fruit_weight,fruit_quality_2022$single_fruit_weight)#>0.05表示方差齐性
t.test(fruit_quality_2021$single_fruit_weight,fruit_quality_2022$single_fruit_weight,var.equal = T)#参数检验
compare_means(single_fruit_weight~time, data = fruit_quality, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验

# 5.1.2 Plot #
single_fruit_weight_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "单果重/g Single fruit weight", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
single_fruit_weight_plot

#### 5 edible_rate####
### 5.1 edible_rate ###
## 4.1.1 statistical analysis##
edible_rate_mean <- aggregate(fruit_quality$edible_rate, by=list(fruit_quality$treatment, fruit_quality$time), FUN=mean)
edible_rate_sd <- aggregate(fruit_quality$edible_rate, by=list(fruit_quality$treatment, fruit_quality$time), FUN=sd)
edible_rate_len <- aggregate(fruit_quality$edible_rate, by=list(fruit_quality$treatment, fruit_quality$time), FUN=length)
df_res <- data.frame(edible_rate_mean, sd=edible_rate_sd$x, len=edible_rate_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)

#2021
leveneTest(edible_rate ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$edible_rate)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2021<-fruit_quality_2021%>%mutate(boxcox_edible_rate =BoxCox(fruit_quality_2021$edible_rate,lambda="auto"))


leveneTest(boxcox_edible_rate ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$boxcox_edible_rate)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,edible_rate~treatment1)
summary(aov_model_fruit_quality_2021)
LSD_model_fruit_quality_2021<-LSD.test(aov_model_fruit_quality_2021,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2021
#DUNCAN
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,edible_rate~treatment1)
summary(aov_model_fruit_quality_2021)
compare_means(data=fruit_quality_2021,edible_rate~treatment1,method = "anova")
duncan_result_soil_edible_rate<- duncan.test(aov_model_fruit_quality_2021,"treatment1")
duncan_result_soil_edible_rate
#Non-parameter test
kruskal.test(edible_rate~treatment1, data = fruit_quality_2021)
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,boxcox_edible_rate~treatment1)
dunnettT3Test(aov_model_fruit_quality_2021,p.adjust.method = "BH")

#2022
leveneTest(edible_rate ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$edible_rate)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2022<-fruit_quality_2022%>%mutate(boxcox_edible_rate =BoxCox(fruit_quality_2022$edible_rate,lambda="auto"))


leveneTest(boxcox_edible_rate ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$boxcox_edible_rate)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,edible_rate~treatment1)
summary(aov_model_fruit_quality_2022)
LSD_model_fruit_quality_2022<-LSD.test(aov_model_fruit_quality_2022,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2022
#DUNCAN
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,edible_rate~treatment1)
summary(aov_model_fruit_quality_2022)
compare_means(data=fruit_quality_2022,edible_rate~treatment1,method = "anova")
duncan_result_soil_edible_rate<- duncan.test(aov_model_fruit_quality_2022,"treatment1")
duncan_result_soil_edible_rate
#Non-parameter test
kruskal.test(edible_rate~treatment1, data = fruit_quality_2022)
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,boxcox_edible_rate~treatment1)
dunnettT3Test(aov_model_fruit_quality_2022,p.adjust.method = "BH")


###t.test 
var.test(fruit_quality_2021$edible_rate,fruit_quality_2022$edible_rate)#>0.05表示方差齐性
t.test(fruit_quality_2021$edible_rate,fruit_quality_2022$edible_rate,var.equal = T)#参数检验
compare_means(edible_rate~time, data = fruit_quality, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验

# 5.1.2 Plot #
edible_rate_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "可食率/% Edible rate", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
edible_rate_plot
#### 6 juice_rate####
### 4.1 juice_rate ###
## 4.1.1 statistical analysis##
juice_rate_mean <- aggregate(fruit_quality$juice_rate, by=list(fruit_quality$treatment, fruit_quality$time), FUN=mean)
juice_rate_sd <- aggregate(fruit_quality$juice_rate, by=list(fruit_quality$treatment, fruit_quality$time), FUN=sd)
juice_rate_len <- aggregate(fruit_quality$juice_rate, by=list(fruit_quality$treatment, fruit_quality$time), FUN=length)
df_res <- data.frame(juice_rate_mean, sd=juice_rate_sd$x, len=juice_rate_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(juice_rate ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$juice_rate)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2021<-fruit_quality_2021%>%mutate(boxcox_juice_rate =BoxCox(fruit_quality_2021$juice_rate,lambda="auto"))


leveneTest(boxcox_juice_rate ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$boxcox_juice_rate)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,juice_rate~treatment1)
summary(aov_model_fruit_quality_2021)
LSD_model_fruit_quality_2021<-LSD.test(aov_model_fruit_quality_2021,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2021
#DUNCAN
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,juice_rate~treatment1)
summary(aov_model_fruit_quality_2021)
compare_means(data=fruit_quality_2021,juice_rate~treatment1,method = "anova")
duncan_result_soil_juice_rate<- duncan.test(aov_model_fruit_quality_2021,"treatment1")
duncan_result_soil_juice_rate
#Non-parameter test
kruskal.test(juice_rate~treatment1, data = fruit_quality_2021)
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,boxcox_juice_rate~treatment1)
dunnettT3Test(aov_model_fruit_quality_2021,p.adjust.method = "BH")

#2022
leveneTest(juice_rate ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$juice_rate)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2022<-fruit_quality_2022%>%mutate(boxcox_juice_rate =BoxCox(fruit_quality_2022$juice_rate,lambda="auto"))


leveneTest(boxcox_juice_rate ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$boxcox_juice_rate)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,juice_rate~treatment1)
summary(aov_model_fruit_quality_2022)
LSD_model_fruit_quality_2022<-LSD.test(aov_model_fruit_quality_2022,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2022
#DUNCAN
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,juice_rate~treatment1)
summary(aov_model_fruit_quality_2022)
compare_means(data=fruit_quality_2022,juice_rate~treatment1,method = "anova")
duncan_result_soil_juice_rate<- duncan.test(aov_model_fruit_quality_2022,"treatment1")
duncan_result_soil_juice_rate
#Non-parameter test
kruskal.test(juice_rate~treatment1, data = fruit_quality_2022)
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,boxcox_juice_rate~treatment1)
dunnettT3Test(aov_model_fruit_quality_2022,p.adjust.method = "BH")
###t.test 
var.test(fruit_quality_2021$juice_rate,fruit_quality_2022$juice_rate)#>0.05表示方差齐性
t.test(fruit_quality_2021$juice_rate,fruit_quality_2022$juice_rate,var.equal = T)#参数检验
compare_means(juice_rate~time, data = fruit_quality, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验
# 5.1.2 Plot #
juice_rate_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "果汁率/% Juice_rate", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
juice_rate_plot


####7 TSS  ####
### 4.1 TSS ###
## 4.1.1 statistical analysis##
TSS_mean <- aggregate(fruit_quality$TSS, by=list(fruit_quality$treatment, fruit_quality$time), FUN=mean)
TSS_sd <- aggregate(fruit_quality$TSS, by=list(fruit_quality$treatment, fruit_quality$time), FUN=sd)
TSS_len <- aggregate(fruit_quality$TSS, by=list(fruit_quality$treatment, fruit_quality$time), FUN=length)
df_res <- data.frame(TSS_mean, sd=TSS_sd$x, len=TSS_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(TSS ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$TSS)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2021<-fruit_quality_2021%>%mutate(boxcox_TSS =BoxCox(fruit_quality_2021$TSS,lambda="auto"))


leveneTest(boxcox_TSS ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$boxcox_TSS)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,TSS~treatment1)
summary(aov_model_fruit_quality_2021)
LSD_model_fruit_quality_2021<-LSD.test(aov_model_fruit_quality_2021,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2021
#DUNCAN
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,TSS~treatment1)
summary(aov_model_fruit_quality_2021)
compare_means(data=fruit_quality_2021,TSS~treatment1,method = "anova")
duncan_result_soil_TSS<- duncan.test(aov_model_fruit_quality_2021,"treatment1")
duncan_result_soil_TSS
#Non-parameter test
kruskal.test(TSS~treatment1, data = fruit_quality_2021)
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,boxcox_TSS~treatment1)
dunnettT3Test(aov_model_fruit_quality_2021,p.adjust.method = "BH")

#2022
leveneTest(TSS ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$TSS)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2022<-fruit_quality_2022%>%mutate(boxcox_TSS =BoxCox(fruit_quality_2022$TSS,lambda="auto"))


leveneTest(boxcox_TSS ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$boxcox_TSS)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,TSS~treatment1)
summary(aov_model_fruit_quality_2022)
LSD_model_fruit_quality_2022<-LSD.test(aov_model_fruit_quality_2022,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2022
#DUNCAN
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,TSS~treatment1)
summary(aov_model_fruit_quality_2022)
compare_means(data=fruit_quality_2022,TSS~treatment1,method = "anova")
duncan_result_soil_TSS<- duncan.test(aov_model_fruit_quality_2022,"treatment1")
duncan_result_soil_TSS
#Non-parameter test
kruskal.test(TSS~treatment1, data = fruit_quality_2022)
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,boxcox_TSS~treatment1)
dunnettT3Test(aov_model_fruit_quality_2022,p.adjust.method = "BH")
###t.test 
var.test(fruit_quality_2021$TSS,fruit_quality_2022$TSS)#>0.05表示方差齐性
t.test(fruit_quality_2021$TSS,fruit_quality_2022$TSS,var.equal = T)#参数检验
compare_means(TSS~time, data = fruit_quality, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验
# 5.1.2 Plot #
TSS_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "可溶性固形物/% Soluble solids content", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
TSS_plot
#### 8 TA ####
### 4.1 TA ###
## 4.1.1 statistical analysis##
TA_mean <- aggregate(fruit_quality$TA, by=list(fruit_quality$treatment, fruit_quality$time), FUN=mean)
TA_sd <- aggregate(fruit_quality$TA, by=list(fruit_quality$treatment, fruit_quality$time), FUN=sd)
TA_len <- aggregate(fruit_quality$TA, by=list(fruit_quality$treatment, fruit_quality$time), FUN=length)
df_res <- data.frame(TA_mean, sd=TA_sd$x, len=TA_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(TA ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$TA)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2021<-fruit_quality_2021%>%mutate(boxcox_TA =BoxCox(fruit_quality_2021$TA,lambda="auto"))


leveneTest(boxcox_TA ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$boxcox_TA)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,TA~treatment1)
summary(aov_model_fruit_quality_2021)
LSD_model_fruit_quality_2021<-LSD.test(aov_model_fruit_quality_2021,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2021
#DUNCAN
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,TA~treatment1)
summary(aov_model_fruit_quality_2021)
compare_means(data=fruit_quality_2021,TA~treatment1,method = "anova")
duncan_result_soil_TA<- duncan.test(aov_model_fruit_quality_2021,"treatment1")
duncan_result_soil_TA
#Non-parameter test
kruskal.test(TA~treatment1, data = fruit_quality_2021)
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,boxcox_TA~treatment1)
dunnettT3Test(aov_model_fruit_quality_2021,p.adjust.method = "BH")

#2022
leveneTest(TA ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$TA)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2022<-fruit_quality_2022%>%mutate(boxcox_TA =BoxCox(fruit_quality_2022$TA,lambda="auto"))


leveneTest(boxcox_TA ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$boxcox_TA)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,TA~treatment1)
summary(aov_model_fruit_quality_2022)
LSD_model_fruit_quality_2022<-LSD.test(aov_model_fruit_quality_2022,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2022
#DUNCAN
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,TA~treatment1)
summary(aov_model_fruit_quality_2022)
compare_means(data=fruit_quality_2022,TA~treatment1,method = "anova")
duncan_result_soil_TA<- duncan.test(aov_model_fruit_quality_2022,"treatment1")
duncan_result_soil_TA
#Non-parameter test
kruskal.test(TA~treatment1, data = fruit_quality_2022)
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,boxcox_TA~treatment1)
dunnettT3Test(aov_model_fruit_quality_2022,p.adjust.method = "BH")
###t.test 
var.test(fruit_quality_2021$TA,fruit_quality_2022$TA)#>0.05表示方差齐性
t.test(fruit_quality_2021$TA,fruit_quality_2022$TA,var.equal = T)#参数检验
compare_means(TA~time, data = fruit_quality, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验
# 5.1.2 Plot #
TA_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "可滴定酸/% Titratable acid content", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
TA_plot
####9 VC  ####
### 4.1 VC ###
## 4.1.1 statistical analysis##
VC_mean <- aggregate(fruit_quality$VC, by=list(fruit_quality$treatment, fruit_quality$time), FUN=mean)
VC_sd <- aggregate(fruit_quality$VC, by=list(fruit_quality$treatment, fruit_quality$time), FUN=sd)
VC_len <- aggregate(fruit_quality$VC, by=list(fruit_quality$treatment, fruit_quality$time), FUN=length)
df_res <- data.frame(VC_mean, sd=VC_sd$x, len=VC_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(VC ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$VC)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2021<-fruit_quality_2021%>%mutate(boxcox_VC =BoxCox(fruit_quality_2021$VC,lambda="auto"))


leveneTest(boxcox_VC ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$boxcox_VC)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,VC~treatment1)
summary(aov_model_fruit_quality_2021)
LSD_model_fruit_quality_2021<-LSD.test(aov_model_fruit_quality_2021,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2021
#DUNCAN
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,VC~treatment1)
summary(aov_model_fruit_quality_2021)
compare_means(data=fruit_quality_2021,VC~treatment1,method = "anova")
duncan_result_soil_VC<- duncan.test(aov_model_fruit_quality_2021,"treatment1")
duncan_result_soil_VC
#Non-parameter test
kruskal.test(VC~treatment1, data = fruit_quality_2021)
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,boxcox_VC~treatment1)
dunnettT3Test(aov_model_fruit_quality_2021,p.adjust.method = "BH")

#2022
leveneTest(VC ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$VC)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2022<-fruit_quality_2022%>%mutate(boxcox_VC =BoxCox(fruit_quality_2022$VC,lambda="auto"))


leveneTest(boxcox_VC ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$boxcox_VC)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,VC~treatment1)
summary(aov_model_fruit_quality_2022)
LSD_model_fruit_quality_2022<-LSD.test(aov_model_fruit_quality_2022,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2022
#DUNCAN
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,VC~treatment1)
summary(aov_model_fruit_quality_2022)
compare_means(data=fruit_quality_2022,VC~treatment1,method = "anova")
duncan_result_soil_VC<- duncan.test(aov_model_fruit_quality_2022,"treatment1")
duncan_result_soil_VC
#Non-parameter test
kruskal.test(VC~treatment1, data = fruit_quality_2022)
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,boxcox_VC~treatment1)
dunnettT3Test(aov_model_fruit_quality_2022,p.adjust.method = "BH")
###t.test 
var.test(fruit_quality_2021$VC,fruit_quality_2022$VC)#>0.05表示方差齐性
t.test(fruit_quality_2021$VC,fruit_quality_2022$VC,var.equal = T)#参数检验
compare_means(VC~time, data = fruit_quality, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验
# 5.1.2 Plot #
VC_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "维生素C/mg·kg-1 Vitamin C content", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
VC_plot
####10 P_L####
### 4.1 P_L ###
## 4.1.1 statistical analysis##
P_L_mean <- aggregate(fruit_quality$P_L, by=list(fruit_quality$treatment, fruit_quality$time), FUN=mean)
P_L_sd <- aggregate(fruit_quality$P_L, by=list(fruit_quality$treatment, fruit_quality$time), FUN=sd)
P_L_len <- aggregate(fruit_quality$P_L, by=list(fruit_quality$treatment, fruit_quality$time), FUN=length)
df_res <- data.frame(P_L_mean, sd=P_L_sd$x, len=P_L_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(P_L ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$P_L)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2021<-fruit_quality_2021%>%mutate(boxcox_P_L =BoxCox(fruit_quality_2021$P_L,lambda="auto"))


leveneTest(boxcox_P_L ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$boxcox_P_L)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,P_L~treatment1)
summary(aov_model_fruit_quality_2021)
LSD_model_fruit_quality_2021<-LSD.test(aov_model_fruit_quality_2021,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2021
#DUNCAN
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,P_L~treatment1)
summary(aov_model_fruit_quality_2021)
compare_means(data=fruit_quality_2021,P_L~treatment1,method = "anova")
duncan_result_soil_P_L<- duncan.test(aov_model_fruit_quality_2021,"treatment1")
duncan_result_soil_P_L
#Non-parameter test
kruskal.test(P_L~treatment1, data = fruit_quality_2021)
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,boxcox_P_L~treatment1)
dunnettT3Test(aov_model_fruit_quality_2021,p.adjust.method = "BH")

#2022
leveneTest(P_L ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$P_L)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2022<-fruit_quality_2022%>%mutate(boxcox_P_L =BoxCox(fruit_quality_2022$P_L,lambda="auto"))


leveneTest(boxcox_P_L ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$boxcox_P_L)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,P_L~treatment1)
summary(aov_model_fruit_quality_2022)
LSD_model_fruit_quality_2022<-LSD.test(aov_model_fruit_quality_2022,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2022
#DUNCAN
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,P_L~treatment1)
summary(aov_model_fruit_quality_2022)
compare_means(data=fruit_quality_2022,P_L~treatment1,method = "anova")
duncan_result_soil_P_L<- duncan.test(aov_model_fruit_quality_2022,"treatment1")
duncan_result_soil_P_L
#Non-parameter test
kruskal.test(P_L~treatment1, data = fruit_quality_2022)
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,boxcox_P_L~treatment1)
dunnettT3Test(aov_model_fruit_quality_2022,p.adjust.method = "BH")
###t.test 
var.test(fruit_quality_2021$P_L,fruit_quality_2022$P_L)#>0.05表示方差齐性
t.test(fruit_quality_2021$P_L,fruit_quality_2022$P_L,var.equal = T)#参数检验
compare_means(P_L~time, data = fruit_quality, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验
# 5.1.2 Plot #
P_L_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
    labs(x = "Year", y = "L（果皮） L(Peel)", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
P_L_plot
#### 11 P_a  ####
### 4.1 P_a ###
## 4.1.1 statistical analysis##
P_a_mean <- aggregate(fruit_quality$P_a, by=list(fruit_quality$treatment, fruit_quality$time), FUN=mean)
P_a_sd <- aggregate(fruit_quality$P_a, by=list(fruit_quality$treatment, fruit_quality$time), FUN=sd)
P_a_len <- aggregate(fruit_quality$P_a, by=list(fruit_quality$treatment, fruit_quality$time), FUN=length)
df_res <- data.frame(P_a_mean, sd=P_a_sd$x, len=P_a_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(P_a ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$P_a)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2021<-fruit_quality_2021%>%mutate(boxcox_P_a =BoxCox(fruit_quality_2021$P_a,lambda="auto"))


leveneTest(boxcox_P_a ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$boxcox_P_a)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,P_a~treatment1)
summary(aov_model_fruit_quality_2021)
LSD_model_fruit_quality_2021<-LSD.test(aov_model_fruit_quality_2021,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2021
#DUNCAN
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,P_a~treatment1)
summary(aov_model_fruit_quality_2021)
compare_means(data=fruit_quality_2021,P_a~treatment1,method = "anova")
duncan_result_soil_P_a<- duncan.test(aov_model_fruit_quality_2021,"treatment1")
duncan_result_soil_P_a
#Non-parameter test
kruskal.test(P_a~treatment1, data = fruit_quality_2021)
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,boxcox_P_a~treatment1)
dunnettT3Test(aov_model_fruit_quality_2021,p.adjust.method = "BH")

#2022
leveneTest(P_a ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$P_a)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2022<-fruit_quality_2022%>%mutate(boxcox_P_a =BoxCox(fruit_quality_2022$P_a,lambda="auto"))


leveneTest(boxcox_P_a ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$boxcox_P_a)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,P_a~treatment1)
summary(aov_model_fruit_quality_2022)
LSD_model_fruit_quality_2022<-LSD.test(aov_model_fruit_quality_2022,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2022
#DUNCAN
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,P_a~treatment1)
summary(aov_model_fruit_quality_2022)
compare_means(data=fruit_quality_2022,P_a~treatment1,method = "anova")
duncan_result_soil_P_a<- duncan.test(aov_model_fruit_quality_2022,"treatment1")
duncan_result_soil_P_a
#Non-parameter test
kruskal.test(P_a~treatment1, data = fruit_quality_2022)
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,boxcox_P_a~treatment1)
dunnettT3Test(aov_model_fruit_quality_2022,p.adjust.method = "BH")
###t.test 
var.test(fruit_quality_2021$P_a,fruit_quality_2022$P_a)#>0.05表示方差齐性
t.test(fruit_quality_2021$P_a,fruit_quality_2022$P_a,var.equal = T)#参数检验
compare_means(P_a~time, data = fruit_quality, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验
# 5.1.2 Plot #
P_a_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "a*（果皮） a*(Peel)", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
P_a_plot
#### 12  P_b ####
### 4.1 P_b ###
## 4.1.1 statistical analysis##
P_b_mean <- aggregate(fruit_quality$P_b, by=list(fruit_quality$treatment, fruit_quality$time), FUN=mean)
P_b_sd <- aggregate(fruit_quality$P_b, by=list(fruit_quality$treatment, fruit_quality$time), FUN=sd)
P_b_len <- aggregate(fruit_quality$P_b, by=list(fruit_quality$treatment, fruit_quality$time), FUN=length)
df_res <- data.frame(P_b_mean, sd=P_b_sd$x, len=P_b_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(P_b ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$P_b)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2021<-fruit_quality_2021%>%mutate(boxcox_P_b =BoxCox(fruit_quality_2021$P_b,lambda="auto"))


leveneTest(boxcox_P_b ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$boxcox_P_b)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,P_b~treatment1)
summary(aov_model_fruit_quality_2021)
LSD_model_fruit_quality_2021<-LSD.test(aov_model_fruit_quality_2021,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2021
#DUNCAN
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,P_b~treatment1)
summary(aov_model_fruit_quality_2021)
compare_means(data=fruit_quality_2021,P_b~treatment1,method = "anova")
duncan_result_soil_P_b<- duncan.test(aov_model_fruit_quality_2021,"treatment1")
duncan_result_soil_P_b
#Non-parameter test
kruskal.test(P_b~treatment1, data = fruit_quality_2021)
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,boxcox_P_b~treatment1)
dunnettT3Test(aov_model_fruit_quality_2021,p.adjust.method = "BH")

#2022
leveneTest(P_b ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$P_b)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2022<-fruit_quality_2022%>%mutate(boxcox_P_b =BoxCox(fruit_quality_2022$P_b,lambda="auto"))


leveneTest(boxcox_P_b ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$boxcox_P_b)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,P_b~treatment1)
summary(aov_model_fruit_quality_2022)
LSD_model_fruit_quality_2022<-LSD.test(aov_model_fruit_quality_2022,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2022
#DUNCAN
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,P_b~treatment1)
summary(aov_model_fruit_quality_2022)
compare_means(data=fruit_quality_2022,P_b~treatment1,method = "anova")
duncan_result_soil_P_b<- duncan.test(aov_model_fruit_quality_2022,"treatment1")
duncan_result_soil_P_b
#Non-parameter test
kruskal.test(P_b~treatment1, data = fruit_quality_2022)
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,boxcox_P_b~treatment1)
dunnettT3Test(aov_model_fruit_quality_2022,p.adjust.method = "BH")
###t.test 
var.test(fruit_quality_2021$P_b,fruit_quality_2022$P_b)#>0.05表示方差齐性
t.test(fruit_quality_2021$P_b,fruit_quality_2022$P_b,var.equal = T)#参数检验
compare_means(P_b~time, data = fruit_quality, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验
# 5.1.2 Plot #
P_b_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "b*（果皮） b*(Peel)", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
P_b_plot

#### 13 fruit_index####
### 4.1 fruit_index ###
## 4.1.1 statistical analysis##
fruit_index_mean <- aggregate(fruit_quality$fruit_index, by=list(fruit_quality$treatment, fruit_quality$time), FUN=mean)
fruit_index_sd <- aggregate(fruit_quality$fruit_index, by=list(fruit_quality$treatment, fruit_quality$time), FUN=sd)
fruit_index_len <- aggregate(fruit_quality$fruit_index, by=list(fruit_quality$treatment, fruit_quality$time), FUN=length)
df_res <- data.frame(fruit_index_mean, sd=fruit_index_sd$x, len=fruit_index_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(fruit_index ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$fruit_index)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2021<-fruit_quality_2021%>%mutate(boxcox_fruit_index =BoxCox(fruit_quality_2021$fruit_index,lambda="auto"))


leveneTest(boxcox_fruit_index ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$boxcox_fruit_index)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,fruit_index~treatment1)
summary(aov_model_fruit_quality_2021)
LSD_model_fruit_quality_2021<-LSD.test(aov_model_fruit_quality_2021,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2021
#DUNCAN
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,fruit_index~treatment1)
summary(aov_model_fruit_quality_2021)
compare_means(data=fruit_quality_2021,fruit_index~treatment1,method = "anova")
duncan_result_soil_fruit_index<- duncan.test(aov_model_fruit_quality_2021,"treatment1")
duncan_result_soil_fruit_index
#Non-parameter test
kruskal.test(fruit_index~treatment1, data = fruit_quality_2021)
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,boxcox_fruit_index~treatment1)
dunnettT3Test(aov_model_fruit_quality_2021,p.adjust.method = "BH")

#2022
leveneTest(fruit_index ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$fruit_index)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2022<-fruit_quality_2022%>%mutate(boxcox_fruit_index =BoxCox(fruit_quality_2022$fruit_index,lambda="auto"))


leveneTest(boxcox_fruit_index ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$boxcox_fruit_index)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,fruit_index~treatment1)
summary(aov_model_fruit_quality_2022)
LSD_model_fruit_quality_2022<-LSD.test(aov_model_fruit_quality_2022,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2022
#DUNCAN
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,fruit_index~treatment1)
summary(aov_model_fruit_quality_2022)
compare_means(data=fruit_quality_2022,fruit_index~treatment1,method = "anova")
duncan_result_soil_fruit_index<- duncan.test(aov_model_fruit_quality_2022,"treatment1")
duncan_result_soil_fruit_index
#Non-parameter test
kruskal.test(fruit_index~treatment1, data = fruit_quality_2022)
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,boxcox_fruit_index~treatment1)
dunnettT3Test(aov_model_fruit_quality_2022,p.adjust.method = "BH")
###t.test 
var.test(fruit_quality_2021$fruit_index,fruit_quality_2022$fruit_index)#>0.05表示方差齐性
t.test(fruit_quality_2021$fruit_index,fruit_quality_2022$fruit_index,var.equal = T)#参数检验
compare_means(fruit_index~time, data = fruit_quality, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验
# 5.1.2 Plot #
fruit_index_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "果形指数 Fruit shape index", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
fruit_index_plot
#### 14 solid_to_acid_ratio####
### 4.1 solid_to_acid_ratio ###
## 4.1.1 statistical analysis##
solid_to_acid_ratio_mean <- aggregate(fruit_quality$solid_to_acid_ratio, by=list(fruit_quality$treatment, fruit_quality$time), FUN=mean)
solid_to_acid_ratio_sd <- aggregate(fruit_quality$solid_to_acid_ratio, by=list(fruit_quality$treatment, fruit_quality$time), FUN=sd)
solid_to_acid_ratio_len <- aggregate(fruit_quality$solid_to_acid_ratio, by=list(fruit_quality$treatment, fruit_quality$time), FUN=length)
df_res <- data.frame(solid_to_acid_ratio_mean, sd=solid_to_acid_ratio_sd$x, len=solid_to_acid_ratio_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(solid_to_acid_ratio ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$solid_to_acid_ratio)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2021<-fruit_quality_2021%>%mutate(boxcox_solid_to_acid_ratio =BoxCox(fruit_quality_2021$solid_to_acid_ratio,lambda="auto"))


leveneTest(boxcox_solid_to_acid_ratio ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$boxcox_solid_to_acid_ratio)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,solid_to_acid_ratio~treatment1)
summary(aov_model_fruit_quality_2021)
LSD_model_fruit_quality_2021<-LSD.test(aov_model_fruit_quality_2021,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2021
#DUNCAN
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,solid_to_acid_ratio~treatment1)
summary(aov_model_fruit_quality_2021)
compare_means(data=fruit_quality_2021,solid_to_acid_ratio~treatment1,method = "anova")
duncan_result_soil_solid_to_acid_ratio<- duncan.test(aov_model_fruit_quality_2021,"treatment1")
duncan_result_soil_solid_to_acid_ratio
#Non-parameter test
kruskal.test(solid_to_acid_ratio~treatment1, data = fruit_quality_2021)
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,boxcox_solid_to_acid_ratio~treatment1)
dunnettT3Test(aov_model_fruit_quality_2021,p.adjust.method = "BH")

#2022
leveneTest(solid_to_acid_ratio ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$solid_to_acid_ratio)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2022<-fruit_quality_2022%>%mutate(boxcox_solid_to_acid_ratio =BoxCox(fruit_quality_2022$solid_to_acid_ratio,lambda="auto"))


leveneTest(boxcox_solid_to_acid_ratio ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$boxcox_solid_to_acid_ratio)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,solid_to_acid_ratio~treatment1)
summary(aov_model_fruit_quality_2022)
LSD_model_fruit_quality_2022<-LSD.test(aov_model_fruit_quality_2022,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2022
#DUNCAN
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,solid_to_acid_ratio~treatment1)
summary(aov_model_fruit_quality_2022)
compare_means(data=fruit_quality_2022,solid_to_acid_ratio~treatment1,method = "anova")
duncan_result_soil_solid_to_acid_ratio<- duncan.test(aov_model_fruit_quality_2022,"treatment1")
duncan_result_soil_solid_to_acid_ratio
#Non-parameter test
kruskal.test(solid_to_acid_ratio~treatment1, data = fruit_quality_2022)
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,boxcox_solid_to_acid_ratio~treatment1)
dunnettT3Test(aov_model_fruit_quality_2022,p.adjust.method = "BH")
###t.test 
var.test(fruit_quality_2021$solid_to_acid_ratio,fruit_quality_2022$solid_to_acid_ratio)#>0.05表示方差齐性
t.test(fruit_quality_2021$solid_to_acid_ratio,fruit_quality_2022$solid_to_acid_ratio,var.equal = T)#参数检验
compare_means(solid_to_acid_ratio~time, data = fruit_quality, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验
# 5.1.2 Plot #
solid_to_acid_ratio_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "固酸比 soluble solids to acidity ratio", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
solid_to_acid_ratio_plot
#### 15 yield####
### 4.1 yield ###
## 4.1.1 statistical analysis##
yield_mean <- aggregate(fruit_quality$yield, by=list(fruit_quality$treatment, fruit_quality$time), FUN=mean)
yield_sd <- aggregate(fruit_quality$yield, by=list(fruit_quality$treatment, fruit_quality$time), FUN=sd)
yield_len <- aggregate(fruit_quality$yield, by=list(fruit_quality$treatment, fruit_quality$time), FUN=length)
df_res <- data.frame(yield_mean, sd=yield_sd$x, len=yield_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(yield ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$yield)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2021<-fruit_quality_2021%>%mutate(boxcox_yield =BoxCox(fruit_quality_2021$yield,lambda="auto"))


leveneTest(boxcox_yield ~ treatment1, data = fruit_quality_2021)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2021$boxcox_yield)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,yield~treatment1)
summary(aov_model_fruit_quality_2021)
LSD_model_fruit_quality_2021<-LSD.test(aov_model_fruit_quality_2021,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2021
#DUNCAN
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,yield~treatment1)
summary(aov_model_fruit_quality_2021)
compare_means(data=fruit_quality_2021,yield~treatment1,method = "anova")
duncan_result_soil_yield<- duncan.test(aov_model_fruit_quality_2021,"treatment1")
duncan_result_soil_yield
#Non-parameter test
kruskal.test(yield~treatment1, data = fruit_quality_2021)
aov_model_fruit_quality_2021<-aov(data=fruit_quality_2021,boxcox_yield~treatment1)
dunnettT3Test(aov_model_fruit_quality_2021,p.adjust.method = "BH")

#2022
leveneTest(yield ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$yield)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
fruit_quality_2022<-fruit_quality_2022%>%mutate(boxcox_yield =BoxCox(fruit_quality_2022$yield,lambda="auto"))


leveneTest(boxcox_yield ~ treatment1, data = fruit_quality_2022)#p>0.05，则满足方差齐性
shapiro.test(fruit_quality_2022$boxcox_yield)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,yield~treatment1)
summary(aov_model_fruit_quality_2022)
LSD_model_fruit_quality_2022<-LSD.test(aov_model_fruit_quality_2022,"treatment1",p.adj = "BH")
LSD_model_fruit_quality_2022
#DUNCAN
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,yield~treatment1)
summary(aov_model_fruit_quality_2022)
compare_means(data=fruit_quality_2022,yield~treatment1,method = "anova")
duncan_result_soil_yield<- duncan.test(aov_model_fruit_quality_2022,"treatment1")
duncan_result_soil_yield
#Non-parameter test
kruskal.test(yield~treatment1, data = fruit_quality_2022)
aov_model_fruit_quality_2022<-aov(data=fruit_quality_2022,boxcox_yield~treatment1)
dunnettT3Test(aov_model_fruit_quality_2022,p.adjust.method = "BH")
###t.test 
var.test(fruit_quality_2021$yield,fruit_quality_2022$yield)#>0.05表示方差齐性
t.test(fruit_quality_2021$yield,fruit_quality_2022$yield,var.equal = T)#参数检验
compare_means(yield~time, data = fruit_quality, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验
# 5.1.2 Plot #
yield_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "单株产量/斤 Yield", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
yield_plot

####all####
fruit_quality_plot <- ggarrange(yield_plot,single_fruit_weight_plot,edible_rate_plot,juice_rate_plot,peel_thickness_plot,
                                TSS_plot,TA_plot,VC_plot,P_moisture_content_plot,F_moisture_content_plot,
                                longitudinal_diameter_plot,diameter_plot,hardness_plot,fruit_index_plot,solid_to_acid_ratio_plot,
                                P_L_plot,P_a_plot,P_b_plot,F_L_plot,F_a_plot,
                                F_b_plot,number_of_seeds_plot,
                                 ncol = 5, nrow = 5,common.legend = T,legend = "right")
fruit_quality_plot
setwd(wdOutput)
getwd()
ggsave(paste("fruit_quality_plot",".pdf",sep=""),
       fruit_quality_plot,device=cairo_pdf,width=270,height=210,dpi = 300,units = "mm")