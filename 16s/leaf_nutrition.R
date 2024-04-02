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
wdOutput <- c("E:/whx_24samples_filtered1/06.env/leaf_nutrition")

#### 3. leaf ####
### 3.1 Import and process data ###
setwd(wdImport)
SampleData <- read.table("24SampleData.txt",header=T,na.strings = c("NA"))
leaf_nutrition<- read.table("env_leaf_nutrition.txt",header=T,na.strings = c("NA"))
leaf_nutrition <- merge(SampleData, leaf_nutrition, by = 'sampleid')
leaf_nutrition$treatment<-factor(leaf_nutrition$treatment,levels=c("NF","CZB","CK"))
leaf_nutrition$time<-factor(leaf_nutrition$time,levels=c("2022","2021"))
leaf_nutrition$treatment1<-factor(leaf_nutrition$treatment1,levels=c("CK_21","CZB_21","NF_21","CK_22","CZB_22","NF_22"))
leaf_nutrition_2021<-filter(leaf_nutrition,time=="2021")
leaf_nutrition_2022<-filter(leaf_nutrition,time=="2022")
#### 4 K####
### 4.1 K ###
## 4.1.1 statistical analysis##
K_mean <- aggregate(leaf_nutrition$K, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=mean)
K_sd <- aggregate(leaf_nutrition$K, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=sd)
K_len <- aggregate(leaf_nutrition$K, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=length)
df_res <- data.frame(K_mean, sd=K_sd$x, len=K_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(K ~ treatment1, data = leaf_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2021$K)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
leaf_nutrition_2021<-leaf_nutrition_2021%>%mutate(boxcox_K =BoxCox(leaf_nutrition_2021$K,lambda="auto"))


leveneTest(boxcox_K ~ treatment1, data = leaf_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2021$boxcox_K)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,K~treatment1)
summary(aov_model_leaf_nutrition_2021)
LSD_model_leaf_nutrition_2021<-LSD.test(aov_model_leaf_nutrition_2021,"treatment1",p.adj = "BH")
LSD_model_leaf_nutrition_2021
#DUNCAN
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,K~treatment1)
summary(aov_model_leaf_nutrition_2021)
compare_means(data=leaf_nutrition_2021,K~treatment1,method = "anova")
duncan_result_soil_K<- duncan.test(aov_model_leaf_nutrition_2021,"treatment1")
duncan_result_soil_K
#Non-parameter test
kruskal.test(K~treatment1, data = leaf_nutrition_2021)
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,boxcox_K~treatment1)
dunnettT3Test(aov_model_leaf_nutrition_2021,p.adjust.method = "BH")

#2022
leveneTest(K ~ treatment1, data = leaf_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2022$K)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
leaf_nutrition_2022<-leaf_nutrition_2022%>%mutate(boxcox_K =BoxCox(leaf_nutrition_2022$K,lambda="auto"))


leveneTest(boxcox_K ~ treatment1, data = leaf_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2022$boxcox_K)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,K~treatment1)
summary(aov_model_leaf_nutrition_2022)
LSD_model_leaf_nutrition_2022<-LSD.test(aov_model_leaf_nutrition_2022,"treatment1",p.adj = "BH")
LSD_model_leaf_nutrition_2022
#DUNCAN
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,K~treatment1)
summary(aov_model_leaf_nutrition_2022)
compare_means(data=leaf_nutrition_2022,K~treatment1,method = "anova")
duncan_result_soil_K<- duncan.test(aov_model_leaf_nutrition_2022,"treatment1")
duncan_result_soil_K
#Non-parameter test
kruskal.test(K~treatment1, data = leaf_nutrition_2022)
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,boxcox_K~treatment1)
dunnettT3Test(aov_model_leaf_nutrition_2022,p.adjust.method = "BH")
###t.test 
var.test(leaf_nutrition_2021$K,leaf_nutrition_2022$K)#>0.05表示方差齐性
t.test(leaf_nutrition_2021$K,leaf_nutrition_2022$K,var.equal = T)#参数检验
compare_means(K~time, data = leaf_nutrition, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验
# 5.1.2 Plot #
K_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "k/%", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
K_plot




#### 5 N####
### 4.1 N ###
## 4.1.1 statistical analysis##
N_mean <- aggregate(leaf_nutrition$N, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=mean)
N_sd <- aggregate(leaf_nutrition$N, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=sd)
N_len <- aggregate(leaf_nutrition$N, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=length)
df_res <- data.frame(N_mean, sd=N_sd$x, len=N_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)

#2021
leveneTest(N ~ treatment1, data = leaf_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2021$N)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
leaf_nutrition_2021<-leaf_nutrition_2021%>%mutate(boxcox_N =BoxCox(leaf_nutrition_2021$N,lambda="auto"))


leveneTest(boxcox_N ~ treatment1, data = leaf_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2021$boxcox_N)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,N~treatment1)
summary(aov_model_leaf_nutrition_2021)
LSD_model_leaf_nutrition_2021<-LSD.test(aov_model_leaf_nutrition_2021,"treatment1",p.adj = "BH")
LSD_model_leaf_nutrition_2021
#DUNCAN
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,N~treatment1)
summary(aov_model_leaf_nutrition_2021)
compare_means(data=leaf_nutrition_2021,N~treatment1,method = "anova")
duncan_result_soil_N<- duncan.test(aov_model_leaf_nutrition_2021,"treatment1")
duncan_result_soil_N
#Non-parameter test
kruskal.test(N~treatment1, data = leaf_nutrition_2021)
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,boxcox_N~treatment1)
dunnettT3Test(aov_model_leaf_nutrition_2021,p.adjust.method = "BH")

#2022
leveneTest(N ~ treatment1, data = leaf_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2022$N)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
leaf_nutrition_2022<-leaf_nutrition_2022%>%mutate(boxcox_N =BoxCox(leaf_nutrition_2022$N,lambda="auto"))


leveneTest(boxcox_N ~ treatment1, data = leaf_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2022$boxcox_N)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,N~treatment1)
summary(aov_model_leaf_nutrition_2022)
LSD_model_leaf_nutrition_2022<-LSD.test(aov_model_leaf_nutrition_2022,"treatment1",p.adj = "BH")
LSD_model_leaf_nutrition_2022
#DUNCAN
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,N~treatment1)
summary(aov_model_leaf_nutrition_2022)
compare_means(data=leaf_nutrition_2022,N~treatment1,method = "anova")
duncan_result_soil_N<- duncan.test(aov_model_leaf_nutrition_2022,"treatment1")
duncan_result_soil_N
#Non-parameter test
kruskal.test(N~treatment1, data = leaf_nutrition_2022)
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,boxcox_N~treatment1)
dunnettT3Test(aov_model_leaf_nutrition_2022,p.adjust.method = "BH")

###t.test 
var.test(leaf_nutrition_2021$N,leaf_nutrition_2022$N)#>0.05表示方差齐性
t.test(leaf_nutrition_2021$N,leaf_nutrition_2022$N,var.equal = T)#参数检验
compare_means(N~time, data = leaf_nutrition, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验
# 5.1.2 Plot #
N_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "N/%", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
N_plot




#### 6 P####
### 4.1 P ###
## 4.1.1 statistical analysis##
P_mean <- aggregate(leaf_nutrition$P, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=mean)
P_sd <- aggregate(leaf_nutrition$P, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=sd)
P_len <- aggregate(leaf_nutrition$P, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=length)
df_res <- data.frame(P_mean, sd=P_sd$x, len=P_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(P ~ treatment1, data = leaf_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2021$P)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
leaf_nutrition_2021<-leaf_nutrition_2021%>%mutate(boxcox_P =BoxCox(leaf_nutrition_2021$P,lambda="auto"))


leveneTest(boxcox_P ~ treatment1, data = leaf_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2021$boxcox_P)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,P~treatment1)
summary(aov_model_leaf_nutrition_2021)
LSD_model_leaf_nutrition_2021<-LSD.test(aov_model_leaf_nutrition_2021,"treatment1",p.adj = "BH")
LSD_model_leaf_nutrition_2021
#DUNCAN
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,P~treatment1)
summary(aov_model_leaf_nutrition_2021)
compare_means(data=leaf_nutrition_2021,P~treatment1,method = "anova")
duncan_result_soil_P<- duncan.test(aov_model_leaf_nutrition_2021,"treatment1")
duncan_result_soil_P
#Non-parameter test
kruskal.test(P~treatment1, data = leaf_nutrition_2021)
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,boxcox_P~treatment1)
dunnettT3Test(aov_model_leaf_nutrition_2021,p.adjust.method = "BH")

#2022
leveneTest(P ~ treatment1, data = leaf_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2022$P)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
leaf_nutrition_2022<-leaf_nutrition_2022%>%mutate(boxcox_P =BoxCox(leaf_nutrition_2022$P,lambda="auto"))


leveneTest(boxcox_P ~ treatment1, data = leaf_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2022$boxcox_P)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,P~treatment1)
summary(aov_model_leaf_nutrition_2022)
LSD_model_leaf_nutrition_2022<-LSD.test(aov_model_leaf_nutrition_2022,"treatment1",p.adj = "BH")
LSD_model_leaf_nutrition_2022
#DUNCAN
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,P~treatment1)
summary(aov_model_leaf_nutrition_2022)
compare_means(data=leaf_nutrition_2022,P~treatment1,method = "anova")
duncan_result_soil_P<- duncan.test(aov_model_leaf_nutrition_2022,"treatment1")
duncan_result_soil_P
#Non-parameter test
kruskal.test(P~treatment1, data = leaf_nutrition_2022)
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,boxcox_P~treatment1)
dunnettT3Test(aov_model_leaf_nutrition_2022,p.adjust.method = "BH")

###t.test 
var.test(leaf_nutrition_2021$P,leaf_nutrition_2022$P)#>0.05表示方差齐性
t.test(leaf_nutrition_2021$P,leaf_nutrition_2022$P,var.equal = T)#参数检验
compare_means(P~time, data = leaf_nutrition, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验
# 5.1.2 Plot #
P_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "P/%", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
P_plot




#### 7 Ca####
### 4.1 Ca ###
## 4.1.1 statistical analysis##
Ca_mean <- aggregate(leaf_nutrition$Ca, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=mean)
Ca_sd <- aggregate(leaf_nutrition$Ca, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=sd)
Ca_len <- aggregate(leaf_nutrition$Ca, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=length)
df_res <- data.frame(Ca_mean, sd=Ca_sd$x, len=Ca_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(Ca ~ treatment1, data = leaf_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2021$Ca)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
leaf_nutrition_2021<-leaf_nutrition_2021%>%mutate(boxcox_Ca =BoxCox(leaf_nutrition_2021$Ca,lambda="auto"))


leveneTest(boxcox_Ca ~ treatment1, data = leaf_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2021$boxcox_Ca)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,Ca~treatment1)
summary(aov_model_leaf_nutrition_2021)
LSD_model_leaf_nutrition_2021<-LSD.test(aov_model_leaf_nutrition_2021,"treatment1",p.adj = "BH")
LSD_model_leaf_nutrition_2021
#DUNCAN
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,Ca~treatment1)
summary(aov_model_leaf_nutrition_2021)
compare_means(data=leaf_nutrition_2021,Ca~treatment1,method = "anova")
duncan_result_soil_Ca<- duncan.test(aov_model_leaf_nutrition_2021,"treatment1")
duncan_result_soil_Ca
#Non-parameter test
kruskal.test(Ca~treatment1, data = leaf_nutrition_2021)
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,boxcox_Ca~treatment1)
dunnettT3Test(aov_model_leaf_nutrition_2021,p.adjust.method = "BH")

#2022
leveneTest(Ca ~ treatment1, data = leaf_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2022$Ca)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
leaf_nutrition_2022<-leaf_nutrition_2022%>%mutate(boxcox_Ca =BoxCox(leaf_nutrition_2022$Ca,lambda="auto"))


leveneTest(boxcox_Ca ~ treatment1, data = leaf_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2022$boxcox_Ca)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,Ca~treatment1)
summary(aov_model_leaf_nutrition_2022)
LSD_model_leaf_nutrition_2022<-LSD.test(aov_model_leaf_nutrition_2022,"treatment1",p.adj = "BH")
LSD_model_leaf_nutrition_2022
#DUNCAN
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,Ca~treatment1)
summary(aov_model_leaf_nutrition_2022)
compare_means(data=leaf_nutrition_2022,Ca~treatment1,method = "anova")
duncan_result_soil_Ca<- duncan.test(aov_model_leaf_nutrition_2022,"treatment1")
duncan_result_soil_Ca
#Non-parameter test
kruskal.test(Ca~treatment1, data = leaf_nutrition_2022)
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,boxcox_Ca~treatment1)
dunnettT3Test(aov_model_leaf_nutrition_2022,p.adjust.method = "BH")
###t.test 
var.test(leaf_nutrition_2021$K,leaf_nutrition_2022$K)#>0.05表示方差齐性
t.test(leaf_nutrition_2021$K,leaf_nutrition_2022$K,var.equal = T)#参数检验
compare_means(K~time, data = leaf_nutrition, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验

# 5.1.2 Plot #
Ca_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "Ca/%", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
Ca_plot




#### 8 Mg####
### 4.1 Mg ###
## 4.1.1 statistical analysis##
Mg_mean <- aggregate(leaf_nutrition$Mg, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=mean)
Mg_sd <- aggregate(leaf_nutrition$Mg, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=sd)
Mg_len <- aggregate(leaf_nutrition$Mg, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=length)
df_res <- data.frame(Mg_mean, sd=Mg_sd$x, len=Mg_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(Mg ~ treatment1, data = leaf_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2021$Mg)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
leaf_nutrition_2021<-leaf_nutrition_2021%>%mutate(boxcox_Mg =BoxCox(leaf_nutrition_2021$Mg,lambda="auto"))


leveneTest(boxcox_Mg ~ treatment1, data = leaf_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2021$boxcox_Mg)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,Mg~treatment1)
summary(aov_model_leaf_nutrition_2021)
LSD_model_leaf_nutrition_2021<-LSD.test(aov_model_leaf_nutrition_2021,"treatment1",p.adj = "BH")
LSD_model_leaf_nutrition_2021
#DUNCAN
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,Mg~treatment1)
summary(aov_model_leaf_nutrition_2021)
compare_means(data=leaf_nutrition_2021,Mg~treatment1,method = "anova")
duncan_result_soil_Mg<- duncan.test(aov_model_leaf_nutrition_2021,"treatment1")
duncan_result_soil_Mg
#Non-parameter test
kruskal.test(Mg~treatment1, data = leaf_nutrition_2021)
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,boxcox_Mg~treatment1)
dunnettT3Test(aov_model_leaf_nutrition_2021,p.adjust.method = "BH")

#2022
leveneTest(Mg ~ treatment1, data = leaf_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2022$Mg)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
leaf_nutrition_2022<-leaf_nutrition_2022%>%mutate(boxcox_Mg =BoxCox(leaf_nutrition_2022$Mg,lambda="auto"))


leveneTest(boxcox_Mg ~ treatment1, data = leaf_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2022$boxcox_Mg)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,Mg~treatment1)
summary(aov_model_leaf_nutrition_2022)
LSD_model_leaf_nutrition_2022<-LSD.test(aov_model_leaf_nutrition_2022,"treatment1",p.adj = "BH")
LSD_model_leaf_nutrition_2022
#DUNCAN
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,Mg~treatment1)
summary(aov_model_leaf_nutrition_2022)
compare_means(data=leaf_nutrition_2022,Mg~treatment1,method = "anova")
duncan_result_soil_Mg<- duncan.test(aov_model_leaf_nutrition_2022,"treatment1")
duncan_result_soil_Mg
#Non-parameter test
kruskal.test(Mg~treatment1, data = leaf_nutrition_2022)
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,boxcox_Mg~treatment1)
dunnettT3Test(aov_model_leaf_nutrition_2022,p.adjust.method = "BH")
###t.test 
var.test(leaf_nutrition_2021$Mg,leaf_nutrition_2022$Mg)#>0.05表示方差齐性
t.test(leaf_nutrition_2021$Mg,leaf_nutrition_2022$Mg,var.equal = T)#参数检验
compare_means(Mg~time, data = leaf_nutrition, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验

# 5.1.2 Plot #
Mg_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "Mg/%", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
Mg_plot












#### 9 Fe####
### 4.1 Fe ###
## 4.1.1 statistical analysis##
Fe_mean <- aggregate(leaf_nutrition$Fe, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=mean)
Fe_sd <- aggregate(leaf_nutrition$Fe, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=sd)
Fe_len <- aggregate(leaf_nutrition$Fe, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=length)
df_res <- data.frame(Fe_mean, sd=Fe_sd$x, len=Fe_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(Fe ~ treatment1, data = leaf_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2021$Fe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
leaf_nutrition_2021<-leaf_nutrition_2021%>%mutate(boxcox_Fe =BoxCox(leaf_nutrition_2021$Fe,lambda="auto"))


leveneTest(boxcox_Fe ~ treatment1, data = leaf_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2021$boxcox_Fe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,Fe~treatment1)
summary(aov_model_leaf_nutrition_2021)
LSD_model_leaf_nutrition_2021<-LSD.test(aov_model_leaf_nutrition_2021,"treatment1",p.adj = "BH")
LSD_model_leaf_nutrition_2021
#DUNCAN
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,Fe~treatment1)
summary(aov_model_leaf_nutrition_2021)
compare_means(data=leaf_nutrition_2021,Fe~treatment1,method = "anova")
duncan_result_soil_Fe<- duncan.test(aov_model_leaf_nutrition_2021,"treatment1")
duncan_result_soil_Fe
#Non-parameter test
kruskal.test(Fe~treatment1, data = leaf_nutrition_2021)
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,boxcox_Fe~treatment1)
dunnettT3Test(aov_model_leaf_nutrition_2021,p.adjust.method = "BH")

#2022
leveneTest(Fe ~ treatment1, data = leaf_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2022$Fe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
leaf_nutrition_2022<-leaf_nutrition_2022%>%mutate(boxcox_Fe =BoxCox(leaf_nutrition_2022$Fe,lambda="auto"))


leveneTest(boxcox_Fe ~ treatment1, data = leaf_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2022$boxcox_Fe)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,Fe~treatment1)
summary(aov_model_leaf_nutrition_2022)
LSD_model_leaf_nutrition_2022<-LSD.test(aov_model_leaf_nutrition_2022,"treatment1",p.adj = "BH")
LSD_model_leaf_nutrition_2022
#DUNCAN
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,Fe~treatment1)
summary(aov_model_leaf_nutrition_2022)
compare_means(data=leaf_nutrition_2022,Fe~treatment1,method = "anova")
duncan_result_soil_Fe<- duncan.test(aov_model_leaf_nutrition_2022,"treatment1")
duncan_result_soil_Fe
#Non-parameter test
kruskal.test(Fe~treatment1, data = leaf_nutrition_2022)
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,boxcox_Fe~treatment1)
dunnettT3Test(aov_model_leaf_nutrition_2022,p.adjust.method = "BH")
###t.test 
var.test(leaf_nutrition_2021$Fe,leaf_nutrition_2022$Fe)#>0.05表示方差齐性
t.test(leaf_nutrition_2021$Fe,leaf_nutrition_2022$Fe,var.equal = T)#参数检验
compare_means(Fe~time, data = leaf_nutrition, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验

# 5.1.2 Plot #
Fe_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "Fe/mg·kg-1", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
Fe_plot




#### 10 Mn####
### 4.1 Mn ###
## 4.1.1 statistical analysis##
Mn_mean <- aggregate(leaf_nutrition$Mn, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=mean)
Mn_sd <- aggregate(leaf_nutrition$Mn, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=sd)
Mn_len <- aggregate(leaf_nutrition$Mn, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=length)
df_res <- data.frame(Mn_mean, sd=Mn_sd$x, len=Mn_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(Mn ~ treatment1, data = leaf_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2021$Mn)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
leaf_nutrition_2021<-leaf_nutrition_2021%>%mutate(boxcox_Mn =BoxCox(leaf_nutrition_2021$Mn,lambda="auto"))


leveneTest(boxcox_Mn ~ treatment1, data = leaf_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2021$boxcox_Mn)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,Mn~treatment1)
summary(aov_model_leaf_nutrition_2021)
LSD_model_leaf_nutrition_2021<-LSD.test(aov_model_leaf_nutrition_2021,"treatment1",p.adj = "BH")
LSD_model_leaf_nutrition_2021
#DUNCAN
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,Mn~treatment1)
summary(aov_model_leaf_nutrition_2021)
compare_means(data=leaf_nutrition_2021,Mn~treatment1,method = "anova")
duncan_result_soil_Mn<- duncan.test(aov_model_leaf_nutrition_2021,"treatment1")
duncan_result_soil_Mn
#Non-parameter test
kruskal.test(Mn~treatment1, data = leaf_nutrition_2021)
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,boxcox_Mn~treatment1)
dunnettT3Test(aov_model_leaf_nutrition_2021,p.adjust.method = "BH")

#2022
leveneTest(Mn ~ treatment1, data = leaf_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2022$Mn)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
leaf_nutrition_2022<-leaf_nutrition_2022%>%mutate(boxcox_Mn =BoxCox(leaf_nutrition_2022$Mn,lambda="auto"))


leveneTest(boxcox_Mn ~ treatment1, data = leaf_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2022$boxcox_Mn)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,Mn~treatment1)
summary(aov_model_leaf_nutrition_2022)
LSD_model_leaf_nutrition_2022<-LSD.test(aov_model_leaf_nutrition_2022,"treatment1",p.adj = "BH")
LSD_model_leaf_nutrition_2022
#DUNCAN
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,Mn~treatment1)
summary(aov_model_leaf_nutrition_2022)
compare_means(data=leaf_nutrition_2022,Mn~treatment1,method = "anova")
duncan_result_soil_Mn<- duncan.test(aov_model_leaf_nutrition_2022,"treatment1")
duncan_result_soil_Mn
#Non-parameter test
kruskal.test(Mn~treatment1, data = leaf_nutrition_2022)
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,boxcox_Mn~treatment1)
dunnettT3Test(aov_model_leaf_nutrition_2022,p.adjust.method = "BH")
###t.test 
var.test(leaf_nutrition_2021$Mn,leaf_nutrition_2022$Mn)#>0.05表示方差齐性
t.test(leaf_nutrition_2021$Mn,leaf_nutrition_2022$Mn,var.equal = T)#参数检验
compare_means(Mn~time, data = leaf_nutrition, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验

# 5.1.2 Plot #
Mn_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "Mn/mg·kg-1", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
Mn_plot




#### 11 Zn####
### 4.1 Zn ###
## 4.1.1 statistical analysis##
Zn_mean <- aggregate(leaf_nutrition$Zn, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=mean)
Zn_sd <- aggregate(leaf_nutrition$Zn, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=sd)
Zn_len <- aggregate(leaf_nutrition$Zn, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=length)
df_res <- data.frame(Zn_mean, sd=Zn_sd$x, len=Zn_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(Zn ~ treatment1, data = leaf_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2021$Zn)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
leaf_nutrition_2021<-leaf_nutrition_2021%>%mutate(boxcox_Zn =BoxCox(leaf_nutrition_2021$Zn,lambda="auto"))


leveneTest(boxcox_Zn ~ treatment1, data = leaf_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2021$boxcox_Zn)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,Zn~treatment1)
summary(aov_model_leaf_nutrition_2021)
LSD_model_leaf_nutrition_2021<-LSD.test(aov_model_leaf_nutrition_2021,"treatment1",p.adj = "BH")
LSD_model_leaf_nutrition_2021
#DUNCAN
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,Zn~treatment1)
summary(aov_model_leaf_nutrition_2021)
compare_means(data=leaf_nutrition_2021,Zn~treatment1,method = "anova")
duncan_result_soil_Zn<- duncan.test(aov_model_leaf_nutrition_2021,"treatment1")
duncan_result_soil_Zn
#Non-parameter test
kruskal.test(Zn~treatment1, data = leaf_nutrition_2021)
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,boxcox_Zn~treatment1)
dunnettT3Test(aov_model_leaf_nutrition_2021,p.adjust.method = "BH")

#2022
leveneTest(Zn ~ treatment1, data = leaf_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2022$Zn)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
leaf_nutrition_2022<-leaf_nutrition_2022%>%mutate(boxcox_Zn =BoxCox(leaf_nutrition_2022$Zn,lambda="auto"))


leveneTest(boxcox_Zn ~ treatment1, data = leaf_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2022$boxcox_Zn)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,Zn~treatment1)
summary(aov_model_leaf_nutrition_2022)
LSD_model_leaf_nutrition_2022<-LSD.test(aov_model_leaf_nutrition_2022,"treatment1",p.adj = "BH")
LSD_model_leaf_nutrition_2022
#DUNCAN
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,Zn~treatment1)
summary(aov_model_leaf_nutrition_2022)
compare_means(data=leaf_nutrition_2022,Zn~treatment1,method = "anova")
duncan_result_soil_Zn<- duncan.test(aov_model_leaf_nutrition_2022,"treatment1")
duncan_result_soil_Zn
#Non-parameter test
kruskal.test(Zn~treatment1, data = leaf_nutrition_2022)
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,boxcox_Zn~treatment1)
dunnettT3Test(aov_model_leaf_nutrition_2022,p.adjust.method = "BH")
###t.test 
var.test(leaf_nutrition_2021$Zn,leaf_nutrition_2022$Zn)#>0.05表示方差齐性
t.test(leaf_nutrition_2021$Zn,leaf_nutrition_2022$Zn,var.equal = T)#参数检验
compare_means(Zn~time, data = leaf_nutrition, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验
# 5.1.2 Plot #
Zn_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "Zn/mg·kg-1", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
Zn_plot








####12 all####
leaf_nutrition_plot <- ggarrange(N_plot,P_plot,K_plot,Ca_plot,Mg_plot,
                                 Fe_plot,Mn_plot,Zn_plot,
                                 ncol = 4, nrow = 2,common.legend = T,legend = "right")
leaf_nutrition_plot
setwd(wdOutput)
getwd()
ggsave(paste("leaf_nutrition_plot",".pdf",sep=""),leaf_nutrition_plot,
       device=cairo_pdf,width=240,height=120,dpi = 600,units = "mm")

####SPAD####
### 4.1 SPAD ###
## 4.1.1 statistical analysis##
leaf_nutrition$treatment<-factor(leaf_nutrition$treatment,levels=c("CK","CZB","NF"))
leaf_nutrition$time<-factor(leaf_nutrition$time,levels=c("2021","2022"))
SPAD_mean <- aggregate(leaf_nutrition$SPAD, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=mean)
SPAD_sd <- aggregate(leaf_nutrition$SPAD, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=sd)
SPAD_len <- aggregate(leaf_nutrition$SPAD, by=list(leaf_nutrition$treatment, leaf_nutrition$time), FUN=length)
df_res <- data.frame(SPAD_mean, sd=SPAD_sd$x, len=SPAD_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(SPAD ~ treatment1, data = leaf_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2021$SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
leaf_nutrition_2021<-leaf_nutrition_2021%>%mutate(boxcox_SPAD =BoxCox(leaf_nutrition_2021$SPAD,lambda="auto"))


leveneTest(boxcox_SPAD ~ treatment1, data = leaf_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2021$boxcox_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,SPAD~treatment1)
summary(aov_model_leaf_nutrition_2021)
LSD_model_leaf_nutrition_2021<-LSD.test(aov_model_leaf_nutrition_2021,"treatment1",p.adj = "BH")
LSD_model_leaf_nutrition_2021

#Non-parameter test
kruskal.test(SPAD~treatment1, data = leaf_nutrition_2021)
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,boxcox_SPAD~treatment1)
dunnettT3Test(aov_model_leaf_nutrition_2021,p.adjust.method = "BH")
#DUNCAN
aov_model_leaf_nutrition_2021<-aov(data=leaf_nutrition_2021,SPAD~treatment1)
summary(aov_model_leaf_nutrition_2021)
compare_means(data=leaf_nutrition_2021,SPAD~treatment1,method = "anova")
duncan_result_soil_SPAD<- duncan.test(aov_model_leaf_nutrition_2021,"treatment1")
duncan_result_soil_SPAD
#2022
leveneTest(SPAD ~ treatment1, data = leaf_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2022$SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
leaf_nutrition_2022<-leaf_nutrition_2022%>%mutate(boxcox_SPAD =BoxCox(leaf_nutrition_2022$SPAD,lambda="auto"))


leveneTest(boxcox_SPAD ~ treatment1, data = leaf_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(leaf_nutrition_2022$boxcox_SPAD)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,SPAD~treatment1)
summary(aov_model_leaf_nutrition_2022)
LSD_model_leaf_nutrition_2022<-LSD.test(aov_model_leaf_nutrition_2022,"treatment1",p.adj = "BH")
LSD_model_leaf_nutrition_2022
#DUNCAN
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,SPAD~treatment1)
summary(aov_model_leaf_nutrition_2022)
compare_means(data=leaf_nutrition_2022,SPAD~treatment1,method = "anova")
duncan_result_soil_SPAD<- duncan.test(aov_model_leaf_nutrition_2022,"treatment1")
duncan_result_soil_SPAD
#Non-parameter test
kruskal.test(SPAD~treatment1, data = leaf_nutrition_2022)
aov_model_leaf_nutrition_2022<-aov(data=leaf_nutrition_2022,boxcox_SPAD~treatment1)
dunnettT3Test(aov_model_leaf_nutrition_2022,p.adjust.method = "BH")

###t.test 
var.test(leaf_nutrition_2021$SPAD,leaf_nutrition_2022$SPAD)#>0.05表示方差齐性
t.test(leaf_nutrition_2021$SPAD,leaf_nutrition_2022$SPAD,var.equal = T)#参数检验
compare_means(SPAD~time, data = leaf_nutrition, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验
# 5.1.2 Plot #
SPAD_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "SPAD", fill = "treatment")+
  scale_fill_manual(values=c("#5C4A46","#DA9464","#888e4a"))+
  FacetTheme
SPAD_plot

setwd(wdOutput)
getwd()
ggsave(paste("SPAD_plot",".pdf",sep=""),
       SPAD_plot,device=cairo_pdf,width=90,height=70,dpi = 300,units = "mm")