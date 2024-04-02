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
wdOutput <- c("E:/whx_24samples_filtered1/06.env/soil_nutrition")

#### 3. soil ####
### 3.1 Import and process data ###
setwd(wdImport)
SampleData <- read.table("24SampleData.txt",header=T,na.strings = c("NA"))
soil_nutrition<- read.table("env_soil_nutrition.txt",header=T,fill=TRUE,na.strings = c("NA"))
soil_nutrition <- merge(SampleData, soil_nutrition, by = 'sampleid')
soil_nutrition$treatment<-factor(soil_nutrition$treatment,levels=c("NF","CZB","CK"))
soil_nutrition$time<-factor(soil_nutrition$time,levels=c("2022","2021"))
soil_nutrition$treatment1<-factor(soil_nutrition$treatment1,levels=c("CK_21","CZB_21","NF_21","CK_22","CZB_22","NF_22"))
soil_nutrition_2021<-filter(soil_nutrition,time=="2021")
soil_nutrition_2022<-filter(soil_nutrition,time=="2022")
#### 4 K####
### 4.1 K ###
## 4.1.1 statistical analysis##
K_mean <- aggregate(soil_nutrition$K, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=mean)
K_sd <- aggregate(soil_nutrition$K, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=sd)
K_len <- aggregate(soil_nutrition$K, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=length)
df_res <- data.frame(K_mean, sd=K_sd$x, len=K_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(K ~ treatment1, data = soil_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2021$K)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
soil_nutrition_2021<-soil_nutrition_2021%>%mutate(boxcox_K =BoxCox(soil_nutrition_2021$K,lambda="auto"))


leveneTest(boxcox_K ~ treatment1, data = soil_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2021$boxcox_K)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,K~treatment1)
summary(aov_model_soil_nutrition_2021)
LSD_model_soil_nutrition_2021<-LSD.test(aov_model_soil_nutrition_2021,"treatment1",p.adj = "BH")
LSD_model_soil_nutrition_2021
#DUNCAN
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,K~treatment1)
summary(aov_model_soil_nutrition_2021)
compare_means(data=soil_nutrition_2021,K~treatment1,method = "anova")
duncan_result_soil_K<- duncan.test(aov_model_soil_nutrition_2021,"treatment1")
duncan_result_soil_K
#Non-parameter test
kruskal.test(K~treatment1, data = soil_nutrition_2021)
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,boxcox_K~treatment1)
dunnettT3Test(aov_model_soil_nutrition_2021,p.adjust.method = "BH")

#2022
leveneTest(K ~ treatment1, data = soil_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2022$K)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
soil_nutrition_2022<-soil_nutrition_2022%>%mutate(boxcox_K =BoxCox(soil_nutrition_2022$K,lambda="auto"))


leveneTest(boxcox_K ~ treatment1, data = soil_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2022$boxcox_K)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,K~treatment1)
summary(aov_model_soil_nutrition_2022)
LSD_model_soil_nutrition_2022<-LSD.test(aov_model_soil_nutrition_2022,"treatment1",p.adj = "BH")
LSD_model_soil_nutrition_2022
#DUNCAN
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,K~treatment1)
summary(aov_model_soil_nutrition_2022)
compare_means(data=soil_nutrition_2022,K~treatment1,method = "anova")
duncan_result_soil_K<- duncan.test(aov_model_soil_nutrition_2022,"treatment1")
duncan_result_soil_K
#Non-parameter test
kruskal.test(K~treatment1, data = soil_nutrition_2022)
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,boxcox_K~treatment1)
dunnettT3Test(aov_model_soil_nutrition_2022,p.adjust.method = "BH")
###t.test 
var.test(soil_nutrition_2021$K,soil_nutrition_2022$K)#>0.05表示方差齐性
t.test(soil_nutrition_2021$K,soil_nutrition_2022$K,var.equal = T)#参数检验
compare_means(K~time, data = soil_nutrition, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验

# 5.1.2 Plot #
K_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "k/%", fill = "treatment")+
  scale_fill_manual(values=c("#809A54","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
K_plot




#### 5 N####
### 4.1 N ###
## 4.1.1 statistical analysis##
N_mean <- aggregate(soil_nutrition$N, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=mean)
N_sd <- aggregate(soil_nutrition$N, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=sd)
N_len <- aggregate(soil_nutrition$N, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=length)
df_res <- data.frame(N_mean, sd=N_sd$x, len=N_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(N ~ treatment1, data = soil_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2021$N)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
soil_nutrition_2021<-soil_nutrition_2021%>%mutate(boxcox_N =BoxCox(soil_nutrition_2021$N,lambda="auto"))


leveneTest(boxcox_N ~ treatment1, data = soil_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2021$boxcox_N)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,N~treatment1)
summary(aov_model_soil_nutrition_2021)
LSD_model_soil_nutrition_2021<-LSD.test(aov_model_soil_nutrition_2021,"treatment1",p.adj = "BH")
LSD_model_soil_nutrition_2021

#Non-parameter test
kruskal.test(N~treatment1, data = soil_nutrition_2021)
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,boxcox_N~treatment1)
dunnettT3Test(aov_model_soil_nutrition_2021,p.adjust.method = "BH")
#DUNCAN
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,N~treatment1)
summary(aov_model_soil_nutrition_2021)
compare_means(data=soil_nutrition_2021,N~treatment1,method = "anova")
duncan_result_soil_N<- duncan.test(aov_model_soil_nutrition_2021,"treatment1")
duncan_result_soil_N
#2022
leveneTest(N ~ treatment1, data = soil_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2022$N)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
soil_nutrition_2022<-soil_nutrition_2022%>%mutate(boxcox_N =BoxCox(soil_nutrition_2022$N,lambda="auto"))


leveneTest(boxcox_N ~ treatment1, data = soil_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2022$boxcox_N)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,N~treatment1)
summary(aov_model_soil_nutrition_2022)
LSD_model_soil_nutrition_2022<-LSD.test(aov_model_soil_nutrition_2022,"treatment1",p.adj = "BH")
LSD_model_soil_nutrition_2022
#DUNCAN
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,N~treatment1)
summary(aov_model_soil_nutrition_2022)
compare_means(data=soil_nutrition_2022,N~treatment1,method = "anova")
duncan_result_soil_N<- duncan.test(aov_model_soil_nutrition_2022,"treatment1")
duncan_result_soil_N
#Non-parameter test
kruskal.test(N~treatment1, data = soil_nutrition_2022)
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,boxcox_N~treatment1)
dunnettT3Test(aov_model_soil_nutrition_2022,p.adjust.method = "BH")
###t.test 
var.test(soil_nutrition_2021$N,soil_nutrition_2022$N)#>0.05表示方差齐性
t.test(soil_nutrition_2021$N,soil_nutrition_2022$N,var.equal = T)#参数检验
compare_means(N~time, data = soil_nutrition, 
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
### 6.1 P2021 ###
## 6.1.1 statistical analysis##
soil_nutrition_2021<- filter(soil_nutrition,time=="2021")
P2021_mean <- aggregate(soil_nutrition_2021$P_bulk, by=list(soil_nutrition_2021$treatment, soil_nutrition_2021$time), FUN=mean)
P2021_sd <- aggregate(soil_nutrition_2021$P_bulk, by=list(soil_nutrition_2021$treatment, soil_nutrition_2021$time), FUN=sd)
P2021_len <- aggregate(soil_nutrition_2021$P_bulk, by=list(soil_nutrition_2021$treatment, soil_nutrition_2021$time), FUN=length)
df_res <- data.frame(P2021_mean, sd=P2021_sd$x, len=P2021_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(P_bulk ~ treatment1, data = soil_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2021$P_bulk)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
soil_nutrition_2021<-soil_nutrition_2021%>%mutate(boxcox_P_bulk =BoxCox(soil_nutrition_2021$P_bulk,lambda="auto"))


leveneTest(boxcox_P_bulk ~ treatment1, data = soil_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2021$boxcox_P_bulk)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,P_bulk~treatment1)
summary(aov_model_soil_nutrition_2021)
LSD_model_soil_nutrition_2021<-LSD.test(aov_model_soil_nutrition_2021,"treatment1",p.adj = "BH")
LSD_model_soil_nutrition_2021
#DUNCAN
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,P_bulk~treatment1)
summary(aov_model_soil_nutrition_2021)
compare_means(data=soil_nutrition_2021,P_bulk~treatment1,method = "anova")
duncan_result_soil_P_bulk<- duncan.test(aov_model_soil_nutrition_2021,"treatment1")
duncan_result_soil_P_bulk
#Non-parameter test
kruskal.test(P_bulk~treatment1, data = soil_nutrition_2021)
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,boxcox_P_bulk~treatment1)
dunnettT3Test(aov_model_soil_nutrition_2021,p.adjust.method = "BH")


# 6.1.2 Plot #
P2021_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "P/%", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
P2021_plot
### 6.2 P2022_bulk ###
## 6.2.1 statistical analysis##
P2022_bulk_mean <- aggregate(soil_nutrition_2022$P_bulk, by=list(soil_nutrition_2022$treatment, soil_nutrition_2022$time), FUN=mean)
P2022_bulk_sd <- aggregate(soil_nutrition_2022$P_bulk, by=list(soil_nutrition_2022$treatment, soil_nutrition_2022$time), FUN=sd)
P2022_bulk_len <- aggregate(soil_nutrition_2022$P_bulk, by=list(soil_nutrition_2022$treatment, soil_nutrition_2022$time), FUN=length)
df_res <- data.frame(P2022_bulk_mean, sd=P2022_bulk_sd$x, len=P2022_bulk_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2022
leveneTest(P_bulk ~ treatment1, data = soil_nutrition_2022)#P_bulk>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2022$P_bulk)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
soil_nutrition_2022<-soil_nutrition_2022%>%mutate(boxcox_P_bulk =BoxCox(soil_nutrition_2022$P_bulk,lambda="auto"))


leveneTest(boxcox_P_bulk ~ treatment1, data = soil_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2022$boxcox_P_bulk)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,P_bulk~treatment1)
summary(aov_model_soil_nutrition_2022)
LSD_model_soil_nutrition_2022<-LSD.test(aov_model_soil_nutrition_2022,"treatment1",p.adj = "BH")
LSD_model_soil_nutrition_2022
#DUNCAN
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,P_bulk~treatment1)
summary(aov_model_soil_nutrition_2022)
compare_means(data=soil_nutrition_2022,P_bulk~treatment1,method = "anova")
duncan_result_soil_P_bulk<- duncan.test(aov_model_soil_nutrition_2022,"treatment1")
duncan_result_soil_P_bulk
#Non-parameter test
kruskal.test(P_bulk~treatment1, data = soil_nutrition_2022)
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,boxcox_P_bulk~treatment1)
dunnettT3Test(aov_model_soil_nutrition_2022,p.adjust.method = "BH")
# 6.1.2 Plot #
P2022_bulk_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "P_bulk/%", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
P2022_bulk_plot
### 6.2 P2022_rhizo ###
## 6.2.1 statistical analysis##
soil_nutrition_2022<- filter(soil_nutrition,time=="2022")
P2022_rhizo_mean <- aggregate(soil_nutrition_2022$P_rhizo, by=list(soil_nutrition_2022$treatment, soil_nutrition_2022$time), FUN=mean)
P2022_rhizo_sd <- aggregate(soil_nutrition_2022$P_rhizo, by=list(soil_nutrition_2022$treatment, soil_nutrition_2022$time), FUN=sd)
P2022_rhizo_len <- aggregate(soil_nutrition_2022$P_rhizo, by=list(soil_nutrition_2022$treatment, soil_nutrition_2022$time), FUN=length)
df_res <- data.frame(P2022_rhizo_mean, sd=P2022_rhizo_sd$x, len=P2022_rhizo_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2022
leveneTest(P_rhizo ~ treatment1, data = soil_nutrition_2022)#P>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2022$P_rhizo)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
soil_nutrition_2022<-soil_nutrition_2022%>%mutate(boxcox_P_rhizo =BoxCox(soil_nutrition_2022$P_rhizo,lambda="auto"))


leveneTest(boxcox_P_rhizo ~ treatment1, data = soil_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2022$boxcox_P_rhizo)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,P_rhizo~treatment1)
summary(aov_model_soil_nutrition_2022)
LSD_model_soil_nutrition_2022<-LSD.test(aov_model_soil_nutrition_2022,"treatment1",p.adj = "BH")
LSD_model_soil_nutrition_2022
#DUNCAN
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,P_rhizo~treatment1)
summary(aov_model_soil_nutrition_2022)
compare_means(data=soil_nutrition_2022,P_rhizo~treatment1,method = "anova")
duncan_result_soil_P_rhizo<- duncan.test(aov_model_soil_nutrition_2022,"treatment1")
duncan_result_soil_P_rhizo
#Non-parameter test
kruskal.test(P_rhizo~treatment1, data = soil_nutrition_2022)
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,boxcox_P_rhizo~treatment1)
dunnettT3Test(aov_model_soil_nutrition_2022,p.adjust.method = "BH")
# 6.1.2 Plot #
P2022_rhizo_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "P_rhizo/%", fill = "treatment")+
  scale_fill_manual(values=c("#888e4a","#DA9464","#5C4A46"))+
  coord_flip()+
  FacetTheme
P2022_rhizo_plot


###t.test 
var.test(soil_nutrition_2022$P_bulk,soil_nutrition_2022$P_rhizo)#>0.05表示方差齐性
t.test(soil_nutrition_2022$P_bulk,soil_nutrition_2022$P_rhizo,var.equal = T)#参数检验
compare_means(soil_nutrition_2022$P_bulk,soil_nutrition_2022$P_rhizo, data = soil_nutrition_2022, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验
#
P_plot <- ggarrange(P2021_plot,ggarrange(P2022_bulk_plot,P2022_rhizo_plot, ncol = 2,common.legend = T,legend = "none"),
                                nrow = 2,common.legend = T,legend = "none")
P_plot

####NPK####
soil_nutrition_NPK_plot <- ggarrange(N_plot,P_plot,K_plot,
                                  nrow = 1,common.legend = T,legend = "right")
soil_nutrition_NPK_plot
setwd(wdOutput)
getwd()
ggsave(paste("soil_nutrition_NPK_plot",".pdf",sep=""),soil_nutrition_NPK_plot,
       device=cairo_pdf,width=240,height=120,dpi = 600,units = "mm")


#### 7 Ca####
### 4.1 Ca ###
## 4.1.1 statistical analysis##
Ca_mean <- aggregate(soil_nutrition$Ca, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=mean)
Ca_sd <- aggregate(soil_nutrition$Ca, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=sd)
Ca_len <- aggregate(soil_nutrition$Ca, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=length)
df_res <- data.frame(Ca_mean, sd=Ca_sd$x, len=Ca_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(Ca ~ treatment1, data = soil_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2021$Ca)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
soil_nutrition_2021<-soil_nutrition_2021%>%mutate(boxcox_Ca =BoxCox(soil_nutrition_2021$Ca,lambda="auto"))


leveneTest(boxcox_Ca ~ treatment1, data = soil_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2021$boxcox_Ca)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,Ca~treatment1)
summary(aov_model_soil_nutrition_2021)
LSD_model_soil_nutrition_2021<-LSD.test(aov_model_soil_nutrition_2021,"treatment1",p.adj = "BH")
LSD_model_soil_nutrition_2021
#DUNCAN
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,Ca~treatment1)
summary(aov_model_soil_nutrition_2021)
compare_means(data=soil_nutrition_2021,Ca~treatment1,method = "anova")
duncan_result_soil_Ca<- duncan.test(aov_model_soil_nutrition_2021,"treatment1")
duncan_result_soil_Ca
#Non-parameter test
kruskal.test(Ca~treatment1, data = soil_nutrition_2021)
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,boxcox_Ca~treatment1)
dunnettT3Test(aov_model_soil_nutrition_2021,p.adjust.method = "BH")

#2022
leveneTest(Ca ~ treatment1, data = soil_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2022$Ca)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
soil_nutrition_2022<-soil_nutrition_2022%>%mutate(boxcox_Ca =BoxCox(soil_nutrition_2022$Ca,lambda="auto"))


leveneTest(boxcox_Ca ~ treatment1, data = soil_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2022$boxcox_Ca)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,Ca~treatment1)
summary(aov_model_soil_nutrition_2022)
LSD_model_soil_nutrition_2022<-LSD.test(aov_model_soil_nutrition_2022,"treatment1",p.adj = "BH")
LSD_model_soil_nutrition_2022
#DUNCAN
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,Ca~treatment1)
summary(aov_model_soil_nutrition_2022)
compare_means(data=soil_nutrition_2022,Ca~treatment1,method = "anova")
duncan_result_soil_Ca<- duncan.test(aov_model_soil_nutrition_2022,"treatment1")
duncan_result_soil_Ca
#Non-parameter test
kruskal.test(Ca~treatment1, data = soil_nutrition_2022)
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,boxcox_Ca~treatment1)
dunnettT3Test(aov_model_soil_nutrition_2022,p.adjust.method = "BH")
###t.test 
var.test(soil_nutrition_2021$Ca,soil_nutrition_2022$Ca)#>0.05表示方差齐性
t.test(soil_nutrition_2021$Ca,soil_nutrition_2022$Ca,var.equal = T)#参数检验
compare_means(Ca~time, data = soil_nutrition, 
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
Mg_mean <- aggregate(soil_nutrition$Mg, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=mean)
Mg_sd <- aggregate(soil_nutrition$Mg, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=sd)
Mg_len <- aggregate(soil_nutrition$Mg, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=length)
df_res <- data.frame(Mg_mean, sd=Mg_sd$x, len=Mg_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(Mg ~ treatment1, data = soil_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2021$Mg)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
soil_nutrition_2021<-soil_nutrition_2021%>%mutate(boxcox_Mg =BoxCox(soil_nutrition_2021$Mg,lambda="auto"))


leveneTest(boxcox_Mg ~ treatment1, data = soil_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2021$boxcox_Mg)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,Mg~treatment1)
summary(aov_model_soil_nutrition_2021)
LSD_model_soil_nutrition_2021<-LSD.test(aov_model_soil_nutrition_2021,"treatment1",p.adj = "BH")
LSD_model_soil_nutrition_2021
#DUNCAN
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,Mg~treatment1)
summary(aov_model_soil_nutrition_2021)
compare_means(data=soil_nutrition_2021,Mg~treatment1,method = "anova")
duncan_result_soil_Mg<- duncan.test(aov_model_soil_nutrition_2021,"treatment1")
duncan_result_soil_Mg
#Non-parameter test
kruskal.test(Mg~treatment1, data = soil_nutrition_2021)
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,boxcox_Mg~treatment1)
dunnettT3Test(aov_model_soil_nutrition_2021,p.adjust.method = "BH")

#2022
leveneTest(Mg ~ treatment1, data = soil_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2022$Mg)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
soil_nutrition_2022<-soil_nutrition_2022%>%mutate(boxcox_Mg =BoxCox(soil_nutrition_2022$Mg,lambda="auto"))


leveneTest(boxcox_Mg ~ treatment1, data = soil_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2022$boxcox_Mg)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,Mg~treatment1)
summary(aov_model_soil_nutrition_2022)
LSD_model_soil_nutrition_2022<-LSD.test(aov_model_soil_nutrition_2022,"treatment1",p.adj = "BH")
LSD_model_soil_nutrition_2022
#DUNCAN
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,Mg~treatment1)
summary(aov_model_soil_nutrition_2022)
compare_means(data=soil_nutrition_2022,Mg~treatment1,method = "anova")
duncan_result_soil_Mg<- duncan.test(aov_model_soil_nutrition_2022,"treatment1")
duncan_result_soil_Mg
#Non-parameter test
kruskal.test(Mg~treatment1, data = soil_nutrition_2022)
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,boxcox_Mg~treatment1)
dunnettT3Test(aov_model_soil_nutrition_2022,p.adjust.method = "BH")
###t.test 
var.test(soil_nutrition_2021$Mg,soil_nutrition_2022$Mg)#>0.05表示方差齐性
t.test(soil_nutrition_2021$Mg,soil_nutrition_2022$Mg,var.equal = T)#参数检验
compare_means(Mg~time, data = soil_nutrition, 
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
### 6.1 Fe2021 ###
## 6.1.1 statistical analysis##
soil_nutrition_2021<- filter(soil_nutrition,time=="2021")
Fe2021_mean <- aggregate(soil_nutrition_2021$Fe_bulk, by=list(soil_nutrition_2021$treatment, soil_nutrition_2021$time), FUN=mean)
Fe2021_sd <- aggregate(soil_nutrition_2021$Fe_bulk, by=list(soil_nutrition_2021$treatment, soil_nutrition_2021$time), FUN=sd)
Fe2021_len <- aggregate(soil_nutrition_2021$Fe_bulk, by=list(soil_nutrition_2021$treatment, soil_nutrition_2021$time), FUN=length)
df_res <- data.frame(Fe2021_mean, sd=Fe2021_sd$x, len=Fe2021_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
#leveneTest(Fe_bulk ~ treatment1, data = soil_nutrition_2021)#p>0.05，则满足方差齐性
#shapiro.test(soil_nutrition_2021$Fe_bulk)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
#soil_nutrition_2021<-soil_nutrition_2021%>%mutate(boxcox_Fe_bulk =BoxCox(soil_nutrition_2021$Fe_bulk,lambda="auto"))


#leveneTest(boxcox_Fe_bulk ~ treatment1, data = soil_nutrition_2021)#p>0.05，则满足方差齐性
#shapiro.test(soil_nutrition_2021$boxcox_Fe_bulk)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
#aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,Fe_bulk~treatment1)
#summary(aov_model_soil_nutrition_2021)
#LSD_model_soil_nutrition_2021<-LSD.test(aov_model_soil_nutrition_2021,"treatment1",p.adj = "BH")
#LSD_model_soil_nutrition_2021
#DUNCAN
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,Fe_bulk~treatment1)
summary(aov_model_soil_nutrition_2021)
compare_means(data=soil_nutrition_2021,Fe_bulk~treatment1,method = "anova")
duncan_result_soil_Fe_bulk<- duncan.test(aov_model_soil_nutrition_2021,"treatment1")
duncan_result_soil_Fe_bulk
#Non-parameter test
#kruskal.test(Fe_bulk~treatment1, data = soil_nutrition_2021)
#aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,boxcox_Fe_bulk~treatment1)
#dunnettT3Test(aov_model_soil_nutrition_2021,p.adjust.method = "BH")


# 6.1.2 Plot #
Fe2021_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "Fe/%", fill = "treatment")+
  scale_fill_manual(values=c("#46C580","#EB9752","#75DAFF"))+
  coord_flip()+
  FacetTheme
Fe2021_plot
### 6.2 Fe2022_bulk ###
## 6.2.1 statistical analysis##
Fe2022_bulk_mean <- aggregate(soil_nutrition_2022$Fe_bulk, by=list(soil_nutrition_2022$treatment, soil_nutrition_2022$time), FUN=mean)
Fe2022_bulk_sd <- aggregate(soil_nutrition_2022$Fe_bulk, by=list(soil_nutrition_2022$treatment, soil_nutrition_2022$time), FUN=sd)
Fe2022_bulk_len <- aggregate(soil_nutrition_2022$Fe_bulk, by=list(soil_nutrition_2022$treatment, soil_nutrition_2022$time), FUN=length)
df_res <- data.frame(Fe2022_bulk_mean, sd=Fe2022_bulk_sd$x, len=Fe2022_bulk_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2022
#leveneTest(Fe_bulk ~ treatment1, data = soil_nutrition_2022)#P>0.05，则满足方差齐性
#shapiro.test(soil_nutrition_2022$Fe_bulk)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
#soil_nutrition_2022<-soil_nutrition_2022%>%mutate(boxcox_Fe_bulk =BoxCox(soil_nutrition_2022$Fe_bulk,lambda="auto"))


#leveneTest(boxcox_Fe_bulk ~ treatment1, data = soil_nutrition_2022)#p>0.05，则满足方差齐性
#shapiro.test(soil_nutrition_2022$boxcox_Fe_bulk)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
#aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,Fe_bulk~treatment1)
#summary(aov_model_soil_nutrition_2022)
#LSD_model_soil_nutrition_2022<-LSD.test(aov_model_soil_nutrition_2022,"treatment1",p.adj = "BH")
#LSD_model_soil_nutrition_2022
#DUNCAN
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,Fe_bulk~treatment1)
summary(aov_model_soil_nutrition_2022)
compare_means(data=soil_nutrition_2022,Fe_bulk~treatment1,method = "anova")
duncan_result_soil_Fe_bulk<- duncan.test(aov_model_soil_nutrition_2022,"treatment1")
duncan_result_soil_Fe_bulk
#Non-parameter test
#kruskal.test(Fe_bulk~treatment1, data = soil_nutrition_2022)
#aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,boxcox_Fe_bulk~treatment1)
#dunnettT3Test(aov_model_soil_nutrition_2022,p.adjust.method = "BH")
# 6.1.2 Plot #
Fe2022_bulk_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "Fe_bulk/%", fill = "treatment")+
  scale_fill_manual(values=c("#46C580","#EB9752","#75DAFF"))+
  coord_flip()+
  FacetTheme
Fe2022_bulk_plot
### 6.2 P2022_rhizo ###
## 6.2.1 statistical analysis##
soil_nutrition_2022<- filter(soil_nutrition,time=="2022")
Fe2022_rhizo_mean <- aggregate(soil_nutrition_2022$Fe_rhizo, by=list(soil_nutrition_2022$treatment, soil_nutrition_2022$time), FUN=mean)
Fe2022_rhizo_sd <- aggregate(soil_nutrition_2022$Fe_rhizo, by=list(soil_nutrition_2022$treatment, soil_nutrition_2022$time), FUN=sd)
Fe2022_rhizo_len <- aggregate(soil_nutrition_2022$Fe_rhizo, by=list(soil_nutrition_2022$treatment, soil_nutrition_2022$time), FUN=length)
df_res <- data.frame(Fe2022_rhizo_mean, sd=Fe2022_rhizo_sd$x, len=Fe2022_rhizo_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2022
#leveneTest(Fe_rhizo ~ treatment1, data = soil_nutrition_2022)#P>0.05，则满足方差齐性
#shapiro.test(soil_nutrition_2022$Fe_rhizo)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
#soil_nutrition_2022<-soil_nutrition_2022%>%mutate(boxcox_Fe_rhizo =BoxCox(soil_nutrition_2022$Fe_rhizo,lambda="auto"))


#leveneTest(boxcox_Fe_rhizo ~ treatment1, data = soil_nutrition_2022)#p>0.05，则满足方差齐性
#shapiro.test(soil_nutrition_2022$boxcox_Fe_rhizo)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
#aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,Fe_rhizo~treatment1)
#summary(aov_model_soil_nutrition_2022)
#LSD_model_soil_nutrition_2022<-LSD.test(aov_model_soil_nutrition_2022,"treatment1",p.adj = "BH")
#LSD_model_soil_nutrition_2022
#DUNCAN
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,Fe_rhizo~treatment1)
summary(aov_model_soil_nutrition_2022)
compare_means(data=soil_nutrition_2022,Fe_rhizo~treatment1,method = "anova")
duncan_result_soil_Fe_rhizo<- duncan.test(aov_model_soil_nutrition_2022,"treatment1")
duncan_result_soil_Fe_rhizo
#Non-parameter test
#kruskal.test(Fe_rhizo~treatment1, data = soil_nutrition_2022)
#aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,boxcox_Fe_rhizo~treatment1)
#dunnettT3Test(aov_model_soil_nutrition_2022,p.adjust.method = "BH")
# 6.1.2 Plot #
Fe2022_rhizo_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "Fe_rhizo/%", fill = "treatment")+
  scale_fill_manual(values=c("#46C580","#EB9752","#75DAFF"))+
  coord_flip()+
  FacetTheme
Fe2022_rhizo_plot


###t.test 
#var.test(soil_nutrition_2022$Fe_bulk,soil_nutrition_2022$Fe_rhizo)#>0.05表示方差齐性
#t.test(soil_nutrition_2022$Fe_bulk,soil_nutrition_2022$Fe_rhizo,var.equal = T)#参数检验
#compare_means(soil_nutrition_2022$Fe_bulk,soil_nutrition_2022$Fe_rhizo, data = soil_nutrition_2022, 
#              ref.group = ".all.", 
#              method = "wilcox.test")#非参数检验
#plot#
Fe_plot <- ggarrange(Fe2021_plot,ggarrange(Fe2022_bulk_plot,Fe2022_rhizo_plot, ncol = 2,common.legend = T,legend = "none"),
                    nrow = 2,common.legend = T,legend = "none")
Fe_plot




#### 10 Mn####
### 4.1 Mn ###
## 4.1.1 statistical analysis##
Mn_mean <- aggregate(soil_nutrition$Mn, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=mean)
Mn_sd <- aggregate(soil_nutrition$Mn, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=sd)
Mn_len <- aggregate(soil_nutrition$Mn, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=length)
df_res <- data.frame(Mn_mean, sd=Mn_sd$x, len=Mn_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(Mn ~ treatment1, data = soil_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2021$Mn)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
soil_nutrition_2021<-soil_nutrition_2021%>%mutate(boxcox_Mn =BoxCox(soil_nutrition_2021$Mn,lambda="auto"))


leveneTest(boxcox_Mn ~ treatment1, data = soil_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2021$boxcox_Mn)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,Mn~treatment1)
summary(aov_model_soil_nutrition_2021)
LSD_model_soil_nutrition_2021<-LSD.test(aov_model_soil_nutrition_2021,"treatment1",p.adj = "BH")
LSD_model_soil_nutrition_2021
#DUNCAN
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,Mn~treatment1)
summary(aov_model_soil_nutrition_2021)
compare_means(data=soil_nutrition_2021,Mn~treatment1,method = "anova")
duncan_result_soil_Mn<- duncan.test(aov_model_soil_nutrition_2021,"treatment1")
duncan_result_soil_Mn
#Non-parameter test
kruskal.test(Mn~treatment1, data = soil_nutrition_2021)
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,boxcox_Mn~treatment1)
dunnettT3Test(aov_model_soil_nutrition_2021,p.adjust.method = "BH")

#2022
leveneTest(Mn ~ treatment1, data = soil_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2022$Mn)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
soil_nutrition_2022<-soil_nutrition_2022%>%mutate(boxcox_Mn =BoxCox(soil_nutrition_2022$Mn,lambda="auto"))


leveneTest(boxcox_Mn ~ treatment1, data = soil_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2022$boxcox_Mn)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,Mn~treatment1)
summary(aov_model_soil_nutrition_2022)
LSD_model_soil_nutrition_2022<-LSD.test(aov_model_soil_nutrition_2022,"treatment1",p.adj = "BH")
LSD_model_soil_nutrition_2022
#DUNCAN
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,Mn~treatment1)
summary(aov_model_soil_nutrition_2022)
compare_means(data=soil_nutrition_2022,Mn~treatment1,method = "anova")
duncan_result_soil_Mn<- duncan.test(aov_model_soil_nutrition_2022,"treatment1")
duncan_result_soil_Mn
#Non-parameter test
kruskal.test(Mn~treatment1, data = soil_nutrition_2022)
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,boxcox_Mn~treatment1)
dunnettT3Test(aov_model_soil_nutrition_2022,p.adjust.method = "BH")
###t.test 
var.test(soil_nutrition_2021$Mn,soil_nutrition_2022$Mn)#>0.05表示方差齐性
t.test(soil_nutrition_2021$Mn,soil_nutrition_2022$Mn,var.equal = T)#参数检验
compare_means(Mn~time, data = soil_nutrition, 
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
Zn_mean <- aggregate(soil_nutrition$Zn, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=mean)
Zn_sd <- aggregate(soil_nutrition$Zn, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=sd)
Zn_len <- aggregate(soil_nutrition$Zn, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=length)
df_res <- data.frame(Zn_mean, sd=Zn_sd$x, len=Zn_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(Zn ~ treatment1, data = soil_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2021$Zn)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
soil_nutrition_2021<-soil_nutrition_2021%>%mutate(boxcox_Zn =BoxCox(soil_nutrition_2021$Zn,lambda="auto"))


leveneTest(boxcox_Zn ~ treatment1, data = soil_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2021$boxcox_Zn)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,Zn~treatment1)
summary(aov_model_soil_nutrition_2021)
LSD_model_soil_nutrition_2021<-LSD.test(aov_model_soil_nutrition_2021,"treatment1",p.adj = "BH")
LSD_model_soil_nutrition_2021
#DUNCAN
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,Zn~treatment1)
summary(aov_model_soil_nutrition_2021)
compare_means(data=soil_nutrition_2021,Zn~treatment1,method = "anova")
duncan_result_soil_Zn<- duncan.test(aov_model_soil_nutrition_2021,"treatment1")
duncan_result_soil_Zn
#Non-parameter test
kruskal.test(Zn~treatment1, data = soil_nutrition_2021)
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,boxcox_Zn~treatment1)
dunnettT3Test(aov_model_soil_nutrition_2021,p.adjust.method = "BH")

#2022
leveneTest(Zn ~ treatment1, data = soil_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2022$Zn)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
soil_nutrition_2022<-soil_nutrition_2022%>%mutate(boxcox_Zn =BoxCox(soil_nutrition_2022$Zn,lambda="auto"))


leveneTest(boxcox_Zn ~ treatment1, data = soil_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2022$boxcox_Zn)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,Zn~treatment1)
summary(aov_model_soil_nutrition_2022)
LSD_model_soil_nutrition_2022<-LSD.test(aov_model_soil_nutrition_2022,"treatment1",p.adj = "BH")
LSD_model_soil_nutrition_2022
#DUNCAN
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,Zn~treatment1)
summary(aov_model_soil_nutrition_2022)
compare_means(data=soil_nutrition_2022,Zn~treatment1,method = "anova")
duncan_result_soil_Zn<- duncan.test(aov_model_soil_nutrition_2022,"treatment1")
duncan_result_soil_Zn
#Non-parameter test
kruskal.test(Zn~treatment1, data = soil_nutrition_2022)
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,boxcox_Zn~treatment1)
dunnettT3Test(aov_model_soil_nutrition_2022,p.adjust.method = "BH")
###t.test 
var.test(soil_nutrition_2021$Zn,soil_nutrition_2022$Zn)#>0.05表示方差齐性
t.test(soil_nutrition_2021$Zn,soil_nutrition_2022$Zn,var.equal = T)#参数检验
compare_means(Zn~time, data = soil_nutrition, 
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





####13 all####
soil_nutrition_low_plot <- ggarrange(Ca_plot,Mg_plot,
                                 Fe_plot,Mn_plot,Zn_plot,
                                 ncol = 3, nrow = 2,common.legend = T,legend = "right")
soil_nutrition_low_plot
setwd(wdOutput)
getwd()
ggsave(paste("soil_nutrition_low_plot",".pdf",sep=""),soil_nutrition_low_plot,
       device=cairo_pdf,width=210,height=120,dpi = 600,units = "mm")

####14 pH####

####pH####
### 4.1 pH ###
## 4.1.1 statistical analysis##
soil_nutrition$treatment<-factor(soil_nutrition$treatment,levels=c("CK","CZB","NF"))
soil_nutrition$time<-factor(soil_nutrition$time,levels=c("2021","2022"))
pH_mean <- aggregate(soil_nutrition$pH, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=mean)
pH_sd <- aggregate(soil_nutrition$pH, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=sd)
pH_len <- aggregate(soil_nutrition$pH, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=length)
df_res <- data.frame(pH_mean, sd=pH_sd$x, len=pH_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(pH ~ treatment1, data = soil_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2021$pH)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
soil_nutrition_2021<-soil_nutrition_2021%>%mutate(boxcox_pH =BoxCox(soil_nutrition_2021$pH,lambda="auto"))

#LSD
leveneTest(boxcox_pH ~ treatment1, data = soil_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2021$boxcox_pH)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,pH~treatment1)
summary(aov_model_soil_nutrition_2021)
LSD_model_soil_nutrition_2021<-LSD.test(aov_model_soil_nutrition_2021,"treatment1",p.adj = "BH")
LSD_model_soil_nutrition_2021
#DUNCAN
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,pH~treatment1)
summary(aov_model_soil_nutrition_2021)
compare_means(data=soil_nutrition_2021,pH~treatment1,method = "anova")
duncan_result_soil_pH<- duncan.test(aov_model_soil_nutrition_2021,"treatment1")
duncan_result_soil_pH
#Non-parameter test
kruskal.test(pH~treatment1, data = soil_nutrition_2021)
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,boxcox_pH~treatment1)
dunnettT3Test(aov_model_soil_nutrition_2021,p.adjust.method = "BH")

#2022
leveneTest(pH ~ treatment1, data = soil_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2022$pH)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
soil_nutrition_2022<-soil_nutrition_2022%>%mutate(boxcox_pH =BoxCox(soil_nutrition_2022$pH,lambda="auto"))


leveneTest(boxcox_pH ~ treatment1, data = soil_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2022$boxcox_pH)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,pH~treatment1)
summary(aov_model_soil_nutrition_2022)
LSD_model_soil_nutrition_2022<-LSD.test(aov_model_soil_nutrition_2022,"treatment1",p.adj = "BH")
LSD_model_soil_nutrition_2022
#DUNCAN
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,pH~treatment1)
summary(aov_model_soil_nutrition_2022)
compare_means(data=soil_nutrition_2022,pH~treatment1,method = "anova")
duncan_result_soil_pH<- duncan.test(aov_model_soil_nutrition_2022,"treatment1")
duncan_result_soil_pH
#Non-parameter test
kruskal.test(pH~treatment1, data = soil_nutrition_2022)
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,boxcox_pH~treatment1)
dunnettT3Test(aov_model_soil_nutrition_2022,p.adjust.method = "BH")
###t.test 
var.test(soil_nutrition_2021$pH,soil_nutrition_2022$pH)#>0.05表示方差齐性
t.test(soil_nutrition_2021$pH,soil_nutrition_2022$pH,var.equal = T)#参数检验
compare_means(pH~time, data = soil_nutrition, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验

# 5.1.2 Plot #
pH_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "pH", fill = "treatment")+
  scale_fill_manual(values=c("#5C4A46","#DA9464","#888e4a"))+
  FacetTheme
pH_plot

####OM####
### 4.1 OM ###
## 4.1.1 statistical analysis##
soil_nutrition$treatment<-factor(soil_nutrition$treatment,levels=c("CK","CZB","NF"))
soil_nutrition$time<-factor(soil_nutrition$time,levels=c("2021","2022"))
OM_mean <- aggregate(soil_nutrition$OM, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=mean)
OM_sd <- aggregate(soil_nutrition$OM, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=sd)
OM_len <- aggregate(soil_nutrition$OM, by=list(soil_nutrition$treatment, soil_nutrition$time), FUN=length)
df_res <- data.frame(OM_mean, sd=OM_sd$x, len=OM_len$x)
colnames(df_res) = c("treatment", "time", "Mean", "Sd", "Count")
df_res
df_res$Se <- df_res$Sd/sqrt(df_res$Count)
#2021
leveneTest(OM ~ treatment1, data = soil_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2021$OM)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
soil_nutrition_2021<-soil_nutrition_2021%>%mutate(boxcox_OM =BoxCox(soil_nutrition_2021$OM,lambda="auto"))


leveneTest(boxcox_OM ~ treatment1, data = soil_nutrition_2021)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2021$boxcox_OM)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,OM~treatment1)
summary(aov_model_soil_nutrition_2021)
LSD_model_soil_nutrition_2021<-LSD.test(aov_model_soil_nutrition_2021,"treatment1",p.adj = "BH")
LSD_model_soil_nutrition_2021
#DUNCAN
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,OM~treatment1)
summary(aov_model_soil_nutrition_2021)
compare_means(data=soil_nutrition_2021,OM~treatment1,method = "anova")
duncan_result_soil_OM<- duncan.test(aov_model_soil_nutrition_2021,"treatment1")
duncan_result_soil_OM
#Non-parameter test
kruskal.test(OM~treatment1, data = soil_nutrition_2021)
aov_model_soil_nutrition_2021<-aov(data=soil_nutrition_2021,boxcox_OM~treatment1)
dunnettT3Test(aov_model_soil_nutrition_2021,p.adjust.method = "BH")

#2022
leveneTest(OM ~ treatment1, data = soil_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2022$OM)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
soil_nutrition_2022<-soil_nutrition_2022%>%mutate(boxcox_OM =BoxCox(soil_nutrition_2022$OM,lambda="auto"))


leveneTest(boxcox_OM ~ treatment1, data = soil_nutrition_2022)#p>0.05，则满足方差齐性
shapiro.test(soil_nutrition_2022$boxcox_OM)#p<0.05 indicates skewed distribution, p>0.05 indicates normal distribution
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,OM~treatment1)
summary(aov_model_soil_nutrition_2022)
LSD_model_soil_nutrition_2022<-LSD.test(aov_model_soil_nutrition_2022,"treatment1",p.adj = "BH")
LSD_model_soil_nutrition_2022
#DUNCAN
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,OM~treatment1)
summary(aov_model_soil_nutrition_2022)
compare_means(data=soil_nutrition_2022,OM~treatment1,method = "anova")
duncan_result_soil_OM<- duncan.test(aov_model_soil_nutrition_2022,"treatment1")
duncan_result_soil_OM
#Non-parameter test
kruskal.test(OM~treatment1, data = soil_nutrition_2022)
aov_model_soil_nutrition_2022<-aov(data=soil_nutrition_2022,boxcox_OM~treatment1)
dunnettT3Test(aov_model_soil_nutrition_2022,p.adjust.method = "BH")
###t.test 
var.test(soil_nutrition_2021$OM,soil_nutrition_2022$OM)#>0.05表示方差齐性
t.test(soil_nutrition_2021$OM,soil_nutrition_2022$OM,var.equal = T)#参数检验
compare_means(OM~time, data = soil_nutrition, 
              ref.group = ".all.", 
              method = "wilcox.test")#非参数检验

# 5.1.2 Plot #
OM_plot <- ggplot(df_res, aes(x=time, y=Mean, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(),
           color="black", width=.8) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean +Se),
                position=position_dodge(.8), width=.2) +
  labs(x = "Year", y = "Organic matter/g·kg-1", fill = "treatment")+
  scale_fill_manual(values=c("#5C4A46","#DA9464","#888e4a"))+
  FacetTheme
OM_plot

soil_pHOM <- ggarrange(pH_plot,OM_plot,ncol = 2, labels = c("A", "B"),common.legend = T,legend = "right")
soil_pHOM


setwd(wdOutput)
getwd()
ggsave(paste("soil_pHOM",".pdf",sep=""),soil_pHOM,
       device=cairo_pdf,width=240,height=120,dpi = 600,units = "mm")

