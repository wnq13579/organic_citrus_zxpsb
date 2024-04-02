#### 1. Loading package####
library(dplyr)
library(linkET)
library(ggplot2)

#### 2 Import and process data ####
wdImport <- c("F:/whx_24samples_filtered1/01.import")
wdOutput <- c("F:/whx_24samples_filtered1/06.env/soil_micro")
data.set.name = '_Filtered1' 


setwd(wdImport)
species <- read.table("24Sample_ASVPercent_filtered_All.txt",comment.char="",check.names=F,stringsAsFactors=F, row.names=1,header = TRUE,sep="")#读取物种数据
fruit_quality <- read.table("fruit_quality_filter.txt",sep="\t",na.strings="",header = T,row.names=1,comment.char = "",check.names = F,stringsAsFactors = F)#读取果实环境因子
leaf_nutrition <- read.table("env_leaf_nutrition.txt",sep="\t",na.strings="",header = T,row.names=1,comment.char = "",check.names = F,stringsAsFactors = F)#读取叶片环境因子
soil_nutrition <- read.table("env_soil2.txt",sep="\t",na.strings="",header = T,row.names=1,comment.char = "",check.names = F,stringsAsFactors = F)#读取土壤环境因子


####3 fruit####

fruit_quality <- fruit_quality[c(10:21),]
### mantel test

fruit_quality_mantel <-mantel_test(species, fruit_quality,
                                   spec_select = list(bulk =1:3078,
                                                      rhizosphere =3079:6156),
                                   spec_dist =  dist_func(.FUN = "vegdist", method = "bray"), # 样本距离使用的vegdist()计算，可以选择适合自己数据的距离指数。
                                   env_dist = dist_func(.FUN = "vegdist", method = "euclidean")
                                   ) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
fruit_quality_mantel
#绘图
fruit_quality_plot<-qcorrplot(correlate(fruit_quality,method = 'spearman'), type = "upper", diag = TRUE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), data = fruit_quality_mantel, curvature = 0.1) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(5, "RdBu")) +#otu与环境因子之间连线的颜色
  scale_size_manual(values = c(1, 2, 3)) +#连线的宽度
  scale_fill_gradientn(colors = c('#0571B0', '#A3CDE2', 'white', '#EC8768', '#D74437'), limits = c(-1, 1)) +  #根据 Spearman 相关指定热图颜色
  scale_colour_manual(values = color_pal(3)) +
  geom_mark(# 添加r值与显著性标记
    sep = '\n', 
    size = 5, 
    only_mark = T,
    sig_level = c(0.05, 0.01), # 显著性水平设置
    sig_thres = 0.05 # 显著性阈值，p值大于阈值的相关性系数不会被绘制。
  ) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey40"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 2),
         fill = guide_colorbar(title = "Spearman's ρ", order = 3)) +
  labs(title = "Mantel test")
fruit_quality_plot
setwd(wdOutput)
getwd()
ggsave(paste("fruit_quality_plot",".pdf",sep=""),fruit_quality_plot,
       device=cairo_pdf,width=270,height=240,dpi = 600,units = "mm")

####4 leaf####

leaf_nutrition <-leaf_nutrition[c(10:21),c(1:11)]
### mantel test
leaf_nutrition_mantel <- mantel_test(species, leaf_nutrition,
                                     spec_select = list(bulk =1:3078,
                                                        rhizosphere =3079:6156),
                                                        spec_dist =  dist_func(.FUN = "vegdist", method = "bray"), # 样本距离使用的vegdist()计算，可以选择适合自己数据的距离指数。
                                                        env_dist = dist_func(.FUN = "vegdist", method = "euclidean")
                                                        ) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
leaf_nutrition_mantel
#绘图
leaf_nutrition_plot<-qcorrplot(correlate(leaf_nutrition,method = 'spearman'), type = "upper", diag = TRUE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), data = leaf_nutrition_mantel, curvature = 0.1) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(5, "RdBu")) +#otu与环境因子之间连线的颜色
  scale_fill_gradientn(colors = c('#0571B0', '#A3CDE2', 'white', '#EC8768', '#D74437'), limits = c(-1, 1)) +  #根据 Spearman 相关指定热图颜色
  scale_size_manual(values = c(1, 2, 3)) +#连线的宽度
  scale_colour_manual(values = color_pal(3)) +
  geom_mark(# 添加r值与显著性标记
    sep = '\n', 
    size = 5, 
    only_mark = T,
    sig_level = c(0.05, 0.01), # 显著性水平设置
    sig_thres = 0.05 # 显著性阈值，p值大于阈值的相关性系数不会被绘制。
  ) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey40"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 2),
         fill = guide_colorbar(title = "Spearman's ρ", order = 3)) +
  labs(title = "Mantel test")
leaf_nutrition_plot
setwd(wdOutput)
getwd()
ggsave(paste("leaf_nutrition_plot",".pdf",sep=""),leaf_nutrition_plot,
       device=cairo_pdf,width=270,height=240,dpi = 600,units = "mm")
####5 soil####

soil_nutrition <-soil_nutrition[,c(1:15)]
### mantel test
soil_nutrition_mantel <-mantel_test(species, soil_nutrition,
                                    spec_select = list(bulk =1:3078,
                                                       rhizosphere =3079:6156),
                                    spec_dist =  dist_func(.FUN = "vegdist", method = "bray"), # 样本距离使用的vegdist()计算，可以选择适合自己数据的距离指数。
                                    env_dist = dist_func(.FUN = "vegdist", method = "euclidean")
) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
soil_nutrition_mantel
#绘图
soil_nutrition_plot<-qcorrplot(correlate(soil_nutrition,method = 'spearman'), type = "upper", diag = TRUE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), data = soil_nutrition_mantel, curvature = 0.1) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(5, "RdBu")) +#otu与环境因子之间连线的颜色
  scale_fill_gradientn(colors = c('#0571B0', '#A3CDE2', 'white', '#EC8768', '#D74437'), limits = c(-1, 1)) +  #根据 Spearman 相关指定热图颜色
  scale_size_manual(values = c(1, 2, 3)) +#连线的宽度
  scale_colour_manual(values = color_pal(3)) +
  geom_mark(# 添加r值与显著性标记
    sep = '\n', 
    size = 5, 
    only_mark = T,
    sig_level = c(0.05, 0.01), # 显著性水平设置
    sig_thres = 0.05 # 显著性阈值，p值大于阈值的相关性系数不会被绘制。
  ) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey40"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 2),
         fill = guide_colorbar(title = "Spearman's ρ", order = 3)) +
  labs(title = "Mantel test")
soil_nutrition_plot
setwd(wdOutput)
getwd()
ggsave(paste("soil_nutrition_plot",".pdf",sep=""),soil_nutrition_plot,
       device=cairo_pdf,width=270,height=240,dpi = 600,units = "mm")

