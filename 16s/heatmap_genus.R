library(pheatmap)
library(psych)
library(stringr)

wdImport <- c("E:/whx_24samples_filtered1/01.import")
wdOutput <- c("E:/whx_24samples_filtered1/06.env/soil_micro")
data.set.name = '_Filtered1' 

####soil####

#bulk
setwd(wdImport)
species_table <- read.table("24Sample_GenusPercent_filtered3.txt",comment.char="",check.names=F,stringsAsFactors=F, row.names=1,header = TRUE,sep="")
map<-read.table("env_soil1.txt",sep="\t",na.strings="",header = T,row.names=1,comment.char = "",check.names = F,stringsAsFactors = F)

map<-map[c(1:12),]
species_table<-species_table[c(1:44),c(1:4,9:12,17:20)]
colname_of_map<-colnames(map)

#将行名设置为属分类名
#rownames(species_table)<-species_table$Fimuly
#去除表格里不需要的注释信息，删除第一列（行名）和最后一列（物种注释信息）
#species_table<-species_table[,-c(1)]

#转置
species_table<-t(species_table)
sum_of_species<-colSums(species_table)
#重排species_table行的顺序，使其和map行名一致
species_table<-species_table[match(rownames(map),rownames(species_table)),]
#将物种丰度和环境因子信息合并
merged_table<-data.frame(map,species_table,check.names = F,check.rows = T)

#计算spearman秩相关
correlation_results<-corr.test(merged_table,method ="spearman",adjust="fdr")
#计算pearson相关（此处我们使用spearman秩相关进行后续的分析，需要计算Pearson相关系数的可以参考此代码，不需要的请忽略）
#correlation_results<-corr.test(merged_table,method ="pearson",adjust="fdr")
#method ="spearman"指明用秩相关的方法
#adjust="fdr"校正错误发现率
#提取相关矩阵和p值矩阵
r<-correlation_results$r
p<-correlation_results$p
#剔除环境因子之间、微生物之间的相关系数，因为此处我们只需要环境因子和微生物的相关系数）
r<-r[-c(1:13),-c(14:113)]
p<-p[-c(1:13),-c(14:113)]
#选择相关系数显著次数最多的20个物种
selected_position_of_species<-head(order(colSums(t(p)<0.05),decreasing =T),30)
#得到筛选后的相关系数，p值矩阵
r<-r[selected_position_of_species,]
p<-p[selected_position_of_species,]

#自定义显著性标记函数，此处显著性标记使用“*”，你也可以根据自己的习惯设置其他的显著性标记。
sig_label<-function(x){ifelse(x<0.001,"**",ifelse(x<0.01,"**",ifelse(x<0.05,"*","")))}
#得到显著性标记矩阵
sig_matrix<-sig_label(p)

soil_pheatmap_genus<-pheatmap(r,fontsize=15,border_color = "black",
                        display_numbers = sig_matrix,fontsize_row =15,fontsize_col = 15,
                        fontsize_number = 22,
                        #显著性标记的符号大小
                        cluster_rows=T,clustering_distance_rows="correlation",
                        #指明行聚类，聚类依据的距离
                        cluster_cols=F,clustering_distance_cols="euclidean",
                        #指明列聚类，聚类依据的距离
                        clustering_method="centroid")
                        #聚类方法


soil_pheatmap_genus

setwd(wdOutput)
getwd()
ggsave(paste("soil_pheatmap_genus_bulk",".pdf",sep=""),soil_pheatmap_genus,
       device=cairo_pdf,width=270,height=270,dpi = 600,units = "mm")
#rhizo
setwd(wdImport)
species_table <- read.table("24Sample_GenusPercent_filtered3.txt",comment.char="",check.names=F,stringsAsFactors=F, row.names=1,header = TRUE,sep="")
map<-read.table("env_soil1.txt",sep="\t",na.strings="",header = T,row.names=1,comment.char = "",check.names = F,stringsAsFactors = F)

map<-map[c(13:24),c(1:3)]
species_table<-species_table[c(1:44),c(5:8,13:16,21:24)]
colname_of_map<-colnames(map)

#将行名设置为属分类名
#rownames(species_table)<-species_table$Fimuly
#去除表格里不需要的注释信息，删除第一列（行名）和最后一列（物种注释信息）
#species_table<-species_table[,-c(1)]

#转置
species_table<-t(species_table)
sum_of_species<-colSums(species_table)
#重排species_table行的顺序，使其和map行名一致
species_table<-species_table[match(rownames(map),rownames(species_table)),]
#将物种丰度和环境因子信息合并
merged_table<-data.frame(map,species_table,check.names = F,check.rows = T)

#计算spearman秩相关
correlation_results<-corr.test(merged_table,method ="spearman",adjust="fdr")
#计算pearson相关（此处我们使用spearman秩相关进行后续的分析，需要计算Pearson相关系数的可以参考此代码，不需要的请忽略）
#correlation_results<-corr.test(merged_table,method ="pearson",adjust="fdr")
#method ="spearman"指明用秩相关的方法
#adjust="fdr"校正错误发现率
#提取相关矩阵和p值矩阵
r<-correlation_results$r
p<-correlation_results$p
#剔除环境因子之间、微生物之间的相关系数，因为此处我们只需要环境因子和微生物的相关系数）
r<-r[-c(1:3),-c(4:103)]
p<-p[-c(1:3),-c(4:103)]
#选择相关系数显著次数最多的20个物种
selected_position_of_species<-head(order(colSums(t(p)<0.05),decreasing =T),11)
#得到筛选后的相关系数，p值矩阵
r<-r[selected_position_of_species,]
p<-p[selected_position_of_species,]

#自定义显著性标记函数，此处显著性标记使用“*”，你也可以根据自己的习惯设置其他的显著性标记。
sig_label<-function(x){ifelse(x<0.001,"**",ifelse(x<0.01,"**",ifelse(x<0.05,"*","")))}
#得到显著性标记矩阵
sig_matrix<-sig_label(p)

soil_pheatmap_genus<-pheatmap(r,fontsize=15,border_color = "black",
                              display_numbers = sig_matrix,fontsize_row =15,fontsize_col = 15,
                              fontsize_number = 22,
                              #显著性标记的符号大小
                              cluster_rows=T,clustering_distance_rows="correlation",
                              #指明行聚类，聚类依据的距离
                              cluster_cols=F,clustering_distance_cols="euclidean",
                              #指明列聚类，聚类依据的距离
                              clustering_method="centroid")
#聚类方法

soil_pheatmap_genus

setwd(wdOutput)
getwd()
ggsave(paste("soil_pheatmap_genus_rhizo",".pdf",sep=""),soil_pheatmap_genus,
       device=cairo_pdf,width=240,height=150,dpi = 600,units = "mm")
####leaf####
#bulk
setwd(wdImport)
species_table <- read.table("24Sample_GenusPercent_filtered3.txt",comment.char="",check.names=F,stringsAsFactors=F, row.names=1,header = TRUE,sep="")
map<-read.table("env_leaf_nutrition.txt",sep="\t",na.strings="",header = T,row.names=1,comment.char = "",check.names = F,stringsAsFactors = F)
map<-map[10:21,1:11]
species_table<-species_table[c(1:44),c(1:4,9:12,17:20)]
colname_of_map<-colnames(map)

#将行名设置为属分类名
#rownames(species_table)<-species_table$Fimuly
#去除表格里不需要的注释信息，删除第一列（行名）和最后一列（物种注释信息）
#species_table<-species_table[,-c(1)]

#转置
species_table<-t(species_table)
sum_of_species<-colSums(species_table)
#重排species_table行的顺序，使其和map行名一致
species_table<-species_table[match(rownames(map),rownames(species_table)),]
#将物种丰度和环境因子信息合并
merged_table<-data.frame(map,species_table,check.names = F,check.rows = T)


#计算spearman秩相关
correlation_results<-corr.test(merged_table,method ="spearman",adjust="fdr")
#计算pearson相关（此处我们使用spearman秩相关进行后续的分析，需要计算Pearson相关系数的可以参考此代码，不需要的请忽略）
#correlation_results<-corr.test(merged_table,method ="pearson",adjust="fdr")
#method ="spearman"指明用秩相关的方法
#adjust="fdr"校正错误发现率
#提取相关矩阵和p值矩阵
r<-correlation_results$r
p<-correlation_results$p




#剔除环境因子之间、微生物之间的相关系数，因为此处我们只需要环境因子和微生物的相关系数）
r<-r[-c(1:11),-c(12:111)]
p<-p[-c(1:11),-c(12:111)]
#选择相关系数显著次数最多的20个物种
selected_position_of_species<-head(order(colSums(t(p)<0.05),decreasing =T),20)
#得到筛选后的相关系数，p值矩阵
r<-r[selected_position_of_species,]
p<-p[selected_position_of_species,]




#自定义显著性标记函数，此处显著性标记使用“*”，你也可以根据自己的习惯设置其他的显著性标记。
sig_label<-function(x){ifelse(x<0.001,"**",ifelse(x<0.01,"**",ifelse(x<0.05,"*","")))}
#得到显著性标记矩阵
sig_matrix<-sig_label(p)

leaf_pheatmap_genus<-pheatmap(r,fontsize=15,border_color = "black",
                              display_numbers = sig_matrix,fontsize_row =15,fontsize_col = 15,
                              fontsize_number = 22,
                              #显著性标记的符号大小
                              cluster_rows=T,clustering_distance_rows="correlation",
                              #指明行聚类，聚类依据的距离
                              cluster_cols=F,clustering_distance_cols="euclidean",
                              #指明列聚类，聚类依据的距离
                              clustering_method="centroid")
#聚类方法

leaf_pheatmap_genus

setwd(wdOutput)
getwd()
ggsave(paste("leaf_pheatmap_genus_bulk",".pdf",sep=""),leaf_pheatmap_genus,
       device=cairo_pdf,width=270,height=240,dpi = 600,units = "mm")
#rhizo
setwd(wdImport)
species_table <- read.table("24Sample_GenusPercent_filtered3.txt",comment.char="",check.names=F,stringsAsFactors=F, row.names=1,header = TRUE,sep="")
map<-read.table("env_leaf_nutrition1.txt",sep="\t",na.strings="",header = T,row.names=1,comment.char = "",check.names = F,stringsAsFactors = F)
map<-map[c(10:21),c(1:11)]
species_table<-species_table[c(1:44),c(5:8,13:16,21:24)]
colname_of_map<-colnames(map)

#将行名设置为属分类名
#rownames(species_table)<-species_table$Fimuly
#去除表格里不需要的注释信息，删除第一列（行名）和最后一列（物种注释信息）
#species_table<-species_table[,-c(1)]

#转置
species_table<-t(species_table)
sum_of_species<-colSums(species_table)
#重排species_table行的顺序，使其和map行名一致
species_table<-species_table[match(rownames(map),rownames(species_table)),]
#将物种丰度和环境因子信息合并
merged_table<-data.frame(map,species_table,check.names = F,check.rows = T)


#计算spearman秩相关
correlation_results<-corr.test(merged_table,method ="spearman",adjust="fdr")
#计算pearson相关（此处我们使用spearman秩相关进行后续的分析，需要计算Pearson相关系数的可以参考此代码，不需要的请忽略）
#correlation_results<-corr.test(merged_table,method ="pearson",adjust="fdr")
#method ="spearman"指明用秩相关的方法
#adjust="fdr"校正错误发现率
#提取相关矩阵和p值矩阵
r<-correlation_results$r
p<-correlation_results$p




#剔除环境因子之间、微生物之间的相关系数，因为此处我们只需要环境因子和微生物的相关系数）
r<-r[-c(1:11),-c(12:111)]
p<-p[-c(1:11),-c(12:111)]
#选择相关系数显著次数最多的20个物种
selected_position_of_species<-head(order(colSums(t(p)<0.05),decreasing =T),18)
#得到筛选后的相关系数，p值矩阵
r<-r[selected_position_of_species,]
p<-p[selected_position_of_species,]




#自定义显著性标记函数，此处显著性标记使用“*”，你也可以根据自己的习惯设置其他的显著性标记。
sig_label<-function(x){ifelse(x<0.001,"**",ifelse(x<0.01,"**",ifelse(x<0.05,"*","")))}
#得到显著性标记矩阵
sig_matrix<-sig_label(p)

leaf_pheatmap_genus<-pheatmap(r,fontsize=15,border_color = "black",
                              display_numbers = sig_matrix,fontsize_row =15,fontsize_col = 15,
                              fontsize_number = 22,
                              #显著性标记的符号大小
                              cluster_rows=T,clustering_distance_rows="correlation",
                              #指明行聚类，聚类依据的距离
                              cluster_cols=F,clustering_distance_cols="euclidean",
                              #指明列聚类，聚类依据的距离
                              clustering_method="centroid")
#聚类方法

leaf_pheatmap_genus

setwd(wdOutput)
getwd()
ggsave(paste("leaf_pheatmap_genus_rhizo",".pdf",sep=""),leaf_pheatmap_genus,
       device=cairo_pdf,width=270,height=240,dpi = 600,units = "mm")
####fruit####
#bulk
setwd(wdImport)
species_table <- read.table("24Sample_GenusPercent_filtered3.txt",comment.char="",check.names=F,stringsAsFactors=F, row.names=1,header = TRUE,sep="")
map<-read.table("fruit_quality_filter.txt",sep="\t",na.strings="",header = T,row.names=1,comment.char = "",check.names = F,stringsAsFactors = F)
map<-map[10:21,]
species_table<-species_table[c(1:44),c(1:4,9:12,17:20)]
colname_of_map<-colnames(map)

#将行名设置为属分类名
#rownames(species_table)<-species_table$Fimuly
#去除表格里不需要的注释信息，删除第一列（行名）和最后一列（物种注释信息）
#species_table<-species_table[,-c(1)]

#转置
species_table<-t(species_table)
sum_of_species<-colSums(species_table)
#重排species_table行的顺序，使其和map行名一致
species_table<-species_table[match(rownames(map),rownames(species_table)),]
#将物种丰度和环境因子信息合并
merged_table<-data.frame(map,species_table,check.names = F,check.rows = T)


#计算spearman秩相关
correlation_results<-corr.test(merged_table,method ="spearman",adjust="fdr")
#计算pearson相关（此处我们使用spearman秩相关进行后续的分析，需要计算Pearson相关系数的可以参考此代码，不需要的请忽略）
#correlation_results<-corr.test(merged_table,method ="pearson",adjust="fdr")
#method ="spearman"指明用秩相关的方法
#adjust="fdr"校正错误发现率
#提取相关矩阵和p值矩阵
r<-correlation_results$r
p<-correlation_results$p




#剔除环境因子之间、微生物之间的相关系数，因为此处我们只需要环境因子和微生物的相关系数）
r<-r[-c(1:12),-c(13:122)]
p<-p[-c(1:12),-c(13:122)]
#选择相关系数显著次数最多的20个物种
selected_position_of_species<-head(order(colSums(t(p)<0.05),decreasing =T),21)
#得到筛选后的相关系数，p值矩阵
r<-r[selected_position_of_species,]
p<-p[selected_position_of_species,]




#自定义显著性标记函数，此处显著性标记使用“*”，你也可以根据自己的习惯设置其他的显著性标记。
sig_label<-function(x){ifelse(x<0.001,"**",ifelse(x<0.01,"**",ifelse(x<0.05,"*","")))}
#得到显著性标记矩阵
sig_matrix<-sig_label(p)

fruit_pheatmap_genus<-pheatmap(r,fontsize=15,border_color = "black",
                              display_numbers = sig_matrix,fontsize_row =15,fontsize_col = 15,
                              fontsize_number = 22,
                              #显著性标记的符号大小
                              cluster_rows=T,clustering_distance_rows="correlation",
                              #指明行聚类，聚类依据的距离
                              cluster_cols=F,clustering_distance_cols="euclidean",
                              #指明列聚类，聚类依据的距离
                              clustering_method="centroid")
#聚类方法

fruit_pheatmap_genus

setwd(wdOutput)
getwd()
ggsave(paste("fruit_pheatmap_genus_bulk",".pdf",sep=""),fruit_pheatmap_genus,
       device=cairo_pdf,width=270,height=240,dpi = 600,units = "mm")

#rhizo
setwd(wdImport)
species_table <- read.table("24Sample_GenusPercent_filtered3.txt",comment.char="",check.names=F,stringsAsFactors=F, row.names=1,header = TRUE,sep="")
map<-read.table("fruit_quality_filter1.txt",sep="\t",na.strings="",header = T,row.names=1,comment.char = "",check.names = F,stringsAsFactors = F)
map<-map[c(10:21),]
species_table<-species_table[c(1:50),c(5:8,13:16,21:24)]
colname_of_map<-colnames(map)

#将行名设置为属分类名
#rownames(species_table)<-species_table$Fimuly
#去除表格里不需要的注释信息，删除第一列（行名）和最后一列（物种注释信息）
#species_table<-species_table[,-c(1)]

#转置
species_table<-t(species_table)
sum_of_species<-colSums(species_table)
#重排species_table行的顺序，使其和map行名一致
species_table<-species_table[match(rownames(map),rownames(species_table)),]
#将物种丰度和环境因子信息合并
merged_table<-data.frame(map,species_table,check.names = F,check.rows = T)


#计算spearman秩相关
correlation_results<-corr.test(merged_table,method ="spearman",adjust="fdr")
#计算pearson相关（此处我们使用spearman秩相关进行后续的分析，需要计算Pearson相关系数的可以参考此代码，不需要的请忽略）
#correlation_results<-corr.test(merged_table,method ="pearson",adjust="fdr")
#method ="spearman"指明用秩相关的方法
#adjust="fdr"校正错误发现率
#提取相关矩阵和p值矩阵
r<-correlation_results$r
p<-correlation_results$p




#剔除环境因子之间、微生物之间的相关系数，因为此处我们只需要环境因子和微生物的相关系数）
r<-r[-c(1:12),-c(13:122)]
p<-p[-c(1:12),-c(13:122)]
#选择相关系数显著次数最多的20个物种
selected_position_of_species<-head(order(colSums(t(p)<0.05),decreasing =T),19)
#得到筛选后的相关系数，p值矩阵
r<-r[selected_position_of_species,]
p<-p[selected_position_of_species,]




#自定义显著性标记函数，此处显著性标记使用“*”，你也可以根据自己的习惯设置其他的显著性标记。
sig_label<-function(x){ifelse(x<0.001,"**",ifelse(x<0.01,"**",ifelse(x<0.05,"*","")))}
#得到显著性标记矩阵
sig_matrix<-sig_label(p)

fruit_pheatmap_genus<-pheatmap(r,fontsize=15,border_color = "black",
                              display_numbers = sig_matrix,fontsize_row =15,fontsize_col = 15,
                              fontsize_number = 22,
                              #显著性标记的符号大小
                              cluster_rows=T,clustering_distance_rows="correlation",
                              #指明行聚类，聚类依据的距离
                              cluster_cols=F,clustering_distance_cols="euclidean",
                              #指明列聚类，聚类依据的距离
                              clustering_method="centroid")
#聚类方法

fruit_pheatmap_genus

setwd(wdOutput)
getwd()
ggsave(paste("fruit_pheatmap_genus_rhizo",".pdf",sep=""),fruit_pheatmap_genus,
       device=cairo_pdf,width=270,height=240,dpi = 600,units = "mm")