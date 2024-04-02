library(UpSetR)
library(ggplot2)#作图 plot
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
wdOutput <- c("E:/whx_24samples_filtered1/04.different_abundance/venn")


#### 3. import ####
### 3.1 Import and process data ###
setwd(wdImport)
sets <- read.csv("ASV_unique_final_average.csv",header = TRUE,row.names=1,sep = ",")
upset_figure<-upset(fromList(sets),order.by = "freq",
      sets = c("NF_R","CZB_R","CK_R","NF_B","CZB_B","CK_B"),
      keep.order = TRUE)
upset_figure
#
setwd(wdOutput)
getwd()
ggsave(paste("upset_figure",".pdf",sep=""),upset_figure,
       device=cairo_pdf,width=210,height=120,dpi = 300,units = "mm")

