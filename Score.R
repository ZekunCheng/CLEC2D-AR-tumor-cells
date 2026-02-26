library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
#有云服务器的，可开启并运算，这里我用4个线程：
library(future)
library(qs)
library(ggsci)
setwd("D:/骨肉瘤/骨肉瘤单细胞")

seurat.data = qread(file = "tumor.qs")
seurat.data

FeaturePlot(seurat.data, features = c("SH3BGRL3"))

options(repr.plot.width = 10, repr.plot.height = 7)
#小提琴图
VlnPlot(seurat.data, features = "SH3BGRL3")
library(dplyr)

seurat.data@meta.data <- seurat.data@meta.data %>%
  mutate(celltype = case_when(
    celltype %in% c("C9", "C2", "C4", "C10", "C8") ~ "SH3-High",
    TRUE ~ "SH3-Low"
  ))
devtools::install_version("multicross", version = "2.1.0", repos = "http://cran.us.r-project.org")
devtools::install_github("jackbibby1/SCPA")
library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(qs)
library(ggsci)
library(ggpubr)
library(SCPA)
library(stringr)
library(ggrepel)

if(T){
  text.size = 12
  text.angle = 45
  text.hjust = 1
  legend.position = "right"
  mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                   axis.ticks = element_line(color = "black"),
                   axis.title = element_text(size = text.size,color ="black"), 
                   axis.text = element_text(size=text.size,color = "black"),
                   axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), #,vjust = 0.5
                   panel.grid=element_blank(), # 去网格线
                   legend.position = legend.position,
                   legend.text = element_text(size= text.size),
                   legend.title= element_text(size= text.size),
                   panel.border = element_rect(fill = NA),
  )
}

options(repr.plot.width = 6, repr.plot.height = 5)
p1 = DimPlot(seurat.data, 
             reduction = "umap",
             group.by = "celltype",
             label = T);p1

Pre.data <- seurat_extract(seurat.data,
                           meta1 = "celltype", value_meta1 = "SH3-Low")

Resp.data <- seurat_extract(seurat.data,
                            meta1 = "celltype", value_meta1 = "SH3-High")

pathways <- "./Rawdata/combined_metabolic_pathways.csv"
Resp.data <- compare_pathways(samples = list(Resp.data, Pre.data),
                              #downsample = 10000,
                              pathways = pathways)

Resp.data <- Resp.data %>%
  mutate(color = case_when(FC > 2 & adjPval < 0.05 ~ '#6dbf88',
                           FC < -2 & adjPval < 0.05 ~ 'mediumseagreen',
                           (FC < 2 | FC > -2) & adjPval < 0.05 ~ '#84b0f0',
                           (FC < 2 | FC > -2) & adjPval > 0.05 ~ 'black'))
head(Resp.data)

marker_path = Resp.data %>% 
  filter(grepl(pattern = "KEGG_OXIDATIVE_PHOSPHORYLATION", ignore.case = T, x = Pathway))
label_path = Resp.data %>% 
  filter(color == 'mediumseagreen') %>%
  top_n(5, qval)
options(repr.plot.width = 5, repr.plot.height = 5)
ggplot(Resp.data, aes(FC, qval)) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", col = 'black', lwd = 0.3) +
  geom_point(cex = 2.6, shape = 21, fill = Resp.data$color, stroke = 0.3) +
  geom_point(data = marker_path, shape = 21, cex = 2.8, fill = "red", color = "black", stroke = 0.3) +
  geom_text_repel(data = label_path, aes(x = FC, y = qval, label = Pathway))+
  xlab("Enrichment") +
  ylab("Qval") + 
  theme_bw()+
  mytheme

label_path$Pathway = str_replace(label_path$Pathway,
                                 pattern = "REACTOME_|KEGG_|HALLMARK_",
                                 replacement = "")

label_path$Pathway = tolower(label_path$Pathway)%>%
  str_to_title()%>%
  str_replace_all(pattern = "\\_",
                  replacement = " ")
ggplot(Resp.data, aes(FC, qval)) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", col = 'black', lwd = 0.3) +
  geom_point(cex = 2.6, shape = 21, fill = Resp.data$color, stroke = 0.3) +
  geom_point(data = marker_path, shape = 21, cex = 2.8,
             fill = "red", color = "black", stroke = 0.3
  ) +
  geom_text_repel(data = label_path, aes(x = FC, y = qval, label = Pathway), size = 3.5,
                  min.segment.length = 0, #始终为标签添加指引线段；若不想添加线段，则改为Inf
                  max.overlaps = Inf #排斥重叠过多标签，设置为Inf则可以保持始终显示所有标签
  )+
  labs(x = "Enrichment", y="Qval", title = "SH3-Low vs. SH3-High") +
  theme_bw()+
  mytheme

ggsave(filename = 'Step1.Scatter_diagram代谢重编程.pdf', width = 5, height = 5)

library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(qs)
library(ggsci)
library(ggpubr)
library(future)
library(stringr)
source("source_R/custom_seurat_functions.R")
if(T){
  text.size = 12
  text.angle = 45
  text.hjust = 1
  legend.position = "right"
  mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                   axis.ticks = element_line(color = "black"),
                   axis.title = element_text(size = text.size,color ="black"), 
                   axis.text = element_text(size=text.size,color = "black"),
                   axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), #,vjust = 0.5
                   panel.grid=element_blank(), # 去网格线
                   legend.position = legend.position,
                   legend.text = element_text(size= text.size),
                   legend.title= element_text(size= text.size)
  )
}
seurat.data <- sc.Pathway.Seurat(obj = seurat.data, 
                                 method = "AUCell", #可选方法："AUCell", "VISION", "ssGSEA","gsva"，单细胞推荐使用AUCell
                                 ncores = 4,
                                 assay.names = "metabolism",
                                 geneList = "./Rawdata/KEGG_metabolism_nc.gmt")

input.path = row.names(seurat.data@assays$metabolism@data)
input.path
devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis)
pdf(file = "Pathway_heatmap.pdf", height = 10, width = 7)
options(repr.plot.width = 10, repr.plot.height = 7)
AverageHeatmap(object = seurat.data, 
               assays = "metabolism",
               htRange = c(-2,0,2),
               markerGene = input.path,
               cluster_rows = T,
               markGenes = input.path[c(1:3,10:12,45,60,65)],
               clusterAnnoName = F,
               showRowNames = F,
               row_title = NULL)
FeaturePlot(seurat.data, features = c('Lipoic acid metabolism'))&
  scale_color_distiller(palette = 'RdBu')
ggsave(filename = '脂肪酸代谢FeaturePlot.pdf', width = 5, height = 4.5)

p2<-DoHeatmap(seurat.data, 
          features = rownames(seurat.data[["metabolism"]]), 
          assay = "metabolism",
          slot = "data")  # 根据数据选择 slot（data/scaled.data）
p2


options(repr.plot.width = 9, repr.plot.height = 3.5)
check_terms = input.path[1:86]
p.dot = DotPlot_2(object = seurat.data,
                  assay = "metabolism",
                  text.size = 12,
                  Combine = F, dot.range.max = 5,
                  dot.range.min = 0, label.size = 3,
                  features = check_terms, legend.key.size = 0.4)&
  scale_color_distiller(palette = 'RdYlBu')&coord_flip();p.dot
ggsave(filename = 'Step2.Dotplot.pdf', plot = p.dot, width = 9, height = 15)


options(repr.plot.width = 9, repr.plot.height = 3.5)
check_terms = input.path[1:86]

# 仅修改 group.by = "celltype"，其他参数原封不动
p.dot = DotPlot_2(
  object = seurat.data,
  assay = "metdata:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAbElEQVR4Xs2RQQrAMAgEfZgf7W9LAguybljJpR3wEse5JOL3ZObDb4x1loDhHbBOFU6i2Ddnw2KNiXcdAXygJlwE8OFVBHDgKrLgSInN4WMe9iXiqIVsTMjH7z/GhNTEibOxQswcYIWYOR/zAjBJfiXh3jZ6AAAAAElFTkSuQmCCabolism",
  group.by = "celltype",  # 唯一修改处：按 SH3-High/SH3-Low 分组
  text.size = 12,
  Combine = F, 
  dot.range.max = 5,
  dot.range.min = 0, 
  label.size = 3,
  features = check_terms, 
  legend.key.size = 0.4
) &
  scale_color_distiller(palette = 'RdYlBu') &
  coord_flip()

# 输出图形
p.dot
ggsave(filename = 'Low_VS_High.Dotplot.pdf', plot = p.dot, width = 9, height = 15)