setwd("D:/骨肉瘤/骨肉瘤单细胞")

library(Seurat)
library(dplyr)
library(monocle)
library(qs)
library(patchwork)
library(ggpubr)

seurat.data = qread(file = "tumor.qs")

current_celltypes <- seurat.data@meta.data$celltype

# 定义需要改为SH3-High的类别
high_celltypes <- c("C9", "C2", "C4", "C10", "C8")

# 如果你想直接替换原celltype列而不是创建新列
seurat.data@meta.data$celltype <- ifelse(current_celltypes %in% high_celltypes, 
                                        "SH3-High", 
                                        "SH3-Low")
sce<- seurat.data

HSMM <- as.CellDataSet(sce)
HSMM

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

HSMM <- detectGenes(HSMM, min_expr = 1)
print(head(fData(HSMM)))

expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))

head(pData(HSMM))

length(VariableFeatures(sce))

HSMM <- setOrderingFilter(HSMM, VariableFeatures(sce))
plot_ordering_genes(HSMM)

disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
length(unsup_clustering_genes$gene_id)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM,
                        max_components = 2,
                        num_dim = 20,
                        #residualModelFormulaStr = "~SampleID", #如果存在批次则指定批次
                        method = 'DDRTree') # DDRTree方式
HSMM <- orderCells(HSMM)
pData(HSMM) %>% head()

GM_state <- function(cds, starting_point, cluster){
  if (length(unique(cds$State)) > 1){
    T0_counts <- table(cds$State, cds@phenoData@data[,cluster])[,starting_point]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

root_start = GM_state(cds = HSMM,starting_point = "SH3-Low",cluster = "celltype")
root_start
HSMM <- monocle::orderCells(HSMM, root_state = root_start)

colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

#主题
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

# 生成三幅轨迹图
a1 <- plot_cell_trajectory(HSMM, color_by = "celltype") + 
  scale_color_manual(values = colour)

a2 <- plot_cell_trajectory(HSMM, color_by = "State") + 
  scale_color_manual(values = colour)

a3 <- plot_cell_trajectory(HSMM, color_by = "Pseudotime") + 
  ggsci::scale_color_gsea()

# 设置图形输出尺寸
options(repr.plot.width = 8, repr.plot.height = 4.5)

# 合并图形
combined_plot <- wrap_plots(a1, a2, a3, ncol = 3)  # 确保使用patchwork包的wrap_plots

# 保存为PDF
pdf("cell_trajectory_plots.pdf", width = 14, height = 4.5)  # 尺寸与options设置一致
print(combined_plot)
dev.off()

input.data = data.frame(celltype = HSMM$celltype,
                        Pseudotime = HSMM$Pseudotime)

options(repr.plot.width = 4.5, repr.plot.height = 4)
ggboxplot(data = input.data,
          fill = "celltype", 
          x = "celltype",  
          y = "Pseudotime")+mytheme

input.data = data.frame(celltype = HSMM$celltype,
                        Pseudotime = HSMM$Pseudotime)
## 基因
input.data$SH3BGRL3 = as.numeric(GetAssayData(object = sce, assay = "RNA",slot = "data")["SH3BGRL3",])
options(repr.plot.width = 6, repr.plot.height = 4)
## Vis
ggplot(input.data, aes(x = Pseudotime, y = SH3BGRL3)) +
  labs(x="Pseudotime",y = "Expression level (log2)")+
  ggpubr::stat_cor(label.sep = "\n",
                   label.y = 3,
                   label.y.npc = "top",
                   size = 4,
                   method = "sp")  +
  geom_smooth(method='loess',size=0.8, color = "black") + theme_bw() +
  ggfun::facet_set(label = "SH3BGRL3")+
  mytheme + theme(legend.position = "none")