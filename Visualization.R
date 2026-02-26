library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(future)
library(qs)
library(ggsci)

seurat.data = qread(file = "D:/共病/Step3.PBMC_annotation.qs")
seurat.data

# 加载函数
source("D:/骨肉瘤/骨肉瘤单细胞/source_R/custom_seurat_functions.R")

options(repr.plot.width = 5, repr.plot.height = 6)
p4 <- TSNE.UMAP.Plot(seurat.data,
                     "celltype",
                     Style=3,
                     plot.title = NA,
                     legend.point.size = 3,
                     reduction = "umap",
                     label = T,
                     label.size = 3,
                     point.size = 1e-1)+labs(title = "Vis3")
p4

options(repr.plot.width = 5, repr.plot.height = 6)
p5 <- TSNE.UMAP.Plot(seurat.data,
                     "celltype",
                     Style=4,
                     plot.title = NA,
                     legend.point.size = 3,
                     reduction = "umap",
                     label = T,
                     label.size = 3,
                     point.size = 1e-1)+labs(title = "Vis4")
p5



options(repr.plot.width = 6, repr.plot.height = 2.5)
check_genes = c(
  # Epithelium (general epithelial markers)
  "EPCAM","KRT8","KRT18",
  
  # B cells
  "CD19","MS4A1","CD79A","CD79B","IGHM","CD74","PAX5","BANK1",
  # T cells
  "CD3D","CD3E","CD3G","TRAC",
  
  # Myeloid cells / Macrophages
  "LYZ","C1QA","C1QB","CD68","CD163",
  
  # Plasma cells
  "MZB1","JCHAIN","IGHA1","SDC1","PRDM1","XBP1","DERL3",
  #GC
  "S100A8","S100A9","FCGR3B","CEACAM8","MPO",
  # Mast cells
  "CPA3","TPSAB1","MS4A2",
  
  # Fibroblasts / CAFs / Stromal
  "FBLN1","ACTA2","TAGLN","COL1A2","COL3A1","COL6A1","ADAMDEC1"
)
p.dot = DotPlot_2(object = seurat.data,
                  Combine = F, dot.range.max = 3,
                  dot.range.min = 0, label.size = 3,
                  features = check_genes, legend.key.size = 0.4)&
  scale_color_distiller(palette = 'RdYlBu');p.dot