library(Seurat)
library(harmony)
library(tidyverse)
library(dplyr)
library(patchwork)
library(tidydr)
library(ggplot2)
library(cowplot)
library(clusterProfiler)  
library(qs)
getwd()
source("./source_R/custom_seurat_functions.R")
options(repr.plot.width = 5, repr.plot.height = 6)
p4 <- TSNE.UMAP.Plot(seurat.tumor,
                     "celltype",
                     Style=3,
                     plot.title = NA,
                     legend.point.size = 3,
                     reduction = "umap",
                     label = T,
                     label.size = 3,
                     point.size = 1e-1)+labs(title = "Vis3")
p4

seurat.os = qread(file = "D:/骨肉瘤/骨肉瘤单细胞/成骨细胞亚群.qs")
# setwd("单细胞打分/")
DimPlot(seurat.os)#umap图
celltype_mapping <- c(
  "C0" = "Malignant",
  "C1" = "non-Malignant",
  "C2" = "Malignant",
  "C3" = "Malignant",
  "C4" = "Malignant",
  "C5" = "Malignant",  
  "C6" = "Malignant",
  "C7" = "Malignant",
  "C8" = "Malignant",
  "C9" = "Malignant",
  "C10" = "Malignant"
)

# 使用 match() 或直接索引映射
seurat.os@meta.data$celltype1 <- celltype_mapping[seurat.os@meta.data$celltype]

# 检查结果
table(seurat.os@meta.data$celltype1, useNA = "always")
library(ggsci)
source("./source_R/custom_seurat_functions.R")
options(repr.plot.width = 5, repr.plot.height = 6)
p4 <- TSNE.UMAP.Plot(seurat.os,
                     "celltype1",
                     Style=3,
                     plot.title = NA,
                     legend.point.size = 3,
                     reduction = "umap",
                     label = T,
                     label.size = 3,
                     point.size = 1e-1)+labs(title = "Vis3")
p4
ggsave(
  filename = "UMAP_celltype.png",  # 文件名
  plot = p4,                      # 要保存的图形对象
  width = 5,                      # 宽度（英寸）
  height = 6,                     # 高度（英寸）
  dpi = 300,                      # 分辨率（默认 300）
  bg = "white"                    # 背景色（避免透明背景）
)


####自定义基因集打分####
proliferation_genes <- list(
  c("TP53", "MAPK8", "MAPK14", "TNF", "DDIT3", "EIF2S1", "MAPK1", "MAPK9", "MAP3K5", "EIF2AK2",
    "CASP3", "CXCL8", "MTOR", "MAPK11", "MAP3K20", "MAP2K6", "CERNA3", "MDM2", "MAP2K4", "ATF3",
    "GDF15", "SRC", "MAPK3", "EIF2AK4", "EGF", "ANG", "PPARG", "FAS", "BCL2L11", "MAP3K7",
    "TRAF2", "MEFV", "CASP8", "EIF4E", "AR", "APOA1", "E2F1", "SERPINE1", "ELAVL1", "HDAC1",
    "EGR1", "NLRP1", "CFLAR", "FASLG", "CASP7", "TNFSF10", "HDAC3", "RPS3", "RARA", "BMP6",
    "NR1H2", "MKNK1", "HCK", "MAP3K2", "CASP6", "MIR181D", "U2AF1", "MAGEB2", "RRN3", "MAGEA2",
    "SF1", "WDR4", "APOA1-AS")
)

# 2. 计算增殖评分
seurat.tumor <- AddModuleScore(
  object = seurat.tumor,
  features = proliferation_genes,
  ctrl = 100,      # 控制基因数量
  name = "Proliferation"  # 基础名称
)

# 3. 检查生成的列名（自动添加数字后缀）
print(grep("Proliferation", names(seurat.tumor@meta.data), value = TRUE))
proliferation_col <- grep("Proliferation", names(seurat.tumor@meta.data), value = TRUE)[1]

# 4. 可视化分析
# 4.1 小提琴图
vln_plot <- VlnPlot(
  seurat.tumor,
  features = proliferation_col,
  pt.size = 0,      # 不显示点
  adjust = 2,       # 调整平滑度
  group.by = "celltype" # 按分组显示
) + 
  labs(title = "Ribotoxic Stress Score Distribution") +
  theme(plot.title = element_text(hjust = 0.5))

# 4.2 点图
dot_plot <- DotPlot(
  seurat.tumor,
  features = proliferation_col,
  group.by = "celltype"
) + 
  coord_flip() +
  labs(y = "Groups", x = "") +
  theme(axis.text.y = element_blank())

# 4.3 UMAP可视化
umap_data <- FetchData(
  seurat.tumor,
  vars = c("UMAP_1", "UMAP_2", proliferation_col, "group")
)

umap_plot <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = .data[[proliferation_col]])) +
  geom_point(size = 0.8, alpha = 0.7) +
  scale_color_gradientn(
    colours = c("#330066", "#3366CC", "#00CCFF", "#66FF99", "#FFFF00", "#FF6600"),
    name = "Proliferation\nScore"
  ) +
  theme_bw() +
  labs(title = "Ribotoxic Stress Score on UMAP") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5)
  )

# 5. 组合图形
combined_plot <- (vln_plot | dot_plot) / 
  (umap_plot | p4) +  # 第二行并列显示umap_plot和p4
  plot_layout(heights = c(1, 2)) +   # 第一行高度:第二行高度 = 1:2
  plot_annotation(tag_levels = 'A')  # 添加子图标签(A, B, C, D)

proliferation_genes <- list(
  c("TP53", "MAPK8", "MAPK14", "TNF", "DDIT3", "EIF2S1", "MAPK1", "MAPK9", "MAP3K5", "EIF2AK2",
    "CASP3", "CXCL8", "MTOR", "MAPK11", "MAP3K20", "MAP2K6", "CERNA3", "MDM2", "MAP2K4", "ATF3",
    "GDF15", "SRC", "MAPK3", "EIF2AK4", "EGF", "ANG", "PPARG", "FAS", "BCL2L11", "MAP3K7",
    "TRAF2", "MEFV", "CASP8", "EIF4E", "AR", "APOA1", "E2F1", "SERPINE1", "ELAVL1", "HDAC1",
    "EGR1", "NLRP1", "CFLAR", "FASLG", "CASP7", "TNFSF10", "HDAC3", "RPS3", "RARA", "BMP6",
    "NR1H2", "MKNK1", "HCK", "MAP3K2", "CASP6", "MIR181D", "U2AF1", "MAGEB2", "RRN3", "MAGEA2",
    "SF1", "WDR4", "APOA1-AS")
)

# 2. 计算增殖评分
seurat.tumor <- AddModuleScore(
  object = seurat.tumor,
  features = proliferation_genes,
  ctrl = 100,      # 控制基因数量
  name = "Proliferation"  # 基础名称
)

# 3. 检查生成的列名（自动添加数字后缀）
print(grep("Proliferation", names(seurat.tumor@meta.data), value = TRUE))
proliferation_col <- grep("Proliferation", names(seurat.tumor@meta.data), value = TRUE)[1]

# 4. 可视化分析
# 4.1 小提琴图
vln_plot <- VlnPlot(
  seurat.tumor,
  features = proliferation_col,
  pt.size = 0,      # 不显示点
  adjust = 2,       # 调整平滑度
  group.by = "celltype" # 按分组显示
) + 
  labs(title = "Ribotoxic Stress Score Distribution") +
  theme(plot.title = element_text(hjust = 0.5))

# 4.2 点图
dot_plot <- DotPlot(
  seurat.tumor,
  features = proliferation_col,
  group.by = "celltype"
) + 
  coord_flip() +
  labs(y = "Groups", x = "") +
  theme(axis.text.y = element_blank())

# 4.3 UMAP可视化
umap_data <- FetchData(
  seurat.tumor,
  vars = c("UMAP_1", "UMAP_2", proliferation_col, "group")
)

umap_plot <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = .data[[proliferation_col]])) +
  geom_point(size = 0.8, alpha = 0.7) +
  scale_color_gradientn(
    colours = c("#330066", "#3366CC", "#00CCFF", "#66FF99", "#FFFF00", "#FF6600"),
    name = "Proliferation\nScore"
  ) +
  theme_bw() +
  labs(title = "Ribotoxic Stress Score on UMAP") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5)
  )

combined_plot <- (vln_plot | dot_plot) / 
  (umap_plot | p4) +  # 第二行并列显示umap_plot和p4
  plot_layout(heights = c(1, 2)) +   # 第一行高度:第二行高度 = 1:2
  plot_annotation(tag_levels = 'A')  # 添加子图标签(A, B, C, D)

ggsave("proliferation_analysis.png", combined_plot, width = 12, height = 10, dpi = 300)