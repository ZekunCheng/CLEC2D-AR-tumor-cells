library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(future)
library(qs)

seurat.data = qread(file = "D:/骨肉瘤/骨肉瘤单细胞/Step2.PBMC_afterQC.qs")
seurat.data

#标准化
seurat.data <- seurat.data %>% NormalizeData(verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F) %>% 
  ScaleData(verbose = F)

#降维和聚类
seurat.data = seurat.data %>% 
  RunPCA(npcs = 30, verbose = F) %>% 
  #RunTSNE(reduction = "pca", dims = 1:30, verbose = F) %>% 
  RunUMAP(reduction = "pca", dims = 1:30, verbose = F)

#检查批次
options(repr.plot.width = 10, repr.plot.height = 4.5)
p1.compare=wrap_plots(ncol = 2,
                      DimPlot(seurat.data, reduction = "pca", group.by = "sampleID")+NoAxes()+ggtitle("Before_PCA"),
                      DimPlot(seurat.data, reduction = "umap", group.by = "sampleID")+NoAxes()+ggtitle("Before_UMAP"),
                      guides = "collect"
)
p1.compare
ggsave(plot=p1.compare, filename="D:/骨肉瘤/骨肉瘤单细胞/Step3.Before_inter_sum.pdf", width = 10 ,height = 4.5)


library(harmony)
#RunHarmony
seurat.data <- seurat.data %>% RunHarmony("sampleID", plot_convergence = T)

seurat.data

#RunUMAP及聚类 
n.pcs = 20
seurat.data <- seurat.data %>% 
  RunUMAP(reduction = "harmony", dims = 1:n.pcs, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:n.pcs)
p2.compare=wrap_plots(ncol = 2,
                      DimPlot(seurat.data, reduction = "harmony", group.by = "sampleID")+NoAxes()+ggtitle("After_PCA (harmony)"),
                      DimPlot(seurat.data, reduction = "umap", group.by = "sampleID")+NoAxes()+ggtitle("After_UMAP"),
                      guides = "collect"
)
p2.compare

options(repr.plot.width = 10, repr.plot.height = 9)
wrap_plots(p1.compare, p2.compare, ncol = 1)
ggsave(plot=p2.compare, filename="D:/骨肉瘤/骨肉瘤单细胞/Step3.After_inter_Harmony.pdf", width = 10 ,height = 4.5)

#
for (res in c(0.05,0.1,0.3,0.5,0.8,1,1.2,1.4,1.5,2)){
print(res)
seurat.data <- FindClusters(seurat.data, resolution = res, algorithm = 1)%>% 
  identity()
}

options(repr.plot.width = 20, repr.plot.height = 8)
#umap可视化
cluster_umap <- wrap_plots(ncol = 5,
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.05", label = T) & NoAxes(),  
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.1", label = T) & NoAxes(),
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.3", label = T)& NoAxes(),
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.5", label = T) & NoAxes(),
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.8", label = T) & NoAxes(), 
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.1", label = T) & NoAxes(),
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.1.2", label = T) & NoAxes(),
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.1.4", label = T)& NoAxes(),
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.1.5", label = T)& NoAxes()
)
cluster_umap
ggsave(cluster_umap,filename = "D:/骨肉瘤/骨肉瘤单细胞/Step3.After_inter.cluster_umap_Harmony.pdf",
       width = 25, height = 9)

# 选择一个合适的分辨率
Idents(object = seurat.data) <- "RNA_snn_res.0.1"
options(repr.plot.width = 6, repr.plot.height = 5)
DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.1", label = T)& NoAxes()

DimPlot(seurat.data, reduction = "umap", label = T)& NoAxes()

#找marker基因
library(COSG)
marker_cosg <- COSG::cosg(
  seurat.data,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  expressed_pct=0.1,
  remove_lowly_expressed = T,
  n_genes_user=200)

write.csv(markers, file = "D:/骨肉瘤/骨肉瘤单细胞/Step3.COSG_res.csv")


#Marker气泡图
options(repr.plot.width = 7.5, repr.plot.height = 7)
check_genes = c("LYZ","S100A9","C1QA","CD68","APOE", #Myeloid cells
                "CD3D","CD3E","CD3G","TRAC", #T cells
                "FBLN1","ACTA2","TAGLN","COL3A1","COL6A1", #Fibroblasts (CAFs)
                "ALPL","RUNX2","CLEC11A", #Osteoblastic cells (OB cells)
                "CD79A","MS4A1","IGHM","CD19", #B cells
                "NKG7","GZMK","GZMA","GZMB", #NKT cells (excluding duplicated T cell markers)
                "VWF","CAV1","CLDN5","EGFL7","PECAM1", #Endothelial cells
                "SDC1", "PRDM1", "XBP1", "MZB1", "DERL3", # Core plasma cell markers
                "CTSK", "ACP5", "TCIRG1",
                "PRF1", "CCL5", "CCL4", "CCL3"
)

DotPlot(object = seurat.data, features = check_genes, assay = "RNA", scale = T) + 
  coord_flip()

#第一次注释
celltype=data.frame(ClusterID=0:10,celltype='NA')

celltype[celltype$ClusterID %in% c(0),2]='Myeloid cells'
celltype[celltype$ClusterID %in% c(5),2]='Myeloid cells'
celltype[celltype$ClusterID %in% c(1),2]='Osteoblasts'
celltype[celltype$ClusterID %in% c(2),2]='NKT'
celltype[celltype$ClusterID %in% c(3),2]='Myeloid cells'
celltype[celltype$ClusterID %in% c(4),2]='CAFs'
celltype[celltype$ClusterID %in% c(6),2]='Plasma cells'
celltype[celltype$ClusterID %in% c(7),2]='Endothelial cells'
celltype[celltype$ClusterID %in% c(8),2]='Osteoclasts'
celltype[celltype$ClusterID %in% c(9),2]='B cells'
celltype[celltype$ClusterID %in% c(10),2]='Treg-CTL'

colnames(celltype) = c("ClusterID","celltype_main")
seurat.data@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  seurat.data@meta.data[which(seurat.data@active.ident == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(seurat.data@meta.data$celltype)

#可视化
options(repr.plot.width = 6, repr.plot.height = 5)
DimPlot(seurat.data, reduction = "umap", group.by = "celltype", label = T)& NoAxes()

# 将celltype设置为默认插槽
Idents(object = seurat.data) <- "celltype"
# 保存
qsave(seurat.data, file = "D:/骨肉瘤/骨肉瘤单细胞/Step3.PBMC_annotation.qs")