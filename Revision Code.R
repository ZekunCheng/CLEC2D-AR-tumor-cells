fileID <- list.files("D:/CLEC2D返修单细胞/GSE274229_RAW/")
path <- "D:/CLEC2D返修单细胞/GSE274229_RAW/"

# 读取
seurat.list <- lapply(fileID, function(sample_id){
  
  data_dir <- file.path(path, sample_id)
  
  seurat_data <- Read10X(data.dir = data_dir)
  
  # 如果是多模态数据（list），默认取 Gene Expression
  if (is.list(seurat_data)) {
    message(sample_id, " contains multiple assays → using Gene Expression")
    seurat_data <- seurat_data[[grep("Gene", names(seurat_data), ignore.case = TRUE)]]
  }
  
  # 检查列名（细胞barcode）
  if (is.null(colnames(seurat_data))) {
    stop("No cell barcodes detected in: ", sample_id)
  }
  
  # 创建 Seurat 对象
  CreateSeuratObject(
    counts = seurat_data,
    min.cells = 3,
    min.features = 200,
    project = sample_id
  )
})

#看样本的信息
names(seurat.list) = fileID
seurat.list

#合并起来
merged_seurat <- merge(x = seurat.list[[1]],
                       y = seurat.list[-1],
                       add.cell.id = fileID)
merged_seurat

#要重新分亚组的时候再用
merged_seurat$sampleID = merged_seurat$orig.ident

#根据样本信息重命名编组
merged_seurat$group <- recode(merged_seurat$sampleID,
                              "CSPC1" = "CSPC",
                              "CSPC2" = "CSPC",
                              "CSPC3" = "CSPC",
                              "CSPC4" = "CSPC",
                              "CSPC5" = "CSPC",
                              "CSPC6" = "CSPC",
                              "CSPC7" = "CSPC",
                              "CSPC8" = "CSPC",
                              "CSPC9" = "CSPC",
                              "CSPC10" = "CSPC",
                              "CSPC11" = "CSPC",
                              "CSPC12" = "CSPC",
                              "CSPC13" = "CSPC",
                              "CSPC14" = "CSPC",
                              "CSPC15" = "CSPC",
                              "CSPC16" = "CSPC",
                              "CSPC17" = "CSPC",
                              "CSPC18" = "CSPC",
                              "CSPC19" = "CSPC",
                              "CSPC20" = "CSPC",
                              "CSPC21" = "CSPC",
                              "CSPC22" = "CSPC",
                              "CSPC23" = "CSPC",
                              "CSPC24" = "CSPC",
                              "CSPC25" = "CSPC",
                              "CRPC1" = "CRPC",
                              "CRPC2" = "CRPC",
                              "CRPC3" = "CRPC",
                              "CRPC4" = "CRPC",
                              "CRPC5" = "CRPC",
                              "CRPC6" = "CRPC")

saveRDS(merged_seurat,file = "D:/CLEC2D返修单细胞/Step1.RawCount_merged_seurat.rds")


library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(future)
library(qs)

seurat.data = read_rds(file = "D:/CLEC2D返修单细胞/Step1.RawCount_merged_seurat.rds")

#有提取亚组的时候改一下
os = subset(seurat.data, group %in% c("CSPC", "CRPC"))
os


#提取线粒体基因
mito_genes=rownames(os)[grep("^MT-", rownames(os))]
mito_genes

os[["percent.mt"]] <- PercentageFeatureSet(os, pattern = "^MT-")
head(os@meta.data, 5)

# 计算核糖体基因
ribo_genes=rownames(os)[grep("^RP[SL]", rownames(os),ignore.case = T)]
os=PercentageFeatureSet(os, "^RP[SL]",col.name = "percent.ribo")

# 计算红细胞基因
hb_genes <- rownames(os)[grep("^HB[^(P)]", rownames(os),ignore.case = T)]
os=PercentageFeatureSet(os, "^HB[^(P)]", col.name = "percent.hb")

##可视化
options(repr.plot.width=10, repr.plot.height=10)
VlnPlot(os,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb"),
        ncol = 3,
        group.by = "group")

#过滤
os.qc <- subset(os, subset = nFeature_RNA > 300 & nFeature_RNA < 7000 & percent.mt < 20 & percent.hb < 1)
os.qc

saveRDS(os.qc, "combined_QC_scRNA.rds")

library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(future)

seurat.data = read_rds(file = "D:/CLEC2D返修单细胞/combined_QC_scRNA.rds")
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
ggsave(plot=p1.compare, filename="D:/CLEC2D返修单细胞/Step3.Before_inter_sum.pdf", width = 10 ,height = 4.5)


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
ggsave(plot=p2.compare, filename="D:/CLEC2D返修单细胞/Step3.After_inter_Harmony.pdf", width = 10 ,height = 4.5)

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

write.csv(marker_cosg, file = "D:/CLEC2D返修单细胞/Step3.COSG_res.csv")


#Marker气泡图
options(repr.plot.width = 7.5, repr.plot.height = 7)
check_genes = c(
  "KRT8","KRT18","CLDN3","CLDN4","CLDN7","AGR2","HOXB13","TSPAN8", # Luminal epithelial
  "KRT5","KRT14","TP63","COL17A1","TRIM29","S100A2", # Basal epithelial
  "CD3D","CD3E","CD3G","CD2","CD7","NKG7","GZMA","GZMM","CD8A","TRAC", # T cells
  "CD79A","MS4A1","CD19","IGHM","IGKC","BANK1","FCRLA","JCHAIN", # B cells
  "CD68","TYROBP","C1QA","C1QB","C1QC","CSF1R","LYZ","TREM2","FCGR2A","LST1", # Myeloid cells
  "VWF","PECAM1","CLDN5","EMCN","CLEC14A","PLVAP","EGFL7","RAMP2", # Endothelial cells
  "DCN","LUM","COL1A1","COL1A2","COL3A1","FBLN1","DPT","SFRP2","MFAP4", # Fibroblasts (CAFs)
  "MYH11","ACTA2","CNN1","TAGLN","RGS5","PDGFRB","MCAM","LMOD1","MYL9","NOTCH3", # Perivascular smooth muscle
  "MKI67","UBE2C","TOP2A","CDK1","AURKB","CCNA2","NUSAP1","CDC20","TPX2","BIRC5", # Cycling cells
  "TPSAB1","TPSB2","CPA3","KIT","HDC","MS4A2","IL1RL1","SIGLEC6","SIGLEC8","LTC4S" # Mast cells
)


DotPlot(object = seurat.data, features = check_genes, assay = "RNA", scale = T) + 
  coord_flip()

#第一次注释
celltype=data.frame(ClusterID=0:9,celltype='NA')

celltype[celltype$ClusterID %in% c(0),2]='Luminal'
celltype[celltype$ClusterID %in% c(1),2]='T'
celltype[celltype$ClusterID %in% c(2),2]='Endothelial'
celltype[celltype$ClusterID %in% c(3),2]='Macrophage'
celltype[celltype$ClusterID %in% c(4),2]='Fibroblast'
celltype[celltype$ClusterID %in% c(5),2]='Basal'
celltype[celltype$ClusterID %in% c(6),2]='B'
celltype[celltype$ClusterID %in% c(7),2]='Smooth muscle'
celltype[celltype$ClusterID %in% c(8),2]='Cycling'
celltype[celltype$ClusterID %in% c(9),2]='Mast'


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
qsave(seurat.data, file = "D:/CLEC2D返修单细胞/Step3.Prostate_annotation.qs")


library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(future)
library(qs)
library(ggsci)

seurat.data = qread(file = "D:/CLEC2D/Step3.Prostate_annotation.qs")
seurat.data

CellDimPlot(seurat.data, group.by = "celltype", reduction = "UMAP")

ht <- GroupHeatmap(
  srt = seurat.data,
  features = c(
    "KRT8","KRT18","CLDN3","CLDN4", # Luminal epithelial
    "CD79A","MS4A1","CD19","IGHM", # B cells
    "CD68","TYROBP","C1QA","C1QB", # Myeloid cells
    "CD3D","CD3E","CD3G","CD2", # T cells
    "KRT5","KRT14","TP63","COL17A1", # Basal epithelial
    "VWF","PECAM1","CLDN5","EMCN", # Endothelial cells
    "TPSAB1","TPSB2","CPA3","KIT", # Mast cells
    "DCN","LUM","COL1A1","COL1A2", # Fibroblasts (CAFs)
    "MKI67","UBE2C","TOP2A","CDK1", # Cycling cells
    "MYH11","ACTA2","CNN1","TAGLN" # Perivascular smooth muscle
  ),
  group.by = "celltype",
  heatmap_palette = "YlOrRd",
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE
)
print(ht$plot)

CellStatPlot(
  seurat.data,
  stat.by ="celltype",
  group.by ="group",
  plot_type ="trend"
)

FeatureStatPlot(
  seurat.data,
  stat.by = "CLEC2D",
  group.by ="celltype",
  comparisons = list(c("Alpha","Beta"), c("Alpha","Delta"))
)

new_seurat <- subset(
  seurat.data,
  subset = celltype %in% c("Luminal", "Basal","Cycling")
)

saveRDS(new_seurat, file = "Epi.rds")



library(copykat)

ncore <- 16
copykat.test <- copykat(rawmat=counts, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names="",output.seg="FLASE", plot.genes="TRUE", genome="hg20",n.cores=16)

pred.test <- data.frame(copykat.test$prediction)
pred.test <- pred.test[which(pred.test$copykat.pred %in% c("aneuploid","diploid")),]  ##keep defined cells
CNA.test <- data.frame(copykat.test$CNAmat)

head(pred.test)


my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

chr <- as.numeric(CNA.test$chrom) %% 2+1
rbPal1 <- colorRampPalette(c('black','grey'))
CHR <- rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR,CHR)

rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
com.preN <- pred.test$copykat.pred
pred <- rbPal5(2)[as.numeric(factor(com.preN))]

cells <- rbind(pred,pred)
col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))

heatmap.3(t(CNA.test[,4:ncol(CNA.test)]),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
          ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
          notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
          keysize=1, density.info="none", trace="none",
          cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
          symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))

legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=0.6, bty="n")

malignant.cells <- rownames(
  copykat.test$prediction[
    copykat.test$prediction$copykat.pred == "aneuploid",
  ]
)

head(malignant.cells)
head(colnames(seurat.data))
malignant.cells <- gsub("\\.1$", "-1", malignant.cells)

seurat.malignant <- subset(
  seurat.data,
  cells = malignant.cells
)

pred <- copykat.test$prediction
rownames(pred) <- gsub("\\.1$", "-1", rownames(pred))

seurat.data$copykat <- pred$copykat.pred[
  match(colnames(seurat.data), rownames(pred))
]
saveRDS(seurat.malignant, file = "seurat_malignant.rds")

library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(future)

seurat.data = read_rds(file = "seurat_malignant.rds")
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
ggsave(plot=p1.compare, filename="恶性细胞去批次前.pdf", width = 10 ,height = 4.5)


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
ggsave(plot=p2.compare, filename="恶性细胞去批次后.pdf", width = 10 ,height = 4.5)

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


library(Seurat)
library(ggplot2)
library(UCell)
library(patchwork)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(GeneNMF)
DefaultAssay(seurat.data) <- "RNA"
seu.list <- SplitObject(seurat.data, split.by = "sampleID")

geneNMF.programs <- multiNMF(seu.list, assay="RNA", k=4:9, min.exp = 0.05)

geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        metric = "cosine",
                                        specificity.weight = 5,
                                        weight.explained = 0.5,
                                        nMP=7)

plotMetaPrograms(geneNMF.metaprograms,
                 similarity.cutoff = c(0.1,1))

geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        metric = "cosine",
                                        weight.explained = 0.8,
                                        nMP=7,
                                        min.confidence = 0.7)

mp.genes <- geneNMF.metaprograms$metaprograms.genes
seurat.data <- AddModuleScore_UCell(seurat.data, features = mp.genes, ncores=4, name = "")


matrix <- seurat.data@meta.data[,names(mp.genes)]

#dimred <- scale(matrix)
dimred <- as.matrix(matrix)

colnames(dimred) <- paste0("MP_",seq(1, ncol(dimred)))
#New dim reduction
seurat.data@reductions[["MPsignatures"]] <- new("DimReduc",
                                                cell.embeddings = dimred,
                                                assay.used = "RNA",
                                                key = "MP_",
                                                global = FALSE)

set.seed(123)
seurat.data <- RunUMAP(seurat.data, reduction="MPsignatures", dims=1:length(seurat.data@reductions[["MPsignatures"]]),
                       metric = "euclidean", reduction.name = "umap_MP")

DimPlot(seurat.data, reduction = "umap_MP", group.by = "sampleID") + theme(aspect.ratio = 1)

library(viridis)
FeaturePlot(seurat.data, features = names(mp.genes), reduction = "umap_MP", ncol=4) &
  scale_color_viridis(option="B") &
  theme(aspect.ratio = 1, axis.text=element_blank(), axis.ticks=element_blank())


expr <- GetAssayData(seurat.data, layer = "data")

CLEC2D <- expr["CLEC2D", ]
AR <- expr["AR", ]

seurat.data$CLEC2D_AR_subtype <- ifelse(
  CLEC2D > 0.25 & AR <= 0,
  "CLEC2D+_AR-",
  "Other"
)

library(dplyr)

meta <- seurat.data@meta.data

sample_summary <- meta %>%
  group_by(sampleID) %>%
  summarise(
    total_cells = n(),
    subtype_cells = sum(CLEC2D_AR_subtype == "CLEC2D+_AR-"),
    proportion = subtype_cells / total_cells,
    detected = subtype_cells > 0
  )

detection_rate <- mean(sample_summary$detected)
detection_rate


DimPlot(
  seurat.data,
  group.by = "CLEC2D_AR_subtype",
  reduction = "umap_MP"
)

library(ggplot2)
sample_summary$disease_type <- substr(sample_summary$sampleID, 1, 4)
ggplot(sample_summary,
       aes(x = reorder(sampleID, proportion),
           y = proportion,
           fill = disease_type)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_classic() +
  ylab("Proportion of CLEC2D+AR- cells") +
  xlab("Sample")

ggplot(sample_summary,
       aes(x = disease_type,
           fill = detected)) +
  geom_bar(position = "fill") +
  ylab("Fraction of samples") +
  theme_classic()

seurat.data <- RunDEtest(
  seurat.data,
  group.by = "CLEC2D_AR_subtype",
  fc.threshold = 1,
  only.pos = FALSE
)

VolcanoPlot(
  seurat.data,
  group.by = "CLEC2D_AR_subtype"
)



saveRDS(seurat.data, file = "seurat_malignant_cells.rds")

seurat.data = read_rds(file = "seurat_malignant_cells.rds")
seurat.data
library(scop)
library(org.Hs.eg.db)

DEGs <- seurat.data@tools$DEtest_CLEC2D_AR_subtype$AllMarkers_wilcox
DEGs <- DEGs[with(DEGs, avg_log2FC > 1 & p_val_adj < 0.05), ]

ht <- FeatureHeatmap(
  seurat.data,
  group.by = "CLEC2D_AR_subtype",
  features = DEGs$gene,
  feature_split = DEGs$group1,
  species = "Homo_sapiens",
  db = c("GO_BP", "KEGG", "WikiPathway"),
  anno_terms = TRUE,
  feature_annotation = c("TF", "CSPA"),
  feature_annotation_palcolor = list(
    c("gold", "steelblue"), c("forestgreen")
  ),
  height = 5, width = 4
)

print(ht$plot)

seurat.data <- RunGSEA(
  seurat.data,
  group.by = "CLEC2D_AR_subtype",
  db = "GO_BP",
  species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05"
)
GSEAPlot(
  seurat.data,
  group.by = "CLEC2D_AR_subtype",
  group_use = "CLEC2D+_AR-",
  id_use = "GO:0016477"
)

GSEAPlot(
  seurat.data,
  group.by = "CLEC2D_AR_subtype",
  group_use = "CLEC2D+_AR-",
  plot_type = "bar",
  direction = "both",
  topTerm = 20
)


gsea <- readRDS("c661ebe589b0067e.rds")

class(gsea)
library(enrichplot)
library(clusterProfiler)
gseaplot2(gsea, geneSetID = "hsa04064")
dotplot(gsea, showCategory=15)
gseaplot2(gsea, geneSetID = "HALLMARK_IL2_STAT5_SIGNALING")


get_gsea_stats <- function(GSEA_obj, pathway_name) {
  stats <- GSEA_obj@result[GSEA_obj@result$Description == pathway_name, ]
  return(list(
    pval = signif(stats$pvalue, 3),
    nes = signif(stats$NES, 2)
  ))
}

# 绘制第一条特定通路
pathway1 <- "HALLMARK_IL2_STAT5_SIGNALING"
stats1 <- get_gsea_stats(gsea, pathway1)
p2_title <- paste0(pathway1, "\n", 
                   "p-value: ", stats1$pval, ", ", 
                   "NES: ", stats1$nes)

{p2 = gseaplot2(gsea, 
                geneSetID = pathway1, 
                title = p2_title,  # 标题包含统计值
                color = "pink")
}
p2
ggsave("p2.pdf", plot = p2, width = 9, height = 6)  # 适当增加宽度容纳标题
dev.off()



library(scTenifoldKnk)
library(Seurat)
library(SeuratObject)
library(scop)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(ggrepel) 
library(tidyverse)
setwd("/mnt/DATA/home/s1352009245/CLEC2D")
seurat.data = read_rds(file = "seurat_malignant_cells.rds")
seurat.data

seurat.data <- FindVariableFeatures(
  seurat.data,
  selection.method = "vst",
  nfeatures = 5000
)

variable_genes<- VariableFeatures(seurat.data)[1:5000]
seurat.data<- subset(seurat.data,
                     features= variable_genes
)
table(seurat.data$celltype)

set.seed(2026)
sampled_cells<- sample(1:ncol(seurat.data),1000, replace = FALSE)
seurat.data<- seurat.data[, sampled_cells]

expression_matrix<- LayerData(seurat.data,
                              assay="RNA",
                              layer= 'counts')

target_gene <-'CLEC2D'
gene_presence <- target_gene %in% rownames(expression_matrix)

set.seed(2026)
perturbation_results <- scTenifoldKnk(
  countMatrix = expression_matrix,
  gKO = target_gene,
  qc_mtThreshold =0.1,
  qc_minLSize =400,
  nc_lambda =0,
  nc_nNet =10,
  nc_nCells =500,
  nc_nComp =3,
  nc_scaleScores = TRUE,
  nc_symmetric = FALSE,
  nc_q =0.9,
  td_K =3,
  td_maxIter =1000,
  td_maxError =1e-05,
  td_nDecimal =3,
  ma_nDim =2
)

significant_threshold<-0.05
differential_genes<- perturbation_results$diffRegulation %>%
  mutate(log_fold_change = log2(FC)) %>%
  filter(gene != target_gene) %>%
  filter(p.value < significant_threshold)
numeric_columns<-2:7
differential_genes[, numeric_columns] <- sapply(differential_genes[, numeric_columns], as.numeric)
output_filename<- paste0('1_', target_gene,"_虚拟扰动结果.xlsx")
write.xlsx(differential_genes, output_filename)

top_n_genes <- 20
differential_genes_top <- differential_genes %>%
  mutate(abs_logfc = abs(log_fold_change)) %>%
  arrange(desc(abs_logfc)) %>%
  slice_head(n = top_n_genes)
ggplot(differential_genes_top,
       aes(x = reorder(gene, log_fold_change),
           y = log_fold_change,
           fill ='pink')) +
  geom_bar(stat='identity', alpha = 0.8) +
  geom_text(aes(label = sprintf("%.2f", log_fold_change)),
            hjust = ifelse(differential_genes_top$log_fold_change> 0, -0.2, 1.2),
            size = 3,
            color ="black") +
  coord_flip() +
  labs(title =paste("Top", top_n_genes,"Differentially Regulated Genes"),
       subtitle =paste("KO gene:", target_gene),
       x ="Gene",
       y ="log2(FC)") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10, color ="black"),
    axis.title = element_text(size = 12, face ="bold"),
    plot.title = element_text(size = 14, face ="bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color ="gray50"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position ="none"
  ) +
  scale_y_continuous(expand= expansion(mult = c(0.05, 0.15))) # 为标签留出空间
ggsave('2_barplot_DR_genes.pdf', width = 10, height = 8)

volcano_data <- perturbation_results$diffRegulation%>%
  mutate(
    significant = p.value < significant_threshold & gene != target_gene,
    gene_label = ifelse(gene %in% head(arrange(differential_genes, desc(log_fold_change))$gene,10), gene, NA),
    neg_log_pval =-log10(p.value),
    perturbation_magnitude = case_when(
      log2(FC) >4& significant ~"High",
      log2(FC) >3& significant ~"Medium",
      significant ~"Low",
      TRUE ~"Not significant"
    ) 
  )
ggplot(volcano_data,
       aes(x = log2(FC), y = neg_log_pval)) +
  geom_point(aes(color = perturbation_magnitude, size = neg_log_pval),
             alpha =0.7, shape =16) +
  scale_color_manual(values = c(
    "High"="#E74C3C",
    "Medium"="#F39C12",
    "Low"="#3498DB",
    "Not significant"="gray70")) +
  geom_hline(yintercept =-log10(significant_threshold),
             linetype ="dashed", color ="blue", alpha =0.7) +
  geom_vline(xintercept =1,
             linetype ="dashed", color ="green", alpha =0.5) +
  geom_text_repel(
    data= subset(volcano_data, !is.na(gene_label)),
    aes(label = gene_label),
    size =3,
    max.overlaps =20,
    box.padding =0.5,
    segment.color ="gray50",
    segment.alpha =0.5
  ) +
  labs( title ="Virtual Gene Perturbation Analysis",
        subtitle = paste("KO gene:", target_gene),
        x ="log2(Perturbation Magnitude)",
        y ="-log10(p-value)",
        caption = paste("Significance threshold: p <", significant_threshold),
        color ="Perturbation\nMagnitude" ) + theme_bw() + theme(
          legend.position ="right",
          axis.text = element_text(size =11),
          axis.title = element_text(size =12, face ="bold"),
          plot.title = element_text(size =14, face ="bold", hjust =0.5),
          plot.subtitle = element_text(size =11, hjust =0.5, color ="gray50"),
          panel.border = element_rect(color ="black", fill = NA, linewidth =1)) +
  guides(size ="none")# 保存图形
ggsave('3_volcano_plot_magnitude.pdf', width =6.4, height =6)

# ===============================
# 3. MA图：扰动程度 vs 平均表达量
# ===============================

## 1️⃣ 计算平均表达量
mean_expression <- rowMeans(expression_matrix)

## 2️⃣ 设置要标记的基因数量
top_label_n <- 5

## 3️⃣ 选择扰动最显著的基因
genes_to_label <- differential_genes %>%
  arrange(desc(log_fold_change)) %>%
  slice_head(n = top_label_n) %>%
  pull(gene)

## 4️⃣ 提取扰动值并匹配顺序
fc_values <- perturbation_results$diffRegulation$FC[
  match(
    rownames(expression_matrix),
    perturbation_results$diffRegulation$gene
  )
]

## 5️⃣ 构建绘图数据框
ma_data <- data.frame(
  gene = rownames(expression_matrix),
  mean_expression = log1p(mean_expression),
  perturbation_magnitude = log2(fc_values),
  is_perturbed = rownames(expression_matrix) %in% differential_genes$gene
) %>%
  mutate(
    gene_label = ifelse(gene %in% genes_to_label, gene, NA),
    
    perturbation_category = case_when(
      is_perturbed & perturbation_magnitude > 2 ~ "High perturbation",
      is_perturbed & perturbation_magnitude > 1 ~ "Medium perturbation",
      is_perturbed ~ "Low perturbation",
      TRUE ~ "Not perturbed"
    )
  )

## 6️⃣ 固定颜色顺序
ma_data$perturbation_category <- factor(
  ma_data$perturbation_category,
  levels = c(
    "High perturbation",
    "Medium perturbation",
    "Low perturbation",
    "Not perturbed"
  )
)

## 7️⃣ 绘图
p <- ggplot(ma_data, aes(x = mean_expression, y = perturbation_magnitude)) +
  
  geom_point(
    aes(
      color = ifelse(is_perturbed, "Significant", "Not significant"),
      size = mean_expression
    ),
    alpha = 0.6,
    shape = 16
  ) +
  
  geom_text_repel(
    data = dplyr::filter(ma_data, !is.na(gene_label)),
    aes(label = gene_label),
    size = 3,
    max.overlaps = 30,
    box.padding = 0.35,
    point.padding = 0.3,
    segment.color = "gray50",
    segment.alpha = 0.5,
    color = "#E74C3C",
    fontface = "bold"
  ) +
  
  geom_smooth(
    method = "loess",
    color = "blue",
    se = FALSE,
    alpha = 0.1,
    linewidth = 0.5
  ) +
  
  scale_color_manual(
    values = c(
      "Significant" = "#E74C3C",
      "Not significant" = "gray70"
    )
  ) +
  
  scale_size_continuous(
    name = "Mean Expression",
    range = c(1, 5)
  ) +
  
  labs(
    title = "Perturbation Magnitude vs Mean Expression",
    subtitle = paste("KO gene:", target_gene),
    x = "log(Mean Expression + 1)",
    y = "log2(Perturbation Magnitude)",
    color = "Significance"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray50"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1)
  )

## 8️⃣ 保存
ggsave(
  "4_perturbation_MA_plot.pdf",
  plot = p,
  width = 6.67,
  height = 5.68
)

# ===============================
# 4. 点图：Top扰动基因可视化
# ===============================

## 1️⃣ 选择前20个扰动最显著基因
top_perturbed <- differential_genes %>%
  arrange(desc(log_fold_change)) %>%
  slice_head(n = 20)

## 2️⃣ 检查是否存在必要列
required_cols <- c("gene", "log_fold_change", "p.value")
missing_cols <- setdiff(required_cols, colnames(top_perturbed))

if (length(missing_cols) > 0) {
  stop(paste("缺少列:", paste(missing_cols, collapse = ", ")))
}

## 3️⃣ 绘图
p_dot <- ggplot(
  top_perturbed,
  aes(
    x = reorder(gene, log_fold_change),
    y = log_fold_change,
    size = -log10(p.value),
    color = log_fold_change
  )
) +
  geom_point(alpha = 0.75) +
  
  scale_color_gradient(
    low = "#3498DB",
    high = "#E74C3C",
    name = "Perturbation\nMagnitude"
  ) +
  
  scale_size_continuous(
    range = c(3, 8),
    name = "-log10(p-value)"
  ) +
  
  coord_flip() +
  
  labs(
    title = "Top Perturbed Genes",
    subtitle = paste("KO:", target_gene),
    x = "Gene",
    y = "Log2(Perturbation Magnitude)"
  ) +
  
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray50"),
    legend.position = "right"
  )

## 4️⃣ 保存
ggsave(
  "5_perturbation_dot_plot.pdf",
  plot = p_dot,
  width = 5.6,
  height = 6.5
)


library(org.Hs.eg.db)
library(clusterProfiler)
library(cowplot)
library(dplyr)
library(ggplot2)
library(enrichplot)

# 设置包含实验数据的文件路径
file_path <- "GSEAinput.txt"

# 读取实验数据表格，以制表符分隔，包含列名，不检查列名是否已存在
data <- read.table(file_path, sep = "\t", header = TRUE, check.names = FALSE)

# 按照 "logFC" 列的值从大到小进行排序
data <- data[order(data$logFC, decreasing = TRUE), ]
##将内容变为向量以便后续分析
geneList <- setNames(data$logFC, data$Symbol)
# 查看排序后的数据框
print(data)
##读取下载好的KEGG基因集
KEGG_GSEA <- read.gmt("c2.cp.kegg_medicus.v2026.1.Hs.symbols.gmt")
##进行分析
{GSEA <- GSEA(geneList, 
              TERM2GENE=KEGG_GSEA,
              pvalueCutoff = 0.5,
              pAdjustMethod = 'BH',
              minGSSize = 20,
              maxGSSize = 500,
              eps = 0
)
}
###导出结果框
GSEA_result<-GSEA@result

GSEA_result$Description <- gsub("^KEGG_MEDICUS_", "", GSEA_result$Description)
# 将更新后的 GSEA_result 保存回 GSEA 对象的 result 文件中
GSEA@result <- GSEA_result

names(GSEA@geneSets) <- gsub("^KEGG_MEDICUS_", "", names(GSEA@geneSets))


if ("ID" %in% colnames(GSEA@result)) {
  GSEA@result$ID <- gsub("^KEGG_MEDICUS_", "", GSEA@result$ID)
  
  # 将ID列设置为行名
  rownames(GSEA@result) <- GSEA@result$ID
  
}


####可视化绘图

get_gsea_stats <- function(GSEA_obj, pathway_name) {
  stats <- GSEA_obj@result[GSEA_obj@result$Description == pathway_name, ]
  return(list(
    pval = signif(stats$pvalue, 3),
    nes = signif(stats$NES, 2)
  ))
}

# 绘制第一条特定通路
pathway1 <- "REFERENCE_PRE_IC_FORMATION"
stats1 <- get_gsea_stats(GSEA, pathway1)
p2_title <- paste0(pathway1, "\n", 
                   "p-value: ", stats1$pval, ", ", 
                   "NES: ", stats1$nes)

{p2 = gseaplot2(GSEA, 
                geneSetID = pathway1, 
                title = p2_title,  # 标题包含统计值
                color = "pink")
}
p2
ggsave("p2.pdf", plot = p2, width = 9, height = 6)  # 适当增加宽度容纳标题
dev.off()

# 绘制第二条特定通路
pathway2 <- "VARIANT_MUTATION_INACTIVATED_PINK1_TO_ELECTRON_TRANSFER_IN_COMPLEX_I"
stats2 <- get_gsea_stats(GSEA, pathway2)
p3_title <- paste0(pathway2, "\n", 
                   "p-value: ", stats2$pval, ", ", 
                   "NES: ", stats2$nes)

{p3 = gseaplot2(GSEA, 
                geneSetID = pathway2, 
                title = p3_title, 
                color = "pink")
}
p3
ggsave("p3.pdf", plot = p3, width = 10, height = 6)  # 通路名称较长，进一步增加宽度
dev.off()
###绘制前5条通路，将下方的geneSetID = 1:5的数字3修改即可改成自己想要的
{p4 = gseaplot2(GSEA, geneSetID = 1:5)
}
p4
ggsave("p4.pdf", plot = p4, width = 8, height = 6)
dev.off()

##绘制气泡图
{p5 = dotplot(GSEA, showCategory=5) + ggtitle("dotplot for GSEA")
}
p5
ggsave("p5.pdf", plot = p5, width = 8, height = 15)
dev.off()


###分开绘制激活和抑制通路
{p6 = dotplot(GSEA,split=".sign")+facet_grid(~.sign)
}
p6
ggsave("p7.pdf", plot = p7, width = 10, height = 10)
dev.off()

############分别绘制top通路############
############分别绘制top通路############
############分别绘制top通路############
############分别绘制top通路############

gsea_res <- GSEA@result
gsea_res$.sign <- ifelse(gsea_res$NES > 0, "up", "down")  # NES>0为激活(up)，<0为抑制(down)

# 按方向分组，每组取NES绝对值最大的top5通路
top5_res <- gsea_res %>%
  dplyr::group_by(.sign) %>%
  dplyr::arrange(dplyr::desc(abs(NES))) %>%  # 按NES绝对值排序
  dplyr::slice_head(n = 3) %>%  # 取每组前5
  dplyr::ungroup()

# 从原始GSEA对象中筛选出这些top5通路（保持GSEA对象结构）
top5_GSEA <- GSEA
top5_GSEA@result <- top5_res
top5_GSEA@geneSets <- top5_GSEA@geneSets[top5_res$Description]  # 同步筛选基因集

# 绘制分组dotplot（此时每组仅含top5）
p8 <- dotplot(top5_GSEA, split = ".sign") + 
  facet_grid(~.sign) +
  ggtitle("Top5 Activated and Inhibited Pathways")

p8
ggsave("p8.pdf", plot = p7, width = 10, height = 10)
dev.off()


library(CellChat)
library(patchwork)
library(Seurat)
library(SeuratData)
library(qs)

seurat.data = qread(file = "/mnt/DATA/home/s1352009245/CLEC2D/Step3.Prostate_annotation.qs")
seurat.data

expr <- GetAssayData(seurat.data, assay = "RNA", layer = "data")
CLEC2D_exp <- expr["CLEC2D", ]
AR_exp     <- expr["AR", ]
CLEC2D_pos <- CLEC2D_exp > 0.5
AR_pos     <- AR_exp > 0
seurat.data$seurat_annotations <- seurat.data$celltype
target_cells <- seurat.data$celltype %in% c("Luminal", "Basal", "Cycling")

# 先全部设为 CLEC2D-AR+
seurat.data$seurat_annotations[target_cells] <- "CLEC2D-AR+"

# 再覆盖 CLEC2D+AR-
seurat.data$seurat_annotations[
  target_cells & CLEC2D_pos & !AR_pos
] <- "CLEC2D+AR-"


table(seurat.data$seurat_annotations)
data.input = seurat.data@assays$RNA@data
meta.data =  seurat.data@meta.data
meta.data = meta.data[!is.na(meta.data$seurat_annotations),]
data.input = data.input[,row.names(meta.data)]

#设置因子水平
meta.data$seurat_annotations = factor(meta.data$seurat_annotations,
                                      levels = c("B", "CLEC2D-AR+", "CLEC2D+AR-", "Endothelial", "Fibroblast", 
                                                 "Macrophage", "Mast", "Smooth muscle", "T"))

### 1.3 Create a CellChat object
cellchat <- createCellChat(object = data.input, 
                           meta = meta.data, 
                           group.by = "seurat_annotations")

cellchat <- addMeta(cellchat, meta = meta.data)

# 设置默认的labels
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use


cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 2) # do parallel
options(future.globals.maxSize = 4 * 1024^3)  # 4GB
cellchat <- identifyOverExpressedGenes(cellchat) #识别高表达基因
cellchat <- identifyOverExpressedInteractions(cellchat) #识别高表达通路

cellchat <- computeCommunProb(cellchat,population.size = F)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)


df.net <- subsetCommunication(cellchat)
head(df.net)
df.pathway = subsetCommunication(cellchat,slot.name = "netP")
head(df.pathway)

levels(cellchat@idents)

df.net.sub <- subsetCommunication(cellchat, sources.use = c(2), targets.use = c(3))
head(df.net.sub)
df.net.sub[order(-df.net.sub$prob), ]
cellchat <- computeCommunProbPathway(cellchat)
head(cellchat@net)

cellchat <- aggregateNet(cellchat)

qsave(cellchat, file = "Step1.CellCha_Res.qs")

library(CellChat)
library(patchwork)
library(Seurat)
library(SeuratData)
library(qs)
library(aplot)
library(ggplotify)

options(repr.plot.width = 4, repr.plot.height = 5)
netVisual_bubble(cellchat, sources.use = 3, targets.use = 2, remove.isolate = FALSE)
options(repr.plot.width = 6, repr.plot.height = 6)

netVisual_chord_gene(cellchat,
                     sources.use = 3,
                     targets.use = 2,
                     lab.cex = 0.5,legend.pos.y = 30)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
gg1 <- netAnalysis_signalingRole_scatter(cellchat)+ggtitle("All pathway")
gg1

netVisual_aggregate(cellchat,
                    signaling = "MIF",
                    layout = "circle")

netAnalysis_contribution(cellchat, signaling = "MK")
