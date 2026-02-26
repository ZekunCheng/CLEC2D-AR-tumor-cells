library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(future)
library(qs)

seurat.data = read_rds(file = "D:/骨肉瘤/骨肉瘤单细胞/Step1.RawCount_merged_seurat.rds")

#有提取亚组的时候改一下
os = subset(seurat.data, group %in% c("OS1", "OS2", 
                                        "OS3", "OS4","OS5","OS6"))
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
os.qc <- subset(os, subset = nFeature_RNA > 500 & nFeature_RNA < 6500 & percent.mt < 10 & percent.hb < 1)
os.qc

qsave(os.qc, file = "D:/骨肉瘤/骨肉瘤单细胞/Step2.PBMC_afterQC.qs")