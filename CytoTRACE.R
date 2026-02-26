library(Seurat)
library(dplyr)
library(readr)
library(future)
library(qs)
library(CytoTRACE)
# 第一次library CytoTRACE 需要 Python 环境来运行一些依赖的功能。你可以选择安装 Miniconda
# 咱选择N就行，没必要安装python环境

# 加载上皮细胞数据
seurat.os = qread(file = "D:/骨肉瘤/骨肉瘤单细胞/成骨细胞亚群.qs")
table(seurat.os$celltype)

# 将 RNA count 数据转换为矩阵
exp1 <- as.matrix(GetAssayData(seurat.os, assay = "RNA", layer = "counts"))
# 过滤掉在少于 5 个细胞中表达的基因
exp1 <- exp1[apply(exp1 > 0, 1, sum) >= 5,]
# 使用 CytoTRACE 进行分析，设置 ncores = 1 表示使用一个 CPU 核心进行计算
results <- CytoTRACE(exp1, ncores = 1)
# 提取细胞类型注释信息，并将其转换为字符向量
phenot <- seurat.os$celltype
phenot <- as.character(phenot)
# 将细胞类型注释的名字设置为元数据中的行名
names(phenot) <- rownames(seurat.os@meta.data)
# 提取 UMAP 降维后的细胞嵌入信息
emb <- seurat.os@reductions[["umap"]]@cell.embeddings
# 使用 CytoTRACE 结果、细胞类型注释和 UMAP 嵌入信息绘制 CytoTRACE 图，并将结果保存到指定目录
plotCytoTRACE(results, phenotype = phenot, emb = emb, outputDir = 'D:/骨肉瘤/骨肉瘤单细胞/CytoTRACE结果/')
# 绘制 CytoTRACE 基因表达图，显示前 30 个基因，并将结果保存到指定目录
plotCytoGenes(results, numOfGenes = 30, outputDir = 'D:/骨肉瘤/骨肉瘤单细胞/CytoTRACE结果/')