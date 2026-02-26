rm(list = ls())
library(infercnv)
library(Seurat)
library(dplyr)
library(readr)
library(future)
library(qs)
outdir = "/public/home/hpc8202221224/InferCNV"
dir.create(outdir)


seurat.data = qread(file = "/public/home/hpc8202221224/Step3.PBMC_annotation.qs")
seurat.data

seurat.endo <- subset(seurat.data, 
                   subset = celltype %in% c("Endothelial cells"))
seurat.os = qread(file = "/public/home/hpc8202221224/成骨细胞亚群.qs")
seurat.data <- merge(seurat.os, seurat.endo)
# 检查合并后的细胞类型
table(seurat.data$celltype)

inferCNV.anno = data.frame(cell.id = rownames(seurat.data@meta.data),
                           group = seurat.data$celltype)
table(inferCNV.anno$group)

count.data <- GetAssayData(seurat.data, slot='counts',assay='RNA') %>% as.data.frame()
count.data = count.data[,inferCNV.anno$cell.id]
dim(count.data)

## 2.3 基因组文件
geneInfor = read.table("/public/home/hpc8202221224/gene_pos.txt")
comm.gene = intersect(geneInfor$V1, rownames(count.data) )
head(comm.gene)

# Write table
write.table(inferCNV.anno,file = paste0(outdir,'groupFiles.txt'),
            sep = '\t',quote = F,col.names = F,row.names = F)

write.table(count.data,file = paste0(outdir,'expFile.txt'),
            sep = '\t',quote = F)

write.table(geneInfor,
            file = paste0(outdir,"geneFile.txt"),
            row.names = F, col.names = F, quote=F, sep="\t")

## 2.4 构建inferCNV对象
expFile = paste0(outdir,'expFile.txt')
groupFiles=paste0(outdir,'groupFiles.txt')
geneFile=paste0(outdir,'geneFile.txt')

group.data = read.table(groupFiles, sep="\t")
dim(group.data)
table(group.data$V2)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix = expFile,
                                    annotations_file = groupFiles,
                                    gene_order_file = geneFile,
                                    ref_group_names = c("Endothelial cells"),
                                    delim="\t")  ## 这个取决于自己的分组信息里面的

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir='Infercnv-results/', 
                             cluster_by_groups = T, # 默认False; 先区分细胞来源，再做层次聚类
                             analysis_mode = 'subclusters',
                             tumor_subcluster_partition_method = 'random_trees',
                             denoise=TRUE,
                             HMM=TRUE,
                             HMM_type = 'i6',
                             write_expr_matrix = F, # 默认False  
                             num_threads=parallel::detectCores())