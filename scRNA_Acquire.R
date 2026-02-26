setwd("D:/骨肉瘤/骨肉瘤单细胞")
library(Seurat)
library(dplyr)

#选择路径
fileID = list.files("D:/骨肉瘤/骨肉瘤单细胞/")
path = "D:/骨肉瘤/骨肉瘤单细胞/"
fileID

#读进去
seurat.list = lapply(fileID, function(file){
  seurat_data <- Read10X(data.dir = paste0(path, file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   min.cells = 0,
                                   min.features = 0,
                                   project = file)
  return(seurat_obj)
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
                              "OS1" = "OS1",
                              "OS2" = "OS2",
                              "OS3" = "OS3",
                              "OS4" = "OS4",
                              "OS5" = "OS5",
                              "OS6" = "OS6")

# 保存一下
saveRDS(merged_seurat,file = "D:/骨肉瘤/骨肉瘤单细胞/Step1.RawCount_merged_seurat.rds")