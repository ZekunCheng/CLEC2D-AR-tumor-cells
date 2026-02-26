library(CellChat)
library(patchwork)
library(Seurat)
library(SeuratData)
library(qs)
pbmc3k = qread(file = "D:/骨肉瘤/骨肉瘤单细胞/tumor.qs")
Seurat.data = qread(file = "D:/骨肉瘤/骨肉瘤单细胞/Step3.PBMC_annotation.qs")

current_celltypes <- pbmc3k@meta.data$celltype
high_celltypes <- c("C9", "C2", "C4", "C10", "C8")
pbmc3k@meta.data$celltype <- ifelse(current_celltypes %in% high_celltypes, 
                                    "SH3-High", 
                                    "SH3-Low")
target_cells1 <- subset(Seurat.data, 
                       subset = celltype %in% c("Myeloid cells", "NKT", "CAFs", 
                                                "Plasma cells", "Endothelial cells", "B cells","Osteoclasts"))
target_cells2<- subset(pbmc3k, 
                       subset = celltype %in% c("SH3-Low"))
# 检查提取的细胞数量
table(target_cells$celltype)

pbmc3k <- merge(target_cells1, target_cells2)
# 检查当前 celltype 列的值分布
table(pbmc3k@meta.data$celltype)  # 假设分组列名为 'celltype'，如果不是请替换

# 将 "SH3-High" 替换为 "tumor"
pbmc3k@meta.data$celltype[pbmc3k@meta.data$celltype == "SH3-Low"] <- "tumor"

# 验证修改结果
table(pbmc3k@meta.data$celltype)  # 应显示 "tumor" 替代了 "SH3-High"

data.input = pbmc3k@assays$RNA@data
meta.data =  pbmc3k@meta.data
meta.data = meta.data[!is.na(meta.data$celltype),]
data.input = data.input[,row.names(meta.data)]

meta.data$celltype = factor(meta.data$celltype,
                            levels = c("Myeloid cells", "NKT", "CAFs", "Plasma cells", "Endothelial cells", 
                                       "B cells", "Osteoclasts","tumor"))

cellchat <- createCellChat(object = data.input, 
                           meta = meta.data, 
                           group.by = "celltype")

cellchat <- addMeta(cellchat, meta = meta.data)

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)
# 在调用 CellChat 函数前设置（单位：字节）
options(future.globals.maxSize = 2 * 1024^3)  # 增加到 2GB
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 2) # do parallel

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

cellchat <- computeCommunProbPathway(cellchat)
head(cellchat@net)

head(cellchat@netP)

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
p1 = netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                      weight.scale = T, label.edge= F,
                      title.name = "Number of interactions")
p2 = netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                      weight.scale = T, label.edge= F,
                      title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  p1 = netVisual_circle(mat2, vertex.weight = groupSize,
                        weight.scale = T, edge.weight.max = max(mat),
                        title.name = rownames(mat)[i])
}
qsave(cellchat, file = "Cellchat_SH3_Low.qs")