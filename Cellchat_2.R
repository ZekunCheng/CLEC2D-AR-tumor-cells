library(CellChat)
library(patchwork)
library(readr)
library(aplot)
library(ggplotify)
getwd()
data.dir <- './Outdata/'
source("./Rawdata/custom_seurat_functions.R")

cellchat.NL <- qread(file = "D:/骨肉瘤/骨肉瘤单细胞/Cellchat_SH3_High.qs")
cellchat.LS <- qread(file = "D:/骨肉瘤/骨肉瘤单细胞/Cellchat_SH3_Low.qs")

object.list <- list("SH3-High" = cellchat.NL, "SH3-Low" = cellchat.LS)
cellchat <- mergeCellChat(
  object.list, 
  add.names = names(object.list),
  cell.prefix = TRUE  # 自动添加数据集名前缀
)

cellchat

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")

options(repr.plot.width = 6, repr.plot.height = 4)
gg1 + gg2

options(repr.plot.width = 12, repr.plot.height = 6)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T,comparison = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",comparison = c(1,2))

gg1 <- netVisual_heatmap(cellchat,comparison = c(1,2))
gg2 <- netVisual_heatmap(cellchat, measure = "weight",comparison = c(1,2))
gg1 + gg2

# 1. 加载必要的包
library(CellChat)
library(patchwork)

# 2. 确保对象列表已命名（例如 "SH3_High" 和 "SH3_Low"）
names(object.list) <- c("SH3_High", "SH3_Low")

# 3. 计算中心性分数
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")
}

#差异性circle plot
options(repr.plot.width = 12, repr.plot.height = 6)
# 总览性circle plot
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

options(repr.plot.width = 10, repr.plot.height = 4)
patchwork::wrap_plots(plots = gg)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional",comparison = c(1,2))
# 确保 reticulate 包已安装
install.packages("reticulate")
library(reticulate)
py_config()
py_install(packages = "umap-learn", method = "pip")

cellchat <- netEmbedding(cellchat, type = "functional",comparison = c(1,2))
cellchat <- netClustering(cellchat, type = "functional",comparison = c(1,2))
options(repr.plot.width = 6, repr.plot.height = 4.5)
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5,comparison = c(1,2))


library(ComplexHeatmap)
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

##outgoing
# 生成热图
ht1 = netAnalysis_signalingRole_heatmap(
  object.list[[1]], 
  pattern = "outgoing", 
  signaling = pathway.union,
  title = names(object.list)[1],
  width = 13, 
  height = 20
)

ht2 = netAnalysis_signalingRole_heatmap(
  object.list[[2]], 
  pattern = "outgoing",
  signaling = pathway.union,
  title = names(object.list)[2],
  width = 13,
  height = 20
)

# 保存为 PDF
pdf("combined_outgoing_heatmaps.pdf", width = 26, height = 20)
draw(
  ht1 + ht2,
  ht_gap = unit(1, "cm"),
  column_title = "Comparative Outgoing Signaling Roles",
  column_title_gp = gpar(fontsize = 18, fontface = "bold")
)
dev.off()

#incoming
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", 
                                        signaling = pathway.union, title = names(object.list)[i],
                                        width = 13, height = 20, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming",
                                        signaling = pathway.union, title = names(object.list)[i+1],
                                        width = 13, height = 20, color.heatmap = "GnBu")
#ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "incoming", 
#                                        signaling = pathway.union, title = names(object.list)[i+2],
#                                        width = 5, height = 6, color.heatmap = "GnBu")
#draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm"))
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# 设置PDF输出（精细调整版）
pdf("combined_incoming_heatmaps.pdf", 
    width = 26.5,  # 13*2 + 0.5(间距)
    height = 20) # 避免中文乱码

# 智能合并热图（自动处理图例和标题）
draw(
  ht1 + ht2,
  ht_gap    = unit(0.8, "cm"),          # 稍大间距防止重叠
  row_title = "Signaling Pathways",     # 行标题
  column_title = paste(names(object.list)[i], "vs", names(object.list)[i+1]), # 动态标题
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_side = "right",        # 图例统一右对齐
  padding = unit(c(2, 15, 2, 2), "mm")  # 增加右侧边距（防止图例截断）
)

# 安全关闭设备
invisible(dev.off())

options(repr.plot.width = 7, repr.plot.height = 4)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)

print(cellchat@meta$datasets)  # 应显示类似 "1" "2" 的标签
print(unique(cellchat@idents))  # 检查所有细胞分组

gg1 <- netVisual_bubble(
  cellchat,
  sources.use = 4, 
  targets.use = 5:11,
  comparison = c("NL", "LS"),  # 改用对象创建时的实际名称
  title.name = "Increased signaling in LS",
  angle.x = 45,
  remove.isolate = TRUE
)

# 方法2：确保 max.dataset 与 comparison 一致
gg1 <- netVisual_bubble(
  cellchat,
  sources.use = 4,
  targets.use = 5:11,
  comparison = c(1, 2),
  max.dataset = 2,  # 必须对应 comparison 的第二个数据集
  title.name = "Increased signaling in LS",
  angle.x = 45,
  remove.isolate = TRUE
)


gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)

options(repr.plot.width = 10, repr.plot.height = 4)
gg1 + gg2