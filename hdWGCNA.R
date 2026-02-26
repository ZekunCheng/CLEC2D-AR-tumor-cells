###设置工作目录
setwd("D:/骨肉瘤/骨肉瘤单细胞")

###读取我们的seurat数据
#load("第一回-设置完软阈值后的膀胱癌hdwgcna所用文件.RData")
#load("设置完软阈值后的膀胱癌hdwgcna所用文件.RData")
seurat.tumor = qread(file = "D:/骨肉瘤/骨肉瘤单细胞/tumor.qs")

##加载包
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(UCell)
library(ggradar)
library(qs)
theme_set(theme_cowplot())
##设置随机种子
set.seed(12345)
##设置8个线程,伙伴们的电脑配置好的话可以尝试调大一点数字，可以加快速度
enableWGCNAThreads(nThreads = 8)
##将我们的单细胞文件mmRNA_harmony文件转换成官方数据所用到的文件名

seurat_obj <- seurat.tumor
target_celltypes <- c("C0", "C2", "C3", "C6", "C8", "C9", "C10")

# 方法1：直接创建子集对象（推荐）
seurat_RS <- subset(seurat_obj, subset = celltype %in% target_celltypes)
seurat_RS$celltype <- "RS"

seurat_obj <- seurat_RS

#以下数据为官方数据，大家可自行下载，自行尝试
#seurat_obj <- readRDS('Zhou_2020_control.rds')
##绘制我们的umap图

p1 <- DimPlot(seurat_obj, group.by='celltype', label=TRUE) +
  umap_theme() + ggtitle('示例数据') + NoLegend()

p1
ggsave("p1.pdf", plot = p1, width = 6, height = 6)
dev.off()
##清空下缓存
gc()

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", 
  #基因选择方法：
  #variable：使用存储在Seurat对象的VariableFeatures
  #fraction：使用在整个数据集或每组细胞中表达的一定比例的细胞中表达的基因
  #custom：使用自定义列表中指定的基因
  fraction = 0.05, # 基因在细胞中表达的比例，超过该比例的基因将被包含
  wgcna_name = "Ribotoxic Stress" # hdWGCNA 实验的名称
)

#构建每个组中的元细胞。单细胞数据中固有的稀疏性和噪声会导致虚假的基因-基因相关性
#从而使共表达网络分析复杂化。此外，单细胞或空间转录组数据的相关结构对于不同的子集
#细胞类型、细胞状态、解剖区域差异很大。scRNA-seq数据中hdWGCNA的工作流程是通过将
#高度相似的细胞折叠成“元细胞”以减少稀疏性，同时保持细胞异质性。
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("celltype", "orig.ident"), # 指定 seurat_obj@meta.data 中分组的列
  reduction = 'harmony', # 选择用于降维的方法
  k = 25, # 最近邻参数
  max_shared = 10, # 两个元细胞之间共享的最大细胞数
  ident.group = 'celltype' # 设置元细胞seurat对象的Idents
)



#接下来对元细胞以及单细胞数据进行规范化、标准化、批次效应校正和降维，这对于单细胞
##和元细胞RNA测序数据分析都是必要的。比如：单细胞RNA测序数据的原始读取数（read counts）
#会受多个因素影响（如捕获效率、测序深度等），这些因素会导致不同细胞之间数据的可比性
#降低。因此，通过规范化，可以减少这些技术偏差，使得后续分析更为可靠。它们的主要目的是：
#减小技术变异：通过规范化和标准化，减少测序深度、捕获效率等技术因素对数据的影响。
#提高数据一致性：确保不同细胞或元细胞数据的可比性，为后续分析提供可靠的基础。
#消除批次效应：整合不同批次的样本，使得数据更具生物学意义。
#降维与可视化：PCA和UMAP等降维技术有助于数据可视化和主要变化趋势的提取，使得分析和解释更加直观。
#通过这些处理步骤，可以更好地揭示数据中的生物学信号，提高分析结果的可靠性和可解释性。
# 规范化元细胞表达矩阵
seurat_obj <- NormalizeMetacells(seurat_obj)
# 获取元细胞对象
metacell_obj <- GetMetacellObject(seurat_obj)
# 再次规范化元细胞
seurat_obj <- NormalizeMetacells(seurat_obj)
# 标准化元细胞
seurat_obj <- ScaleMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
# 运行PCA分析
seurat_obj <- RunPCAMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
# 输出变量特征的数量
length(VariableFeatures(seurat_obj))
# 运行Harmony分析
seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars='orig.ident')
# 运行UMAP降维
seurat_obj <- RunUMAPMetacells(seurat_obj, reduction='harmony', dims=1:15)
# 绘制元细胞的UMAP图，按细胞类型分组
p2 <- DimPlotMetacells(seurat_obj, group.by='celltype') + umap_theme() + ggtitle("Cell Type")
p2
# 绘制元细胞的UMAP图，按原始样本分组
p3 <- DimPlotMetacells(seurat_obj, group.by='orig.ident') + umap_theme() + ggtitle("orig.ident")
p3
# 并排显示两个图
p4 = p2 | p3
p4
ggsave("p4.pdf", plot = p4, width = 12, height = 6)
dev.off()

#设置表达矩阵，同时挑选分析的细胞群，
#使用hdWGNCA对挑选的细胞群Epithelial cells进行共表达网络分析
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c("RS"),  # 所有需要的细胞类型
  group.by = 'celltype',  # 指定分组列
  assay = 'RNA',          # 使用RNA assay
  slot = 'data'           # 使用标准化数据
)
###清空下缓存
gc()
#或者也可以将目标设置为多个细胞类群,如"Epithelial cells", "Monocyte cell"。
{#seurat_obj <- SetDatExpr(
  #seurat_obj,
  #group_name = c("Epithelial cells", "Monocyte cell"),
  #group.by='Celltype'
  )
}

#选择软功率阈值，这是极其重要的一步。hdWGCNA构建基因-基因相关邻接矩阵来推断基因
#之间的共表达关系。目的是为确定多组基因中的相关性，以减少相关矩阵中存在的噪声量，
#从而保留强连接并删除弱连接。因此，确软功率阈值的适当值至关重要。
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' #你也可以使用"unsigned"或"signed hybrid"
  #有符号网络Signed：
  #有符号网络不仅考虑基因表达相关性的强度，还考虑其方向。正相关和负相关会被区分开来。
  #无符号网络Unsigned：
  #无符号网络不考虑基因表达相关性的方向。也就是说，正相关和负相关都被视为相同的强度，只是方向不同。
  #签名混合网络Hybrid：
  #结合了有符号和无符号网络的特性。它对正相关的基因对使用标准相关系数，而对负相关的基因对使用零权重。
)
gc()

# 可视化软阈值结果
plot_list <- PlotSoftPowers(seurat_obj)
# 使用patchwork进行拼装
p5 <- wrap_plots(plot_list, ncol=2)
p5
# 保存为PDF
ggsave("p5.pdf", plot = p5, width = 9, height = 6) 
dev.off()

# 获取功率表
power_table <- GetPowerTable(seurat_obj)
head(power_table)
#更改软阈值，构建共表达网络
  
###清空缓存
gc()
##保存工作环境所有文件
save.image(file = "第一回-设置完软阈值后的膀胱癌hdwgcna所用文件.RData")
gc()

# 绘图
pdf(file = "p6.pdf", width = 9, height = 5)
#绘制树状图，“灰色”模块由未归入任何共表达模块的基因组成。在所有下游分析和解释中，应忽略灰色模块。
PlotDendrogram(seurat_obj, main='RS cells hdWGCNA Dendrogram')
# 关闭 PDF 设备
dev.off()
# 获取拓扑重叠矩阵
TOM <- GetTOM(seurat_obj)
# 需要先运行 ScaleData，否则 Harmony 会报错
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

###########计算基因模块相关性############
#通过调用 ModuleConnectivity 函数来计算基于特征基因module eigengene的
#连通性（kME）。具体来说，该函数根据基因表达数据和模块划分，计算每个
#基因与特征基因的皮尔逊相关系数。这些相关系数也被称为模块成员值
#（Module Membership values, kME），用来衡量每个基因在其所属模块中的
#重要性或连接度。计算基于特征基因的连通性 (kME)。
#在这里，选择了 "Epithelial cells"（上皮细胞）。这意味着仅计算上皮细胞
#中的基因表达数据，进而计算这些细胞类型中的基因模块特征和连通性。

# 在完整的单细胞数据集中计算所有模块特征基因
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="orig.ident"
)
gc()

seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'celltype', group_name = 'RS'
)

# 根据Epithelial cells-M重命名模块
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "RS-M"
)

# 绘制每个模块中按kME排序的基因图
### 下面的ncol可调整
p7 <- PlotKMEs(seurat_obj, ncol=3)#ncol=3意味着每行放3个图，可自行调整
p7
ggsave("p7.pdf", plot = p7, width = 10, height = 10)
dev.off()
#从 seurat_obj 对象中提取模块信息，并过滤掉灰色模块。灰色模块通常代表
#那些未分类或没有显著共表达模式的基因，因此在很多分析中会被排除。
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')
# 显示前6列
head(modules[,1:6])
# 获取中心基因
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
view(hub_df)
gc()
##保存一下计算模块及特征基因后的seurat_obj
saveRDS(seurat_obj,"计算模块后的seurat_obj.rds")
##加载一下
#seurat_obj=readRDS("计算模块后的seurat_obj.rds")

# 使用UCell方法计算每个模块中前25个中心基因的基因评分
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='Seurat'
)
gc()

######绘制特征基因和中心基因的图
#特征基因（Module Eigengene, ME）和中心基因（Hub Gene）是基因共表达网络分析中的两个
#不同概念，尽管它们在某些情况下可能具有相似的作用，但它们的定义和计算方式不同，因此不
#能完全等同。

# 绘制每个模块的特征基因图
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # 绘制hMEs
  order=TRUE # 按hMEs从高到低排序
)

# 使用patchwork将多个图拼接在一起，ncol=3可修改
p8 = wrap_plots(plot_list, ncol=5)
p8##这里绘图可能需要一点时间，出图后再继续往下
ggsave("p8特征基因.pdf", plot = p8, width = 9, height = 5)
dev.off()
# 绘制每个模块的中心基因得分图
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # 绘制中心基因得分
  order='shuffle', # 随机排列细胞
  ucell = TRUE # 根据Seurat或UCell方法进行基因评分
)

# 使用patchwork将多个图拼接在一起，ncol=3可修改
p9 = wrap_plots(plot_list, ncol=5)
p9##这里绘图可能需要一点时间，出图后再继续往下
ggsave("p9中心基因.pdf", plot = p9, width = 9, height = 5)
dev.off()

# 绘制模块相关图
pdf(file = "p10.pdf", width = 8, height = 5)
ModuleCorrelogram(seurat_obj)
dev.off()

# Seurat官方推荐方式
seurat_obj$annotation <- as.character(Idents(seurat_obj))

###雷达图呈现,但你的细胞类群需要进行亚分群，其中亚分群的列名为annotation，
###如不是，则需要修改下方的“annotation”
seurat_obj$cluster <- do.call(rbind, strsplit(as.character(seurat_obj$annotation), ' '))[,1]

p20 = ModuleRadarPlot(
  seurat_obj,
  group.by = 'cluster',
  barcodes = seurat_obj@meta.data %>% subset(celltype == 'RS') %>% rownames(),
  axis.label.size=4,
  grid.label.size=4
)
ggsave("雷达图.pdf", plot = p20, width = 10, height = 10)


# 从Seurat对象中获取hMEs，同样排除grey模块
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# 将hMEs添加到Seurat元数据中
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# 绘制所有细胞类群的气泡图
p11 <- DotPlot(seurat_obj, features=mods, group.by = 'celltype')
p11
ggsave("p11.pdf", plot = p11, width = 13, height = 5)
dev.off()

# 翻转x/y轴，旋转轴标签，并改变颜色方案
p12 <- p11 +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# 绘制输出
p12
ggsave("p12.pdf", plot = p12, width = 8, height = 5)
dev.off()


###小提琴图绘制某个模块与其他细胞类群的相关性
p13 <- VlnPlot(
  seurat_obj,
  features = 'RS',
  group.by = 'celltype',
  pt.size = 0 # don't show actual data points
)
p13
ggsave("p13.pdf", plot = p13, width = 8, height = 5)
dev.off()

#添加箱线
p14 <- p + geom_boxplot(width=.25, fill='white')
p14 <- p + xlab('') + ylab('hME') + NoLegend()
p14
ggsave("p14.pdf", plot = p14, width = 8, height = 5)
dev.off()
gc()

#############绘制枢纽基因网络图################
ModuleNetworkPlot(
  seurat_obj,
  outdir = 'ModuleNetworks'#创建一个名字为ModuleNetworks的文件夹进行存放图
)

ModuleNetworkPlot(
  seurat_obj, 
  outdir='ModuleNetworks2', # 新文件夹名称
  n_inner = 20, # 内环中的基因数量
  n_outer = 30, # 外环中的基因数量
  n_conns = Inf, # 显示所有连接
  plot_size=c(10,10), # 较大的绘图区域
  vertex.label.cex=1 # 字体大小
)


####导出基因，并进行功能富集分析#####
####导出基因，并进行功能富集分析#####
####导出基因，并进行功能富集分析#####
# 获取中心基因
gokegg_hub_df <- GetHubGenes(seurat_obj, n_hubs = 50)
view(gokegg_hub_df)

filtered_gokegg_hub_df <- gokegg_hub_df[grepl("M(1[3]?|[24569])$", gokegg_hub_df$module), ]
write.table(gokegg_hub_df, 'hdWGCNA_hubgene.txt', sep="\t", quote=FALSE, row.names = FALSE) 
