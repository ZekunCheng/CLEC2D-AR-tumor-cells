library(rjags)
library(infercnv)
library(Seurat)
library(harmony)
library(tidyverse)
library(dplyr)
library(patchwork)
library(tidydr)
library(ggplot2)
library(cowplot)
library(qs)
getwd()


seurat.data = qread(file = "D:/骨肉瘤/骨肉瘤单细胞/Step3.PBMC_annotation.qs")
seurat.os = qread(file = "D:/骨肉瘤/骨肉瘤单细胞/成骨细胞亚群.qs")
seurat.endo <- subset(seurat.data, 
                      subset = celltype %in% c("Endothelial cells"))
seurat.data <- merge(seurat.os, seurat.endo)
table(seurat.data$celltype)

matrix_counts <- as.matrix(GetAssayData(seurat.data, assay="RNA", layer='counts'))

#保存分组
write.table(seurat.data$celltype, "celltype.label.txt", sep = "\t", quote = F, col.names = F)
# Anno=read.table("celltype.label.txt",row.names = 1)

#gencode和单细胞基因取交集
library(readr)

#加载基因组注释文件
gencode <- read_tsv("hg38_gencode_v27.txt", col_names = c("gene", "chr", "start", "end"))

#gencode和单细胞基因取交集
gencode <- gencode[!duplicated(gencode$gene),]
common_genes <- intersect(gencode$gene, rownames(matrix_counts))#基因取交集


sort(unique(seurat.data$celltype))
#创建inferfercnvObject
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = matrix_counts[common_genes,],#counts矩阵
                                    annotations_file="./celltype.label.txt",#celltype信息
                                    delim="\t",
                                    gene_order_file="./hg38_gencode_v27.txt",#基因组注释文件
                                    ref_group_names=c("Endothelial cells" ),#组正常组织的上皮的细胞作为reference
                                    chr_exclude=c("chrY", "chrM"))#选择不需要的染色体，查看帮助函数去除



# infercnv分析，这一步很慢，比较消耗cpu和内存
#具体参数可参照帮助函数，一般使用默认或者看看参考文献
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="infercnv/",  #分析结果out文件夹名，会自动创建到当前路径
                             no_prelim_plot = T,
                             cluster_by_groups=T,##是否根据患者类型
                             denoise=TRUE,# 去噪音
                             HMM=F,#是否要去背景噪音，选择F节省时间且不影响结果
                             min_cells_per_gene = 10,
                             num_threads=10,#线程数
                             # num_threads = parallel::detectCores(),  # 使用的线程数，这里设置为使用所有可用的CPU核心
                             write_expr_matrix = T
                             #这里的默认参数是F
                             #这里要选择T，分析后将结果exp矩阵导出，出现infercnv.observations.txt结果文件
)
#保存infercnv_obj
save(infercnv_obj,file = "infercnv_obj.Rdata")
load("infercnv_obj.Rdata")

#分析完成后在infercnv文件夹中会有很多的文件，包括文献中一直见到的热图
#大都数文件是没有用的，不用保存，需要看或者之后用的文件有：
#infercnv.png
#infercnv.observations.txt
#infercnv.references.txt
#run.final.infercnv_obj



#infercnv.png是默认的出图，infercnv也提供了作图函数，我们可以对热图颜色进行修改
#然后将文件保存为pdf
library(RColorBrewer)
infercnv::plot_cnv(infercnv_obj, 
                   output_filename = "inferCNV_heatmap",
                   output_format = "pdf",
                   custom_color_pal =  color.palette(c("#2067AE","white","#B11F2B"))) #改颜色


####cnvscore方法一（其实就是对cnv进行定量，方法很多，但往往比较局限）####

#加载infercnv.observation_groupings.txt 文件：包含细胞群体的信息，即每个细胞被分配到的组。
# infercnv.observations.txt 文件：包含实验组的每个细胞的cnv值，用于CNV分析。
grp=read.table("./infercnv/infercnv.observation_groupings.txt",sep = "",header = T)
obs = data.table::fread("./infercnv/infercnv.observations.txt",
                        data.table = F) %>%
  column_to_rownames(var = 'V1')

#查看拷贝数最大值和最小值
max(obs)
# [1] 1.111734
min(obs)
# [1] 0.8882661


library(scales)  # 加载scales包，用于数据缩放
# 定义cnvScore函数，用于计算CNV评分
cnvScore <- function(data) {
  data <- data %>% 
    as.matrix() %>%           # 将数据转换为矩阵
    t() %>%                   # 转置矩阵
    scale() %>%               # 标准化数据，使其均值为0，标准差为1
    rescale(to = c(-1, 1)) %>%  # 重新缩放数据，使其范围在-1到1之间
    t()                       # 再次转置矩阵，使其恢复原始形状
  
  cnv_score <- as.data.frame(colSums(data * data))  # 计算每列的平方和，并转换为数据框
  return(cnv_score)  # 返回CNV评分数据框
}
# 使用cnvScore函数计算CNV评分
cnv_score <- cnvScore(obs)


# 提取Dendrogram.Group列中第二个下划线前的字符并创建新列cluster
grp$cluster <- substr(grp$Dendrogram.Group, 1, 3)  # 提取第1到第3个字符
sort(unique(grp$cluster))
grp$cluster <- factor(grp$cluster)
# 检查细胞名称顺序能否对得上
identical(row.names(cnv_score),row.names(grp))
# grp=grp[row.names(cnv_score)
cnv_score$cluster=grp$cluster
colnames(cnv_score)=c("score","cluster")
# cnvscore作图
library(ggpubr)
ggboxplot(cnv_score,"cluster","score",fill = "cluster")+ coord_cartesian(ylim = c(0, 1000))
table(cnv_score$cluster)
library(ggplot2)
identical(row.names(cnv_score), row.names(grp))