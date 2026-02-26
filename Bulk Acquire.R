# 设置工作目录，设置当前工作目录为TCGA下载与整理的文件夹
setwd("D:/骨肉瘤") 

# 加载必要的包
library(jsonlite)
library(tidyverse)
library(dplyr)

# 读取JSON文件中的元数据
meta_data <- fromJSON("metadata.cart.2025-07-05.json") # 读取JSON格式的元数据文件

# 构建样本信息数据框
# 使用map2_df函数结合元数据中的关联实体信息与文件名，创建一个样本信息的数据框
samples_df <- map2_df(meta_data$associated_entities, seq_along(meta_data$associated_entities),
                      ~ tibble(sample_id = .x[,1], file_name = meta_data$file_name[.y]))

# 获取计数数据文件列表
count_file_paths <- list.files('gdc_download_20250705_074012.501028/', pattern = '*.tsv', recursive = TRUE) # 获取当前目录下所有tsv文件的路径

# 提取纯文件名
file_names_only <- sapply(strsplit(count_file_paths, split='/'), function(x) x[2]) # 提取文件名，不包含路径信息

# 创建空的表达矩阵
expr_matrix <- data.frame() # 初始化空的数据框，用于存储表达数据

# 循环读取每个文件，并生成TPM类型文件
for (file_index in seq_along(count_file_paths)) {
  full_path <- paste0('gdc_download_20250705_074012.501028/', count_file_paths[file_index]) # 拼接完整的文件路径
  count_data <- read.delim(full_path, fill = TRUE, header = FALSE, row.names = 1) # 读取计数数据，使用行名作为列名
  colnames(count_data) <- count_data[2,] # 设置列名，跳过第一行
  count_data <- count_data[-(1:6),] # 删除前六行，这些行包含 header 信息
  
  # 提取感兴趣的TPM数据列
  tpm_data <- count_data[6]
  
  # 设置TPM数据的列名
  sample_id <- samples_df$sample_id[which(samples_df$file_name == file_names_only[file_index])][[1]]
  colnames(tpm_data) <- sample_id
  
  if (nrow(expr_matrix) == 0) {
    # 如果expr_matrix还是空的，直接赋值
    expr_matrix <- tpm_data
  } else {
    # 检查行名是否一致，确保数据的一致性
    if (!all(rownames(expr_matrix) == rownames(tpm_data))) {
      stop("行名不一致，无法合并数据")
    }
    expr_matrix <- cbind(expr_matrix, tpm_data)
  }
}


# 循环读取每个文件，并生成count类型文件
for (file_index in seq_along(count_file_paths)) {
  full_path <- paste0('gdc_download_20250705_074012.501028/', count_file_paths[file_index]) # 拼接完整的文件路径
  count_data <- read.delim(full_path, fill = TRUE, header = FALSE, row.names = 1) # 读取计数数据，使用行名作为列名
  colnames(count_data) <- count_data[2,] # 设置列名，跳过第一行
  count_data <- count_data[-(1:6),] # 删除前六行，这些行包含 header 信息
  
  # 提取感兴趣的TPM数据列
  tpm_data <- count_data[3]
  
  # 设置TPM数据的列名
  sample_id <- samples_df$sample_id[which(samples_df$file_name == file_names_only[file_index])][[1]]
  colnames(tpm_data) <- sample_id
  
  if (nrow(expr_matrix) == 0) {
    # 如果expr_matrix还是空的，直接赋值
    expr_matrix <- tpm_data
  } else {
    # 检查行名是否一致，确保数据的一致性
    if (!all(rownames(expr_matrix) == rownames(tpm_data))) {
      stop("行名不一致，无法合并数据")
    }
    expr_matrix <- cbind(expr_matrix, tpm_data)
  }
}
# 循环读取每个文件，并生成fpkm类型文件
for (file_index in seq_along(count_file_paths)) {
  full_path <- paste0('gdc_download_20250705_074012.501028/', count_file_paths[file_index]) # 拼接完整的文件路径
  count_data <- read.delim(full_path, fill = TRUE, header = FALSE, row.names = 1) # 读取计数数据，使用行名作为列名
  colnames(count_data) <- count_data[2,] # 设置列名，跳过第一行
  count_data <- count_data[-(1:6),] # 删除前六行，这些行包含 header 信息
  
  # 提取感兴趣的TPM数据列
  tpm_data <- count_data[7]
  
  # 设置TPM数据的列名
  sample_id <- samples_df$sample_id[which(samples_df$file_name == file_names_only[file_index])][[1]]
  colnames(tpm_data) <- sample_id
  
  if (nrow(expr_matrix) == 0) {
    # 如果expr_matrix还是空的，直接赋值
    expr_matrix <- tpm_data
  } else {
    # 检查行名是否一致，确保数据的一致性
    if (!all(rownames(expr_matrix) == rownames(tpm_data))) {
      stop("行名不一致，无法合并数据")
    }
    expr_matrix <- cbind(expr_matrix, tpm_data)
  }
} 
# 读取基因信息
gene_info <- as.matrix(read.delim(paste0('gdc_download_20250705_074012.501028/', count_file_paths[1]), fill = TRUE, header = FALSE, row.names = 1))
genes <- gene_info[-(1:6), 1] # 获取基因名称
gene_types <- gene_info[-(1:6), 2] # 获取基因类型

# 添加基因符号和类型
expr_matrix <- cbind(gene_type = gene_types, gene_symbol = genes, expr_matrix) # 将基因类型和符号添加到表达矩阵中

# 聚合数据，保留最大表达量
Bexpr_matrix <- aggregate(. ~ gene_symbol, data = expr_matrix, max) # 对每个基因符号的最大表达量进行聚合
expr_matrix <- expr_matrix[, 1:89]  # 保留第1到89列

# 整合
expr_matrix <- semi_join(expr_matrix, Bexpr_matrix, by = "gene_symbol") # 合并原始表达矩阵和聚合后的表达矩阵
expr_matrix <- Bexpr_matrix # 将聚合后的矩阵赋值给expr_matrix


# 将gene_symbol列设为行名,并转化为导出格式
rownames(expr_matrix) <- expr_matrix[, "gene_symbol"] # 将基因符号列设为行名
expr_matrix <- expr_matrix[, -1] # 移除基因符号列
expr_matrix <- data.frame(ID = rownames(expr_matrix), expr_matrix) # 创建ID列，并将表达数据转化为data.frame格式
colnames(expr_matrix) <- gsub('[.]', '-', colnames(expr_matrix)) # 将列名中的点替换为连字符，以符合导出文件的要求
expr_matrix <- subset(x = expr_matrix, gene_type == "protein_coding") # 筛选出蛋白质编码基因的表达数据
expr_matrix <- expr_matrix[, -1] # 移除多余的列
expr_matrix <- expr_matrix[, -1] # 再次移除多余的列（可能是由于上一步筛选导致的额外列）
# 导出最终表格
write.table(expr_matrix, 'TARGET-OS.txt', sep="\t", quote=FALSE, col.names = NA) # 将最终的表达矩阵导出为tsv格式文件