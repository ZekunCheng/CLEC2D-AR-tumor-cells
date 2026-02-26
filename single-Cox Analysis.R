#设置工作目录
setwd("D:/骨肉瘤")

library(survival)
library(survminer)


TPM="TARGET-OS.txt"   
cli="TARGET-OS-Time.txt"
Deg="hdWGCNA_hubgene.txt"
GROUP="TARGET-OS-Group.txt"
outFile="单因素forest.pdf"    

#读取输入文件
data=read.table(TPM,header=T,sep="\t",row.names=1,check.names=F)
data2=read.table(Deg,header=T,sep="\t",row.names=1,check.names=F)
sample_info=read.table(GROUP,header=T,sep="\t",row.names=1,check.names=F)
cli=read.table(cli,header=T,sep="\t",row.names=1,check.names=F)

####数据处理#########
#转化为matrix
data_numeric <- as.data.frame(lapply(data, as.numeric))
rownames(data_numeric) <- rownames(data)
colnames(data_numeric) <- colnames(data)
data = data_numeric
#合并表达矩阵
data = data[rownames(data2),]
#临床数据处理
cli$time=cli$time/365
#设置筛选标准
p.value = 0.05

# 将样本的组别信息转换为因子
factors <- factor(sample_info$Group)
# 提取唯一的组别
groups <- unique(sample_info$Group)
# 显示组别
groups
# 更新样本信息中的组别为因子形式
sample_info$Group <- factors
# 显示更新后的样本组别信息
sample_info$Group

#仅保留肿瘤样本
tumor_samples <- rownames(sample_info)[sample_info$Group == "Tumor"]

# 仅保留 data 中这些样本的列
data <- data[, tumor_samples]


#转置
data=t(data)

#获取共同样本
sameSample <- intersect(rownames(data), rownames(cli)) 
rt <- cbind(cli[sameSample, ], data[sameSample, ]) 


# 初始化空列表来存储结果
results <- list()
unique_counts <- apply(rt[, 3:ncol(rt)], 2, function(x) length(unique(x)))

# 筛选出需要保留的列（唯一值 > 1 的列）
valid_cols <- names(unique_counts)[unique_counts > 1]  # 保留非单一值列
rt <- rt[, c("time", "state", valid_cols)]      # 保留时间、状态和有效基因列


# 循环进行Cox回归分析
for (i in colnames(rt)[3:ncol(rt)]) {
  # 进行Cox回归分析
  cox <- coxph(Surv(time, state) ~ rt[, i], data = rt)
  coxSummary <- summary(cox)
  
  # 提取p值
  coxP <- coxSummary$coefficients[, "Pr(>|z|)"]
  
  # 检查p值是否小于阈值
  if (coxP < p.value) {
    results[[length(results) + 1]] <- with(coxSummary, {
      data.frame(
        id = i,
        HR = conf.int[, "exp(coef)"],
        HR.95L = conf.int[, "lower .95"],
        HR.95H = conf.int[, "upper .95"],
        pvalue = coefficients[, "Pr(>|z|)"]
      )
    })
  }
}

# 将结果列表转换为数据框
unicox <- do.call(rbind, results)

#保存结果
write.table(unicox, file="unicox.txt", sep="\t", row.names=F, quote=F)
##画图
dup_values <- unicox[, 1][duplicated(unicox[, 1])]
table(unicox[, 1])[table(unicox[, 1]) > 1]
unicox <- unicox[!duplicated(unicox[, 1]), ]  # 保留首次出现的行
rownames(unicox) <- unicox[, 1]
#将unicox文件的第一列设置为行名
rownames(unicox) <- unicox[, 1]
#移除第一列，因为它已经被用作行名
unicox <- unicox[, -1]
###转化为数值类型
unicox$HR <- as.numeric(unicox$HR)
unicox$HR.95L <- as.numeric(unicox$HR.95L)
unicox$HR.95H <- as.numeric(unicox$HR.95H)
unicox$pvalue <- as.numeric(unicox$pvalue)
gene <- rownames(unicox)
hr <- sprintf("%.3f",unicox$"HR")
hrLow  <- sprintf("%.3f",unicox$"HR.95L")
hrHigh <- sprintf("%.3f",unicox$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(unicox$pvalue<0.001, "<0.001", sprintf("%.3f", unicox$pvalue))

#绘制森林图
# 打开PDF输出文件，并根据行数设置图形高度
n <- nrow(unicox)
pdf(outFile, width = 7, height = 15)

# 设定图形整体布局和边界
layout(matrix(c(1, 2), nc = 2), widths = c(3, 2.5))
par(mar = c(4, 2.5, 2, 1))

# 绘制基因名和P值
plot(1, type = "n", xlim = c(0, 3), ylim = c(1, n + 1), axes = FALSE, xlab = "", ylab = "")
text.cex <- 0.8
text(0, n:1, gene, adj = 0, cex = text.cex)
text(1.5, n:1, pVal, adj = 1, cex = text.cex)
text(1.5, n + 1, 'pvalue', cex = text.cex, font = 2, adj = 1)

# 绘制风险比和标题
text(3, n:1, Hazard.ratio, adj = 1, cex = text.cex)
text(3, n + 1, 'Hazard ratio', cex = text.cex, font = 2, adj = 1)

# 切换到第二个图形的绘制
par(mar = c(4, 1, 2, 1), mgp = c(2, 0.5, 0))
max_hr <- max(as.numeric(hrLow), as.numeric(hrHigh))
plot(1, type = "n", xlim = c(0, max_hr), ylim = c(1, n + 1), axes = FALSE, xlab = "Hazard ratio", ylab = "")
arrows(as.numeric(hrLow), n:1, as.numeric(hrHigh), n:1, angle = 90, code = 3, length = 0.05, col = "darkblue", lwd = 2.5)
abline(v = 1, col = "black", lty = 2, lwd = 2)

# 绘制中间点和根据风险比设定点颜色
boxcolor <- ifelse(as.numeric(hr) > 1, "pink", "blue")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex = 1.6)
axis(1)

# 关闭PDF设备
dev.off()


##挑选一部分画图################
##挑选一部分画图################
##挑选一部分画图################
##挑选一部分画图################
unicox <- do.call(rbind, results)
unicox = unicox[1:15,]


#将unicox文件的第一列设置为行名
rownames(unicox) <- unicox[, 1]
#移除第一列，因为它已经被用作行名
unicox <- unicox[, -1]
###转化为数值类型
unicox$HR <- as.numeric(unicox$HR)
unicox$HR.95L <- as.numeric(unicox$HR.95L)
unicox$HR.95H <- as.numeric(unicox$HR.95H)
unicox$pvalue <- as.numeric(unicox$pvalue)
gene <- rownames(unicox)
hr <- sprintf("%.3f",unicox$"HR")
hrLow  <- sprintf("%.3f",unicox$"HR.95L")
hrHigh <- sprintf("%.3f",unicox$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(unicox$pvalue<0.001, "<0.001", sprintf("%.3f", unicox$pvalue))

#绘制森林图
# 打开PDF输出文件，并根据行数设置图形高度
n <- nrow(unicox)
pdf(file = "15单因素cox图.pdf", width = 7, height = n / 13 + 5)

# 设定图形整体布局和边界
layout(matrix(c(1, 2), nc = 2), widths = c(3, 2.5))
par(mar = c(4, 2.5, 2, 1))

# 绘制基因名和P值
plot(1, type = "n", xlim = c(0, 3), ylim = c(1, n + 1), axes = FALSE, xlab = "", ylab = "")
text.cex <- 0.8
text(0, n:1, gene, adj = 0, cex = text.cex)
text(1.5, n:1, pVal, adj = 1, cex = text.cex)
text(1.5, n + 1, 'pvalue', cex = text.cex, font = 2, adj = 1)

# 绘制风险比和标题
text(3, n:1, Hazard.ratio, adj = 1, cex = text.cex)
text(3, n + 1, 'Hazard ratio', cex = text.cex, font = 2, adj = 1)

# 切换到第二个图形的绘制
par(mar = c(4, 1, 2, 1), mgp = c(2, 0.5, 0))
max_hr <- max(as.numeric(hrLow), as.numeric(hrHigh))
plot(1, type = "n", xlim = c(0, max_hr), ylim = c(1, n + 1), axes = FALSE, xlab = "Hazard ratio", ylab = "")
arrows(as.numeric(hrLow), n:1, as.numeric(hrHigh), n:1, angle = 90, code = 3, length = 0.05, col = "darkblue", lwd = 2.5)
abline(v = 1, col = "black", lty = 2, lwd = 2)

# 绘制中间点和根据风险比设定点颜色
boxcolor <- ifelse(as.numeric(hr) > 1, "pink", "blue")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex = 1.6)
axis(1)

# 关闭PDF设备
dev.off()



###按条件筛选HR大于等于2的数据
unicox <- do.call(rbind, results)
unicox <- unicox[unicox$HR >= 1.1 | unicox$HR <= 0.9, ]

#将unicox文件的第一列设置为行名
rownames(unicox) <- unicox[, 1]
#移除第一列，因为它已经被用作行名
unicox <- unicox[, -1]
###转化为数值类型
unicox$HR <- as.numeric(unicox$HR)
unicox$HR.95L <- as.numeric(unicox$HR.95L)
unicox$HR.95H <- as.numeric(unicox$HR.95H)
unicox$pvalue <- as.numeric(unicox$pvalue)
gene <- rownames(unicox)
hr <- sprintf("%.3f",unicox$"HR")
hrLow  <- sprintf("%.3f",unicox$"HR.95L")
hrHigh <- sprintf("%.3f",unicox$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(unicox$pvalue<0.001, "<0.001", sprintf("%.3f", unicox$pvalue))

#绘制森林图
# 打开PDF输出文件，并根据行数设置图形高度
n <- nrow(unicox)
pdf(file = "符合条件的单因素cox图.pdf", width = 7, height = n / 13 + 5)

# 设定图形整体布局和边界
layout(matrix(c(1, 2), nc = 2), widths = c(3, 2.5))
par(mar = c(4, 2.5, 2, 1))

# 绘制基因名和P值
plot(1, type = "n", xlim = c(0, 3), ylim = c(1, n + 1), axes = FALSE, xlab = "", ylab = "")
text.cex <- 0.8
text(0, n:1, gene, adj = 0, cex = text.cex)
text(1.5, n:1, pVal, adj = 1, cex = text.cex)
text(1.5, n + 1, 'pvalue', cex = text.cex, font = 2, adj = 1)

# 绘制风险比和标题
text(3, n:1, Hazard.ratio, adj = 1, cex = text.cex)
text(3, n + 1, 'Hazard ratio', cex = text.cex, font = 2, adj = 1)

# 切换到第二个图形的绘制
par(mar = c(4, 1, 2, 1), mgp = c(2, 0.5, 0))
max_hr <- max(as.numeric(hrLow), as.numeric(hrHigh))
plot(1, type = "n", xlim = c(0, max_hr), ylim = c(1, n + 1), axes = FALSE, xlab = "Hazard ratio", ylab = "")
arrows(as.numeric(hrLow), n:1, as.numeric(hrHigh), n:1, angle = 90, code = 3, length = 0.05, col = "darkblue", lwd = 2.5)
abline(v = 1, col = "black", lty = 2, lwd = 2)

# 绘制中间点和根据风险比设定点颜色
boxcolor <- ifelse(as.numeric(hr) > 1, "red", "blue")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex = 1.6)
axis(1)

# 关闭PDF设备
dev.off()
