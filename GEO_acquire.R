#设置工作目录
setwd("D:/骨肉瘤")
#加载
library(GEOquery)
library(dplyr)
#读取文件
gset <- getGEO("GSE39055",destdir = "D:/骨肉瘤",AnnotGPL = F,getGPL = F) 

#获取表达矩阵，下方数字需要按照实际情况修改
GSE=exprs(gset[[1]])
GSE=as.data.frame(GSE)

###判断是否需要处理
{
  ex <- GSE
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  
  if (LogC) { ex[which(ex <= 0)] <- NaN
  GSE <- log2(ex)
  print("log2 transform finished")}else{print("log2 transform not needed")}
}
#读入txt注释文件
gpl=data.table::fread("GPL14951-11332.txt",
                      header = TRUE,sep = "\t")
###处理
gpl = gpl[5:nrow(gpl), ]
colnames =as.character(gpl[1, ])
gpl = gpl[-1, ]
colnames(gpl) = colnames

#提取探针ID及基因symbol，同时删掉不必要的东西
ids=gpl[,c("ID","Symbol")]

##删除//后的所有内容，比如VSS//SF，基因则为VSS
##删除//后的所有内容，比如VSS//SF，基因则为VSS
##删除//后的所有内容，比如VSS//SF，基因则为VSS
ids$`Symbol` = gsub("//.*", "", ids$`Symbol`)

#以//为分割，保留第二个元素。比如ABC//CDF//EFF，会保留CDF
#以//为分割，保留第二个元素。比如ABC//CDF//EFF，会保留CDF
#以//为分割，保留第二个元素。比如ABC//CDF//EFF，会保留CDF
#以//为分割，保留第二个元素。比如ABC//CDF//EFF，会保留CDF
#ids$`Gene Symbol` <- sapply(strsplit(ids$`Gene Symbol`, "//"), function(x) x[2])

###合并
ids=ids[ids$ID %in% rownames(GSE),]
GSE=GSE[ids$ID,]
table(rownames(GSE) == ids$ID)
colnames(ids)=c('ID','GeneSymbol')
GSE <- bind_cols(ids, GSE)

###多个探针对应一个基因，检查是否有重复，如无，则返回false，有则返回True
any(duplicated(GSE$GeneSymbol))
####去重方法1.保留第一个出现的基因
####去重方法1.保留第一个出现的基因
####去重方法1.保留第一个出现的基因
GSE <- subset(GSE, !duplicated(GSE$GeneSymbol))
##再次检查
any(duplicated(GSE$GeneSymbol))

####去重方法2，保留平均值，且样本开头以GSM开头
####去重方法2，保留平均值，且样本开头以GSM开头
####去重方法2，保留平均值，且样本开头以GSM开头
GSE <- GSE %>%
  group_by(GeneSymbol) %>%
  summarise(across(starts_with("GSM"), mean, na.rm = TRUE)) %>%
  ungroup()
any(duplicated(GSE$GeneSymbol))

####去重方法3，保留表达量最大，且样本开头以GSM开头
####去重方法3，保留表达量最大，且样本开头以GSM开头
####去重方法3，保留表达量最大，且样本开头以GSM开头
####去重方法3，保留表达量最大，且样本开头以GSM开头
GSE <- GSE %>%
  group_by(GeneSymbol) %>%
  slice_max(order_by = rowMeans(across(starts_with("GSM"))), n = 1) %>%
  ungroup()
any(duplicated(GSE$GeneSymbol))
####去重方法4，保留表达量最小，且样本开头以GSM开头
####去重方法4，保留表达量最小，且样本开头以GSM开头
####去重方法4，保留表达量最小，且样本开头以GSM开头
{
  GSE <- GSE %>%
    group_by(GeneSymbol) %>%
    slice_min(order_by = rowMeans(across(starts_with("GSM"))), n = 1) %>%
    ungroup()
}


##删除不必要的列
GSE = GSE[,-1]
###再次检查GSE的GeneSymbol列是否有异常
View(GSE)
##删掉异常值
install.packages()
#GSE = GSE[-515,]

#导出
write.table(GSE, file = "GSE39055.txt", sep = "\t", row.names = F)

#获取临床信息
clinical=pData(gset[[1]])

#输出临床信息
write.csv(clinical,'clinical_GSE39055.csv',row.names = TRUE)
