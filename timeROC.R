library(timeROC)
library(survival)

data1 <- read.delim("timeROC.txt",
                    row.names = 1,
                    header = TRUE,      # 第一行为列名
                    sep = "\t",         # 分隔符为制表符
                    stringsAsFactors = FALSE)  # 不自动转换字符为因子

## 构建timeROC
ROC <- timeROC(T=data1$os.time, #生存时间
               delta=data1$os_status,   #生存状态
               marker=data1$riskScore, #计算timeROC的变量
               cause=1,                #阳性结局指标数值(1表示死亡)
               weighting="marginal",   #计算方法，默认为marginal
               times=c(1, 3, 5),       #时间点，选取1年，3年和5年的生存率
               iid=TRUE)
ROC

pdf("timeROC曲线.pdf",width = 12,height = 8)
plot(ROC,
     time=1, col="red", lty=1,lwd=2, title = "")   #time是时间点，col是线条颜色、lty为图例线条类型、lwd为图例线条宽度
plot(ROC,
     time=3, col="blue", add=TRUE, lty=1,lwd=2)    #add指是否添加在上一张图中
plot(ROC,
     time=5, col="orange", add=TRUE, lty=1,lwd=2)
## 添加图例
legend("bottomright",#图例画在右下角
       c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],2)), #提取1年AUC构建图例标签
         paste0("AUC at 3 year: ",round(ROC[["AUC"]][2],2)), #提取3年AUC构建图例标签
         paste0("AUC at 5 year: ",round(ROC[["AUC"]][3],2))),#提取5年AUC构建图例标签
       col=c("red",
             "blue",
             "orange"), #设置1，3，5年AUC图例标签的图例颜色，注意与曲线保持对应
       lty=1,  
       lwd=2,  
       bty = "n" #o表示用框框把图例部分框起来，为默认。n表示不画框框
)
dev.off()