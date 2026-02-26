 
library(survival)
library(survminer)
setwd("D:/骨肉瘤")
###前期准备
Gene="SH3BGRL3"  
bioCol=c("Firebrick3","MediumSeaGreen","#D5652C","#F00646")
data=read.table("TARGET生存分析.txt", header=T, sep="\t", check.names=F, row.names=1,quote = "")

## 选作：该部分代码为更新代码，与视频不一致
## 直接运行即可，前提是你的数据是TCGA数据且样本名是完整的样本名，比如TCGA-EJ-SDFA-01A-32R-3242-05,而不是TCGA-EJ-SDFA就没了

## 选作：该部分代码为更新代码，与视频不一致
## 直接运行即可，前提是你的数据是TCGA数据且样本名是完整的样本名，比如TCGA-EJ-SDFA-01A-32R-3242-05,而不是TCGA-EJ-SDFA就没了

## 选作：该部分代码为更新代码，与视频不一致
## 直接运行即可，前提是你的数据是TCGA数据且样本名是完整的样本名，比如TCGA-EJ-SDFA-01A-32R-3242-05,而不是TCGA-EJ-SDFA就没了
## 提取行名的第14个字符
sample_names <- rownames(data)
fourteenth_char <- substr(sample_names, 14, 14)  # 提取第14个字符
data <- data[fourteenth_char %in% c("0", "2"), , drop = FALSE]

## 以下为必作部分，与原视频一致
## 以下为必作部分，与原视频一致
## 以下为必作部分，与原视频一致
## 以下为必作部分，与原视频一致

##转化年份
data <- mutate(data, os.time = os.time / 365)

####按照中位数分组，MIR22HG需要替换成你的基因名
group = ifelse(data[,"CRIP1"]>quantile(data[,"CRIP1"], seq(0,1,1/2))[2],"High","Low")
data <- mutate(data, group = group)
length = length(levels(factor(group)))

##数据处理
data=data[,c("os.time","os_status",Gene,"group")]
colnames(data)=c("os.time","os_status","Gene")

#比较高低表达生存差异
diffsurvival=survdiff(Surv(os.time, os_status) ~group,data =data)
pValue=1-pchisq(diffsurvival$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.04f",pValue))
}
rt <- survfit(Surv(os.time, os_status) ~group,data = data)

#绘制生存曲线
surPlot=ggsurvplot(rt, 
                   data=data,
                   conf.int=TRUE,
                   pval.size=4,
                   pval=pValue,
                   legend.labs=c("High", "Low"),
                   legend.title="MIR22HG exp",
                   xlab="Time(years)",
                   risk.table=T,
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "green"),
                   risk.table.height=.24)
surPlot
dev.off()
####换一种形式
surPlot <- ggsurvplot(
  fit = rt,                # survfit 对象
  data = data,             # 包含分组信息的临床数据
  conf.int = TRUE,         # 显示置信区间
  pval = pValue,           # 显示 p 值
  pval.size = 4,           # p 值字体大小
  legend.labs = c("High", "Low"),   # 图例标签
  legend.title = "SH3BGRL3 exp",       # 图例标题
  xlab = "Time (years)",            # x轴标签
  break.time.by = 1,                # x轴间隔
  risk.table = TRUE,                # 显示风险表
  risk.table.title = "",            # 风险表标题为空
  risk.table.height = 0.24,         # 风险表高度
  palette = c("red", "blue")        # 颜色
)

# 保存为 PDF 文件
pdf("survival_plot.pdf", width = 7, height = 6)
print(surPlot)
dev.off()
####换一种形式
bioCol=bioCol[1:length]
surPlot=ggsurvplot(rt, 
                   data=data,
                   conf.int=F,
                   pval=pValue,
                   pval.size=5,
                   legend.title="MIR22HG exp",
                   legend.labs=c("High", "Low"),
                   legend = c(0.9, 0.9),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by = 2,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.24)
surPlot
