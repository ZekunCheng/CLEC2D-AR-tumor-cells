setwd("D:/骨肉瘤")
##加载r包
library(survival)
library(survminer)
library(dplyr)

# 读取输入文件并处理数据
data <- read.table("风险因子评分.txt", header=TRUE, sep="\t", check.names=FALSE, row.names=1)
data <- data[data$time > 30, ]
data <- mutate(data, time = time / 365)
# 查看分组数量
high_risk_count <- sum(data$Risk == "High") 
low_risk_count <- sum(data$Risk == "Low") 
# 绘制风险评分图
plot_risk_distribution <- function(risk_data, risk_col, score_col, 
                                   low_color = "#0096FF", high_color = "pink") {
  # 排序数据
  risk_data <- risk_data[order(risk_data[[score_col]]), ]
  
  # 计算分组信息
  risk_groups <- risk_data[[risk_col]]
  low_count <- sum(risk_groups == "Low")
  high_count <- sum(risk_groups == "High")
  cutoff_score <- max(risk_data[[score_col]][risk_groups == "Low"])
  
  # 设置颜色
  point_colors <- ifelse(risk_groups == "Low", low_color, high_color)
  
  # 绘制图形
  plot(risk_data[[score_col]], type = "p", pch = 20,
       ylab = "Risk score",
       col = point_colors,
       main = "")
  abline(h = cutoff_score, v = low_count, lty = 2)
  legend("topright", c("High risk", "Low risk"), bty = "m", pch = 18, 
         col = c(high_color, low_color), cex = 1.3)
}

# 辅助函数：绘制生存状态图
plot_survival_status <- function(risk_data, time_col, status_col, 
                                 event_color = "pink", censor_color = "#0096FF") {
  # 排序数据（保持与风险评分图一致）
  risk_data <- risk_data[order(risk_data$riskScore), ]
  
  # 计算分组信息
  low_count <- sum(risk_data$Risk == "Low")
  
  # 设置颜色
  status_colors <- ifelse(risk_data[[status_col]] == 1, event_color, censor_color)
  
  # 绘制图形
  plot(risk_data[[time_col]], pch = 18,
       ylab = "Survival time (years)",
       col = status_colors,
       main = "")
  legend("topright", c("Dead", "Alive"), bty = "m", pch = 18, 
         col = c(event_color, censor_color), cex = 1.3)
  abline(v = low_count, lty = 2)
}

# 生成风险评分分布图
pdf("risk_distribution.pdf", width = 6, height = 5)
plot_risk_distribution(data, "Risk", "riskScore")
dev.off()

# 生成生存状态图
pdf("survival_status.pdf", width = 6, height = 5)
plot_survival_status(data, "time", "situation")
dev.off()
