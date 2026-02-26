setwd("D:/骨肉瘤") #设置工作路径到上图的路径下

# 载入必要的R包,用于数据处理
library(jsonlite)
library(dplyr) 

# 从JSON文件中读取数据
json_data <- fromJSON("metadata.cart.2025-07-05.json")

# 提取ID和case_id
ID <- sapply(json_data$associated_entities,
             function(x){x[,1]})
case_id <- sapply(json_data$associated_entities,
                  function(x){x[,3]})

# 将提取的ID和case_id合并成一个数据框
sample_case <- data.frame(ID, case_id, stringsAsFactors = FALSE)

# 读取临床信息数据
clinical_data1 <- read.delim("clinical.cart.2025-07-05\\clinical.tsv", header = TRUE, sep = "\t")
clinical_data2 <- read.delim("clinical.cart.2025-07-05\\follow_up.tsv", header = TRUE, sep = "\t")
clinical_data3 <- read.delim("clinical_PANCAN_patient_with_followup.tsv", header = TRUE, sep = "\t")

# 以下代码是针对clinical_data2文件进行处理
## 由于数据库再次更新，下面这段代码我视频里提到的列名需要以自己文件的实际列名为主
## 比如case_id要看看是否是cases.case_id还是其他，days_to_follow_up是否是follow_ups.days_to_follow_up
## follow_ups.timepoint_category是否是对的
clinical_data2_cleaned <- clinical_data2 %>%
  group_by(cases.case_id) %>%
  arrange(desc(follow_ups.timepoint_category == "Last Contact"), desc(follow_ups.days_to_follow_up)) %>%
  slice(1) %>%  # 每组保留第一行
  ungroup()

# 以下代码是针对clinical_data1文件进行处理
## 由于数据库再次更新，下面这段代码我视频里提到的列名需要以自己文件的实际列名为主
## 比如case_id要看看是否是cases.case_id还是其他，demographic.days_to_death是否是demographic.days_to_death还是其他
clinical_data1_cleaned <- clinical_data1 %>%
  group_by(cases.case_id) %>%
  filter(
    # 如果 demographic.days_to_death 相同，保留第一行
    (n() > 1 & duplicated(demographic.days_to_death)) |
      # 如果 demographic.days_to_death 是缺失值，保留第一行
      (n() > 1 & is.na(demographic.days_to_death)) |
      # 如果没有重复，保留所有行
      (n() == 1)
  ) %>%
  slice(1) %>%  # 每组保留第一行
  ungroup()

## 注意！！！
## 注意！！！
# 下面代码是通过case_id合并两个数据集，那么你就要问自己，你的是否是case_id还是cases.case_id
#  要以文件里实际的列名为准！
clinical_data <- full_join(clinical_data1_cleaned, clinical_data2_cleaned, by = "cases.case_id")

## 注意！！！
## 注意！！！
## 由于数据库再次更新，下面这段代码我视频里提到的列名，如cases.case_id、demographic.days_to_death、follow_ups.days_to_follow_up
## follow_ups.timepoint_category等列名需要以自己文件的实际列名为准！
clinical_data <- clinical_data %>%
  select(cases.case_id, demographic.days_to_death, follow_ups.days_to_follow_up,follow_ups.timepoint_category, everything())

#合并days_to_last_followup
## 下面代码意思是把clinical_data3文件中的bcr_patient_barcode, days_to_last_followup按照clinical_data文件中的cases.submitter_id.x合并
##如bcr_patient_barcode、days_to_last_followup以及cases.submitter_id.x列要以实际的列名为准！
clinical_data <- clinical_data %>%
  left_join(clinical_data3 %>% select(bcr_patient_barcode, days_to_last_followup),
            by = c("cases.submitter_id.x" = "bcr_patient_barcode")) %>%
  select(cases.submitter_id.x, days_to_last_followup, everything())


# 下面代码是通过case_id将sample_case和clinical_data合并，那你就要问问你自己，你这两个数据框
# 是否都有一个相同列名为case_id的东西，如果不是，以实际的列名为准！
# 如果你的sample_case文件没有case_id，则通过以下代码将case_id改成你想要的列名
# 运行前记得将代码前面的符号####去掉！！！！！！！！！！！
colnames(sample_case)[colnames(sample_case) == "case_id"] <- "cases.case_id"

# 修改完成列名，或者不需要修改之后，运行下面的代码，代码中的cases.case_id要以实际为准
clinical <- left_join(sample_case, clinical_data, by = "cases.case_id")

# 去除不必要的列
clinical <- clinical[,-2]

# 将处理后的临床数据保存为CSV文件
write.csv(clinical, "clinical.csv", row.names = FALSE)
