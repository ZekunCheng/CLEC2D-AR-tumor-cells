setwd("D:/骨肉瘤")
#加载
library(GEOquery)
library(dplyr)

gset <- getGEO("GSE154540",destdir = "D:/骨肉瘤",AnnotGPL = F,getGPL = F) 

GSE=exprs(gset[[1]])
GSE=as.data.frame(GSE)

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

gpl=data.table::fread("GPL96-57554.txt",
                      header = TRUE,sep = "\t")

gpl = gpl[5:nrow(gpl), ]
colnames =as.character(gpl[1, ])
gpl = gpl[-1, ]
colnames(gpl) = colnames


ids=gpl[,c("ID","GENE_SYMBOL")]

ids$`GENE_SYMBOL` = gsub("//.*", "", ids$`GENE_SYMBOL`)


ids=ids[ids$ID %in% rownames(GSE),]
GSE=GSE[ids$ID,]
table(rownames(GSE) == ids$ID)
colnames(ids)=c('ID','GENE_SYMBOL')
GSE <- bind_cols(ids, GSE)

any(duplicated(GSE$GENE_SYMBOL))

GSE <- subset(GSE, !duplicated(GSE$GENE_SYMBOL))
any(duplicated(GSE$GENE_SYMBOL))

GSE = GSE[,-1]

View(GSE)

install.packages()

write.table(GSE, file = "GSE16088.txt", sep = "\t", row.names = F)

clinical=pData(gset[[1]])

write.csv(clinical,'clinical_GSE16088.csv',row.names = TRUE)
