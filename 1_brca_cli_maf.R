
rm(list=ls())
setwd("E:/work/240624_brca/01_page_rank/data/")

## 加载临床数据
load("msk_met_2021_clinical.rda")
brca_msk_clinical <- msk_met_2021_clinical[which(msk_met_2021_clinical$cancer_type=="breast cancer"), ]
sta_tim_msk <- brca_msk_clinical[, c("sample_id", "os_days", "os_status")]
sta_tim_msk[which(sta_tim_msk$os_status == "dead"), 3] <- 1
sta_tim_msk[which(sta_tim_msk$os_status == "alive"), 3] <- 0

## 样本筛选(根据生存时间是否为NA)
sta_tim_msk <- na.omit(sta_tim_msk)
brca_msk_clinical <- brca_msk_clinical[rownames(sta_tim_msk), ]


## 加载突变数据（4605个样本）
library(maftools)
load("msk_met_2021_maf.rda")
msk_met_2021_maf@data$Tumor_Sample_Barcode <- as.character(msk_met_2021_maf@data$Tumor_Sample_Barcode)
brca_maf <- subsetMaf(maf = msk_met_2021_maf, tsb = brca_msk_clinical$sample_id)
final_sample <- unique(as.character(brca_maf@data$Tumor_Sample_Barcode))
length(final_sample)


## 样本筛选（根据是否存在突变数据）
brca_msk_clinical <- brca_msk_clinical[final_sample, ]
sta_tim_msk <- sta_tim_msk[final_sample, ]


## 保存数据：临床、生存、突变
save(brca_msk_clinical, file = "brca_msk_clinical.rda")
save(sta_tim_msk, file = "brca_sta_tim_msk.rda")
save(brca_maf, file = "brca_maf_msk.rda")


