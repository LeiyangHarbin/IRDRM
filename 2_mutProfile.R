rm(list = ls())
library(maftools)
setwd("E:/work/240624_brca/01_page_rank/data/")

load('brca_msk_clinical.rda')
load('msk_cancer_type_list.rda')
class(msk_cancer_type_list)

# 查看数据
data <- msk_cancer_type_list[["breast cancer"]]@data; data$Tumor_Sample_Barcode <- as.character(data$Tumor_Sample_Barcode)
variantsPerSample <- msk_cancer_type_list[["breast cancer"]]@variants.per.sample
variantsTypeSummary <- msk_cancer_type_list[["breast cancer"]]@variant.type.summary
vcs <- msk_cancer_type_list[["breast cancer"]]@variant.classification.summary
geneSummary <- msk_cancer_type_list[["breast cancer"]]@gene.summary
summary <- msk_cancer_type_list[["breast cancer"]]@summary
mafSilent <- msk_cancer_type_list[["breast cancer"]]@maf.silent
clinicalData <- msk_cancer_type_list[["breast cancer"]]@clinical.data

# 构建空矩阵
mutProfile <- matrix(0,
                     nrow = length(unique(data$Hugo_Symbol)), 
                     ncol = length(unique(data$Tumor_Sample_Barcode))
)
row.names(mutProfile) <- as.character(unique(data$Hugo_Symbol))
colnames(mutProfile) <- as.character(unique(data$Tumor_Sample_Barcode))


# 根据突变信息填充矩阵
data <- as.data.frame(data)
for(i in 1:ncol(mutProfile)){
  sample <- colnames(mutProfile)[i]
  genes <- data[which(data$Tumor_Sample_Barcode==sample), 1]
  mutProfile[genes, sample] = 1
}
mutProfile <- mutProfile[, rownames(brca_msk_clinical)]

# 保存mutProfile
save(mutProfile, file = "brca_mutProfile.rda")
