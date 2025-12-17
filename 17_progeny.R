rm(list = ls())
setwd('E:/work/240624_brca/37_progeny/')

# 加载包
library(progeny)
# 使用ssgsea结果中training set对应的样本
load('E:/work/240624_brca/02_ssgsea/result/brca_pagerank_ssgsva_scale.rda')
load('E:/work/240624_brca/21_cluster_diff/data/label_train.rda')
trainLab$cluster <- ifelse(trainLab$cluster==0, 'Cluster 1', 'Cluster 2')
trainLab$sample <- rownames(trainLab)
trainLab <- trainLab[order(trainLab$cluster, decreasing = F), ]
gene_expression <- pagerank_NES[, rownames(trainLab)]

# 计算通路活性
pathways <- progeny(gene_expression, scale=TRUE,
                    organism="Human",
                    top = 100, perm = 1)
save(pathways, file = 'data/pathways_14.rda')
save(trainLab, file = 'data/trainLab.rda')

