rm(list = ls())
setwd('E:/work/240624_brca/04_variance_reduce_cox/')
load('E:/work/240624_brca/03_delta_rank/result/brca_pagerank_drm_zscore.rda')
load('E:/work/240624_brca/03_delta_rank/result/train_test_samples.rda')

## 加载包
library(matrixStats)

## DRM矩阵，行为基因，列为样本
## 计算每行（每个基因）的变异性
geneVariability <- rowVars(as.matrix(drmzs))
## 初始化一个列表来存储每次迭代的结果
results <- list()
## 从20%开始到2%，每次减少3%
for(percent in seq(0.20, 0.01, by = -0.01)){
  ## 计算要提取的基因数
  topPercentIndex <- order(geneVariability, decreasing = T)[1:floor(length(geneVariability)*percent)]
  ## 获取这些基因的表达数据
  topPercentGenes <- drmzs[topPercentIndex, ]
  ## 将结果保存到列表中

  traindrmzs <- topPercentGenes[, trainsamp]
  testdrmzs <- topPercentGenes[, testsamp]
  results[[paste0('top', percent*100, '%')]] <- list(train = traindrmzs, test = testdrmzs)
}

save(results, file = 'result/variance_reduce_train_test.rda')
