setwd('E:/work/240624_brca/03_delta_rank/')
source('代码/3_DRM.R')

## 表达谱数据
load('E:/work/240624_brca/02_ssgsea/result/brca_pagerank_ssgsva_scale.rda')
load('data/immune_net.rda')
string <- immNet[, 1:2]
##

## 构建DRM
rm <- rank.matrix(pagerank_NES) # rank函数将表达量相同基因的秩取平均值
drm <- delta.rank(rm, string)
drmzs <- scale(drm)
##

save(drm, file = 'result/brca_pagerank_drm.rda')
save(drmzs, file = 'result/brca_pagerank_drm_zscore.rda')







'''
#-----外部评价指标
## adjusted rand index
## python实现
from sklearn import metrics
labels_true = [0, 0, 0, 1, 1, 1]
labels_pred = [0, 0, 1, 1, 2, 2]
 
metrics.adjusted_rand_score(labels_true, labels_pred)  

## R实现
labels_true = c(0, 0, 0, 1, 1, 1)
labels_pred = c(0, 0, 1, 1, 2, 2)
adjustedRandIndex(labels_true, labels_pred)


## normalized mutual information
## python实现
labels1 <- c(1, 1, 2, 2, 3, 3)
labels2 <- c(1, 1, 1, 2, 3, 3)
result_NMI=metrics.normalized_mutual_info_score(A, B)

## R实现
library(NMI)
NMI(c1, c2)
library(aricode)
aricode::NMI(labels1, labels2)


## purity
'''