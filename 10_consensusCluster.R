setwd("E:/work/240624_brca/07_consensus_clustering/3984/train/")
#读取数据，行是基因，列是样本
load("E:/work/240624_brca/05_sdcn/3000/train/001/data/train_2353_1631_3984.rda")
data1 <- t(train3984)

#做一致性聚类，聚两类
library(ConsensusClusterPlus)

results <- ConsensusClusterPlus(as.matrix(data1), maxK = 10, 
                                reps = 50, 
                                pItem = 0.8,
                                pFeature = 0.8, 
                                clusterAlg = "hc", 
                                distance = "pearson", 
                                tmyPal = c('white', '#17489c'), 
                                title = 'Consensus_cluster_for_BRCA', 
                                plot = "png", 
                                seed = 123)

# # 输出K=3时的一致性矩阵
# consensusTree <- results[[3]][["consensusTree"]]
# consensusTree
# 
# # hclust选项
# consensusMatrix <- results[[3]][["consensusMatrix"]]
# consensusMatrix[1:5, 1:5]

# 样本分类
consensusClass <- results[[2]][["consensusClass"]]
consensusCluster <- data.frame(cluster = as.numeric(consensusClass))
rownames(consensusCluster) <- names(consensusClass)
save(consensusCluster, file = "result/consensus_2_cluster.rda")

# #计算聚类一致性 (cluster-consensus) 和样品一致性 (item-consensus)
# icl <- calcICL(results, 
#                #title = title,
#                plot = "png")
# ## 返回了具有两个元素的list，然后分别查看一下
# dim(icl[["clusterConsensus"]])
# 
# icl[["clusterConsensus"]] 
# 
# dim(icl[["itemConsensus"]])
# 
# icl[["itemConsensus"]][1:5,] 




#相关参数：
# pItem (item resampling, proportion of items to sample) : 80%
# pFeature (gene resampling, proportion of features to sample) : 80%
# maxK (a maximum evalulated k, maximum cluster number to evaluate) : 6
# reps (resamplings, number of subsamples) : 50
# clusterAlg (agglomerative heirarchical clustering algorithm) : 'hc' (hclust)
# distance : 'pearson' (1 - Pearson correlation)