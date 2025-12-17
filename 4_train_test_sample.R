rm(list = ls())
setwd('E:/work/240624_brca/03_delta_rank/')

## 根据生存数据对样本分层抽样，划分training set和testing set
load('E:/work/240624_brca/02_ssgsea/result/brca_pagerank_ssgsva_scale.rda')
load('E:/work/240624_brca/01_page_rank/data/brca_sta_tim_msk.rda')

data <- cbind(sta_tim_msk[colnames(pagerank_NES), 2:3], t(pagerank_NES))

set.seed(123)
samples0 <- rownames(data)[which(data$os_status==0)]
samples1 <- rownames(data)[which(data$os_status==1)]
samples07 <- sample(samples0, round(length(samples0)*0.7), replace = FALSE)
samples17 <- sample(samples1, round(length(samples1)*0.7), replace = FALSE)
samples03 <- setdiff(samples0, samples07)
samples13 <- setdiff(samples1, samples17)
trainsamp <- c(samples07, samples17)
testsamp <- c(samples03, samples13)
save(trainsamp, testsamp, file = "result/train_test_samples.rda")

###
# 1756个train样本，753个test样本