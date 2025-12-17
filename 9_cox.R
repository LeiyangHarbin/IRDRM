rm(list = ls())
setwd('E:/work/240624_brca/04_variance_reduce_cox')
load('result/variance_reduce_train_test.rda')
load('E:/work/240624_brca/01_page_rank/data/brca_sta_tim_msk.rda')
train20 <- results[['top20%']]$train
train20 <- cbind(sta_tim_msk[colnames(train20), 2:3], t(train20))
train20$os_days <- as.numeric(train20$os_days)
train20$os_status <- as.numeric(train20$os_status)

library(survival)
coxoutput <- NULL
for (i in 3:ncol(train20)) {
  cox<-coxph(Surv(os_days, os_status)~as.numeric(train20[, i]) > median(as.numeric(train20[, i])), data = train20)
  coxSummary = summary(cox)
  coxoutput = rbind.data.frame(coxoutput,
  							               data.frame(Term = colnames(train20)[i],
                                          HR = as.numeric(coxSummary$coefficients[,"exp(coef)"]),
                                          z = as.numeric(coxSummary$coefficients[,"z"]),
                                          pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"]),
                                          lower = as.numeric(coxSummary$conf.int[,3]),
                                          upper = as.numeric(coxSummary$conf.int[,4]),
                                          stringsAsFactors = F))
}
head(coxoutput)
rownames(coxoutput) <- coxoutput[,1]
save(coxoutput, file = "result/brca_var20edge_coxoutput.rda")


rm(list = ls())
setwd('E:/work/240624_brca/04_variance_reduce_cox')
load('result/brca_var20edge_coxoutput.rda')
load('result/variance_reduce_train_test.rda')

varReducCox <- list()
for(percent in seq(0.20, 0.01, by = -0.01)){
  varReducCox[[paste0('top', percent*100, '%')]] <- coxoutput[rownames(results[[paste0('top', percent*100, '%')]]$train), ]
}
save(varReducCox, file = 'result/brca_varreduce_coxoutput.rda')