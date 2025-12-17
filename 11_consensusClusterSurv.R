
#####################################################
rm(list=ls())
library(survminer)
library(survival)
library(ggplot2)
library(RColorBrewer)
setwd("E:/work/240624_brca/05_sdcn/3000/train/001")
load("data/train_surv_2353_1631_3984.rda")
load("E:/work/240624_brca/07_consensus_clustering/3984/train/result/consensus_2_cluster.rda")

consensusData <- trainSurv3984[, 1:2]
consensusData$cc <- consensusCluster[rownames(trainSurv3984), 1]
consensusData$os_days <- as.numeric(round(consensusData$os_days/30, 1))
consensusData$os_status <- as.numeric(consensusData$os_status)
consensusData$cc <- factor(as.character(consensusData$cc), levels = c("1", "2"))


fit<-survfit(Surv(`os_days`, `os_status`)~cc, data = consensusData)
Bcox<-coxph(Surv(`os_days`, `os_status`)~cc, data = consensusData)
summcph<-summary(Bcox)
summcph
summcph$sctest[3]

p1<-ggsurvplot(fit,legend.title="Subtype",
               main="times",
               legend.labs = paste0(paste0('Cluster ', 1:2, sep = ""), '   ', sep = ""),
               pval=F, conf.int = T,
               risk.table =T, # Add risk table
               #risk.table.col = "strata", # Change risk table color by groups
               linetype = "strata", # Change line type by groups
               palette = c('#F18C8D', '#90EE90'), 
               surv.median.line = "hv", # Specify median survival
               ggtheme=theme_survminer(base_size=13,base_family="serif",
                                       font.main=c(15, "bold","black"),
                                       font.submain= c(14, "bold", "black"),
                                       font.caption = c(14, "bold", "black"),
                                       font.x=c(15, "bold", "black"),font.legend = c(16,"bold"),
                                       font.y=c(15, "bold", "black"),
                                       font.tickslab=c(14, "bold")))

p2<-p1+labs(x="Month")
p3<-p2$plot+ggtitle("Consensus Clustering")+
  theme(plot.title=element_text(size=15,family="serif",color="black",face="bold",hjust=0.5))
p4<-p3+annotate("text", x=26, y=0.25, label="P=6.99e-11",color="black",cex=6)
p5<-p2$table+theme(legend.position="none")
p6<-ggarrange(p4,p5,heights=c(5,1.8),ncol=1,nrow=2, align = "v")
p6
survival_cc<-p6
ggsave("E:/work/240624_brca/07_consensus_clustering/3984/train/result/brca_c1_c2_survival_consensus_clustering.pdf",width = 6,height = 6)
save(survival_cc,file="E:/work/240624_brca/07_consensus_clustering/3984/train/result/brca_c1_c2_survival_consensus_clustering.rda")
#####################################

