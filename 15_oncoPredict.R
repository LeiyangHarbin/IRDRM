
######################
rm(list = ls())
library(GEOquery)
library(oncoPredict)
setwd('E:/work/240624_brca/28_oncopredict/')
load('E:/work/240624_brca/02_ssgsea/result/brca_pagerank_ssgsva_scale.rda')
GDSC2_Expr<-readRDS(file = "data/calcPhenotype_Output/GDSC2_Expr (RMA Normalized and Log Transformed).rds") #oncoPredict包给出该文件
GDSC2_Res<-readRDS(file = "data/calcPhenotype_Output/GDSC2_Res.rds")  #用药反应数据
GDSC2_Res<-exp(GDSC2_Res)
testExpr<-pagerank_NES #自己的NES矩阵

calcPhenotype(trainingExprData = GDSC2_Expr,  #使用大规模的基因表达和药物筛选数据建立岭回归模型并应用于新的基因表达数据，实现预测临床化疗反应。
              trainingPtype = GDSC2_Res,  #需要输入的数据有训练集表达数据、训练集表型数据以及测试的表达谱数据
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )
#################################################
