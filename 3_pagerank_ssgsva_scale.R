
##############################################
rm(list=ls())
library(GSVA)
setwd("E:/work/240624_brca/02_ssgsea/")
# brca_pagerank <- read.table('E:/work/240624_brca/01_page_rank/result/brca_pagerank.txt', 
#                              sep = "\t", 
#                              header = T, 
#                              check.names = F, 
#                              row.names = 1)
# brca_pagerank <- brca_pagerank[, -which(apply(brca_pagerank, 2, function(x){length(which(x==0))})==16134)]
# save(brca_pagerank, file = 'data/brca_pagerank.rda')
load('E:/work/240624_brca/02_ssgsea/data/brca_pagerank.rda')

#############################################
distCor <- function(x) as.dist(1-cor(x))
zClust <- function(x, scale="row", zlim=c(-3,3), method="average") {
  if (scale=="row") z <- t(scale(t(x))) else z<-x
  if (scale=="col") z <- scale(x)
  z <- pmin(pmax(z, zlim[1]), zlim[2])
  # hcl_row <- hclust(distCor(t(z)), method=method)
  # hcl_col <- hclust(distCor(z), method=method)
  return(list(data=z))
}
#------------------------------------------------------------

uni_matrix <- brca_pagerank
genes <- rownames(uni_matrix)
genes_type=list()
for (i in 1:length(genes)){
  b=as.character(genes[i])
  genes_type[i]=list(b)
  names(genes_type)[i]=genes[i]
}

matrix1 <- as.matrix(uni_matrix)
gsva_es <- gsva(matrix1, genes_type, method="ssgsea", abs.ranking=F)
id <- which(apply(gsva_es, 1, sd) == 0) #不存在基因对应的行的标准差为0，因此不计算gsva_es1和random_NES1
# gsva_es1 <- gsva_es[-which(apply(gsva_es, 1, sd) == 0),] #删除标准差为0的基因
z <- zClust(gsva_es) #归一化
# z1 <- zClust(gsva_es1)
pagerank_NES <- z$data
# random_NES1 <- z1$data

save(gsva_es, pagerank_NES, file = "result/brca_pagerank_ssgsva_scale.rda")

#gsva_es和random_NES分别是没有删除标准差为0的行的ssgsva结果和归一化后的结果
#gsva_es1和random_NES1分别是删除了标准差为0的行的ssgsva结果和归一化后的结果

