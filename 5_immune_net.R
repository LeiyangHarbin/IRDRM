rm(list = ls())
setwd('E:/work/240624_brca/03_delta_rank')

imm <- read.csv('data/Geneappend3.csv', 
                header = T)

## 由于imm中存在symbol和synonyms，综合这两列得到免疫相关基因
## symbol列1811个无重复基因，纳入synonyms共有7681个基因
immGene <- unique(c(imm$Symbol, 
                    unlist(strsplit(imm[which(imm$Synonyms!='-'), "Synonyms"], split = '\\|'))))


## 加载string网络，两个基因都属于免疫相关基因的边被纳入
load('E:/work/240624_brca/01_page_rank/data/string0.7.rda')
## 以下脚本证明string0.7为有向图
x1$edge1 <- paste(x1$protein1, x1$protein2, sep = "|")
x1$edge2 <- paste(x1$protein2, x1$protein1, sep = "|")
length(unique(x1$edge1)); length(unique(x1$edge2))
length(unique(c(x1$edge1, x1$edge2)))
table(table(c(x1$edge1)))
table(table(c(x1$edge1, x1$edge2)))

## 先将x1转变成无向图
# x2 <- NULL
# for(i in 1:nrow(x1)){
#   print(i)
#   if(i == 1){
#     x2 = rbind(x2, x1[i, ])
#     next
#   }
#   if(!(x1[i, 4]%in%x2[, 4]) & !(x1[i, 5]%in%x2[, 4])){
#     x2 = rbind(x2, x1[i, ])
#   }
# }
x2 <- x1[which(x1$protein1 >= x1$protein2), ] #需要去重操作
edge <- unique(x2$edge1)
x2 <- x2[match(edge, x2$edge1), ]

id1 <- which(x2$protein1%in%immGene)
id2 <- which(x2$protein2%in%immGene)
immNet <- x2[intersect(id1, id2), ]
save(immNet, file = 'E:/work/240624_brca/03_delta_rank/data/immune_net.rda')
