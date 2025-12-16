
####################################################
rm(list=ls())
library(HGNChelper)
library(BioNet)
library(dnet)
library(igraph)
library(diffusr)
setwd("D:\\OneDrive\\yanglei004\\BRCA-mutation\\network diffusion\\string0.7\\")
load("string0.7.rda")
load("BRCA.msk.met.2021.MM.rda")
setwd("D:\\OneDrive\\yanglei004\\BRCA-mutation\\network diffusion\\string0.7\\page rank\\")


mutCount_BRCA<-BRCA.msk.met.2021.MM
mutCount_BRCA[mutCount_BRCA > 1] <- 1

string0.7<-x1[,c(1,2)]
string0.7<-as.matrix(string0.7)

network<-graph.edgelist(string0.7, directed=FALSE)
network2<-igraph.to.graphNEL(network)
length(nodes(network2))

nodes<-matrix(0,length(nodes(network2)),1)
rownames(nodes)<-nodes(network2)

matrix<-matrix(0,ncol=ncol(mutCount_BRCA),nrow=length(nodes(network2)))
rownames(matrix)<-nodes(network2)
colnames(matrix)<-colnames(mutCount_BRCA)

i=1
try(for (i in 1:ncol(mutCount_BRCA)){
  A1<-as.matrix(mutCount_BRCA[,i])
  rownames(A1)<-rownames(mutCount_BRCA)
  mut_genes<-rownames(A1)[which(A1==1)]
  A2<-nodes
  A2[rownames(A2)%in%mut_genes,]<-1
  A3<-page_rank(network,personalized = A2)
  A4<-matrix(A3$vector)
  matrix[,i]<-A4
},silent = TRUE)

length(which(matrix[1,]==0))

matrix1<-matrix[,-which(matrix[1,]==0)]
save(matrix1,file = "page_rank_matrix.rda")

#######################################