##################################3
#tSNE二维：
#行为样本，列为特征的矩阵
rm(list = ls())
load("E:/work/240624_brca/05_sdcn/3000/train/001/data/train_2353_1631_3984.rda")
cluster <- read.table("E:/work/240624_brca/05_sdcn/3000/train/001/result/3984/result/sdcn/pred2.txt", header = F, sep = "\t")

allcolour <- c("#E15759","#59A14F")



setwd("E:/work/240624_brca/26_tsne/train/")
library(Rtsne)
library(ggplot2)
library(magrittr)
library(dplyr)
library(tidydr)
set.seed(123)
tsne_out = Rtsne(train3984, 
                 perplexity = 30,
                 max_iter = 1000, 
                 theta = 0.3, 
                 check_duplicates = FALSE, 
                 eta = 50)
pdat = data.frame(tsne_out$Y)
pdat$Subtype = cluster[, 1]
colnames(pdat)[1:2] = c("tSNE_1","tSNE_2")
head(pdat)
pdat$Subtype <- ifelse(pdat$Subtype==0, 'Cluster 1', 'Cluster 2')
rownames(pdat) <- rownames(train3984)

label <- pdat %>%
  group_by(Subtype) %>% 
  summarise(tSNE_1 = median(tSNE_1), tSNE_2 = median(tSNE_2))



#ggplot2绘图：
mytheme <- theme_void() + theme(plot.margin = margin(5.5,15,5.5,5.5))

p <- ggplot(data = pdat, aes(x = tSNE_1, y = tSNE_2)) +
  # stat_ellipse(aes(fill = Subtype), level = 0.95, linetype = 1,
  #              show.legend = F, geom = 'polygon', alpha = 0.1) +
  geom_point(aes(color = Subtype), size = 2, alpha = 0.5) +
  mytheme +
  theme_dr(xlength = 0.2, ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"),
                               ends = 'last', type = "closed")) +
  theme(panel.grid = element_blank()) +
  geom_text(data = label, aes(x = tSNE_1, y = tSNE_2, label = Subtype), fontface = "bold", family = 'serif',  color = 'black', size = 5.5) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_fill_manual(values = allcolour) +
  scale_color_manual(values = allcolour) + 
  theme(axis.title.x = element_text(size=15,family="serif",color="black",face="bold", hjust = 0.05), 
        axis.title.y = element_text(size=15,family="serif",color="black",face="bold", hjust = 0.05),
        legend.title=element_text(size=16,family="serif",color="black",face="bold"),
        legend.text=element_text(size=13,family="serif",color="black",face="bold"))
p
ggsave('result/tsne_2cluster_train.pdf', width = 7.5, height = 7.5)
save(p, file = 'result/tsne_2cluste_train.rda')




##################################3
#tSNE三维：
#行为样本，列为特征的矩阵
rm(list = ls())
load("E:/work/240624_brca/05_sdcn/3000/train/001/data/train_2353_1631_3984.rda")
cluster <- read.table("E:/work/240624_brca/05_sdcn/3000/train/001/result/3984/result/sdcn/pred2.txt", header = F, sep = "\t")

allcolour <- c("#E15759","#59A14F")



setwd("E:/work/240624_brca/26_tsne/train/")
library(Rtsne)
library(ggplot2)
library(magrittr)
set.seed(123)
tsne_out = Rtsne(train3984, 
                 perplexity = 30,
                 max_iter = 1000, 
                 theta = 0.3, 
                 check_duplicates = FALSE, 
                 eta = 50, 
                 dims = 3)
pdat = data.frame(tsne_out$Y)
pdat$Subtype = cluster[, 1]
colnames(pdat)[1:3] = c("tSNE_1","tSNE_2", "tSNE_3")
head(pdat)
pdat$Subtype <- ifelse(pdat$Subtype==0, 'Cluster 1', 'Cluster 2')
rownames(pdat) <- rownames(train3984)

label <- pdat %>%
  group_by(Subtype) %>% 
  summarise(tSNE_1 = median(tSNE_1), tSNE_2 = median(tSNE_2), tSNE_3 = median(tSNE_3))



#ggplot2绘图：
#3DtSNE图绘制：
p <- plot_ly(
  data = pdat,
  x = ~tSNE_1,
  y = ~tSNE_2,
  z = ~tSNE_3,
  type = 'scatter3d',
  mode = 'markers',
  marker = list(size = 2),
  color = Subtype,
  colors = allcolour,
  text = Subtype,
  hoverinfo = 'text'
) 
p1 <- p %>% layout(
  xaxis = list(title = list(text = "Sepal Length", font = list(family = "serif", size = 18, color = "black")),
               tickfont = list(family = "serif", size = 12, color = "black")),
  yaxis = list(title = list(text = "Sepal Width", font = list(family = "serif", size = 18, color = "black")),
               tickfont = list(family = "serif", size = 12, color = "black")), 
  legend = list(font = list(family = "serif", face = 'bold', size = 14, color = "black"))
)
