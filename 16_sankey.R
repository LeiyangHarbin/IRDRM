rm(list = ls())
setwd('E:/work/240624_brca/41_SKC_sankey/')

load('data/sdcn_km_cc.rda')
km <- km[rownames(sdcn), ]
cc <- cc[rownames(sdcn), ]
label <- cbind(sdcn, km, cc)[, c(1, 3, 5)]
label1 <- label[, c(1, 2)]
label2 <- label[, c(1, 3)]

## sdcn-km、sdcn-cc
library(ggplot2)
library(ggsankey)
df <- label1 %>%
  make_long(`SDCN`, `Kmeans Clustering`)
col <- c('#ff7f50', '#1f77b4', '#5fb404', '#db7093', '#37ade7', '#db8f1e', '#48d1cc', '#9400d3', '#dda7c0')
sk <- ggplot(df, aes(x = x, 
                     next_x = next_x, 
                     node = node, 
                     next_node = next_node,
                     fill = factor(node),
                     label = node)) +
  geom_sankey(flow.alpha = 0.5, node.color = 1) +
  geom_sankey_label(size = 5, color = 1, fill = "white", family = 'serif') +
  scale_fill_manual(values = col) +
  labs(y = '', x = '') + 
  theme_sankey(base_size = 17) + 
  theme(text=element_text(size=17,family="serif",color="black",face= "bold"),
        axis.text.x=element_text(size=14,family="serif",color="black",face="bold"), 
        legend.position = "none")
sk
ggsave("result/brca_SK_sankey_train.pdf", width = 6, height = 6, units = "in")
save(sk, file = "result/brca_SK_sankey_train.rda")




df <- label2 %>%
  make_long(`SDCN`, `Consensus Clustering`)
col <- c('#ff7f50', '#1f77b4', '#5fb404', '#db7093', '#37ade7', '#db8f1e', '#48d1cc', '#9400d3', '#dda7c0')
sk <- ggplot(df, aes(x = x, 
                     next_x = next_x, 
                     node = node, 
                     next_node = next_node,
                     fill = factor(node),
                     label = node)) +
  geom_sankey(flow.alpha = 0.5, node.color = 1) +
  geom_sankey_label(size = 5, color = 1, fill = "white", family = 'serif') +
  scale_fill_manual(values = col) +
  labs(y = '', x = '') + 
  theme_sankey(base_size = 17) + 
  theme(text=element_text(size=17,family="serif",color="black",face= "bold"),
        axis.text.x=element_text(size=14,family="serif",color="black",face="bold"), 
        legend.position = "none")
sk
ggsave("result/brca_SC_sankey_train.pdf", width = 6, height = 6, units = "in")
save(sk, file = "result/brca_SC_sankey_train.rda")


## km-sdcn-cc
df <- label %>%
  make_long(`Kmeans Clustering`, `SDCN`, `Consensus Clustering`)
col <- c('#ff7f50', '#1f77b4', '#9400d3', '#5fb404', '#37ade7', '#db8f1e', '#db7093', '#48d1cc', '#dda7c0')
sk <- ggplot(df, aes(x = x, 
                     next_x = next_x, 
                     node = node, 
                     next_node = next_node,
                     fill = factor(node),
                     label = node)) +
  geom_sankey(flow.alpha = 0.5, node.color = 1) +
  geom_sankey_label(size = 5, color = 1, fill = "white", family = 'serif') +
  scale_fill_manual(values = col) +
  labs(y = '', x = '') + 
  theme_sankey(base_size = 17) + 
  theme(text=element_text(size=17,family="serif",color="black",face= "bold"),
        axis.text.x=element_text(size=14,family="serif",color="black",face="bold"), 
        legend.position = "none")
sk
ggsave("result/brca_SKC_sankey_train.pdf", width = 6, height = 6, units = "in")
save(sk, file = "result/brca_SKC_sankey_train.rda")
