library(tidyverse)
library(RColorBrewer)
library(gplots)
library(reshape2)

all_averaged <- read_csv("results/all_averaged.csv")

all_averaged_cross <- all_averaged %>% 
  dcast(perturbagen ~ treatment, value.var = "mean_similarity") %>% 
  column_to_rownames("perturbagen") %>% 
  as.matrix()

colors <- brewer.pal(8, "RdBu")

png(filename = "figures/concordance-heatmap.png", width = 1920, height = 1384)
heatmap(all_averaged_cross[,c(1,4:6,8,9:14)], Colv = NA, Rowv = NA, scale = "column", col = colors)
dev.off()
