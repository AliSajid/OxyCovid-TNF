library(tidyverse)
library(RColorBrewer)
library(gplots)
library(reshape2)

all_averaged <- read_csv("results/all_averaged.csv")

all_averaged_cross_inf <- all_averaged %>% 
  select(IL1A, IL1B, IL6, IL6R, IL6ST, TNF) %>% 
  dcast(perturbagen ~ treatment, value.var = "mean_similarity") %>% 
  column_to_rownames("perturbagen") %>% 
  as.matrix()

all_averaged_cross_imm <- all_averaged %>% 
  select(ARG1, CD19, TLR7, TLR9, CD40, CD44, CD46, CTLA4) %>% 
  dcast(perturbagen ~ treatment, value.var = "mean_similarity") %>% 
  column_to_rownames("perturbagen") %>% 
  as.matrix()

colors <- rev(brewer.pal(8, "RdBu"))


png(filename = "figures/average-concordance-heatmap-inflammation.png", width = 1920, height = 1384)
heatmap(all_averaged_cross_inf, Colv = NA, Rowv = NA, scale = "column", col = colors)
dev.off()

png(filename = "figures/average-concordance-heatmap-immune.png", width = 1920, height = 1384)
heatmap(all_averaged_cross_imm, Colv = NA, Rowv = NA, scale = "column", col = colors)
dev.off()
