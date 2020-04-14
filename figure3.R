library(tidyverse)
library(RColorBrewer)
library(gplots)
library(reshape2)

all_results <- read_csv("results/all_results.csv")

all_results_cross_inf <- all_results %>% 
  select(IL1, IL6, IL6R, IL6ST, TNF, -cellline) %>% 
  dcast(perturbagen ~ treatment, value.var = "similarity") %>% 
  column_to_rownames("perturbagen") %>% 
  as.matrix()

all_results_cross_imm <- all_results %>% 
  select(ARG1, CD19, TLR7, TLR9, CD40, CD44, CD46, CTLA4, -cellline) %>% 
  dcast(perturbagen ~ treatment, value.var = "similarity") %>% 
  column_to_rownames("perturbagen") %>% 
  as.matrix()

colors <- rev(brewer.pal(8, "RdBu"))


png(filename = "figures/max-concordance-heatmap-inflammation.png", width = 1920, height = 1384)
heatmap(all_results_cross_inf, Colv = NA, Rowv = NA, scale = "column", col = colors)
dev.off()

png(filename = "figures/max-concordance-heatmap-immune.png", width = 1920, height = 1384)
heatmap(all_results_cross_imm, Colv = NA, Rowv = NA, scale = "column", col = colors)
dev.off()

