library(tidyverse)
library(RColorBrewer)
library(reshape2)

all_results <- read_csv("results/all_results.csv")

all_results_cross <- all_results %>% 
  select(-cellline) %>% 
  dcast(perturbagen ~ treatment, value.var = "similarity") %>% 
  column_to_rownames("perturbagen") %>% 
  as.matrix()

all_results_cross_inf <- all_results %>% 
  select(-cellline) %>% 
  dcast(perturbagen ~ treatment, value.var = "similarity") %>% 
  select_if(names(.) %in% c("IL1A", "IL1B", "IL1R1", "IL1R2", "IL1RN",
                            "IL1RAP", "IL1RL1", "IL1RL2", "IL6", "IL6R", "IL6ST", "TNF", "perturbagen")) %>% 
  column_to_rownames("perturbagen") %>% 
  as.matrix()

all_results_cross_imm <- all_results %>% 
  select(-cellline) %>% 
  dcast(perturbagen ~ treatment, value.var = "similarity") %>% 
  select_if(names(.) %in% c("ARG1", "CD19", "TLR7", "TLR9", "CD40", "CD44", "CD46", "CTLA4", "perturbagen")) %>% 
  column_to_rownames("perturbagen") %>% 
  as.matrix()

colors <- colorRampPalette(c("red", "black", "green"))(n=11)

png(filename = "figures/max-concordance-heatmap-inflammation.png", width = 1920, height = 1384)
heatmap(all_results_cross_inf, Colv = NA, Rowv = NA, scale = "column", col = colors)
dev.off()

png(filename = "figures/max-concordance-heatmap-immune.png", width = 1920, height = 1384)
heatmap(all_results_cross_imm, Colv = NA, Rowv = NA, scale = "column", col = colors)
dev.off()

png(filename = "figures/max-concordance-heatmap-all.png", width = 1920, height = 1384)
heatmap(all_results_cross, Colv = NA, Rowv = NA, scale = "column", col = colors)
dev.off()