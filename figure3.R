library(tidyverse)
library(RColorBrewer)
library(gplots)
library(reshape2)

order <- rev(c("Carbetocin", "Chloroquine", "Hydroxychloroquine", "Lopinavir", "Bupropion"))

all_results <- read_csv("results/all_results.csv")

all_results_cross <- all_results %>% 
  select(-cellline) %>% 
  dcast(perturbagen ~ treatment, value.var = "similarity") %>% 
  filter(perturbagen != "Desmopressin") %>% 
  column_to_rownames("perturbagen") %>% 
  as.matrix()
all_results_cross <- all_results_cross[order,]

all_results_cross_inf <- all_results %>% 
  select(-cellline) %>% 
  dcast(perturbagen ~ treatment, value.var = "similarity") %>%  
  filter(perturbagen != "Desmopressin") %>% 
  select_if(names(.) %in% c("IL1A", "IL1B", "IL1R1", "IL1R2", "IL1RN",
                            "IL1RAP", "IL1RL1", "IL1RL2", "IL6", "IL6R", "IL6ST", "TNF", "perturbagen")) %>% 
  column_to_rownames("perturbagen") %>% 
  as.matrix()
all_results_cross_inf <- all_results_cross_inf[order,]

all_results_cross_imm <- all_results %>% 
  select(-cellline) %>% 
  dcast(perturbagen ~ treatment, value.var = "similarity") %>%  
  filter(perturbagen != "Desmopressin") %>% 
#  select_if(names(.) %in% c("ARG1", "CD19", "TLR9", "CD40", "CD44", "CD46", "CTLA4", "perturbagen")) %>% 
  select_if(names(.) %in% c("ARG1", "TLR9", "CD40", "CD46", "CTLA4", "perturbagen")) %>% 
  column_to_rownames("perturbagen") %>% 
  as.matrix()
all_results_cross_imm <- all_results_cross_imm[order,]

colors <- colorRampPalette(c("red", "black", "green"))(n=11)

png(filename = "figures/max-concordance-heatmap-inflammation.png", width = 1920, height = 1384)
heatmap.2(all_results_cross_inf, dendrogram = "none", col = colors, colsep = 1:3, rowsep = 1:5,
          trace = "none", cexRow = 3, cexCol = 3, Rowv = FALSE, Colv = FALSE,
          density.info = "none", keysize = 1, margins = c(15, 25), notecex=3.0,
          main = "Concordance of the Perturbagens with Inflammation-driving genes",
          key.title = NA, key.xlab = NA, key.par = list(cex = 1.3))
dev.off()

png(filename = "figures/max-concordance-heatmap-immune.png", width = 1920, height = 1384)
heatmap.2(all_results_cross_imm, dendrogram = "none", col = colors, colsep = 1:6, rowsep = 1:5,
          trace = "none", cexRow = 3, cexCol = 3, Rowv = FALSE, Colv = FALSE,
          density.info = "none", keysize = 1, margins = c(15, 25), notecex=3.0,
          main = "Concordance of the Perturbagens with Immune genes",
          key.title = NA, key.xlab = NA, key.par = list(cex = 1.3))
dev.off()

png(filename = "figures/max-concordance-heatmap-all.png", width = 1920, height = 1384)
heatmap.2(all_results_cross, dendrogram = "none", col = colors, colsep = 1:13, rowsep = 1:5,
          trace = "none", cexRow = 3, cexCol = 3, Rowv = FALSE, Colv = FALSE,
          density.info = "none", keysize = 1, margins = c(15, 25), notecex=1.5,
          main = "Concordance of the Perturbagens with all genes",
          key.title = NA, key.xlab = NA, key.par = list(cex = 1.3))
dev.off()


pdf(file = "figures/max-concordance-heatmap-inflammation.pdf", width = 1920, height = 1384)
heatmap.2(all_results_cross_inf, dendrogram = "none", col = colors, colsep = 1:3, rowsep = 1:5,
          cellnote = all_results_cross_inf, notecol = "white", trace = "none", cexRow = 3, cexCol = 3,
          density.info = "none", keysize = 1, margins = c(15, 25), notecex=3.0, Rowv = FALSE, Colv = FALSE,
          main = "Concordance of the Perturbagens with Inflammation-driving genes",
          key.title = NA, key.xlab = NA, key.par = list(cex = 1.3))
dev.off()

pdf(file = "figures/max-concordance-heatmap-immune.pdf", width = 1920, height = 1384)
heatmap.2(all_results_cross_imm, dendrogram = "none", col = colors, colsep = 1:6, rowsep = 1:5,
          cellnote = all_results_cross_imm, notecol = "white", trace = "none", cexRow = 3, cexCol = 3,
          density.info = "none", keysize = 1, margins = c(15, 25), notecex=3.0, Rowv = FALSE, Colv = FALSE,
          main = "Concordance of the Perturbagens with Immune genes",
          key.title = NA, key.xlab = NA, key.par = list(cex = 1.3))
dev.off()

pdf(file = "figures/max-concordance-heatmap-all.pdf", width = 1920, height = 1384)
heatmap.2(all_results_cross_imm, dendrogram = "none", col = colors, colsep = 1:6, rowsep = 1:5,
          cellnote = all_results_cross_imm, notecol = "white", trace = "none", cexRow = 3, cexCol = 3,
          density.info = "none", keysize = 1, margins = c(15, 25), notecex=3.0, Rowv = FALSE, Colv = FALSE,
          main = "Concordance of the Perturbagens with Immune genes",
          key.title = NA, key.xlab = NA, key.par = list(cex = 1.3))
dev.off()
