library(tidyverse)
library(RColorBrewer)
library(gplots)
library(reshape2)

order <-
  rev(c(
    "Carbetocin",
    "Hydroxychloroquine",
    "Desmopressin",
    "Lopinavir",
    "Bupropion"
  ))

inf_genes <- c(
  "IL6",
  "TNF",
  "IL2",
  "IL2RA",
  "IL2RB"
)

imm_genes <- c("ARG1",
               "TLR9",
               "CD40",
               "CD46",
               "CEACAM1",
               "CD83",
               "CCL20",
               "TGFB1",
               "TGFBR1",
               "TGFBR2")

all_genes <- c(inf_genes, imm_genes)

all_results_cross <- read_csv("results/all_results_cross.csv") %>% 
  column_to_rownames("Perturbagen") %>% 
  as.matrix() %>% 
  t()

all_results_cross <- all_results_cross[all_genes,]
side_color <- c(rep("Blue", length(inf_genes)), rep("Red", length(imm_genes)))

colors <- colorRampPalette(c("red", "black", "green"))(n = 100)

png(filename = "figures/max-concordance-heatmap-all-annotated.png",
    width = 1920,
    height = 1384)
heatmap.2(
  all_results_cross,
  dendrogram = "none",
  col = colors,
  colsep = 1:5,
  rowsep = 1:15,
  trace = "none",
  cexRow = 3,
  cexCol = 3,
  Rowv = FALSE,
  Colv = FALSE,
  density.info = "none",
  keysize = 1,
  margins = c(15, 25),
  notecex = 1.5,
  key.title = NA,
  key.xlab = NA,
  key.par = list(cex = 1.3),
  cellnote = all_results_cross,
  notecol = "white",
  RowSideColors = side_color,
  srtCol = 45,
  labCol = c("BUP", "LOP", "DES", "HCQ", "CAB")
)
dev.off()
