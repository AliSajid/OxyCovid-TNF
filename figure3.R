library(tidyverse)
library(RColorBrewer)
library(gplots)
library(reshape2)

# order <- rev(c("Carbetocin", "Desmopressin", "Hydroxychloroquine", "Chloroquine",
#                "Bupropion", "Ritonavir", "Lopinavir", "Benazepril", "Captopril", "Enalapril",
#                         "Fosinopril", "Lisinopril", "Moexipril","Perindopril", "Quinapril",
#                         "Ramipril", "Telmisartan", "Valsartan", "Olmesartan"))

order <-
  rev(c(
    "Carbetocin",
    "Hydroxychloroquine",
    "Desmopressin",
    #      "Chloroquine",
    "Lopinavir",
    "Bupropion"
  ))

all_all_genes <-  c("IL1A", "IL1B", "IL2", "IL2RA", "IL2RB", "IL2RG", "IL6", "IL6R",  "TNF",
                        "CD8A", "CD8B", "CD4", "CTLA4", "CD19",
                        "CD20", "CD3G", "CD11B", "TLR9", "TLR7", "ARG1",
                        "CD40", "CD46", "CD44", "CD81", "CD83", "AGT", "AGTR1", "ACE", "ACE2",
                        "DPP4", "ANPEP", "CEACAM1", "LAP3", "MMEL1", "CXCL1", "CXCL2", "CXCL3",
                        "CXCL5", "CXCL8", "CCL20", "HEY1", "MUC21", "CCL2", "TGFB1", "TGFBR1", "TGFBR2", "TGFBR3")

inf_genes <- c(
  "IL1A",
  "IL1B",
  "IL1R1",
  "IL1R2",
  "IL1RN",
  "IL1RAP",
  "IL1RL1",
  "IL1RL2",
  "IL6",
  "IL6ST",
  "TNF",
  "IL2",
  "IL2RA",
  "IL2RB",
  "IL2RG"
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

all_results <- read_csv("results/all_results.csv")

all_results_cross <- all_results %>%
  select(-cellline) %>%
  filter(perturbagen %in% order) %>% 
  dcast(perturbagen ~ treatment, value.var = "similarity") %>%
  # filter(perturbagen != "Desmopressin") %>%
  select_if(
    names(.) %in% c(
all_genes,
      "perturbagen"
    )
  ) %>%
  column_to_rownames("perturbagen") %>%
  as.matrix()


order <- order[order %in% rownames(all_results_cross)]
all_results_cross <- all_results_cross[order, ]

all_results_cross %>% as.data.frame() %>% rownames_to_column("Perturbagen") %>% filter(
  Perturbagen %in% c(
    "Carbetocin",
    "Chloroquine",
    "Hydroxychloroquine",
    "Desmopressin",
    "Lopinavir",
    "Bupropion"
  )
) %>%  write_csv("results/all_results_cross.csv")


all_results_cross_inf <- all_results %>%
  select(-cellline) %>%
  filter(perturbagen %in% order) %>% 
  dcast(perturbagen ~ treatment, value.var = "similarity") %>%
  select_if(
    names(.) %in% c(inf_genes,
      "perturbagen"
    )
  ) %>%
  column_to_rownames("perturbagen") %>%
  as.matrix()
all_results_cross_inf <- all_results_cross_inf[order, ]

all_results_cross_imm <- all_results %>%
  select(-cellline) %>%
  filter(perturbagen %in% order) %>% 
  dcast(perturbagen ~ treatment, value.var = "similarity") %>%
  select_if(names(.) %in% c(imm_genes,
                            "perturbagen")) %>%
  column_to_rownames("perturbagen") %>%
  as.matrix()
all_results_cross_imm <- all_results_cross_imm[order,]

colors <- colorRampPalette(c("red", "black", "green"))(n = 100)

png(filename = "figures/max-concordance-heatmap-inflammation.png",
    width = 1920,
    height = 1384)
heatmap.2(
  all_results_cross_inf,
  dendrogram = "none",
  col = colors,
  colsep = 1:7,
  rowsep = 1:5,
  trace = "none",
  cexRow = 3,
  cexCol = 3,
  Rowv = FALSE,
  Colv = FALSE,
  density.info = "none",
  keysize = 1,
  margins = c(15, 25),
  notecex = 3.0,
#  main = "Concordance of the Perturbagens with Inflammation-driving genes",
  key.title = NA,
  key.xlab = NA,
  key.par = list(cex = 1.3),
  cellnote = all_results_cross_inf,
  notecol = "white"
)
dev.off()

png(filename = "figures/max-concordance-heatmap-immune.png",
    width = 1920,
    height = 1384)
heatmap.2(
  all_results_cross_imm,
  dendrogram = "none",
  col = colors,
  colsep = 1:4,
  rowsep = 1:5,
  trace = "none",
  cexRow = 3,
  cexCol = 3,
  Rowv = FALSE,
  Colv = FALSE,
  density.info = "none",
  keysize = 1,
  margins = c(15, 25),
  notecex = 3.0,
  # main = "Concordance of the Perturbagens with Immune genes",
  key.title = NA,
  key.xlab = NA,
  key.par = list(cex = 1.3),
  cellnote = all_results_cross_imm,
  notecol = "white"
)
dev.off()

png(filename = "figures/max-concordance-heatmap-all.png",
    width = 1920,
    height = 1384)
heatmap.2(
  all_results_cross,
  dendrogram = "none",
  col = colors,
  colsep = 1:10,
  rowsep = 1:5,
  trace = "none",
  cexRow = 3,
  cexCol = 3,
  Rowv = FALSE,
  Colv = FALSE,
  density.info = "none",
  keysize = 1,
  margins = c(15, 25),
  notecex = 1.5,
  #          main = "Concordance of the Perturbagens with all genes",
  key.title = NA,
  key.xlab = NA,
  key.par = list(cex = 1.3),
  cellnote = all_results_cross,
  notecol = "white"
)
dev.off()