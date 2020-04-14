library(tidyverse)

carbetocin <- read_tsv("signature_data/Carbetocin-Signatures.tsv")

desmopressin <- read_tsv("signature_data/Desmopressin-Signatures.tsv")

hydroxychloroquine <- read_tsv("signature_data/Hydroxychloroquine-Signatures.tsv")

chloroquine <- read_tsv("signature_data/Chloroquine-Signatures.tsv")

bupropion <- read_tsv("signature_data/Bupropion-Signatures.tsv")

lopinavir <- read_tsv("signature_data/Lopinavir-Signatures.tsv")

ritonavir <- read_tsv("signature_data/Ritonavir-Signatures.tsv")

drugs_list <- list(bupropion, carbetocin, chloroquine, desmopressin, lopinavir, ritonavir)
selected_lines <- list()

for (i in 1:length(drugs_list)) {
  selected_lines[[i]] <- drugs_list[[i]]$CellLine
}

cell_lines <- reduce(selected_lines, intersect)

IL6 <- read_tsv("signature_data/IL6-Signatures.tsv")

IL6R <- read_tsv("signature_data/IL6-Signatures.tsv")

IL6ST <- read_tsv("signature_data/IL6-Signatures.tsv")

drugs <- bind_rows(carbetocin, desmopressin, hydroxychloroquine, chloroquine, bupropion) %>% 
  filter(CellLine %in% cell_lines)

write_file(paste(drugs$SignatureId, collapse = "\n"), "drugs_signature_ids")

targets <- bind_rows(IL6, IL6R, IL6ST)  %>% 
  filter(CellLine %in% cell_lines)

write_file(paste(targets$SignatureId, collapse = "\n"), "targets_signature_ids")

drugs %>% select(SignatureId, Perturbagen, CellLine) %>% 
  write_csv("signature_data/id-name-cellline_mapping.csv")
