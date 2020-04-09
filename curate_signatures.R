library(tidyverse)

cell_lines <- c("A375", "HA1E", "MCF7", "PC3")

carbetocin <- read_tsv("signature_data/Carbetocin-Signature.tsv") %>% 
  filter(CellLine %in% cell_lines)

desmopressin <- read_tsv("signature_data/Desmopressin-Signature.tsv") %>% 
  filter(CellLine %in% cell_lines)

hydroxychloroquine <- read_tsv("signature_data/Hydroxychloroquine-Signatures.tsv") %>% 
  filter(CellLine %in% cell_lines)

chloroquine <- read_tsv("signature_data/Chloroquine-Signatures.tsv") %>% 
  filter(CellLine %in% cell_lines)

bupropion <- read_tsv("signature_data/Bupropion-Signatures.tsv") %>% 
  filter(CellLine %in% cell_lines)

IL6 <- read_tsv("signature_data/IL6-Signatures.tsv") %>% 
  filter(TargetGene == "IL6", CellLine %in% cell_lines)

IL6R <- read_tsv("signature_data/IL6-Signatures.tsv") %>% 
  filter(TargetGene == "IL6R", CellLine %in% cell_lines)

IL6ST <- read_tsv("signature_data/IL6-Signatures.tsv") %>% 
  filter(TargetGene == "IL6ST", CellLine %in% cell_lines)

drugs <- bind_rows(carbetocin, desmopressin, hydroxychloroquine, chloroquine, bupropion)
write_file(paste(drugs$SignatureId, collapse = "\n"), "drugs_signature_ids")
targets <- bind_rows(IL6, IL6R, IL6ST)
write_file(paste(targets$SignatureId, collapse = "\n"), "targets_signature_ids")
drugs %>% select(SignatureId, Perturbagen, CellLine) %>% 
  write_csv("signature_data/id-name-cellline_mapping.csv")
