library(tidyverse)

carbetocin <- read_tsv("signature_data/Carbetocin-Signatures.tsv")

desmopressin <- read_tsv("signature_data/Desmopressin-Signatures.tsv")

hydroxychloroquine <- read_tsv("signature_data/Hydroxychloroquine-Signatures.tsv")

chloroquine <- read_tsv("signature_data/Chloroquine-Signatures.tsv")

bupropion <- read_tsv("signature_data/Bupropion-Signatures.tsv")

lopinavir <- read_tsv("signature_data/Lopinavir-Signatures.tsv")

ritonavir <- read_tsv("signature_data/Ritonavir-Signatures.tsv")

benazepril <- read_tsv("signature_data/Benazepril-Signatures.tsv")

captopril <- read_tsv("signature_data/Captopril-Signatures.tsv")

enalapril <- read_tsv("signature_data/Enalapril-Signatures.tsv")

fosinopril <- read_tsv("signature_data/Fosinopril-Signatures.tsv")

lisinopril <- read_tsv("signature_data/Lisinopril-Signatures.tsv")

moexipril <- read_tsv("signature_data/Moexipril-Signatures.tsv")

olmesartan <- read_tsv("signature_data/Olmesartan-Signatures.tsv")

perindopril <- read_tsv("signature_data/Perindopril-Signatures.tsv")

quinapril <- read_tsv("signature_data/Quinapril-Signatures.tsv")

ramipril <- read_tsv("signature_data/Ramipril-Signatures.tsv")

telmisartan <- read_tsv("signature_data/Telmisartan-Signatures.tsv")

valsartan <- read_tsv("signature_data/Valsartan-Signatures.tsv")

drugs_list <- list(bupropion, carbetocin, chloroquine, desmopressin, lopinavir, ritonavir,
                   benazepril, captopril, enalapril, fosinopril, lisinopril, moexipril,
                   olmesartan, perindopril, quinapril, ramipril, telmisartan, valsartan, hydroxychloroquine)
selected_lines <- list()

for (i in 1:length(drugs_list)) {
  selected_lines[[i]] <- drugs_list[[i]]$CellLine
}

cell_lines <- reduce(selected_lines, union)

drugs <- bind_rows(bupropion, carbetocin, chloroquine, desmopressin, lopinavir, ritonavir,
                   benazepril, captopril, enalapril, fosinopril, lisinopril, moexipril,
                   olmesartan, perindopril, quinapril, ramipril, telmisartan, valsartan, hydroxychloroquine) %>% 
  filter(CellLine %in% cell_lines)

write_file(paste(drugs$SignatureId, collapse = "\n"), "drugs_signature_ids")

drugs %>% select(SignatureId, Perturbagen, CellLine) %>% 
  write_csv("signature_data/id-name-cellline_mapping.csv")
