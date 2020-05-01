library(tidyverse)
source("pipeline.R")

drugs_signatures <- unlist(strsplit(read_file("drugs_signature_ids"), split = "\n"))

for (drug in drugs_signatures) {
  prefix <- paste("data", "signatures", sep = "/")
  filename <- paste(paste(drug, "Signature", sep = "-"), "tsv", sep = ".")
  fullname <- paste(prefix, filename, sep = "/")
  if (!file.exists(fullname)) {
    print(paste("Now Processing Signature:", drug))
    sig <- get_l1000_signature(drug)
    write_tsv(sig, path = fullname)
  } else {
    print(paste(drug, "Signature already exists."))
  }
}

for (drug in drugs_signatures) {
  prefix <- paste("data", "signatures", sep = "/")
  filename <- paste(paste(drug, "Signature", sep = "-"), "tsv", sep = ".")
  fullname <- paste(prefix, filename, sep = "/")
  up <- generate_filtered_signature(signature_file = fullname, direction = "up")
  down <- generate_filtered_signature(signature_file = fullname, direction = "down")
  consensus <- bind_rows(up, down)
  
  consensus_prefix <- paste("data", "filtered", sep = "/")
  consensus_filename <- paste(paste(drug, "Consensus", "Signature", sep = "-"), "tsv", sep = ".")
  consensus_fullname <- paste(consensus_prefix, consensus_filename, sep = "/")

  write_tsv(consensus, consensus_fullname)
  
}

drug_results <- list.files("data/filtered/", pattern = "LINCSCP", full.names = TRUE)

genes <- c("IL1A", "IL1B", "IL2", "IL2RA", "IL2RB", "IL2RG", "IL6", "IL6R",  "TNF",
          "CD8A", "CD8B", "CD4", "CTLA4", "CD19",
          "CD20", "CD3G", "CD11B", "TLR9", "TLR7", "ARG1",
          "CD40", "CD46", "CD44", "CD81", "CD83", "AGT", "AGTR1", "ACE", "ACE2")

for (drug in drug_results) {
  name <- str_split(drug, "/")[[1]][4]
  name <- str_split(name, "-")[[1]][1]
  prefix <- paste("results", "drugs", sep = "/")
  filename <- paste(paste(name, "Concordant", sep = "-"), "tsv", sep = ".")
  fullname <- paste(prefix, filename, sep = "/")
  if (!file.exists(fullname)) {
    print(paste("Processing concordance for signature:", name))
    concordance <- get_concordant_signatures(drug, library = "LIB_6") %>% 
      filter(treatment %in% genes) %>% 
      mutate(Source_Signature = name)
    write_tsv(concordance, fullname)
  } else {
    print(paste(name, "concordant signature already exists. Skipping."))
  }

}