library(tidyverse)
library(httr)
library(jsonlite)
source("pipeline.R")

drugs_signatures <- unlist(strsplit(read_file("drugs_signature_ids"), split = "\n"))
targets_signatures <- unlist(strsplit(read_file("targets_signature_ids"), split = "\n"))

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

for (target in targets_signatures) {
  prefix <- paste("data", "signatures", sep = "/")
  filename <- paste(paste(target, "Signature", sep = "-"), "tsv", sep = ".")
  fullname <- paste(prefix, filename, sep = "/")
  if (!file.exists(fullname)) {
    print(paste("Now Processing Signature:", drug))
    sig <- get_l1000_signature(target)
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

# for (target in targets_signatures) {
#   prefix <- paste("data", "signatures", sep = "/")
#   filename <- paste(paste(target, "Signature", sep = "-"), "tsv", sep = ".")
#   fullname <- paste(prefix, filename, sep = "/")
#   up <- generate_filtered_signature(signature_file = fullname, direction = "up")
#   down <- generate_filtered_signature(signature_file = fullname, direction = "down")
#   consensus <- bind_rows(up, down)
#   
#   consensus_prefix <- paste("data", "filtered", sep = "/")
#   consensus_filename <- paste(paste(target, "Consensus", "Signature", sep = "-"), "tsv", sep = ".")
#   consensus_fullname <- paste(consensus_prefix, consensus_filename, sep = "/")
#   
#   write_tsv(consensus, consensus_fullname)
# }

drug_results <- list.files("data/filtered/", pattern = "LINCSCP", full.names = TRUE)

genes <- c("IL6", "IL6R", "IL6ST", "TNF", "TNF", "TNFAIP1", "TNFAIP3", "TNFRSF10A",
          "TNFRSF10B", "TNFRSF11A", "TNFRSF12A", "TNFRSF13C", "TNFRSF14", "TNFRSF18",
          "TNFRSF19", "TNFRSF1A", "TNFRSF21", "TNFRSF6B", "TNFSF10", "TNFSF13", "TNFSF13B",
          "TNFSF15", "TNFSF4", "TNFSF8", "TNFRSF13B", "TNFRSF17", "C1QTNF6")

for (drug in drug_results) {
  name <- str_split(drug, "/")[[1]][4]
  name <- str_split(name, "-")[[1]][1]
  print(paste("Processing concordance for signature:", name))
  concordance <- get_concordant_signatures(drug, library = "LIB_6") %>% 
  filter(treatment %in% genes) %>% 
  mutate(Source_Signature = name)
  prefix <- paste("results", "drugs", sep = "/")
  filename <- paste(paste(name, "Concordant", sep = "-"), "tsv", sep = ".")
  fullname <- paste(prefix, filename, sep = "/")
  write_tsv(concordance, fullname)
}

# drugs <- c("Fluoxetine", "Paroxetine", "Bupropion")
# 
# for (target in target_results) {
#   name <- str_split(target, "/")[[1]][4]
#   concordance <- get_concordant_signatures(target) %>% 
#     filter(compound %in% drugs)
#   prefix <- paste("results", "targets", sep = "/")
#   filename <- paste(paste(name, "Concordant", sep = "-"), "tsv", sep = ".")
#   fullname <- paste(prefix, filename, sep = "/")
#   write_tsv(concordance, fullname)
# }