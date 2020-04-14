library(tidyverse)
library(reshape2)

drugs_signatures <- unlist(strsplit(read_file("drugs_signature_ids"), split = "\n"))

prefix <- paste("results", "drugs", sep = "/")
filenames <- paste(paste(drugs_signatures, "Concordant", sep = "-"), "tsv", sep = ".")

files <- paste(prefix, filenames, sep = "/")

metadata <- read_csv("signature_data/id-name-cellline_mapping.csv",
                     col_types = cols(
                       SignatureId = col_character(),
                       Perturbagen = col_character(),
                       CellLine = col_character()
                     ))

col_spec <- cols(
  similarity = col_double(),
  pValue = col_double(),
  nGenes = col_double(),
  treatment = col_character(),
  perturbagenID = col_character(),
  time = col_character(),
  signatureid = col_character(),
  cellline = col_character(),
  Source_Signature = col_character()
)

dfs <- list()

for (i in 1:length(files)) {
  df <- read_tsv(files[i], col_types = col_spec)
  dfs[[i]] <- df
}

df <- reduce(dfs, bind_rows)

drugs <- c("Carbetocin", "Desmopressin", "Hydroxychloroquine", "Chloroquine",
           "Bupropion", "Ritonavit", "Lopinavir")

complete <- inner_join(df, metadata, by = c("Source_Signature" = "SignatureId")) %>% 
  mutate(perturbagen = str_to_title(Perturbagen)) %>% 
  filter(perturbagen %in% drugs)

filter_data <- function(data) {
    dataframe <- data
    output <- dataframe %>%
    group_by(treatment, perturbagen) %>%
    filter(abs(similarity) == max(abs(similarity))) %>%
    ungroup() %>% 
    select(signatureid, treatment, perturbagen, similarity, pValue, cellline)
    return(output)
  }

write_csv(complete, 
          paste("results", paste(paste("complete", "result", sep = "-"), "csv", sep = "."), sep = "/"))

analysed <- filter_data(complete)

#cell_lines <- c("A375", "HA1E", "MCF7", "PC3")

for (cell in unique(analysed$cellline)) {
  outfile <- paste("results", paste(paste(cell, "result", sep = "-"), "csv", sep = "."), sep = "/")

  analysed %>%
    filter(cellline == cell) %>%
    select(perturbagen, treatment, cellline, similarity) %>%
    write_csv(outfile)
}

result_files <- list.files("results/", pattern = "result")

all_results <- analysed %>% 
  select(perturbagen, treatment, cellline, similarity)

write_csv(all_results, "results/all_results.csv")

all_averaged <- all_results %>% 
  group_by(perturbagen, treatment) %>% 
  summarise(mean_similarity = mean(similarity))

write_csv(all_averaged, "results/all_averaged.csv")

process_gene <- function(dataset, gene) {
  g <- dataset %>% 
    filter(treatment == gene) %>% 
    ungroup() %>% 
    select(-treatment) %>% 
    arrange(cellline, similarity)
  
  gcross <- g %>% 
    dcast(cellline ~ perturbagen)
  
  prefix <- "results/"
  file <- paste(gene, "csv", sep = ".")
  crossfile <- paste(paste(gene, "crosstab", sep = "-"), "csv", sep = ".")
  
  write_csv(g, paste(prefix, file, sep = "/"))
  write_csv(gcross, paste(prefix, crossfile, sep = "/"))
  invisible(list(g, gcross))
}

process_gene(all_results, "TNF") # Selected and subset to HA1E
process_gene(all_results, "TLR7") # Selecting HA1E
process_gene(all_results, "TLR9") # Selecting HA1E
process_gene(all_results, "ARG1") # Selecting HA1E
process_gene(all_results, "CD40") # Does not have Carbetocin in result
process_gene(all_results, "CD46") # Does not have a direct comparison with Carbetocin
process_gene(all_results, "CD83") # Selecting HA1E
process_gene(all_results, "CD44") # Does not have a direct comparison with Carbetocin

carbetocin <- all_results %>% 
  filter(perturbagen == "Carbetocin") %>% 
  ungroup() %>% 
  select(-perturbagen, -cellline) %>% 
  arrange(similarity)

write_csv(carbetocin, "results/carbetocin.csv")