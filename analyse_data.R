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

drugs <- c("Carbetocin", "Desmopressin", "Hydroxychloroquine", "Chloroquine", "Bupropion")

complete <- inner_join(df, metadata, by = c("Source_Signature" = "SignatureId", "cellline" = "CellLine")) %>% 
  mutate(perturbagen = str_to_title(Perturbagen)) %>% 
  filter(perturbagen %in% drugs)

filter_data <- function(data, cell_line, cutoff) {
    dataframe <- data
    output <- dataframe %>%
    filter(cellline == cell_line) %>%
    group_by(cellline, treatment, Perturbagen) %>%
    filter(similarity == max(similarity)) %>%
    ungroup() %>% 
    select(signatureid, treatment, Perturbagen, similarity, pValue, cellline)
    return(output)
  }

write_csv(complete, 
          paste("results", paste(paste("complete", "result", sep = "-"), "csv", sep = "."), sep = "/"))

analysed <- complete %>% 
  group_by(cellline, treatment, perturbagen) %>% 
  filter(similarity == max(similarity))

cell_lines <- c("A375", "HA1E", "MCF7", "PC3")

for (cell in cell_lines) {
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

common_cell_lines <- c("MCF7, A375, HA1E")

tnf <- all_results %>% 
  filter(treatment == "TNF") %>% 
  ungroup() %>% 
  select(-treatment) %>% 
  arrange(cellline, similarity)

write_csv(tnf, "results/tnf.csv")

tnfcrosstab <- dcast(tnf, cellline ~ perturbagen)

write.csv(tnfcrosstab, "results/tnfcrosstab.csv")

carbetocin <- all_results %>% 
  filter(perturbagen == "Carbetocin") %>% 
  ungroup() %>% 
  select(-perturbagen) %>% 
  arrange(cellline, similarity)

write_csv(carbetocin, "results/carbetocin.csv")

carbetocincrosstab <- dcast(carbetocin, cellline ~ treatment)

write.csv(carbetocincrosstab, "results/carbetocincrosstab.csv")

