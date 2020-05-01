library(tidyverse)
library(reshape2)
library(gplots)
library(RColorBrewer)

carbetocin <- read_tsv("signature_data/Carbetocin-Signatures.tsv")
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

drugs_signatures <- c(carbetocin$SignatureId, benazepril$SignatureId, captopril$SignatureId, enalapril$SignatureId,
                    fosinopril$SignatureId, lisinopril$SignatureId, moexipril$SignatureId, olmesartan$SignatureId, perindopril$SignatureId,
                    quinapril$SignatureId, ramipril$SignatureId, telmisartan$SignatureId, valsartan$SignatureId)

drugs_list <- list(carbetocin,
                   benazepril, captopril, enalapril, fosinopril, lisinopril, moexipril,
                   olmesartan, perindopril, quinapril, ramipril, telmisartan, valsartan)

metadata <- bind_rows(drugs_list) %>% 
  select(SignatureId, Perturbagen) %>% 
  group_by(SignatureId, Perturbagen) %>% 
  summarise_all(n) %>% 
  ungroup() %>% 
  column_to_rownames("SignatureId")

data <- list()

for (drug in drugs_signatures) {
  prefix <- paste("data", "signatures", sep = "/")
  filename <- paste(paste(drug, "Signature", sep = "-"), "tsv", sep = ".")
  fullname <- paste(prefix, filename, sep = "/")
  if (file.exists(fullname)) {
    print(paste("Now Processing Signature:", drug))
    sig <- read_tsv(fullname)
    data[[drug]] <- sig
  }
}

drugA <- c()
drugA_Name <- c()
drugB <- c()
drugB_Name <- c()
similarity <- c()

for (drug1 in drugs_signatures) {
  for (drug2 in drugs_signatures) {
    d1 <- data[[drug1]] %>% arrange(Name_GeneSymbol)
    d2 <- data[[drug2]] %>% arrange(Name_GeneSymbol)
    c <- cor(d1$Value_LogDiffExp, d2$Value_LogDiffExp)
    drugA <- append(drugA, drug1)
    drugA_Name <- append(drugA_Name, metadata[drug1, "Perturbagen"])
    drugB <- append(drugB, drug2)
    drugB_Name <- append(drugB_Name, metadata[drug2, "Perturbagen"])
    similarity <- append(similarity, c)
  }
}


similarity_df <- tibble(drug_A_id = drugA,
                        drug_A = drugA_Name,
                        drug_B_id = drugB,
                        drug_B = drugB_Name,
                        similarity = similarity)

drugs <- c("Carbetocin", "Benazepril", "Captopril", "Enalapril",
           "Fosinopril", "Lisinopril", "Moexipril", "Olmesartan", "Perindopril", "Quinapril",
           "Ramipril", "Telmisartan", "Valsartan")

filtered <- similarity_df %>% 
  select(drug_A, drug_B, similarity) %>% 
  group_by(drug_A, drug_B) %>% 
  filter(abs(similarity) == max(abs(similarity))) %>% 
  summarise(similarity = mean(similarity))

cross <- filtered %>% 
  filter(drug_A %in% drugs & drug_B %in% drugs) %>% 
  dcast(drug_A ~ drug_B) %>% 
  column_to_rownames("drug_A") %>% 
  as.matrix()

order <- rev(c("Carbetocin", "Benazepril", "Captopril", "Enalapril",
           "Fosinopril", "Lisinopril", "Moexipril", "Olmesartan", "Perindopril", "Quinapril",
           "Ramipril", "Telmisartan", "Valsartan"))

cross <- cross[order,order]

colors <- colorRampPalette(c("red", "black", "green"))(n=11)

png("figures/drug_correlations.png",  width = 1920, height = 1384)
heatmap.2(cross, dendrogram = "none", col = colors, colsep = 1:12, rowsep = 1:12,
          notecol = "white", trace = "none", cexRow = 3, cexCol = 3,
          density.info = "none", keysize = 1, margins = c(15, 25), notecex=3.0, Rowv = FALSE, Colv = FALSE,
          main = "Concordance of the Perturbagens with Immune genes",
          key.title = NA, key.xlab = NA, key.par = list(cex = 1.3))
dev.off()
