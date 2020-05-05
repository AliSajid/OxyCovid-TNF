library(tidyverse)
source("pipeline.R")

drugs <- c("Carbetocin", "Losartan", "Candesartan")

genes <- c("MAS1", "TNF", "AGT")

mapping <- read_csv("signature_data/id-name-cellline_mapping.csv") %>% 
  filter(Perturbagen %in% drugs)

return_data <- function(signature_id) {
  signature <- get_l1000_signature(signature_id)
  filter_up <- generate_filtered_signature(signature, direction = "up")
  filter_dn <- generate_filtered_signature(signature, direction = "down")
  filter_con <- bind_rows(filter_up, filter_dn)
  
  results_kd <- get_concordant_signatures(signature_df = filter_con, library = "LIB_6") %>% 
    filter(treatment %in% genes) %>% 
    mutate(Source_Signature = signature_id)
  results_oe <- get_concordant_signatures(signature_df = filter_con, library = "LIB_11") %>% 
    filter(treatment %in% genes) %>% 
    mutate(Source_Signature = signature_id)
  
  output <- list(
    knockdown = results_kd,
    overexpress = results_oe
  )
  
  return(output)
}

sigs <- lapply(mapping$SignatureId, return_data)

saveRDS(sigs, "ACE2Data.Rds")
comp <- reduce(sigs, bind_rows)

complete <- inner_join(comp, mapping, by = c("Source_Signature" = "SignatureId")) %>% 
  mutate(Class = if_else(str_detect(signatureid, "LINCSKD"), "Knockdown", "Overexpress")) %>% 
  select(treatment, Perturbagen, Class, similarity)

summary_data <- complete %>% 
  group_by(Perturbagen, treatment, Class) %>% 
  filter(abs(similarity) == max(abs(similarity))) %>% 
  rename()

oex <- summary_data %>% 
  filter(Class == "Overexpress")

p <- ggplot(oex, aes(Perturbagen, similarity, fill = treatment, width = 0.5))
p + geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_x_discrete(limits = drugs) +
  scale_y_continuous(limits = c(-1,1)) +
  theme_minimal() +
  labs(
    title = "Comparison of Concordance for MAS1 gene",
    x = "Drug",
    y = "Concordance"
  ) +
  guides(fill = guide_legend(title = "Gene"))
