# 06_pathway_analysis.R
# Map top metabolite features to KEGG/HMDB, summarize brain-gut relevant pathways.

#Mapping Tabelle

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

root_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

mtb_map <- read_tsv(file.path(root_dir, "data", "raw", "FRANZOSA_IBD_2019", "mtb.map.tsv"), show_col_types = FALSE)

names(mtb_map)
head(mtb_map, 5)
dim(mtb_map)

#Pathway Analyse ausführen

top_features <- read_csv(file.path(root_dir, "output", "tables", "feature_importance_top20_per_model.csv"), show_col_types = FALSE)

mapped <- top_features %>%
  inner_join(mtb_map, by = c("feature" = "Compound"))

message(sprintf("Features gematcht: %d von %d", nrow(mapped), nrow(top_features))) #--> Das zeigt wie viele deiner Top-20-Features in der Mapping-Tabelle gefunden wurden. Dann schauen wir was rauskommt:

mapped %>%
  select(model, feature, Compound.Name, HMDB, KEGG, importance_mean) %>%
  arrange(desc(importance_mean)) %>%
  print(n = 40)


#Brain-Gut-Achse Relevanz prüfen

brain_gut_keywords <- c(
  "butyrate", "propionate", "acetate", "short chain",
  "tryptophan", "kynurenine", "serotonin", "indole",
  "bile acid", "gaba", "dopamine", "glutamate", "glutamic",
  "sphingo", "carnitine", "tyrosine", "phenylalanine",
  "urobilin", "bilirubin"
)

mapped_relevance <- mapped %>%
  mutate(
    text_blob = str_to_lower(paste(
      feature, Compound.Name, HMDB, KEGG, 
      Putative.Chemical.Class, sep = " "
    )),
    brain_gut_flag = str_detect(text_blob, paste(brain_gut_keywords, collapse = "|"))
  ) %>%
  select(-text_blob)

brain_gut_hits <- mapped_relevance %>%
  filter(brain_gut_flag == TRUE) %>%
  select(model, feature, Compound.Name, KEGG, importance_mean)

message(sprintf("Brain-Gut relevante Features: %d", nrow(brain_gut_hits)))
print(brain_gut_hits, n = 30)


#Speichern
tab_dir <- file.path(root_dir, "output", "tables")

write_csv(mapped, file.path(tab_dir, "pathway_mapped_top_metabolites.csv"))
write_csv(mapped_relevance, file.path(tab_dir, "pathway_brain_gut_relevance_flags.csv"))

summary_tbl <- mapped_relevance %>%
  count(model, brain_gut_flag, name = "n_features") %>%
  arrange(model, desc(brain_gut_flag))

write_csv(summary_tbl, file.path(tab_dir, "pathway_brain_gut_summary.csv"))
print(summary_tbl)

message("Pathway analysis completed!")
