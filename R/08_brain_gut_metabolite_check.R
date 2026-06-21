library(readr)
library(dplyr)
library(stringr)

root_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

mtb_map <- read_tsv(
  file.path(root_dir, "data", "raw", "FRANZOSA_IBD_2019", "mtb.map.tsv"),
  show_col_types = FALSE
)

all_importance <- read_csv(
  file.path(root_dir, "output", "tables", "feature_importance_all_models.csv"),
  show_col_types = FALSE
)

brain_gut_terms <- c(
  "butyrate", "butyric", "propionate", "propionic", "acetate", "acetic",
  "tryptophan", "kynurenine", "kynurenic", "serotonin", "indole",
  "bile acid", "cholic", "deoxycholic", "lithocholic", "ursodeoxycholic",
  "gaba", "gamma-aminobutyric",
  "dopamine", "tyrosine", "phenylalanine",
  "glutamate", "glutamic", "glutamine",
  "sphingosine", "sphingomyelin", "ceramide",
  "carnitine",
  "urobilin", "bilirubin",
  "trimethylamine", "TMAO",
  "histamine", "histidine"
)

pattern <- paste(brain_gut_terms, collapse = "|")

bg_in_map <- mtb_map %>%
  filter(str_detect(str_to_lower(paste(Compound, Compound.Name, Putative.Chemical.Class)), pattern)) %>%
  select(Compound, Compound.Name, Putative.Chemical.Class, HMDB, KEGG)

message(sprintf("Brain-Gut-relevante Metaboliten im Datensatz: %d", nrow(bg_in_map)))
print(bg_in_map, n = 50)

bg_importance <- bg_in_map %>%
  inner_join(all_importance, by = c("Compound" = "feature")) %>%
  select(Compound, Compound.Name, model, importance_mean) %>%
  arrange(desc(importance_mean))

message(sprintf("Davon mit Feature-Importance-Wert: %d", nrow(bg_importance)))
print(bg_importance, n = 60)

bg_summary <- bg_importance %>%
  group_by(Compound, Compound.Name) %>%
  summarise(
    max_importance = max(importance_mean, na.rm = TRUE),
    n_models = n_distinct(model),
    .groups = "drop"
  ) %>%
  arrange(desc(max_importance))

message("Brain-Gut Metaboliten - Ranking:")
print(bg_summary, n = 40)

tab_dir <- file.path(root_dir, "output", "tables")
write_csv(bg_in_map, file.path(tab_dir, "brain_gut_metabolites_in_dataset.csv"))
write_csv(bg_summary, file.path(tab_dir, "brain_gut_metabolites_importance_ranking.csv"))

message("Fertig! Ergebnisse gespeichert.")

