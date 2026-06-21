# 07_identifizierung_unbekannte_metaboliten.R
# Document top metabolite features without compound annotation.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

root_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
tab_dir <- file.path(root_dir, "output", "tables")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

top_features <- read_csv(
  file.path(tab_dir, "feature_importance_top20_per_model.csv"),
  show_col_types = FALSE
)

mtb_map <- read_tsv(
  file.path(root_dir, "data", "raw", "FRANZOSA_IBD_2019", "mtb.map.tsv"),
  show_col_types = FALSE
)

top_met <- top_features %>%
  filter(str_detect(feature, "^(C18|HILIC|C8)"))

unk <- top_met %>%
  inner_join(mtb_map, by = c("feature" = "Compound")) %>%
  filter(
    is.na(Compound.Name) | Compound.Name == "" | Compound.Name == "NA"
  ) %>%
  select(
    feature, model, importance_mean,
    m.z, Retention.Time, Adduct,
    Putative.Chemical.Class, HMDB, KEGG, Compound.Name
  ) %>%
  arrange(desc(importance_mean))

unk_unique <- unk %>%
  group_by(feature) %>%
  slice_max(order_by = importance_mean, n = 1, with_ties = FALSE) %>%
  ungroup()

write_csv(unk_unique, file.path(tab_dir, "unknown_top_metabolites_for_review.csv"))
write_csv(unk, file.path(tab_dir, "unknown_top_metabolites_all_models.csv"))

message(sprintf("Unknown top metabolites (unique): %d", nrow(unk_unique)))
message("Saved unknown_top_metabolites_for_review.csv and unknown_top_metabolites_all_models.csv")
