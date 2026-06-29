# 01_load_data.R
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

dataset_name <- "FRANZOSA_IBD_2019"
base_url <- "https://raw.githubusercontent.com/borenstein-lab/microbiome-metabolome-curated-data/main/data/processed_data"

root_dir <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE)
if (!dir.exists(file.path(root_dir, "R"))) {
  root_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

data_dir <- file.path(root_dir, "data", "raw", dataset_name)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

files_to_download <- c(
  "genera.tsv",
  "mtb.tsv",
  "metadata.tsv",
  "mtb.map.tsv"
)

download_one <- function(file_name) {
  src <- sprintf("%s/%s/%s", base_url, dataset_name, file_name)
  dest <- file.path(data_dir, file_name)
  if (!file.exists(dest)) {
    download.file(src, dest, mode = "wb", quiet = FALSE)
  }
  dest
}

local_paths <- vapply(files_to_download, download_one, FUN.VALUE = character(1))

genus <- read_tsv(local_paths[1], show_col_types = FALSE)
metabolites <- read_tsv(local_paths[2], show_col_types = FALSE)
metadata <- read_tsv(local_paths[3], show_col_types = FALSE)
mtb_map <- read_tsv(local_paths[4], show_col_types = FALSE)

sample_col <- "Sample"

common_samples <- Reduce(intersect, list(
  genus[[sample_col]],
  metabolites[[sample_col]],
  metadata[[sample_col]]
))

genus <- genus %>% filter(.data[[sample_col]] %in% common_samples)
metabolites <- metabolites %>% filter(.data[[sample_col]] %in% common_samples)
metadata <- metadata %>% filter(.data[[sample_col]] %in% common_samples)

message("Dataset loaded and harmonized.")
message(sprintf("Genus table: %s samples x %s features", nrow(genus), ncol(genus) - 1))
message(sprintf("Metabolite table: %s samples x %s features", nrow(metabolites), ncol(metabolites) - 1))
message(sprintf("Metadata table: %s samples x %s variables", nrow(metadata), ncol(metadata) - 1))

group_col <- "Study.Group"
if (!is.na(group_col)) {
  print(table(metadata[[group_col]], useNA = "ifany"))
}

output_rds <- file.path(root_dir, "data", "raw", dataset_name, "franzosa2019_loaded.rds")
saveRDS(
  list(
    genus = genus,
    metabolites = metabolites,
    metadata = metadata,
    mtb_map = mtb_map,
    sample_col = sample_col,
    group_col = group_col
  ),
  output_rds
)

message(sprintf("Saved merged object to: %s", output_rds))
source("R/01_load_data.R")

