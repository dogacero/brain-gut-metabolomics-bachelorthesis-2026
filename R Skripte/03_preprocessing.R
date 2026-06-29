# 03_preprocessing.R
# Missing value handling, feature filtering, transformations, and train/test split.

#Pakete + Daten laden

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(caret)
})

root_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
input_rds <- file.path(root_dir, "data", "raw", "FRANZOSA_IBD_2019", "franzosa2019_loaded.rds")
obj <- readRDS(input_rds)

genus       <- obj$genus
metabolites <- obj$metabolites
metadata    <- obj$metadata
sample_col  <- "Sample"
group_col   <- "Study.Group"


#Seltene Genera filtern

genus_x <- genus[, -1, drop = FALSE]
met_x   <- metabolites[, -1, drop = FALSE]

nonzero_prop <- colMeans(genus_x > 0, na.rm = TRUE)
genus_keep <- names(nonzero_prop)[nonzero_prop >= 0.25]
genus_f <- genus_x[, genus_keep, drop = FALSE]

message(sprintf("Genera vor Filterung: %d -> nach Filterung: %d", ncol(genus_x), ncol(genus_f)))


#Imputation

impute_median <- function(df) {
  for (j in seq_along(df)) {
    v <- df[[j]]
    if (anyNA(v)) {
      med <- median(v, na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      v[is.na(v)] <- med
      df[[j]] <- v
    }
  }
  df
}

genus_f <- impute_median(genus_f)
met_x   <- impute_median(met_x)

message(sprintf("NAs in Genera: %d", sum(is.na(genus_f))))
message(sprintf("NAs in Metabolite: %d", sum(is.na(met_x))))


#Transformation auch Normalisierung

genus_t <- asin(sqrt(pmax(as.matrix(genus_f), 0)))
met_t   <- log(as.matrix(met_x) + 1)

#Samples zusammenführen + Labels

rownames(genus_t) <- genus[[sample_col]]
rownames(met_t)   <- metabolites[[sample_col]]

common_samples <- intersect(rownames(genus_t), rownames(met_t))
common_samples <- intersect(common_samples, metadata[[sample_col]])

genus_t <- genus_t[common_samples, , drop = FALSE]
met_t   <- met_t[common_samples, , drop = FALSE]
meta_f  <- metadata %>% filter(.data[[sample_col]] %in% common_samples)

y <- as.factor(meta_f[[group_col]])
names(y) <- meta_f[[sample_col]]

message(sprintf("Samples: %d", length(common_samples)))
print(table(y))

#Kombinieren + Near-Zero-Variance entfernen


combined_x <- cbind(
  as.data.frame(genus_t[meta_f[[sample_col]], , drop = FALSE]),
  as.data.frame(met_t[meta_f[[sample_col]], , drop = FALSE])
)

near_zero <- nearZeroVar(combined_x)
message(sprintf("Near-Zero-Variance Features entfernt: %d", length(near_zero)))

if (length(near_zero) > 0) {
  combined_x <- combined_x[, -near_zero, drop = FALSE]
}

message(sprintf("Finale Feature-Anzahl: %d", ncol(combined_x)))


#Train/Test-Split

set.seed(42)
idx <- createDataPartition(y, p = 0.8, list = FALSE)

train_data <- combined_x[idx, , drop = FALSE]
test_data  <- combined_x[-idx, , drop = FALSE]
train_y    <- y[idx]
test_y     <- y[-idx]

message(sprintf("Training: %d Samples | Test: %d Samples", nrow(train_data), nrow(test_data)))
message("Training-Verteilung:")
print(table(train_y))
message("Test-Verteilung:")
print(table(test_y))

#Speichern

out <- list(
  train_data = train_data,
  test_data  = test_data,
  train_y    = train_y,
  test_y     = test_y,
  metadata   = meta_f,
  group_col  = group_col,
  sample_id_col = sample_col
)

saveRDS(out, file.path(root_dir, "data", "raw", "FRANZOSA_IBD_2019", "franzosa2019_preprocessed.rds"))

write_csv(
  tibble(
    metric = c("n_train", "n_test", "n_features"),
    value  = c(nrow(train_data), nrow(test_data), ncol(train_data))
  ),
  file.path(root_dir, "output", "tables", "preprocessing_summary.csv")
)

message("Preprocessing completed and saved!")

