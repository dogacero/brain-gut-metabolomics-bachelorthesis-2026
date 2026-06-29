# 00_setup.R
# Install and load required packages for the thesis pipeline.

required_cran <- c(
  "caret",
  "randomForest",
  "glmnet",
  "e1071",
  "pROC",
  "ggplot2",
  "pheatmap",
  "FactoMineR",
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "tibble",
  "kernlab"
)

install_if_missing <- function(pkgs) {
  missing_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(missing_pkgs) > 0) {
    install.packages(missing_pkgs, dependencies = TRUE)
  }
}

install_if_missing(required_cran)

invisible(lapply(required_cran, require, character.only = TRUE))

message("Setup finished. All required packages are installed and loaded.")
