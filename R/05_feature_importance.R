# 05_feature_importance.R
# Extract and visualize top features from trained models.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(caret)
  library(tidyr)
})

root_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
fig_dir <- file.path(root_dir, "output", "figures")
tab_dir <- file.path(root_dir, "output", "tables")

models <- readRDS(
  file.path(root_dir, "data", "raw", "FRANZOSA_IBD_2019", "franzosa2019_models.rds")
)

extract_importance <- function(model, model_name) {
  imp <- varImp(model)
  if (is.null(imp$importance)) return(NULL)
  
  imp_df <- imp$importance %>%
    tibble::rownames_to_column("feature")
  
  imp_df <- imp_df %>%
    mutate(importance_mean = rowMeans(select(., -feature), na.rm = TRUE)) %>%
    select(feature, importance_mean) %>%
    arrange(desc(importance_mean)) %>%
    mutate(model = model_name)
  
  imp_df
}

all_importances <- bind_rows(lapply(names(models), function(nm) extract_importance(models[[nm]], nm)))
write_csv(all_importances, file.path(tab_dir, "feature_importance_all_models.csv"))

top20 <- all_importances %>%
  group_by(model) %>%
  slice_max(order_by = importance_mean, n = 20, with_ties = FALSE) %>%
  ungroup()

write_csv(top20, file.path(tab_dir, "feature_importance_top20_per_model.csv"))

p <- ggplot(top20, aes(x = reorder(feature, importance_mean), y = importance_mean, fill = model)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~model, scales = "free_y") +
  theme_minimal(base_size = 11) +
  labs(
    title = "Top 20 Features per Model",
    x = "Feature",
    y = "Importance"
  )

print(p)
ggsave(file.path(fig_dir, "feature_importance_top20_per_model.png"), p, width = 12, height = 10, dpi = 300)

message("Feature importance tables and plot generated!")


#Lieber als Tabelle lesen als Grafik
top20 %>%
  group_by(model) %>%
  slice_max(order_by = importance_mean, n = 10, with_ties = FALSE) %>%
  select(model, feature, importance_mean) %>%
  print(n = 40)

#Overlap Code
overlap <- top20 %>%
  count(feature, sort = TRUE) %>%
  filter(n >= 2)

print(overlap, n = 30)


