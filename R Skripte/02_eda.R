# 02_eda.R
# Explorative analysis: diversity, PCA, heatmap, boxplots, and univariate tests.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(pheatmap)
  library(FactoMineR)
  library(factoextra)
  library(tidyr)
})



root_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

input_rds <- file.path(root_dir, "data", "raw", "FRANZOSA_IBD_2019", "franzosa2019_loaded.rds")
obj <- readRDS(input_rds)

genus       <- obj$genus
metabolites <- obj$metabolites
metadata    <- obj$metadata
sample_col  <- "Sample"
group_col   <- "Study.Group"

#nrow(genus) number of rows also die proben
#ncol(genus) anzahl spalten features+1
#dim(genus) beides auf einmal
#head(genus[,1:3]) die ersten paar zeilen und spalten anschauen

#ordner anlegen

fig_dir <- file.path(root_dir, "output", "figures")
tab_dir <- file.path(root_dir, "output", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)



#shannon diversität

calc_shannon <- function(x) {
  x <- as.numeric(x)
  x <- x[x > 0]
  if (length(x) == 0) return(NA_real_)
  p <- x / sum(x)
  -sum(p * log(p))
}

genus_matrix <- as.matrix(genus[, -1])
rownames(genus_matrix) <- genus[[sample_col]]
shannon <- apply(genus_matrix, 1, calc_shannon)

shannon_df <- metadata %>%
  mutate(Shannon = shannon[match(.data[[sample_col]], names(shannon))])

p_shannon <- ggplot(shannon_df, aes(x = .data[[group_col]], y = Shannon, fill = .data[[group_col]])) +
  geom_boxplot(alpha = 0.8) +
  theme_minimal(base_size = 12) +
  guides(fill = "none") +
  labs(title = "Shannon Diversity by Study Group", x = "Group", y = "Shannon Index")

print(p_shannon)
ggsave(file.path(fig_dir, "eda_shannon_diversity.png"), p_shannon, width = 8, height = 5, dpi = 300)

#PCA auf Genus Ebene

genus_log <- log10(genus_matrix + 1e-6)
pca_genus <- PCA(genus_log, graph = FALSE, scale.unit = TRUE)

pca_genus_df <- data.frame(
  SampleID = rownames(pca_genus$ind$coord),
  PC1 = pca_genus$ind$coord[, 1],
  PC2 = pca_genus$ind$coord[, 2]
) %>%
  left_join(metadata, by = c("SampleID" = sample_col))

p_genus <- ggplot(pca_genus_df, aes(x = PC1, y = PC2, color = .data[[group_col]])) +
  geom_point(alpha = 0.8, size = 2.5) +
  theme_minimal(base_size = 12) +
  labs(title = "PCA - Genus Level", color = "Group")

print(p_genus)
ggsave(file.path(fig_dir, "eda_pca_genus.png"), p_genus, width = 8, height = 5, dpi = 300)

#PCA auf Metabolit Ebene

met_matrix <- as.matrix(metabolites[, -1])
rownames(met_matrix) <- metabolites[[sample_col]]
met_log <- log10(met_matrix + 1)

pca_met <- PCA(met_log, graph = FALSE, scale.unit = TRUE)

pca_met_df <- data.frame(
  SampleID = rownames(pca_met$ind$coord),
  PC1 = pca_met$ind$coord[, 1],
  PC2 = pca_met$ind$coord[, 2]
) %>%
  left_join(metadata, by = c("SampleID" = sample_col))

p_met <- ggplot(pca_met_df, aes(x = PC1, y = PC2, color = .data[[group_col]])) +
  geom_point(alpha = 0.8, size = 2.5) +
  theme_minimal(base_size = 12) +
  labs(title = "PCA - Metabolite Level", color = "Group")

print(p_met)
ggsave(file.path(fig_dir, "eda_pca_metabolites.png"), p_met, width = 8, height = 5, dpi = 300)



#Heatmap der Top-30 Genera

genus_var <- apply(genus_log, 2, var, na.rm = TRUE)
top_genera <- names(sort(genus_var, decreasing = TRUE))[1:min(30, length(genus_var))]
heat_mat <- genus_log[, top_genera, drop = FALSE]

ann <- data.frame(Group = metadata[[group_col]])
rownames(ann) <- metadata[[sample_col]]
ann_for_heatmap <- ann[rownames(heat_mat), , drop = FALSE]

pheatmap(
  heat_mat,
  scale = "row",
  show_rownames = FALSE,
  annotation_row = ann_for_heatmap,
  filename = file.path(fig_dir, "eda_heatmap_top30_genera.png"),
  width = 8,
  height = 10
)

message("Heatmap gespeichert!")

#Boxplots Top Metabolite

met_var <- apply(met_log, 2, var, na.rm = TRUE)
top_mets <- names(sort(met_var, decreasing = TRUE))[1:min(6, length(met_var))]

box_df <- metabolites %>%
  select(all_of(c(sample_col, top_mets))) %>%
  pivot_longer(cols = -all_of(sample_col), names_to = "Metabolite", values_to = "Value") %>%
  left_join(metadata, by = setNames(sample_col, sample_col))

p_box <- ggplot(box_df, aes(x = .data[[group_col]], y = log10(Value + 1), fill = .data[[group_col]])) +
  geom_boxplot(outlier.alpha = 0.3) +
  facet_wrap(~Metabolite, scales = "free_y") +
  theme_minimal(base_size = 11) +
  guides(fill = "none") +
  labs(title = "Top Variable Metabolites by Group", x = "Group", y = "log10(Value + 1)")

print(p_box)
ggsave(file.path(fig_dir, "eda_boxplots_top_metabolites.png"), p_box, width = 10, height = 8, dpi = 300)



#Univariate Tests

group_vals <- metadata[[group_col]]

run_test <- function(feature_vec) {
  df <- data.frame(v = feature_vec, g = group_vals)
  df <- df[complete.cases(df), ]
  if (nrow(df) < 10) return(NA_real_)
  suppressWarnings(kruskal.test(v ~ g, data = df)$p.value)
}

genus_p <- apply(genus_log, 2, run_test)
met_p <- apply(met_log, 2, run_test)

genus_res <- tibble(
  feature = names(genus_p),
  p_value = genus_p,
  p_adj = p.adjust(genus_p, method = "fdr")
) %>% arrange(p_adj)

met_res <- tibble(
  feature = names(met_p),
  p_value = met_p,
  p_adj = p.adjust(met_p, method = "fdr")
) %>% arrange(p_adj)


#log2-Fold-Change (Effektstärke) gegenüber nonIBD

# Gruppen-Vektor passend zur Zeilenreihenfolge der jeweiligen Matrix zuordnen
grp_genus <- metadata[[group_col]][match(rownames(genus_matrix), metadata[[sample_col]])]
grp_met   <- metadata[[group_col]][match(rownames(met_matrix),   metadata[[sample_col]])]

# Funktion: log2FC einer Gruppe gegen nonIBD (Median-basiert)
log2fc_vs_control <- function(mat, group, case, control = "nonIBD", pseudo) {
  med_case    <- apply(mat[group == case,    , drop = FALSE], 2, median, na.rm = TRUE)
  med_control <- apply(mat[group == control, , drop = FALSE], 2, median, na.rm = TRUE)
  log2((med_case + pseudo) / (med_control + pseudo))
}

# Genus: relative Abundanzen -> kleiner Pseudocount
lfc_genus_CD <- log2fc_vs_control(genus_matrix, grp_genus, "CD", pseudo = 1e-6)
lfc_genus_UC <- log2fc_vs_control(genus_matrix, grp_genus, "UC", pseudo = 1e-6)

genus_res <- genus_res %>%
  mutate(
    log2FC_CD_vs_nonIBD = lfc_genus_CD[feature],
    log2FC_UC_vs_nonIBD = lfc_genus_UC[feature]
  )

# Metabolite: Intensitäten -> Pseudocount 1
lfc_met_CD <- log2fc_vs_control(met_matrix, grp_met, case = "CD", control = "Control", pseudo = 1)
lfc_met_UC <- log2fc_vs_control(met_matrix, grp_met, case = "UC", control = "Control", pseudo = 1)

met_res <- met_res %>%
  mutate(
    log2FC_CD_vs_nonIBD = lfc_met_CD[feature],
    log2FC_UC_vs_nonIBD = lfc_met_UC[feature]
  )


#Speichern (jetzt inkl. log2FC-Spalten)

write_csv(genus_res, file.path(tab_dir, "eda_univariate_genus.csv"))
write_csv(met_res, file.path(tab_dir, "eda_univariate_metabolites.csv"))

message(sprintf("Signifikante Genera (FDR < 0.05): %d von %d",
                sum(genus_res$p_adj < 0.05, na.rm = TRUE), nrow(genus_res)))
message(sprintf("Signifikante Metabolite (FDR < 0.05): %d von %d",
                sum(met_res$p_adj < 0.05, na.rm = TRUE), nrow(met_res)))

head(genus_res, 10)
head(met_res, 10)



# 1) Stimmen die Gruppen-Labels?
table(grp_met, useNA = "always")     # sollte CD / UC / nonIBD zeigen, keine reinen NAs

# 2) Passen die Namen zusammen?
head(names(lfc_met_CD))              # wie heißen die Einträge im FC-Vektor?
head(met_res$feature)                # wie heißen sie in met_res?
mean(met_res$feature %in% names(lfc_met_CD))   # wenn 0 -> Namens-Mismatch = Ursache


# Beispiel-Feature: stimmen die Mediane wirklich überein?
f <- "C18-neg_Cluster_2083: NA"
tapply(met_matrix[, f], grp_met, median)   # Mediane je Gruppe CD / Control / UC


# Funktion (falls noch nicht definiert)
log2fc_vs_control <- function(mat, group, case, control, pseudo) {
  med_case    <- apply(mat[group == case,    , drop = FALSE], 2, median, na.rm = TRUE)
  med_control <- apply(mat[group == control, , drop = FALSE], 2, median, na.rm = TRUE)
  log2((med_case + pseudo) / (med_control + pseudo))
}

# Gruppenvektor passend zu den Zeilen der Genus-Matrix
grp_gen <- metadata$Study.Group[match(rownames(genus_matrix), metadata$Sample)]

# kleiner Pseudocount fuer relative Abundanzen
ps_gen <- 1e-6

lfc_gen_CD <- log2fc_vs_control(genus_matrix, grp_gen, "CD", "Control", ps_gen)
lfc_gen_UC <- log2fc_vs_control(genus_matrix, grp_gen, "UC", "Control", ps_gen)

# an genus_res anhaengen (Zuordnung ueber den vollen Feature-Namen)
genus_res$log2FC_CD_vs_nonIBD <- lfc_gen_CD[genus_res$feature]
genus_res$log2FC_UC_vs_nonIBD <- lfc_gen_UC[genus_res$feature]

head(genus_res, 10)

