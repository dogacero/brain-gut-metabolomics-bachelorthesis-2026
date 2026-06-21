# Brain-Gut Metabolomics — Bachelor Thesis Pipeline

Reproduzierbare R-Analysepipeline zur Bachelorarbeit  
**Brain-Gut Axis in Neurodegenerative Diseases Based on Metabolomics Data Using Machine Learning Methods**

**Autorin:** Doga Ceren Bozkurt  
**Jahr:** 2026  
**Repository:** https://github.com/dogacero/brain-gut-metabolomics-bachelorthesis-2026

---

## Überblick

Dieses Repository enthält den **Quellcode** (`R/`) und die **Ergebnisse** (`output/`) einer kombinierten Mikrobiom-Metabolom-Analyse auf Basis des Franzosa-IBD-Datensatzes. Untersucht werden Unterschiede zwischen Morbus Crohn (CD), Colitis ulcerosa (UC) und nonIBD-Kontrollen. Ziel ist die Identifikation krankheitsassoziierter Signaturen und die explorative Einordnung Brain-Gut-relevanter Metabolite.

**Externe Rohdatenquelle:** [microbiome-metabolome-curated-data](https://github.com/borenstein-lab/microbiome-metabolome-curated-data) (Müller et al., 2022)  
**Datensatz:** Franzosa et al. (2019), *Gut microbiome structure and metabolic activity in inflammatory bowel disease*

---

## Repository-Struktur

```
├── R/                          # Analysepipeline (Skripte 00–08)
├── output/
│   ├── tables/                 # Tabellarische Ergebnisse (CSV, TXT)
│   └── figures/                # Abbildungen (PNG)
├── renv.lock                   # Feste Paketversionen
├── .Rprofile                   # renv-Aktivierung
└── README.md
```

Rohdaten werden beim Ausführen von `01_load_data.R` automatisch heruntergeladen und lokal unter `data/raw/FRANZOSA_IBD_2019/` abgelegt (nicht im Repository enthalten).

---

## Pipeline (Ausführungsreihenfolge)

| Skript | Zweck |
|--------|--------|
| `R/00_setup.R` | Pakete installieren und laden |
| `R/01_load_data.R` | Download, Import, Harmonisierung der Proben-IDs (220 Proben) |
| `R/02_eda.R` | Shannon-Diversität, PCA, Heatmap, Boxplots, Kruskal-Wallis |
| `R/03_preprocessing.R` | Filter, Imputation, Transformation, NZV-Filter, Train-Test-Split |
| `R/04_ml_models.R` | Random Forest, Ridge-Regression, SVM, Elastic Net |
| `R/05_feature_importance.R` | Feature Importance je Modell |
| `R/06_pathway_analysis.R` | Metabolit-Mapping, Brain-Gut-Screening (Top 20) |
| `R/07_identifizierung_unbekannte_metaboliten.R` | Dokumentation unannotierter Top-Metabolite |
| `R/08_brain_gut_metabolite_check.R` | Erweitertes Brain-Gut-Screening über Gesamtdatensatz |

---

## Reproduktion

**Voraussetzungen:** R 4.5.x, RStudio optional, Internetzugang für Datendownload

1. Repository klonen  
   `git clone https://github.com/dogacero/brain-gut-metabolomics-bachelorthesis-2026.git`

2. Projektordner in RStudio öffnen

3. Pakete wiederherstellen  
   `renv::restore()`

4. Skripte der Reihe nach ausführen (`R/00_setup.R` bis `R/08_brain_gut_metabolite_check.R`)

5. Neue Ergebnisse werden unter `output/tables/` und `output/figures/` geschrieben

---

## Ergebnistabellen (`output/tables/`)

| Datei | Inhalt |
|-------|--------|
| `preprocessing_summary.csv` | Stichprobe und Featurezahlen nach Vorverarbeitung |
| `eda_univariate_genus.csv` | Univariate Genus-Tests mit FDR-Korrektur |
| `eda_univariate_metabolites.csv` | Univariate Metabolit-Tests mit FDR-Korrektur |
| `abbildung_4_6_fdr_signifikanz.csv` | Zusammenfassung FDR-signifikanter Features |
| `ml_model_performance.csv` | Accuracy, Kappa, AUC aller Modelle |
| `ml_confusion_RandomForest.csv` | Konfusionsmatrix Random Forest |
| `ml_confusion_LogisticRegression.csv` | Konfusionsmatrix Ridge-Regression |
| `ml_confusion_SVM.csv` | Konfusionsmatrix SVM |
| `ml_confusion_ElasticNet.csv` | Konfusionsmatrix Elastic Net |
| `feature_importance_all_models.csv` | Vollständige Importance-Werte |
| `feature_importance_top20_per_model.csv` | Top-20-Features je Modell |
| `pathway_mapped_top_metabolites.csv` | Metabolit-Mapping der Top-Features |
| `pathway_brain_gut_relevance_flags.csv` | Brain-Gut-Flags (Schlagwortscreening) |
| `pathway_brain_gut_summary.csv` | Zusammenfassung Brain-Gut-Treffer |
| `unknown_top_metabolites_for_review.csv` | Unannotierte Top-Metabolite (Review) |
| `unknown_top_metabolites_all_models.csv` | Unannotierte Top-Metabolite (alle Modelle) |
| `brain_gut_metabolites_in_dataset.csv` | Brain-Gut-Treffer im Gesamtdatensatz |
| `brain_gut_metabolites_importance_ranking.csv` | Brain-Gut-Ranking mit Importance |
| `metaboanalyst_input_*.csv/txt` | Export für externe MetaboAnalyst-Nutzung |

---

## Abbildungen (`output/figures/`)

| Datei | Inhalt |
|-------|--------|
| `eda_shannon_diversity.png` | Shannon-Diversität nach Studiengruppe |
| `eda_pca_genus.png` | PCA Genus-Ebene |
| `eda_pca_metabolites.png` | PCA Metabolit-Ebene |
| `eda_heatmap_top30_genera.png` | Heatmap Top-30-Genera |
| `eda_boxplots_top_metabolites.png` | Boxplots variabelste Metabolite |
| `feature_importance_top20_per_model.png` | Top-20-Features je Modell |
| `ml_roc_*.png` | ROC-Kurven je Modell |

---

## Hinweise

- **IBD als Modellsystem:** Es werden keine Parkinson- oder Alzheimer-Patienten untersucht. Brain-Gut-Bezüge sind hypothesengenerierend.
- **Untargeted Metabolomics:** Viele Features sind unannotiert (m/z, RT, Cluster-ID).
- **Kein Pathway-Enrichment:** Brain-Gut-Prüfung erfolgt schlagwort- und annotationsbasiert.

---

## Zitierung

> Bozkurt, D. C. (2026). *Brain-Gut Metabolomics Bachelor Thesis Pipeline*. GitHub.  
> https://github.com/dogacero/brain-gut-metabolomics-bachelorthesis-2026
