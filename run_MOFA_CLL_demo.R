# run_MOFA_CLL_demo.R
# Demo: apply MOFA2 to CLL (Chronic Lymphocytic Leukemia) cohort  
# Based on the public demo: https://usda-ree-ars.github.io/SEAStats/multiomics_demo/Demo2_MOFA_chronicleukemia.html

# 0. (Optional) set a reproducible seed  
set.seed(12345)

# 1. Load libraries
library(data.table)
library(ggplot2)
library(tidyverse)
library(MOFA2)
library(randomForest)
library(survival)
library(survminer)

# 2. Load data (from web)
CLL_data <- readRDS(url("https://usda-ree-ars.github.io/SEAStatsData/omics/CLL_data.rds"))

#  Optionally inspect dimensions
print(lapply(CLL_data, dim))

# 3. Load sample metadata
CLL_metadata <- fread("https://usda-ree-ars.github.io/SEAStatsData/omics/sample_metadata.txt")
head(CLL_metadata)

# 4. Create the MOFA object
MOFAobject <- create_mofa(CLL_data)
MOFAobject  # examine object

# 5. (Optional) plot data overview
plot_data_overview(MOFAobject)

# 6. Define MOFA model options  
model_opts <- get_default_model_options(MOFAobject)
model_opts$likelihoods["Mutations"] <- "bernoulli"  # because mutation data is binary
model_opts$num_factors <- 15
print(model_opts)

# 7. Prepare MOFA object with options
MOFAobject <- prepare_mofa(MOFAobject, model_options = model_opts)

# 8. Train the MOFA model  
#   This may take several minutes
#  MOFAobject <- run_mofa(MOFAobject, use_basilisk = TRUE)

## If you don’t wish to re-train, you can instead load a precomputed model:
MOFAobject <- readRDS(url("https://usda-ree-ars.github.io/SEAStatsData/omics/MOFA2_CLL.rds"))

# 9. Add sample metadata  
stopifnot(all(CLL_metadata$sample %in% samples_metadata(MOFAobject)$sample))
samples_metadata(MOFAobject) <- CLL_metadata

# 10. (Optional) Rename features for interpretability  
drug_metadata <- fread("https://usda-ree-ars.github.io/SEAStatsData/omics/drugs.txt.gz")
gene_metadata <- fread("https://usda-ree-ars.github.io/SEAStatsData/omics/Hsapiens_genes_BioMart.87.txt.gz")

# Copy object before renaming if you want to keep original
MOFAobject.ensembl <- MOFAobject

# Retrieve current feature names
old_names <- features_names(MOFAobject)

# For Drugs: map drug IDs to drug names
tmp <- drug_metadata$name
names(tmp) <- drug_metadata$drug_id
new_drugs <- stringr::str_replace_all(old_names[["Drugs"]], tmp)

# For mRNA: map ENSG IDs to gene symbols (or fallback to ENSG)
gene_metadata[, symbol := ifelse(symbol == "", ens_id, symbol)]
tmp2 <- gene_metadata$symbol
names(tmp2) <- gene_metadata$ens_id
new_mrna <- stringr::str_replace_all(old_names[["mRNA"]], tmp2)

# Avoid name clashes with mutation view if needed
dupes <- new_mrna %in% old_names[["Mutations"]]
new_mrna[dupes] <- paste0(new_mrna[dupes], "_mRNA")

# Assign new names
new_names <- old_names
new_names[["Drugs"]] <- new_drugs
new_names[["mRNA"]] <- new_mrna
features_names(MOFAobject) <- new_names

# 11. Variance decomposition: how much variance each factor explains per view
plot_variance_explained(MOFAobject, max_r2 = 10)

# 12. Inspect factor values — e.g. factors 1 and 3
plot_factors(MOFAobject, factors = c(1, 3), dot_size = 2.5)

# 13. Inspect feature weights — for example: mutations for Factor 1
plot_weights(MOFAobject, view = "Mutations", factor = 1, nfeatures = 10, scale = TRUE)
plot_top_weights(MOFAobject, view = "Mutations", factor = 1, nfeatures = 10, scale = TRUE)

# 14. Plot factor 1 by a covariate — e.g. IGHV mutation status
plot_factor(MOFAobject,
            factors = 1,
            color_by = "IGHV",
            add_violin = TRUE,
            dodge = TRUE,
            show_missing = FALSE)

# 15. For mRNA — top weighted genes for Factor 1
plot_weights(MOFAobject, view = "mRNA", factor = 1, nfeatures = 10)

# 16. Scatter plots: factor values vs expression of top genes  
plot_data_scatter(MOFAobject,
                  view = "mRNA",
                  factor = 1,
                  features = 4,
                  sign = "positive",
                  color_by = "IGHV"
) + labs(y = "RNA expression")

plot_data_scatter(MOFAobject,
                  view = "mRNA",
                  factor = 1,
                  features = 4,
                  sign = "negative",
                  color_by = "IGHV"
) + labs(y = "RNA expression")

# 17. Heatmap: expression of top 25 genes (mRNA) vs Factor 1
plot_data_heatmap(MOFAobject,
                  view = "mRNA",
                  factor = 1,
                  features = 25,
                  cluster_rows = FALSE,
                  cluster_cols = FALSE,
                  show_rownames = TRUE,
                  show_colnames = FALSE,
                  scale = "row")

# 18. Gene set enrichment analysis (GSEA) for mRNA view  
reactomeGS <- readRDS(url("https://usda-ree-ars.github.io/SEAStatsData/omics/reactomeGS.rds"))

enrichment.results <- run_enrichment(
  object = MOFAobject.ensembl,
  view = "mRNA",
  feature.sets = reactomeGS,
  set.statistic = "mean.diff",   # or 'rank.sum'
  statistical.test = "parametric" # or 'permutation'
)

# Basic summary:
print(names(enrichment.results))
print(head(enrichment.results$pval.adj))

# Plot enrichment heatmap
plot_enrichment_heatmap(enrichment.results)

# 19. (Optional) Extract the raw factors & weights for custom / further analyses  
weights_df <- get_weights(MOFAobject, views = "all", factors = "all", as.data.frame = TRUE)
factors_df <- get_factors(MOFAobject, factors = "all", as.data.frame = TRUE)

head(weights_df)
head(factors_df)

# # 20. (Optional) Build predictive models — e.g. predict clinical markers (IGHV, trisomy12) from latent factors  
# df_factors <- as.data.frame(get_factors(MOFAobject, factors = c(1, 2, 3))[[1]])
# df_factors$IGHV <- as.factor(samples_metadata(MOFAobject)$IGHV)
# df_factors$trisomy12 <- as.factor(samples_metadata(MOFAobject)$trisomy12)

# # Example: random forest to predict IGHV status  
# rf_IGHV <- randomForest(IGHV ~ ., data = df_factors[!is.na(df_factors$IGHV), ])
# df_factors$IGHV_pred <- predict(rf_IGHV, df_factors)

# # Example: random forest to predict trisomy12
# rf_tri <- randomForest(trisomy12 ~ ., data = df_factors[!is.na(df_factors$trisomy12), ])
# df_factors$trisomy12_pred <- predict(rf_tri, df_factors)

# # (You may continue with plotting predictions or evaluating performance.)

# # 21. Print session info
# sessionInfo()
