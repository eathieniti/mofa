###############################
# Minimal MOFA CLL Demo Script
# FAST, clean, ready for teaching
###############################

# Output folder
OUTDIR <- "mofa_cll_demo_output"
dir.create(OUTDIR, showWarnings = FALSE)

# ---- 1. Packages ----
pkgs <- c("MOFA2", "data.table", "ggplot2")
need <- pkgs[!pkgs %in% installed.packages()[,1]]
if (length(need)) install.packages(need, repos="https://cloud.r-project.org")
lapply(pkgs, library, character.only = TRUE)

# ---- 2. Load data from web ----
CLL_data <- readRDS(url(
  "https://usda-ree-ars.github.io/SEAStatsData/omics/CLL_data.rds"
))
CLL_meta <- data.table::fread(
  "https://usda-ree-ars.github.io/SEAStatsData/omics/sample_metadata.txt"
)

# ---- 3. Load precomputed MOFA model (FAST) ----
MOFAobject <- readRDS(url(
  "https://usda-ree-ars.github.io/SEAStatsData/omics/MOFA2_CLL.rds"
))

# ---- 4. Add metadata ----
samples_metadata(MOFAobject) <- as.data.frame(CLL_meta)

# ---- 5. Basic MOFA plots ----

# Variance explained plot
png(file.path(OUTDIR, "variance_explained.png"), 900, 600)
plot_variance_explained(MOFAobject)
dev.off()

# Factor 1 vs 2, colored by IGHV status
png(file.path(OUTDIR, "factor1_vs_factor2_IGHV.png"), 900, 600)
plot_factors(MOFAobject,
             factors = c(1, 2),
             color_by = "IGHV",
             dot_size = 2.5)
dev.off()

# Top mutation features for Factor 1
png(file.path(OUTDIR, "top_mutations_factor1.png"), 900, 600)
plot_top_weights(MOFAobject,
                 view = "Mutations",
                 factor = 1,
                 nfeatures = 10,
                 scale = TRUE)
dev.off()

# Top mRNA features for Factor 1
png(file.path(OUTDIR, "top_mrna_factor1.png"), 900, 600)
plot_top_weights(MOFAobject,
                 view = "mRNA",
                 factor = 1,
                 nfeatures = 10,
                 scale = TRUE)
dev.off()

# ---- 6. Save factor matrix ----
factors <- as.data.frame(get_factors(MOFAobject, factors = "all")$group1)
data.table::fwrite(factors,
                   file = file.path(OUTDIR, "MOFA_factors.tsv"),
                   sep = "\t")

message("DONE! Results saved to: ", normalizePath(OUTDIR))
