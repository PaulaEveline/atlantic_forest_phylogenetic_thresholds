# ============================================================
# 01_VIF COLLINEARITY ANALYSIS
# Reduces environmental predictors from 9 to 7 for each taxon
#
# Input:
#   predictors_all_wide.csv — matrix of 9 land-use predictors
#   across all sampling units (birds n=32 and anurans n=15).
#   Columns 1-3: hex_id, latitude, longitude.
#   Columns 4-12: Forest, Eucalyptus, Pioneer.vegetation,
#   Urban.areas, Agriculture, Bodies.of.water, Wetlands,
#   Bare.soil, Heterogeneity (Shannon diversity index).
#
# Procedure (two-step):
#   Step 1 — vifcor: removes predictors with pairwise Pearson
#             correlation |r| > 0.7
#   Step 2 — vifstep: removes predictors with VIF > 5
#   Forest cover is retained throughout regardless of
#   collinearity diagnostics (ecological centrality).
#
# Output:
#   predictors_selected_with_coords.csv — final predictor set
#   with coordinates, ready for GradientForest and GDM analyses
#   correlation_plot.png — pairwise Pearson correlation matrix
#   of retained predictors (visual check)
#
# Note: the script is run separately for birds and anurans by
#   subsetting the input to the relevant sampling units before
#   running the VIF procedure. The excluded predictors differ
#   between taxa: bare soil and agriculture excluded for
#   anurans; bare soil and heterogeneity excluded for birds.
#
# Reference: Zuur et al. (2010) Methods Ecol. Evol. 1, 3-14.
# ============================================================

library(usdm)
library(readr)
library(dplyr)
library(corrplot)

# ---- Load data ----
env <- read_csv("predictors_all_wide.csv")

# Separate coordinates from predictors
hex_coords <- env[, 1:3]   # hex_id, latitude, longitude
predictors <- env[, -c(1:3)]

# ---- Clean predictor matrix ----
# Convert all columns to numeric
predictors_numeric <- data.frame(lapply(predictors, function(x) {
  if (!is.numeric(x)) x <- as.numeric(as.character(x))
  return(x)
}))

# Remove columns with all NAs or zero variance (constant values)
predictors_clean <- predictors_numeric[, sapply(predictors_numeric, function(x) {
  is.numeric(x) && !all(is.na(x)) && length(unique(na.omit(x))) > 1
})]

cat("Predictors before VIF procedure:", ncol(predictors_clean), "\n")
cat("Predictor names:", paste(colnames(predictors_clean), collapse = ", "), "\n")

# ---- Check Forest cover is present ----
forest_col <- "Forest"   # update if column name differs
if (!forest_col %in% colnames(predictors_clean)) {
  stop(paste("Forest cover column '", forest_col,
             "' not found. Check column names."))
}

# ---- Step 1: vifcor — remove pairwise Pearson |r| > 0.7 ----
cat("\nStep 1: vifcor (Pearson |r| > 0.7)\n")
cor_filtered <- vifcor(predictors_clean, th = 0.7)
print(cor_filtered)

vars_after_cor <- cor_filtered@results$Variables
cat("Predictors retained after vifcor:", length(vars_after_cor), "\n")
cat(paste(vars_after_cor, collapse = ", "), "\n")

# Force Forest cover retention even if flagged by vifcor
if (!forest_col %in% vars_after_cor) {
  cat("Note: Forest cover was excluded by vifcor but retained",
      "due to ecological centrality.\n")
  vars_after_cor <- c(forest_col, vars_after_cor)
}

# ---- Step 2: vifstep — remove VIF > 5 ----
cat("\nStep 2: vifstep (VIF > 5)\n")
vif_selected <- vifstep(predictors_clean[, vars_after_cor, drop = FALSE],
                        th = 5)
print(vif_selected)

selected_vars <- unique(c(forest_col, vif_selected@results$Variables))
cat("\nFinal predictors retained:", length(selected_vars), "\n")
cat(paste(selected_vars, collapse = ", "), "\n")

cat("\nPredictors excluded:\n")
excluded <- setdiff(colnames(predictors_clean), selected_vars)
cat(paste(excluded, collapse = ", "), "\n")

# ---- Combine with coordinates ----
final_dataset <- cbind(hex_coords,
                       predictors_clean[, selected_vars, drop = FALSE])

# ---- Save output ----
write_csv(final_dataset, "predictors_selected_with_coords.csv")
cat("\nFinal dataset saved to predictors_selected_with_coords.csv\n")
cat("Dimensions:", nrow(final_dataset), "rows x",
    ncol(final_dataset), "columns\n")

# ---- Correlation matrix of retained predictors ----
cor_matrix <- cor(predictors_clean[, selected_vars],
                  use = "pairwise.complete.obs")

cat("\nPearson correlation matrix of retained predictors:\n")
print(round(cor_matrix, 2))

png("correlation_plot.png",
    width = 7, height = 6, units = "in", res = 300, bg = "white")
corrplot(cor_matrix,
         method     = "color",
         type       = "upper",
         order      = "hclust",
         tl.cex     = 0.8,
         tl.col     = "black",
         addCoef.col = "black",
         number.cex = 0.7,
         mar        = c(0, 0, 1, 0),
         title      = "Pearson correlations — retained predictors")
dev.off()
cat("Correlation plot saved to correlation_plot.png\n")