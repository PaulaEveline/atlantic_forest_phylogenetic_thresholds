# ============================================================
# SUPPLEMENTARY ANALYSIS — GDM WITH PHYLOGENETIC BETA DIVERSITY
# Robustness check for GradientForest genus-level PG approach
# Uses PhyloSor and betaMNTD as response variables
# Accepts pre-matched objects from tree matching scripts
#
# KEY DESIGN DECISIONS:
# geo = FALSE — geographic distance excluded to ensure direct
#   comparability with GradientForest which uses only land-use
#   predictors. Geographic distance dominated anuran GDM when
#   included (geo=TRUE), reflecting dispersal limitation rather
#   than environmental effects, which is not the target inference.
#
# Two metrics used:
#   PhyloSor — proportion of unique branch length between sites
#     (analogous to Sorensen; sensitive to whole-tree turnover)
#   betaMNTD — mean nearest-taxon distance between communities
#     (sensitive to phylogenetic clustering/overdispersion;
#     captures clade-level exclusion that PhyloSor may miss)
# Using both is standard practice; agreement strengthens
# conclusions, disagreement reveals turnover structure.
#
# RERUN CONTROL:
# Set rerun_distances = FALSE to skip PhyloSor/betaMNTD
# computation and load from saved RDS files instead (faster)
# ============================================================

library(gdm)
library(picante)
library(ape)
library(vegan)
if (requireNamespace("phytools", quietly = TRUE)) library(phytools)

rerun_distances <- FALSE   # set TRUE only for a fresh run from scratch

# ============================================================
# STEP 1 — LOAD PRE-MATCHED OBJECTS
# Run anuran_tree_matching.R and bird_tree_matching.R first
# ============================================================

tree_birds   <- readRDS("bird_tree_pruned.rds")
comm_birds   <- readRDS("bird_comm_matched.rds")
tree_anurans <- readRDS("anuran_tree_pruned.rds")
comm_anurans <- readRDS("anuran_comm_matched.rds")

# Environmental data — same files used in GradientForest
envir_birds_raw   <- read.table("sites_envir_birds_coord.txt",
                                header = TRUE, row.names = 1, sep = "\t")
envir_anurans_raw <- read.table("sites_envir_anura_coord.txt",
                                header = TRUE, row.names = 1, sep = "\t")

# Separate coordinates from predictors
# Coordinates kept for site-pair table formatting but geo = FALSE
# means GDM will NOT use them as a predictor
coord_cols <- c(1, 2)

coords_birds   <- envir_birds_raw[,  coord_cols]
coords_anurans <- envir_anurans_raw[, coord_cols]
envir_birds    <- envir_birds_raw[,  -coord_cols]
envir_anurans  <- envir_anurans_raw[, -coord_cols]

colnames(coords_birds)   <- c("xCoord", "yCoord")
colnames(coords_anurans) <- c("xCoord", "yCoord")

cat("Birds:   ", nrow(comm_birds),   "sites,", ncol(comm_birds),
    "species,", ncol(envir_birds),   "predictors\n")
cat("Anurans: ", nrow(comm_anurans), "sites,", ncol(comm_anurans),
    "species,", ncol(envir_anurans), "predictors\n")

# ============================================================
# COLOUR PALETTE — Wong (2011) colourblind-friendly
# ============================================================

col_birds   <- "#0072B2"   # blue
col_anurans <- "#E69F00"   # orange

# ============================================================
# STEP 2 — COMPUTE PHYLOGENETIC BETA DIVERSITY
# PhyloSor: proportion of unique branch length (0-1, GDM-ready)
# betaMNTD: raw phylogenetic distance in Ma, rescaled to 0-1
#   by dividing by tree root age (maximum possible distance)
# ============================================================

if (rerun_distances) {
  
  cat("\nComputing phylogenetic beta diversity — Birds\n")
  phylosor_birds     <- as.dist(phylosor(comm_birds, tree_birds))
  betamntd_birds_raw <- as.dist(comdistnt(
    comm_birds, cophenetic(tree_birds), abundance.weighted = FALSE))
  
  cat("PhyloSor range (birds):",
      round(range(phylosor_birds, na.rm = TRUE), 3), "\n")
  cat("betaMNTD range before rescaling (birds):",
      round(range(betamntd_birds_raw, na.rm = TRUE), 3), "\n")
  
  max_dist_birds   <- max(cophenetic(tree_birds)) / 2
  cat("Bird tree root age (Ma):", round(max_dist_birds, 1), "\n")
  betamntd_birds   <- betamntd_birds_raw / max_dist_birds
  betamntd_birds[betamntd_birds > 1] <- 1
  cat("betaMNTD range after rescaling (birds):",
      round(range(betamntd_birds, na.rm = TRUE), 3), "\n")
  
  cat("\nComputing phylogenetic beta diversity — Anurans\n")
  phylosor_anurans     <- as.dist(phylosor(comm_anurans, tree_anurans))
  betamntd_anurans_raw <- as.dist(comdistnt(
    comm_anurans, cophenetic(tree_anurans), abundance.weighted = FALSE))
  
  cat("PhyloSor range (anurans):",
      round(range(phylosor_anurans, na.rm = TRUE), 3), "\n")
  cat("betaMNTD range before rescaling (anurans):",
      round(range(betamntd_anurans_raw, na.rm = TRUE), 3), "\n")
  
  max_dist_anurans <- max(cophenetic(tree_anurans)) / 2
  cat("Anuran tree root age (Ma):", round(max_dist_anurans, 1), "\n")
  betamntd_anurans <- betamntd_anurans_raw / max_dist_anurans
  betamntd_anurans[betamntd_anurans > 1] <- 1
  cat("betaMNTD range after rescaling (anurans):",
      round(range(betamntd_anurans, na.rm = TRUE), 3), "\n")
  
  saveRDS(phylosor_birds,   "gdm_phylosor_birds.rds")
  saveRDS(betamntd_birds,   "gdm_betamntd_birds.rds")
  saveRDS(phylosor_anurans, "gdm_phylosor_anurans.rds")
  saveRDS(betamntd_anurans, "gdm_betamntd_anurans.rds")
  cat("\nDistance matrices saved\n")
  
} else {
  
  cat("Loading saved distance matrices (rerun_distances = FALSE)\n")
  phylosor_birds   <- readRDS("gdm_phylosor_birds.rds")
  betamntd_birds   <- readRDS("gdm_betamntd_birds.rds")
  phylosor_anurans <- readRDS("gdm_phylosor_anurans.rds")
  betamntd_anurans <- readRDS("gdm_betamntd_anurans.rds")
  cat("Loaded\n")
  
}

# ============================================================
# STEP 3 — MANTEL TESTS
# Compare phylogenetic vs taxonomic beta diversity
# Negative correlation = sites taxonomically similar are
# phylogenetically dissimilar = phylogenetic overdispersion
# at disturbed sites (generalists from distant lineages co-occur)
# ============================================================

cat("\n--- Mantel tests ---\n")

taxon_birds   <- vegdist(comm_birds,   method = "bray")
taxon_anurans <- vegdist(comm_anurans, method = "bray")

mantel_birds   <- mantel(phylosor_birds,   taxon_birds,
                         method = "spearman", permutations = 999)
mantel_anurans <- mantel(phylosor_anurans, taxon_anurans,
                         method = "spearman", permutations = 999)

cat("Birds   — PhyloSor vs Bray-Curtis: r =",
    round(mantel_birds$statistic, 3),
    "  p =", mantel_birds$signif, "\n")
cat("Anurans — PhyloSor vs Bray-Curtis: r =",
    round(mantel_anurans$statistic, 3),
    "  p =", mantel_anurans$signif, "\n")

# ============================================================
# STEP 4 — PREPARE GDM SITE-PAIR TABLES
# geo = FALSE: geographic distance excluded from all models
# Coordinates still required by formatsitepair for formatting
# GDM bioFormat=3 requires numeric site IDs
# ============================================================

prepare_gdm_table <- function(dist_mat, envir, coords) {
  
  site_names  <- rownames(as.matrix(dist_mat))
  site_id_map <- setNames(seq_along(site_names), site_names)
  
  dist_mat_num <- as.matrix(dist_mat)
  rownames(dist_mat_num) <- site_id_map[rownames(dist_mat_num)]
  colnames(dist_mat_num) <- site_id_map[colnames(dist_mat_num)]
  dist_df <- as.data.frame(dist_mat_num)
  dist_df <- cbind(site = as.integer(rownames(dist_df)), dist_df)
  
  site_data        <- cbind(coords, envir)
  site_data$site   <- site_id_map[rownames(site_data)]
  rownames(site_data) <- site_data$site
  
  gdm::formatsitepair(
    bioData    = dist_df,
    bioFormat  = 3,
    XColumn    = "xCoord",
    YColumn    = "yCoord",
    siteColumn = "site",
    predData   = site_data
  )
}

# ============================================================
# STEP 5 — FIT GDM MODELS (geo = FALSE throughout)
# ============================================================

cat("\n--- Fitting GDM models (geo = FALSE) ---\n")

cat("Birds PhyloSor...\n")
st_ps_birds  <- prepare_gdm_table(phylosor_birds,   envir_birds,   coords_birds)
gdm_ps_birds <- gdm(st_ps_birds, geo = FALSE)
cat("  Deviance explained:", round(gdm_ps_birds$explained, 3), "\n")

cat("Birds betaMNTD...\n")
st_bm_birds  <- prepare_gdm_table(betamntd_birds,   envir_birds,   coords_birds)
gdm_bm_birds <- gdm(st_bm_birds, geo = FALSE)
cat("  Deviance explained:", round(gdm_bm_birds$explained, 3), "\n")

cat("Anurans PhyloSor...\n")
st_ps_anurans  <- prepare_gdm_table(phylosor_anurans, envir_anurans, coords_anurans)
gdm_ps_anurans <- gdm(st_ps_anurans, geo = FALSE)
cat("  Deviance explained:", round(gdm_ps_anurans$explained, 3), "\n")

cat("Anurans betaMNTD...\n")
st_bm_anurans  <- prepare_gdm_table(betamntd_anurans, envir_anurans, coords_anurans)
gdm_bm_anurans <- gdm(st_bm_anurans, geo = FALSE)
cat("  Deviance explained:", round(gdm_bm_anurans$explained, 3), "\n")

# ============================================================
# STEP 6 — PERMUTATION IMPORTANCE TESTS (geo = FALSE)
# ============================================================

n_perm <- 999   # reduce to 99 for a quick test run

cat("\n--- Permutation tests (", n_perm, "permutations, geo = FALSE) ---\n")

cat("Birds PhyloSor...\n")
perm_ps_birds   <- gdm.varImp(st_ps_birds,   geo = FALSE, nPerm = n_perm,
                              parallel = FALSE)

cat("Birds betaMNTD...\n")
perm_bm_birds   <- gdm.varImp(st_bm_birds,   geo = FALSE, nPerm = n_perm,
                              parallel = FALSE)

cat("Anurans PhyloSor...\n")
perm_ps_anurans <- gdm.varImp(st_ps_anurans, geo = FALSE, nPerm = n_perm,
                              parallel = FALSE)

cat("Anurans betaMNTD...\n")
perm_bm_anurans <- gdm.varImp(st_bm_anurans, geo = FALSE, nPerm = n_perm,
                              parallel = FALSE)

cat("\nVariable importance — Birds PhyloSor:\n");   print(perm_ps_birds[[3]])
cat("\nVariable importance — Birds betaMNTD:\n");   print(perm_bm_birds[[3]])
cat("\nVariable importance — Anurans PhyloSor:\n"); print(perm_ps_anurans[[3]])
cat("\nVariable importance — Anurans betaMNTD:\n"); print(perm_bm_anurans[[3]])

# ============================================================
# STEP 7 — EXTRACT I-SPLINES AND APPROXIMATE THRESHOLDS
# ============================================================

isplines_ps_birds   <- isplineExtract(gdm_ps_birds)
isplines_bm_birds   <- isplineExtract(gdm_bm_birds)
isplines_ps_anurans <- isplineExtract(gdm_ps_anurans)
isplines_bm_anurans <- isplineExtract(gdm_bm_anurans)

get_ispline_peak <- function(isplines, pred) {
  x    <- isplines$x[, pred]
  y    <- isplines$y[, pred]
  dy   <- diff(y) / diff(x)
  xmid <- (x[-1] + x[-length(x)]) / 2
  xmid[which.max(dy)]
}

make_imp_table <- function(isplines_ps, isplines_bm, label) {
  preds <- colnames(isplines_ps$x)
  # use only predictors present in both (geo=FALSE may drop Geographic)
  preds_bm <- colnames(isplines_bm$x)
  preds_shared <- intersect(preds, preds_bm)
  data.frame(
    taxon          = label,
    predictor      = preds_shared,
    ispline_sum_ps = colSums(isplines_ps$y[, preds_shared, drop = FALSE]),
    ispline_sum_bm = colSums(isplines_bm$y[, preds_shared, drop = FALSE]),
    threshold_ps   = sapply(preds_shared, get_ispline_peak,
                            isplines = isplines_ps),
    threshold_bm   = sapply(preds_shared, get_ispline_peak,
                            isplines = isplines_bm)
  )
}

imp_birds   <- make_imp_table(isplines_ps_birds,   isplines_bm_birds,   "Birds")
imp_anurans <- make_imp_table(isplines_ps_anurans, isplines_bm_anurans, "Anurans")
imp_birds   <- imp_birds[order(-imp_birds$ispline_sum_ps), ]
imp_anurans <- imp_anurans[order(-imp_anurans$ispline_sum_ps), ]

cat("\nTop 3 predictors — Birds (PhyloSor):\n")
print(imp_birds[1:3,   c("predictor", "ispline_sum_ps", "threshold_ps")])
cat("\nTop 3 predictors — Anurans (PhyloSor):\n")
print(imp_anurans[1:3, c("predictor", "ispline_sum_ps", "threshold_ps")])

write.csv(imp_birds,   "GDM_summary_Birds.csv",   row.names = FALSE)
write.csv(imp_anurans, "GDM_summary_Anurans.csv", row.names = FALSE)

# ============================================================
# STEP 8 — I-SPLINE PLOTS PER TAXON
# Row 1: PhyloSor | Row 2: betaMNTD
# Top 3 predictors by PhyloSor I-spline sum
# ============================================================

plot_isplines <- function(isplines_ps, isplines_bm, imp_table,
                          taxon_label, col_line) {
  
  top3     <- imp_table$predictor[1:min(3, nrow(imp_table))]
  out_file <- paste0("GDM_isplines_", taxon_label, ".png")
  
  png(out_file, width = 13, height = 8,
      units = "in", res = 300, bg = "white")
  
  par(mfrow  = c(2, 3), mar = c(4, 4.5, 2.5, 0.8),
      oma = c(0, 0, 3, 0), mgp = c(2.8, 0.7, 0), family = "serif")
  
  for (k in seq_along(top3)) {
    pred <- top3[k]
    ymax <- max(isplines_ps$y[, pred]) * 1.2
    plot(isplines_ps$x[, pred], isplines_ps$y[, pred],
         type = "l", lwd = 2.5, col = col_line,
         ylim = c(0, ymax),
         xlab = gsub("\\.", " ", pred),
         ylab = "Partial ecological distance",
         bty = "l", cex.axis = 1.1, cex.lab = 1.0)
    mtext(gsub("\\.", " ", pred), side = 3, line = 0.3, font = 2, cex = 0.95)
    abline(v = get_ispline_peak(isplines_ps, pred),
           lty = 2, col = "grey50", lwd = 1)
    if (k == 1) mtext("PhyloSor", side = 2, line = 3.5, font = 2, cex = 0.9)
  }
  
  for (k in seq_along(top3)) {
    pred <- top3[k]
    ymax <- max(isplines_bm$y[, pred]) * 1.2
    plot(isplines_bm$x[, pred], isplines_bm$y[, pred],
         type = "l", lwd = 2.5, col = col_line,
         ylim = c(0, ymax),
         xlab = gsub("\\.", " ", pred),
         ylab = "Partial ecological distance",
         bty = "l", cex.axis = 1.1, cex.lab = 1.0)
    abline(v = get_ispline_peak(isplines_bm, pred),
           lty = 2, col = "grey50", lwd = 1)
    if (k == 1) mtext("betaMNTD", side = 2, line = 3.5, font = 2, cex = 0.9)
  }
  
  mtext(paste("GDM I-splines —", taxon_label),
        outer = TRUE, side = 3, line = 0.8,
        font = 2, cex = 1.3, family = "serif")
  dev.off()
  cat("Saved", out_file, "\n")
}

plot_isplines(isplines_ps_birds,   isplines_bm_birds,
              imp_birds,   "Birds",   col_birds)
plot_isplines(isplines_ps_anurans, isplines_bm_anurans,
              imp_anurans, "Anurans", col_anurans)

# ============================================================
# STEP 9 — CROSS-TAXA COMPARISON FIGURE
# Top 3 shared predictors — PhyloSor row 1, betaMNTD row 2
# ============================================================

shared_preds <- intersect(imp_birds$predictor, imp_anurans$predictor)
top_shared   <- shared_preds[1:min(3, length(shared_preds))]

if (length(top_shared) > 0) {
  
  png("GDM_comparison_taxa.png",
      width = 13, height = 8, units = "in", res = 300, bg = "white")
  
  par(mfrow = c(2, length(top_shared)),
      mar   = c(4, 4.5, 2.5, 0.8), oma = c(0, 0, 3, 0),
      mgp   = c(2.8, 0.7, 0), family = "serif")
  
  for (row in 1:2) {
    isp_b <- if (row == 1) isplines_ps_birds   else isplines_bm_birds
    isp_a <- if (row == 1) isplines_ps_anurans else isplines_bm_anurans
    metric <- if (row == 1) "PhyloSor" else "betaMNTD"
    
    for (k in seq_along(top_shared)) {
      pred <- top_shared[k]
      xr   <- range(c(isp_b$x[, pred], isp_a$x[, pred]))
      yr   <- c(0, max(c(isp_b$y[, pred], isp_a$y[, pred])) * 1.15)
      
      plot(xr, yr, type = "n",
           xlab = gsub("\\.", " ", pred),
           ylab = "Partial ecological distance",
           bty = "l", cex.axis = 1.1, cex.lab = 1.0)
      
      lines(isp_b$x[, pred], isp_b$y[, pred], lwd = 2.5, col = col_birds)
      lines(isp_a$x[, pred], isp_a$y[, pred], lwd = 2.5, col = col_anurans)
      mtext(gsub("\\.", " ", pred), side = 3, line = 0.3, font = 2, cex = 0.95)
      
      if (k == 1) {
        mtext(metric, side = 2, line = 3.5, font = 2, cex = 0.9)
        if (row == 1)
          legend("topleft", legend = c("Birds", "Anurans"),
                 col = c(col_birds, col_anurans),
                 lwd = 2.5, bty = "n", cex = 0.9)
      }
    }
  }
  
  mtext("GDM — phylogenetic beta diversity across taxa",
        outer = TRUE, side = 3, line = 0.8,
        font = 2, cex = 1.3, family = "serif")
  dev.off()
  cat("Saved GDM_comparison_taxa.png\n")
  
} else {
  cat("No shared predictors between taxa — cross-taxa figure skipped\n")
}

# ============================================================
# FINAL SUMMARY
# ============================================================

cat("\n", rep("=", 60), "\n", sep = "")
cat("SUMMARY (geo = FALSE)\n")
cat(rep("=", 60), "\n", sep = "")
cat("\nDeviance explained:\n")
cat("  Birds   — PhyloSor:", round(gdm_ps_birds$explained,   3), "\n")
cat("  Birds   — betaMNTD:", round(gdm_bm_birds$explained,   3), "\n")
cat("  Anurans — PhyloSor:", round(gdm_ps_anurans$explained, 3), "\n")
cat("  Anurans — betaMNTD:", round(gdm_bm_anurans$explained, 3), "\n")
cat("\nMantel correlations (PhyloSor vs Bray-Curtis):\n")
cat("  Birds:   r =", round(mantel_birds$statistic,   3),
    "  p =", mantel_birds$signif,   "\n")
cat("  Anurans: r =", round(mantel_anurans$statistic, 3),
    "  p =", mantel_anurans$signif, "\n")
cat("\nOutputs:\n")
cat("  GDM_isplines_Birds.png\n")
cat("  GDM_isplines_Anurans.png\n")
cat("  GDM_comparison_taxa.png\n")
cat("  GDM_summary_Birds.csv\n")
cat("  GDM_summary_Anurans.csv\n")
cat(rep("=", 60), "\n", sep = "")
