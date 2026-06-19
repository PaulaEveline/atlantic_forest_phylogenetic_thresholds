# ============================================================
# 02_GRADIENTFOREST ANALYSIS — BIRDS & ANURANS
#
# Produces:
#   Figure 2  — Bootstrap median variable importance (both taxa)
#   Figure 3a — Bird panel plot (split density + cumulative curves)
#   Figure 3b — Anuran panel plot (split density + cumulative curves)
#   importance_birds.csv / importance_anurans.csv
#   thresholds_birds.csv / thresholds_anurans.csv
#
# Inputs:
#   site_PG_birds_new.txt  — bird phylogenetic group occurrence matrix
#   site_PG_Anura.txt      — anuran phylogenetic group occurrence matrix
#   sites_envir_birds_coord.txt — bird environmental predictors + coords
#   sites_envir_anura_coord.txt — anuran environmental predictors + coords
#   (output of 01_vif_analysis.R — 7 predictors per taxon)
#
# Reproducibility:
#   set.seed(42) is applied before the analysis loop.
#   Set rerun_analysis = TRUE to run bootstraps from scratch (~2-3 hrs).
#   Set rerun_analysis = FALSE to load saved CSVs and regenerate
#   figures only (requires gf_model_birds and gf_model_anurans
#   objects to be in the R environment from a previous session).
#
# Reference: Ellis et al. (2012) Ecology 93, 156-168.
# ============================================================

library(gradientForest)

# ============================================================
# OUTPUT FOLDER
# Change this path to match your working directory.
# All figures and CSV summaries will be saved here.
# ============================================================

output_dir <- getwd()
# output_dir <- "C:/your/path/here"   # uncomment to set explicitly
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat("Outputs will be saved to:", output_dir, "\n")

# ============================================================
# RERUN CONTROL
# TRUE  — runs full analysis including bootstraps (~2-3 hours)
# FALSE — loads saved CSVs and regenerates figures only
#         (requires gf_model objects in R environment)
# ============================================================

rerun_analysis <- TRUE   # set FALSE to skip bootstraps
set.seed(42)             # ensures reproducibility across runs

# ============================================================
# TAXA CONFIGURATION
# ============================================================

taxa_config <- list(
  birds = list(
    comm_file  = "site_PG_birds_new.txt",
    envir_file = "sites_envir_birds_coord.txt",
    coord_cols = c(1, 2),
    ntree_full = 5000,   # ntree selected based on stability analysis
    ntree_boot = 1000,
    output_fig = file.path(output_dir, "Figure3a_birds_GF.png"),
    output_imp = file.path(output_dir, "importance_birds.csv"),
    output_thr = file.path(output_dir, "thresholds_birds.csv"),
    label      = "Birds"
  ),
  anurans = list(
    comm_file  = "site_PG_Anura.txt",
    envir_file = "sites_envir_anura_coord.txt",
    coord_cols = c(1, 2),
    ntree_full = 2000,   # stable from ntree=1000; 2000 used as conservative choice
    ntree_boot = 1000,
    output_fig = file.path(output_dir, "Figure3b_anurans_GF.png"),
    output_imp = file.path(output_dir, "importance_anurans.csv"),
    output_thr = file.path(output_dir, "thresholds_anurans.csv"),
    label      = "Anurans"
  )
)

# ============================================================
# PHYLOGENETIC GROUP TO GENUS LOOKUP TABLES
# ============================================================

pg_to_genus_birds <- c(
  PG_1="Amazona", PG_2="Aphantochroa", PG_3="Aramides", PG_4="Aratinga",
  PG_5="Arremon", PG_6="Asemospiza", PG_7="Attila", PG_8="Automolus",
  PG_9="Basileuterus", PG_10="Brotogeris", PG_11="Campephilus",
  PG_12="Camptostoma", PG_13="Caracara", PG_14="Celeus",
  PG_15="Certhiaxis", PG_16="Chionomesa", PG_17="Chiroxiphia",
  PG_18="Chlorostilbon", PG_19="Chondrohierax", PG_20="Coereba",
  PG_21="Colaptes", PG_22="Colonia", PG_23="Columbina",
  PG_24="Conirostrum", PG_25="Conopophaga", PG_26="Coryphospingus",
  PG_27="Corythopis", PG_28="Crypturellus", PG_29="Cyanocorax",
  PG_30="Cyclarhis", PG_31="Dacnis", PG_32="Dendroma",
  PG_33="Drymophila", PG_34="Dryocopus", PG_35="Dysithamnus",
  PG_36="Elaenia", PG_37="Empidonomus", PG_38="Euphonia",
  PG_39="Eupsittula", PG_40="Florisuga", PG_41="Forpus",
  PG_42="Furnarius", PG_43="Geotrygon", PG_44="Habia",
  PG_45="Hemithraupis", PG_46="Hemitriccus", PG_47="Herpetotheres",
  PG_48="Herpsilochmus", PG_49="Hylophilus", PG_50="Hypoedaleus",
  PG_51="Icterus", PG_52="Ictinia", PG_53="Lathrotriccus",
  PG_54="Legatus", PG_55="Leptodon", PG_56="Leptopogon",
  PG_57="Leptotila", PG_58="Leucochloris", PG_59="Lochmias",
  PG_60="Malacoptila", PG_61="Manacus", PG_62="Megarynchus",
  PG_63="Melanerpes", PG_64="Micrastur", PG_65="Milvago",
  PG_66="Mimus", PG_67="Mionectes", PG_68="Myiarchus",
  PG_69="Myiodynastes", PG_70="Myiopagis", PG_71="Myiophobus",
  PG_72="Myiornis", PG_73="Myiothlypis", PG_74="Myiozetetes",
  PG_75="Myrmoderus", PG_76="Nemosia", PG_77="Neopelma",
  PG_78="Pachyramphus", PG_79="Patagioenas", PG_80="Penelope",
  PG_81="Phacellodomus", PG_82="Phaethornis", PG_83="Phyllomyias",
  PG_84="Phylloscartes", PG_85="Piaya", PG_86="Picumnus",
  PG_87="Pionus", PG_88="Pitangus", PG_89="Platyrinchus",
  PG_90="Poecilotriccus", PG_91="Procnias", PG_92="Psarocolius",
  PG_93="Psittacara", PG_94="Pyriglena", PG_95="Pyroderus",
  PG_96="Ramphastos", PG_97="Rhopias", PG_98="Rupornis",
  PG_99="Saltator", PG_100="Schiffornis", PG_101="Setophaga",
  PG_102="Sittasomus", PG_103="Sporophila", PG_104="Stelgidopteryx",
  PG_105="Stilpnia", PG_106="Synallaxis", PG_107="Syndactyla",
  PG_108="Tachyphonus", PG_109="Tangara", PG_110="Tapera",
  PG_111="Tersina", PG_112="Thalurania", PG_113="Thamnophilus",
  PG_114="Thlypopsis", PG_115="Thraupis", PG_116="Tityra",
  PG_117="Todirostrum", PG_118="Tolmomyias", PG_119="Trichothraupis",
  PG_120="Troglodytes", PG_121="Trogon", PG_122="Turdus",
  PG_123="Tyrannus", PG_124="Veniliornis", PG_125="Vireo",
  PG_126="Volatinia", PG_127="Xenops", PG_128="Xiphorhynchus",
  PG_129="Zenaida", PG_130="Zonotrichia"
)

pg_to_genus_anurans <- c(
  PG_1="Adenomera", PG_2="Aplastodiscus", PG_3="Aquarana",
  PG_4="Boana", PG_5="Bokermannohyla", PG_6="Brachycephalus",
  PG_7="Dendropsophus", PG_8="Elachistocleis", PG_9="Haddadus",
  PG_10="Hylodes", PG_11="Ischnocnema", PG_12="Leptodactylus",
  PG_13="Odontophrynus", PG_14="Phyllomedusa", PG_15="Physalaemus",
  PG_16="Proceratophrys", PG_17="Rhinella", PG_18="Scinax",
  PG_19="Vitreorana"
)

# ============================================================
# COLOUR PALETTE — Wong (2011) colourblind-friendly
# ============================================================

col_splits <- "#000000"    # black      — split density
col_data   <- "#E69F00"    # orange     — data density + rug
col_ratio  <- "#0072B2"    # blue       — threshold line + IQR band
col_ratio2 <- "#D55E00"    # vermillion — normalised ratio curve
col_shade  <- "#0072B230"  # blue 19% opacity — bootstrap IQR band

pg_cols_base <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999",
  "#44AA99", "#882255", "#DDCC77", "#332288"
)

predictor_cols <- c(
  "Forest"             = "#2D6A4F",
  "Pioneer vegetation" = "#95D5B2",
  "Eucalyptus"         = "#74C69D",
  "Agriculture"        = "#F4A261",
  "Urban areas"        = "#E63946",
  "Bodies of water"    = "#4CC9F0",
  "Wetlands"           = "#4361EE",
  "Heterogeneity"      = "#9B5DE5"
)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

clean_pred_name <- function(pred) gsub("\\.", " ", pred)

pg_label <- function(pg_codes, lookup) {
  labels <- lookup[pg_codes]
  labels[is.na(labels)] <- pg_codes[is.na(labels)]
  labels
}

get_driving_PGs <- function(obj, pred, threshold, window = 5, top_n = 7) {
  pred_data      <- obj$res[obj$res$var == pred, ]
  near_threshold <- pred_data[abs(pred_data$split - threshold) <= window, ]
  if (nrow(near_threshold) == 0) {
    cat("No splits found near threshold for", pred, "\n")
    return(NULL)
  }
  pg_imp <- aggregate(improve.norm ~ spec, data = near_threshold, FUN = sum)
  pg_imp <- pg_imp[order(-pg_imp$improve.norm), ]
  head(pg_imp, top_n)
}

get_driving_PGs_auto <- function(obj, pred, threshold, top_n = 7) {
  pred_range <- range(obj$X[[pred]], na.rm = TRUE)
  window     <- (pred_range[2] - pred_range[1]) * 0.15
  pg_df      <- get_driving_PGs(obj, pred, threshold,
                                window = window, top_n = top_n)
  if (is.null(pg_df)) return(character(0))
  pg_codes <- as.character(pg_df$spec)
  rsq_vals <- obj$imp.rsq[pred, pg_codes]
  pg_codes[!is.na(rsq_vals)]
}

getCumulativeSplits <- function(res, pred, dens) {
  resA <- res[res$var == pred, ]
  tapply(1:nrow(resA), resA$spec, function(i) {
    split  <- resA$split[i]
    height <- pmax(0, resA$improve.norm[i]) *
      approx(dens[[pred]]$x, 1/dens[[pred]]$y, split, rule = 2)$y
    height <- height[order(split)]
    list(x = sort(split), y = height)
  })
}

plot_splits_density <- function(obj, pred, thresh_row) {
  res_pred <- obj$res[obj$res$var == pred, ]
  w        <- pmax(0, res_pred$improve.norm)
  if (sum(w) == 0) w <- rep(1, length(w))
  w        <- w / sum(w)
  d_splits <- density(res_pred$split, weights = w, n = 512)
  obs_vals <- obj$X[[pred]]
  d_data   <- density(obs_vals, n = 512,
                      from = min(d_splits$x), to = max(d_splits$x))
  ratio_y  <- d_splits$y / (d_data$y + 1e-10)
  ratio_y  <- ratio_y / max(ratio_y, na.rm = TRUE)
  xrange   <- range(d_splits$x)
  ymax     <- max(d_splits$y) * 1.2
  
  plot(xrange, c(0, ymax), type = "n",
       xlab = "", ylab = "", cex.axis = 1.3, las = 1, bty = "l")
  
  if (!is.null(thresh_row) && nrow(thresh_row) > 0) {
    is_unstable <- isTRUE(thresh_row$unstable)
    rect(thresh_row$IQR_lower, 0, thresh_row$IQR_upper, ymax,
         col = col_shade, border = NA)
    abline(v = thresh_row$median, col = col_ratio, lwd = 2,
           lty = if (is_unstable) 2 else 1)
  }
  
  lines(d_splits$x, ratio_y * ymax * 0.55, col = col_ratio2, lwd = 1.5, lty = 2)
  lines(d_data$x, d_data$y / max(d_data$y) * ymax * 0.35, col = col_data, lwd = 1.5)
  lines(d_splits$x, d_splits$y, col = col_splits, lwd = 2)
  rug(obs_vals, col = col_data, ticksize = 0.03, lwd = 0.5)
  mtext(clean_pred_name(pred), side = 3, line = 0.4, cex = 1.1, font = 2)
}

plot_cumulative_curves <- function(obj, pred, pg_list, lookup) {
  cumsplits <- getCumulativeSplits(obj$res, pred, obj$dens)
  CU        <- cumimp(obj, pred)
  ymax      <- max(CU$y) * 1.2
  
  available <- unique(obj$res$spec[obj$res$var == pred])
  pg_list   <- pg_list[pg_list %in% available]
  rsq_check <- obj$imp.rsq[pred, pg_list]
  pg_list   <- pg_list[!is.na(rsq_check)]
  
  genus_labels <- pg_label(pg_list, lookup)
  n_pg         <- length(pg_list)
  cols         <- setNames(pg_cols_base[seq_len(n_pg)], pg_list)
  
  plot(range(CU$x), c(0, ymax), type = "n",
       xlab = "", ylab = "", cex.axis = 1.3, las = 1, bty = "l")
  
  for (j in seq_along(pg_list)) {
    pg <- pg_list[j]
    xy <- cumsplits[[pg]]
    if (is.null(xy) || length(xy$x) == 0) next
    isub <- seq(1, length(xy$x), len = pmin(500, length(xy$x)))
    lines(xy$x[isub],
          (cumsum(xy$y) / sum(xy$y, na.rm = TRUE) *
             obj$imp.rsq[pred, pg])[isub],
          type = "s", col = cols[pg], lwd = 1.8)
  }
  
  isub <- seq(1, length(CU$x), len = pmin(500, length(CU$x)))
  lines(CU$x[isub], CU$y[isub], type = "s", lwd = 2.5, col = "black")
  
  left_half <- CU$x <= median(CU$x)
  left_ymax <- if (any(left_half)) max(CU$y[left_half]) else 0
  leg_pos   <- if (left_ymax < ymax * 0.45) "topleft" else "bottomright"
  
  legend(leg_pos,
         legend = c(genus_labels, "Overall"),
         col    = c(cols, "black"),
         lwd    = c(rep(1.8, n_pg), 2.5),
         bty    = "n", cex = 1.0, seg.len = 1.2)
}

# ============================================================
# MAIN ANALYSIS LOOP
# ============================================================

imp_all <- list()   # stores importance summaries for Figure 2

for (taxon in names(taxa_config)) {
  
  cfg <- taxa_config[[taxon]]
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("RUNNING:", cfg$label, "\n")
  cat(rep("=", 60), "\n", sep = "")
  
  if (rerun_analysis) {
    
    # ---- Data loading ----------------------------------------
    comm  <- read.table(cfg$comm_file,  header = TRUE,
                        row.names = 1, sep = "\t")
    envir <- read.table(cfg$envir_file, header = TRUE,
                        row.names = 1, sep = "\t")
    envir <- envir[, -cfg$coord_cols]
    rownames(envir) <- rownames(comm)
    
    # ---- Remove constant PGs and convert to factors ----------
    comm_filt <- comm[, apply(comm, 2, function(x) length(unique(x)) > 1)]
    cat("PGs after removing constants:", ncol(comm_filt),
        "(removed:", ncol(comm) - ncol(comm_filt), ")\n")
    cat("Removed:",
        paste(colnames(comm)[!colnames(comm) %in% colnames(comm_filt)],
              collapse = ", "), "\n")
    comm_filt[] <- lapply(comm_filt, factor)
    comm        <- comm_filt
    
    # ---- Ntree stability check -------------------------------
    cat("\nNtree stability check:\n")
    ntree_vals    <- c(500, 1000, 2000, 5000)
    rsq_stability <- sapply(ntree_vals, function(nt) {
      suppressWarnings(
        gf_test <- gradientForest(
          cbind(envir, comm),
          predictor.vars = colnames(envir),
          response.vars  = colnames(comm),
          ntree = nt, maxLevel = 5, corr.threshold = 0.5
        )
      )
      c(mean_R2  = mean(gf_test$result, na.rm = TRUE),
        top_imp  = importance(gf_test)[1],
        n_pos_R2 = length(gf_test$result))
    })
    colnames(rsq_stability) <- ntree_vals
    print(round(rsq_stability, 4))
    
    # ---- Full model ------------------------------------------
    cat("\nFitting full model (ntree =", cfg$ntree_full, ")...\n")
    gf_model <- gradientForest(
      cbind(envir, comm),
      predictor.vars = colnames(envir),
      response.vars  = colnames(comm),
      ntree          = cfg$ntree_full,
      maxLevel       = 5,
      corr.threshold = 0.5
    )
    assign(paste0("gf_model_", taxon), gf_model)
    
    # ---- Bootstrap 1: threshold locations -------------------
    cat("\nBootstrap 1 — threshold locations (n = 100)\n")
    n_boot          <- 100
    thresholds_boot <- list()
    
    for (i in 1:n_boot) {
      boot_idx  <- sample(1:nrow(envir), replace = TRUE)
      boot_pred <- envir[boot_idx, ]
      boot_comm <- comm[boot_idx, ]
      rownames(boot_comm) <- rownames(boot_pred)
      
      valid_pred <- colnames(boot_pred)[
        apply(boot_pred, 2, function(x) var(x, na.rm = TRUE) > 0)]
      valid_comm <- colnames(boot_comm)[
        sapply(boot_comm, function(x) length(unique(x)) > 1)]
      
      if (length(valid_pred) > 1 && length(valid_comm) > 1) {
        suppressWarnings(
          boot_gf <- gradientForest(
            cbind(boot_pred[, valid_pred, drop = FALSE],
                  boot_comm[, valid_comm, drop = FALSE]),
            predictor.vars = valid_pred,
            response.vars  = valid_comm,
            ntree = cfg$ntree_boot, maxLevel = 5
          )
        )
        for (pred in valid_pred) {
          ps <- boot_gf$res[boot_gf$res$var == pred, ]
          if (nrow(ps) > 0)
            thresholds_boot[[pred]][i] <- ps$split[which.max(ps$improve.norm)]
        }
      }
      if (i %% 10 == 0) cat("  Threshold replicate", i, "done\n")
    }
    
    threshold_summary <- do.call(rbind, lapply(names(thresholds_boot), function(pred) {
      vals <- thresholds_boot[[pred]][!is.na(thresholds_boot[[pred]])]
      data.frame(predictor = pred,
                 median    = median(vals),
                 IQR_lower = quantile(vals, 0.25),
                 IQR_upper = quantile(vals, 0.75),
                 IQR_width = IQR(vals),
                 unstable  = IQR(vals) > 10)
    }))
    cat("\nThreshold summary:\n")
    print(threshold_summary)
    
    # ---- Bootstrap 2: variable importance -------------------
    cat("\nBootstrap 2 — variable importance (n = 100)\n")
    imp_boot <- matrix(NA, nrow = n_boot, ncol = ncol(envir),
                       dimnames = list(NULL, colnames(envir)))
    
    for (i in 1:n_boot) {
      boot_idx  <- sample(1:nrow(envir), replace = TRUE)
      boot_pred <- envir[boot_idx, ]
      boot_comm <- comm[boot_idx, ]
      rownames(boot_comm) <- rownames(boot_pred)
      
      valid_pred <- colnames(boot_pred)[
        apply(boot_pred, 2, function(x) var(x, na.rm = TRUE) > 0)]
      valid_comm <- colnames(boot_comm)[
        sapply(boot_comm, function(x) length(unique(x)) > 1)]
      
      if (length(valid_pred) > 1 && length(valid_comm) > 1) {
        suppressWarnings(
          boot_gf <- gradientForest(
            cbind(boot_pred[, valid_pred, drop = FALSE],
                  boot_comm[, valid_comm, drop = FALSE]),
            predictor.vars = valid_pred,
            response.vars  = valid_comm,
            ntree = cfg$ntree_boot, maxLevel = 5
          )
        )
        imp_vals <- importance(boot_gf)
        imp_boot[i, names(imp_vals)] <- imp_vals
      }
      if (i %% 10 == 0) cat("  Importance replicate", i, "done\n")
    }
    
    imp_summary <- data.frame(
      predictor = colnames(imp_boot),
      median    = apply(imp_boot, 2, median,   na.rm = TRUE),
      IQR_lower = apply(imp_boot, 2, quantile, 0.25, na.rm = TRUE),
      IQR_upper = apply(imp_boot, 2, quantile, 0.75, na.rm = TRUE)
    )
    imp_summary <- imp_summary[order(-imp_summary$median), ]
    cat("\nImportance summary:\n")
    print(imp_summary)
    
    imp_summary$taxon <- cfg$label
    imp_all[[taxon]]  <- imp_summary
    
    write.csv(imp_summary,       cfg$output_imp, row.names = FALSE)
    write.csv(threshold_summary, cfg$output_thr, row.names = FALSE)
    cat("Saved:", cfg$output_imp, "\n")
    cat("Saved:", cfg$output_thr, "\n")
    
    # ---- R² per PG ------------------------------------------
    rsq_df <- data.frame(PG = names(gf_model$result),
                         R2 = gf_model$result)
    rsq_df <- rsq_df[order(-rsq_df$R2), ]
    cat("\nR\u00B2 per PG:\n")
    print(rsq_df)
    cat("Total PGs retained:", nrow(rsq_df), "\n")
    cat("Mean R\u00B2:", round(mean(rsq_df$R2), 3), "\n")
    cat("Range:", round(min(rsq_df$R2), 3),
        "to", round(max(rsq_df$R2), 3), "\n")
    
  } else {
    
    # ---- Load saved results (skip bootstraps) ----------------
    cat("Loading saved results (rerun_analysis = FALSE)\n")
    gf_model          <- get(paste0("gf_model_", taxon))
    imp_summary       <- read.csv(cfg$output_imp)
    threshold_summary <- read.csv(cfg$output_thr)
    imp_summary$taxon <- cfg$label
    imp_all[[taxon]]  <- imp_summary
    cat("Importance summary:\n")
    print(imp_summary[, c("predictor", "median")])
    cat("Threshold summary:\n")
    print(threshold_summary[, c("predictor", "median", "unstable")])
    
  }
  
  # ---- Driving PGs per threshold ----------------------------
  cat("\nDriving PGs per threshold:\n")
  for (pred in threshold_summary$predictor) {
    thresh_val <- threshold_summary[
      threshold_summary$predictor == pred, "median"]
    cat("\n--- Top PGs for", pred,
        "at", round(thresh_val, 1), "---\n")
    print(get_driving_PGs(gf_model, pred, threshold = thresh_val))
  }
  
  # ---- Select PG lookup table and top 3 predictors ----------
  pg_lookup      <- if (taxon == "birds") pg_to_genus_birds else pg_to_genus_anurans
  most_important <- imp_summary$predictor[1:3]
  
  imp_species <- setNames(
    lapply(seq_len(3), function(k) {
      pred       <- most_important[k]
      thresh_val <- threshold_summary[
        threshold_summary$predictor == pred, "median"]
      get_driving_PGs_auto(gf_model, pred,
                           threshold = thresh_val, top_n = 7)
    }),
    c("pred1", "pred2", "pred3")
  )
  cat("\nDriving PGs for top 3 predictors:\n")
  print(imp_species)
  
  assign(paste0("threshold_summary_", taxon), threshold_summary)
  assign(paste0("imp_species_",        taxon), imp_species)
  assign(paste0("most_important_",     taxon), most_important)
  assign(paste0("pg_lookup_",          taxon), pg_lookup)
  
  # ---- Panel plot (Figure 3a / 3b) --------------------------
  png(cfg$output_fig, width = 13, height = 9,
      units = "in", res = 300, bg = "white")
  
  par(mfrow = c(2, 3),
      mgp   = c(3.2, 0.8, 0),
      mar   = c(2.5, 5.5, 3.2, 1.2),
      oma   = c(0, 0, 2.5, 0))
  
  for (k in 1:3) {
    pred       <- most_important[k]
    thresh_row <- threshold_summary[threshold_summary$predictor == pred, ]
    plot_splits_density(gf_model, pred, thresh_row)
    if (k == 1) {
      mtext("Splits & data density", side = 2,
            line = 4.0, cex = 1.0, font = 2)
      legend("topright",
             legend = c("Split density", "Data density (scaled)",
                        "Ratio (norm.)", "Stable threshold",
                        "Unstable threshold", "Bootstrap IQR"),
             col    = c(col_splits, col_data, col_ratio2,
                        col_ratio, col_ratio, col_shade),
             lwd    = c(2, 1.5, 1.5, 2, 2, 8),
             lty    = c(1, 1, 2, 1, 2, 1),
             bty    = "n", cex = 0.95)
    }
  }
  
  for (k in 1:3) {
    pred       <- most_important[k]
    thresh_row <- threshold_summary[threshold_summary$predictor == pred, ]
    pg_list    <- imp_species[[paste0("pred", k)]]
    plot_cumulative_curves(gf_model, pred, pg_list, pg_lookup)
    if (k == 1)
      mtext("Cumulative importance (R\u00B2)", side = 2,
            line = 4.0, cex = 1.0, font = 2)
  }
  
  mtext(cfg$label, outer = TRUE, side = 3, line = 0.8,
        font = 2, cex = 1.4, family = "serif")
  dev.off()
  cat("Panel plot saved to", cfg$output_fig, "\n")
}

# ============================================================
# FIGURE 2 — VARIABLE IMPORTANCE
# Separate panels per taxon; Heterogeneity excluded for birds
# by VIF procedure (see 01_vif_analysis.R)
# ============================================================

if (!exists("imp_all") || length(imp_all) < 2) {
  stop("imp_all not found. Reload:\n",
       "  imp_birds         <- read.csv('importance_birds.csv')\n",
       "  imp_anurans       <- read.csv('importance_anurans.csv')\n",
       "  imp_birds$taxon   <- 'Birds'\n",
       "  imp_anurans$taxon <- 'Anurans'\n",
       "  imp_all <- list(birds = imp_birds, anurans = imp_anurans)")
}

cat("\nGenerating Figure 2\n")

imp_birds_df             <- imp_all[["birds"]]
imp_anurans_df           <- imp_all[["anurans"]]
imp_birds_df$predictor   <- gsub("\\.", " ", imp_birds_df$predictor)
imp_anurans_df$predictor <- gsub("\\.", " ", imp_anurans_df$predictor)

imp_birds_df   <- imp_birds_df[order(imp_birds_df$median), ]
imp_anurans_df <- imp_anurans_df[order(imp_anurans_df$median), ]

get_cols  <- function(preds) {
  cols <- predictor_cols[preds]
  cols[is.na(cols)] <- "#888888"
  cols
}
get_ranks <- function(vals) rank(-vals, ties.method = "min")

xmax <- max(imp_birds_df$median, imp_anurans_df$median,
            na.rm = TRUE) * 1.22

plot_importance_panel <- function(df, title, show_yaxis) {
  preds <- df$predictor
  vals  <- df$median
  n     <- nrow(df)
  cols  <- get_cols(preds)
  ranks <- get_ranks(vals)
  y_pos <- seq_len(n)
  bar_h <- 0.38
  
  plot(NULL, xlim = c(0, xmax), ylim = c(0.5, n + 0.5),
       xlab = "", ylab = "", axes = FALSE, xaxs = "i", yaxs = "i")
  
  for (i in seq_len(n)) {
    rect(0, y_pos[i] - bar_h, vals[i], y_pos[i] + bar_h,
         col = cols[i], border = NA)
    text(vals[i] + xmax * 0.015, y_pos[i],
         labels = paste0("#", ranks[i]),
         adj = 0, cex = 0.78, col = "grey30", family = "serif")
  }
  
  axis(1, cex.axis = 0.88, family = "serif")
  mtext("Variable importance", side = 1, line = 3.0,
        cex = 0.92, family = "serif")
  
  if (show_yaxis)
    axis(2, at = y_pos, labels = preds,
         las = 1, cex.axis = 0.92, tick = FALSE,
         family = "serif")
  
  mtext(title, side = 3, line = 0.5,
        font = 2, cex = 1.15, family = "serif")
  
  abline(v = mean(vals, na.rm = TRUE),
         col = "grey60", lty = 2, lwd = 0.9)
  box(bty = "l")
}

left_in <- 1.60
top_in  <- 0.50
bot_in  <- 0.90
gap_in  <- 0.20
dev_w   <- 13.0
dev_h   <- 5.8

panel_w <- (dev_w - 2 * left_in - gap_in) / 2
l1 <- 0
r1 <- (left_in + panel_w) / dev_w
l2 <- r1 + gap_in / dev_w
r2 <- 1.0

fig2_path <- file.path(output_dir, "Figure2_importance.png")
png(fig2_path, width = dev_w, height = dev_h,
    units = "in", res = 300, bg = "white")

par(fig = c(l1, r1, 0, 1), new = FALSE,
    mai = c(bot_in, left_in, top_in, 0.05),
    mgp = c(2.8, 0.65, 0), family = "serif")
plot_importance_panel(imp_birds_df,   "Birds",   show_yaxis = TRUE)

par(fig = c(l2, r2, 0, 1), new = TRUE,
    mai = c(bot_in, left_in, top_in, 0.05),
    mgp = c(2.8, 0.65, 0), family = "serif")
plot_importance_panel(imp_anurans_df, "Anurans", show_yaxis = TRUE)

dev.off()
cat("Figure 2 saved to", fig2_path, "\n")

# ============================================================
# SUMMARY
# ============================================================

cat("\n", rep("=", 60), "\n", sep = "")
cat("All outputs saved to:", output_dir, "\n")
cat("  Figure2_importance.png\n")
cat("  Figure3a_birds_GF.png\n")
cat("  Figure3b_anurans_GF.png\n")
cat("  importance_birds.csv\n")
cat("  importance_anurans.csv\n")
cat("  thresholds_birds.csv\n")
cat("  thresholds_anurans.csv\n")
cat(rep("=", 60), "\n", sep = "")