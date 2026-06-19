# ============================================================
# 03_ANURAN TREE MATCHING AND PRUNING
#
# Matches anuran species names between the community matrix
# and the Jetz & Pyron (2018) amphibian phylogeny, prunes the
# tree to study species, and makes it ultrametric via chronos.
# Output is used as input for 05_gdm_phylogenetic_analysis.R.
#
# Inputs:
#   sites_species_anura_names.txt — species-level presence/absence
#     matrix (not the PG matrix; used for GDM only)
#   amph_shl_new_Consensus_7238.tre — Jetz & Pyron (2018) consensus
#     tree. Download from Dryad:
#     https://doi.org/10.5061/dryad.db265
#
# Outputs:
#   anuran_tree_pruned.rds  — pruned ultrametric phylogeny (36 tips)
#   anuran_comm_matched.rds — species matrix matched to tree tips
#   anuran_tree_pruned.nex  — pruned tree in Nexus format
#   anuran_pruned_tree.png  — visual check of tree topology
#
# Key decisions:
#   Haddadus binotatus excluded — absent from Jetz & Pyron tree
#     under any synonym including historical Eleutherodactylus
#   9 species synonymised to match tree nomenclature (see Step 2)
#   Tree pruned BEFORE chronos for speed (36 tips vs 7238)
#   chronos() with relaxed clock (lambda = 0) to make ultrametric
#
# Reference: Jetz, W. & Pyron, R.A. (2018) Nat. Ecol. Evol. 2, 850-858.
# ============================================================

library(ape)
if (requireNamespace("phytools", quietly = TRUE)) library(phytools)
library(picante)

# ============================================================
# STEP 1 — LOAD YOUR DATA
# ============================================================

comm <- read.table("sites_species_anura_names.txt",
                   header = TRUE, row.names = 1,
                   sep = "\t", dec = ",")

tree <- read.tree("amph_shl_new_Consensus_7238.tre")

cat("Species in community matrix:", ncol(comm), "\n")
cat("Tips in tree:", length(tree$tip.label), "\n")

# Remove outgroup
tree <- drop.tip(tree, "Homo_sapiens")
cat("Tips after removing outgroup:", length(tree$tip.label), "\n")

# Standardise tree tip labels
tree$tip.label <- gsub(" ", "_", tree$tip.label)

# ============================================================
# STEP 2 — SYNONYM TABLE
# Current name in community matrix -> name in Jetz & Pyron tree
# ============================================================

synonym_table <- c(
  "Aquarana_catesbeiana"        = "Rana_catesbeiana",
  "Boana_albopunctata"          = "Hypsiboas_albopunctatus",
  "Boana_bischoffi"             = "Hypsiboas_bischoffi",
  "Boana_faber"                 = "Hypsiboas_faber",
  "Boana_polytaenia"            = "Hypsiboas_polytaenius",
  "Boana_prasina"               = "Hypsiboas_prasinus",
  "Boana_semilineata"           = "Hypsiboas_semilineatus",
  "Ischnocnema_aff._guentheri"  = "Ischnocnema_guentheri",
  "Vitreorana_parvula"          = "Vitreorana_uranoscopa"
  # Add more here if any species remain unmatched after Step 3
)

# Species genuinely absent from Jetz & Pyron 2018 tree
# Excluded from GDM only — still present in GradientForest analysis
# Haddadus binotatus was historically placed in Eleutherodactylus
# but is not represented in the tree under any synonym
species_absent_from_tree <- c("Haddadus_binotatus")

# ============================================================
# STEP 3 — STANDARDISE NAMES, APPLY SYNONYMS, CHECK MATCHES
# ============================================================

# Standardise community matrix names (spaces to underscores)
colnames(comm) <- gsub(" ", "_", colnames(comm))

# Apply synonym table to community matrix names
comm_names_mapped <- colnames(comm)
for (i in seq_along(comm_names_mapped)) {
  nm <- comm_names_mapped[i]
  if (nm %in% names(synonym_table)) {
    comm_names_mapped[i] <- synonym_table[nm]
  }
}

# Remove species absent from tree BEFORE updating column names
# This keeps comm and comm_names_mapped the same length
keep_idx          <- !comm_names_mapped %in% species_absent_from_tree
comm_names_mapped <- comm_names_mapped[keep_idx]
comm              <- comm[, keep_idx, drop = FALSE]

cat("\nSpecies excluded from GDM (absent from tree):",
    species_absent_from_tree, "\n")
cat("Species remaining:", ncol(comm), "\n")

# Update column names to tree names
colnames(comm) <- comm_names_mapped

# Check matches
matched   <- intersect(comm_names_mapped, tree$tip.label)
unmatched <- setdiff(comm_names_mapped, tree$tip.label)

cat("\n--- Matching summary ---\n")
cat("Species in community:", length(comm_names_mapped), "\n")
cat("Matched to tree:", length(matched), "\n")
cat("Unmatched (need synonym):", length(unmatched), "\n")

if (length(unmatched) > 0) {
  cat("\nUnmatched species — search for synonyms:\n")
  print(unmatched)
  
  cat("\nPossible matches in tree (partial string search):\n")
  for (sp in unmatched) {
    genus   <- strsplit(sp, "_")[[1]][1]
    matches <- grep(paste0("^", genus, "_"), tree$tip.label, value = TRUE)
    cat(sp, "->",
        if (length(matches) > 0) paste(matches, collapse = ", ")
        else "NOT FOUND IN TREE", "\n")
  }
  cat("\nAdd missing synonyms to synonym_table above and rerun Steps 2-3\n")
  stop("Resolve unmatched species before proceeding to Step 4.")
}

cat("All species matched — proceeding to pruning\n")

# ============================================================
# STEP 4 — PRUNE TREE TO STUDY SPECIES
# Pruning before chronos is much faster:
# chronos on 36-tip tree takes seconds vs 20+ min on 7238 tips
# ============================================================

# Keep only matched species in community matrix
comm_matched <- comm[, matched, drop = FALSE]

# Prune tree to matched species only
tips_to_drop <- setdiff(tree$tip.label, matched)
tree_pruned  <- drop.tip(tree, tips_to_drop)

cat("\n--- After pruning ---\n")
cat("Species retained:", ncol(comm_matched), "\n")
cat("Tree tips:", length(tree_pruned$tip.label), "\n")
cat("Is ultrametric:", is.ultrametric(tree_pruned), "\n")

# ============================================================
# STEP 5 — MAKE ULTRAMETRIC VIA CHRONOS (fast on small tree)
# ============================================================

if (!is.ultrametric(tree_pruned)) {
  cat("\nMaking tree ultrametric via chronos (relaxed clock)...\n")
  tree_pruned <- chronos(tree_pruned,
                         lambda = 0,
                         model  = "relaxed",
                         quiet  = FALSE)
  cat("Ultrametric after chronos:", is.ultrametric(tree_pruned), "\n")
} else {
  cat("Tree already ultrametric — skipping chronos\n")
}

# ============================================================
# STEP 6 — VERIFY WITH PHYLOGENETIC PLOT
# ============================================================

png("anuran_pruned_tree.png", width = 8, height = 10,
    units = "in", res = 200, bg = "white")
plot(tree_pruned, cex = 0.85,
     main = "Pruned anuran phylogeny (Jetz & Pyron 2018)")
dev.off()
cat("Pruned tree saved to anuran_pruned_tree.png\n")

# ============================================================
# STEP 7 — SANITY CHECKS
# ============================================================

cat("\nFaith's PD per site:\n")
pd_vals <- pd(comm_matched, tree_pruned, include.root = TRUE)
print(round(pd_vals, 2))

cat("\nTest PhyloSor computation...\n")
ps_test <- phylosor(comm_matched, tree_pruned)
cat("PhyloSor range:", round(range(ps_test, na.rm = TRUE), 3), "\n")
cat("NA values:", sum(is.na(ps_test)), "\n")

cat("\nTest betaMNTD computation...\n")
bm_test <- comdistnt(comm_matched, cophenetic(tree_pruned),
                     abundance.weighted = FALSE)
cat("betaMNTD range:", round(range(bm_test, na.rm = TRUE), 3), "\n")

cat("\nAll sanity checks passed — ready for GDM analysis\n")

# ============================================================
# STEP 8 — SAVE CLEANED OBJECTS
# ============================================================

saveRDS(tree_pruned,  "anuran_tree_pruned.rds")
saveRDS(comm_matched, "anuran_comm_matched.rds")
write.nexus(tree_pruned, file = "anuran_tree_pruned.nex")

cat("\nSaved:\n")
cat("  anuran_tree_pruned.rds  — pruned ultrametric tree\n")
cat("  anuran_comm_matched.rds — community matrix (matched species only)\n")
cat("  anuran_tree_pruned.nex  — pruned tree in Nexus format\n")
cat("  anuran_pruned_tree.png  — visual check of tree topology\n")
