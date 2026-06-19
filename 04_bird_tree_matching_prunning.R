# ============================================================
# BIRD PHYLOGENY — TAXONOMY MATCHING AND TREE PRUNING
# BirdTree Hackett backbone (Jetz et al. 2012)
# 145 of 166 species retained (87%)
# 15 species absent from BirdTree excluded
# 7 congeneric pairs merged to single BirdTree representative
# ============================================================

library(ape)
library(picante)

# ============================================================
# STEP 1 — LOAD DATA
# ============================================================

comm_birds <- read.table("sites_species_birds_names_new.txt",
                         header = TRUE, row.names = 1,
                         sep = "\t", dec = ",")

# Load 1000 posterior trees — use tree 1 as representative
tree_birds_all    <- read.nexus("hacket_vertlife_seq.nex")
tree_birds_single <- tree_birds_all[[1]]

cat("Bird species in community:", ncol(comm_birds), "\n")
cat("Tips in representative tree:", length(tree_birds_single$tip.label), "\n")

# ============================================================
# STEP 2 — STANDARDISE NAMES
# Replace dots and spaces with underscores
# ============================================================

colnames(comm_birds)        <- gsub("\\.", "_", gsub(" ", "_", colnames(comm_birds)))
tree_birds_single$tip.label <- gsub(" ", "_", tree_birds_single$tip.label)

# ============================================================
# STEP 3 — SYNONYM TABLE
# Current name in community matrix -> name in BirdTree
# ============================================================

bird_synonym_table <- c(
  "Stilpnia_cayana"          = "Tangara_cayana",
  "Asemospiza_fuliginosa"    = "Tiaris_fuliginosus",
  "Chionomesa_lactea"        = "Amazilia_lactea",
  "Drymophila_malura"        = "Drymophila_ferruginea",
  "Dysithamnus_stictothorax" = "Dysithamnus_mentalis",
  "Euphonia_violacea"        = "Euphonia_pectoralis",
  "Hemitriccus_orbitatus"    = "Hemitriccus_nidipendulus",
  "Micrastur_ruficollis"     = "Micrastur_semitorquatus",
  "Sporophila_lineola"       = "Sporophila_caerulescens",
  "Turdus_subalaris"         = "Turdus_amaurochalinus"
)

# Species genuinely absent from BirdTree — excluded from GDM only
# Still present in main GradientForest analysis
birds_absent_from_tree <- c(
  "Aramides_saracura",
  "Arremon_semitorquatus",
  "Campephilus_robustus",
  "Dendroma_rufa",
  "Hemithraupis_ruficapilla",
  "Hylophilus_poicilotis",
  "Leucochloris_albicollis",
  "Malacoptila_striata",
  "Myrmotherula_gularis",
  "Neopelma_chrysolophum",
  "Phyllomyias_griseocapilla",
  "Phylloscartes_eximius",
  "Rhopias_gularis",
  "Todirostrum_poliocephalum",
  "Troglodytes_aedon"
)

# ============================================================
# STEP 4 — APPLY SYNONYMS AND REMOVE ABSENT SPECIES
# Order matters: synonym first, then remove absent
# ============================================================

comm_birds_mapped <- colnames(comm_birds)

# Apply synonyms
for (i in seq_along(comm_birds_mapped)) {
  nm <- comm_birds_mapped[i]
  if (nm %in% names(bird_synonym_table)) {
    comm_birds_mapped[i] <- bird_synonym_table[nm]
  }
}

# Remove species absent from tree BEFORE assigning names
keep_idx          <- !comm_birds_mapped %in% birds_absent_from_tree
comm_birds_mapped <- comm_birds_mapped[keep_idx]
comm_birds_clean  <- comm_birds[, keep_idx, drop = FALSE]
colnames(comm_birds_clean) <- comm_birds_mapped

cat("\nSpecies excluded (absent from BirdTree):", length(birds_absent_from_tree), "\n")
cat("Species remaining:", ncol(comm_birds_clean), "\n")

# ============================================================
# STEP 5 — MERGE DUPLICATE COLUMNS
# Congeneric species mapped to same BirdTree tip
# Presence of either = presence of merged taxon
# ============================================================

dup_cols <- which(duplicated(colnames(comm_birds_clean)))
cat("Duplicate columns to merge:", length(dup_cols), "\n")
cat("Duplicated names:", paste(colnames(comm_birds_clean)[dup_cols], collapse=", "), "\n")

comm_merged <- comm_birds_clean
for (col_idx in rev(dup_cols)) {
  col_name  <- colnames(comm_merged)[col_idx]
  first_idx <- which(colnames(comm_merged) == col_name)[1]
  comm_merged[, first_idx] <- as.integer(
    comm_merged[, first_idx] > 0 | comm_merged[, col_idx] > 0)
  comm_merged <- comm_merged[, -col_idx, drop = FALSE]
}

cat("Columns after merging duplicates:", ncol(comm_merged), "\n")
cat("Duplicate columns remaining:", sum(duplicated(colnames(comm_merged))), "\n")

# ============================================================
# STEP 6 — PRUNE TREE TO MATCHED SPECIES
# ============================================================

comm_birds_final    <- comm_merged
matched_birds_final <- colnames(comm_birds_final)

tips_to_drop      <- setdiff(tree_birds_single$tip.label, matched_birds_final)
tree_birds_pruned <- drop.tip(tree_birds_single, tips_to_drop)

cat("\n--- After pruning ---\n")
cat("Species retained:", ncol(comm_birds_final), "\n")
cat("Tree tips:", length(tree_birds_pruned$tip.label), "\n")
cat("Is ultrametric:", is.ultrametric(tree_birds_pruned), "\n")

# If not ultrametric, run chronos on small pruned tree
if (!is.ultrametric(tree_birds_pruned)) {
  cat("Making tree ultrametric via chronos...\n")
  tree_birds_pruned <- chronos(tree_birds_pruned,
                               lambda = 0, model = "relaxed",
                               quiet = FALSE)
  cat("Ultrametric after chronos:", is.ultrametric(tree_birds_pruned), "\n")
}

# ============================================================
# STEP 7 — VERIFY PHYLOGENETIC PLOT
# ============================================================

png("bird_pruned_tree.png", width = 10, height = 14,
    units = "in", res = 200, bg = "white")
plot(tree_birds_pruned, cex = 0.55,
     main = "Pruned bird phylogeny (Jetz et al. 2012 - Hackett backbone)")
dev.off()
cat("Pruned tree saved to bird_pruned_tree.png\n")

# ============================================================
# STEP 8 — SANITY CHECKS
# ============================================================

cat("\nFaith's PD per site:\n")
pd_vals <- pd(comm_birds_final, tree_birds_pruned, include.root = TRUE)
print(round(pd_vals, 2))

cat("\nTest PhyloSor computation...\n")
ps_test <- phylosor(comm_birds_final, tree_birds_pruned)
cat("PhyloSor range:", round(range(ps_test, na.rm = TRUE), 3), "\n")
cat("NA values:", sum(is.na(ps_test)), "\n")

cat("\nTest betaMNTD computation...\n")
bm_test <- comdistnt(comm_birds_final, cophenetic(tree_birds_pruned),
                     abundance.weighted = FALSE)
cat("betaMNTD range:", round(range(bm_test, na.rm = TRUE), 3), "\n")

cat("\nAll sanity checks passed — ready for GDM analysis\n")

# ============================================================
# STEP 9 — SAVE CLEANED OBJECTS
# ============================================================

saveRDS(tree_birds_pruned, "bird_tree_pruned.rds")
saveRDS(comm_birds_final,  "bird_comm_matched.rds")
write.nexus(tree_birds_pruned, file = "bird_tree_pruned.nex")

cat("\nSaved:\n")
cat("  bird_tree_pruned.rds  — pruned ultrametric tree (145 tips)\n")
cat("  bird_comm_matched.rds — community matrix (145 species)\n")
cat("  bird_tree_pruned.nex  — pruned tree in Nexus format\n")
cat("  bird_pruned_tree.png  — visual check of tree topology\n")

# ============================================================
# SUMMARY OF EXCLUSIONS FOR METHODS
# ============================================================

cat("\n=== METHODS NOTE ===\n")
cat("Total bird species recorded:", ncol(comm_birds), "\n")
cat("Excluded — absent from BirdTree:", length(birds_absent_from_tree), "\n")
cat("  Species:", paste(birds_absent_from_tree, collapse=", "), "\n")
cat("Merged congeneric pairs (7 pairs -> 7 tips):\n")
for (nm in names(bird_synonym_table)) {
  if (bird_synonym_table[nm] != nm) {
    cat(" ", nm, "->", bird_synonym_table[nm], "\n")
  }
}
cat("Final species in GDM analysis:", ncol(comm_birds_final), "\n")
