#!/usr/bin/env Rscript
# ==============================================================================
# run_pipeline.R - DNA Typewriter Lineage Tree Reconstruction Pipeline
# ==============================================================================
#
# Usage:
#   1. Edit R/00_config.R to set your data paths and parameters
#   2. Run this script: Rscript run_pipeline.R
#      Or source it interactively in RStudio
#
# Prerequisites:
#   - The included reference tree and merged annotations for the default run.
#   - Optional raw TAPE v2 inputs are needed only when REBUILD_TREE is enabled.
#   - Install R packages: ape, dplyr, readr, ggtree, ggnewscale, phangorn,
#     tibble, here, scales
#
# ==============================================================================

# ---- Source all modules ----
source(file.path("R", "00_config.R"))
source(file.path("R", "01_filter_cells.R"))
source(file.path("R", "02_build_tree.R"))
source(file.path("R", "03_assign_clades.R"))
source(file.path("R", "04_plot_tree.R"))
source(file.path("R", "05_export_tree.R"))

# ---- (Optional) Source external distance-matrix script ----
# Uncomment to load the accompanying analysis helpers:
# source(CHOI_SCRIPT_PATH)

# Use the included merged E8.5 noExVE tree by default. Rebuilding from
# raw TAPE data is optional and requires the matching one-library-per-cell v2 inputs.
if (REBUILD_TREE) {
  message("=== Filtering merged cells and rebuilding the tree ===")
  filter_result <- filter_cells(
    tape_file = TAPE_FILE, long_paths = TAPE_LONG_FILES,
    annot_file = ANNOT_FILE, min_tapes = MIN_TAPES,
    excluded_types = EXCLUDED_CELL_TYPES
  )
  cell_annot_filtered <- filter_result$cell_annot_filtered
  tree <- build_tree(
    dm_file = DM_FILE, cell_list = filter_result$cell_list,
    cell_list_filtered = filter_result$cell_list_filtered,
    method = HCLUST_METHOD
  )
} else {
  message("=== Loading the included merged E8.5 noExVE reference tree ===")
  tree <- ape::read.tree(TREE_FILE)
  annotations <- read.csv(ANNOT_FILE, stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(tree$tip.label), !anyDuplicated(annotations$Cell),
            all(tree$tip.label %in% annotations$Cell))
  cell_annot_filtered <- annotations[match(tree$tip.label, annotations$Cell), , drop = FALSE]
  stopifnot(!any(is.na(cell_annot_filtered$subcluster)),
            !any(cell_annot_filtered$subcluster %in% EXCLUDED_CELL_TYPES))
}

# Save Newick format
nwk_path <- file.path(OUTPUT_DIR, "E8_tree.nwk")
ape::write.tree(tree, file = nwk_path, digits = 17)
message(sprintf("Newick tree saved to: %s", nwk_path))

# ==============================================================================
# Step 3: Assign clades
# ==============================================================================
message("\n=== Step 3: Assigning clades ===")
clade_result <- assign_clades(tree, k = K_CLADES)
tr_grp       <- clade_result$tr_grp
clade_map    <- clade_result$clade_map

# Write clade IDs into cell annotation
cell_annot_filtered <- write_clade_into_meta(cell_annot_filtered, clade_map)
annot_out <- file.path(OUTPUT_DIR, "cell_annot_with_clades.csv")
readr::write_csv(cell_annot_filtered, annot_out)
message(sprintf("Annotated metadata saved to: %s", annot_out))

# ==============================================================================
# Step 4: Plot tree
# ==============================================================================
message("\n=== Step 4: Plotting tree ===")
p <- plot_fan_tree(tr_grp, clade_palette = CLADE_PALETTE)
if (interactive()) print(p)

# Save plot
plot_path <- file.path(OUTPUT_DIR, "E8_tree_fan.pdf")
ggplot2::ggsave(plot_path, p, width = 12, height = 12)
message(sprintf("Tree plot saved to: %s", plot_path))

# ==============================================================================
# Step 5: Export tree structure to CSV
# ==============================================================================
message("\n=== Step 5: Exporting tree structure ===")
export_tree_csv(tree, output_dir = OUTPUT_DIR)

message("\n=== Pipeline complete! ===")
