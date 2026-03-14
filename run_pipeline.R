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
#   - extra_script_choi_brief.R (Choi et al.) for computing normalized distance
#     matrices. This script is NOT included; see README.md for details.
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
# Uncomment if you need to recompute the distance matrix from raw TAPE data:
# source(CHOI_SCRIPT_PATH)

# ==============================================================================
# Step 1: Filter cells by tape recovery
# ==============================================================================
message("=== Step 1: Filtering cells ===")
filter_result <- filter_cells(
  tape_file  = TAPE_FILE,
  long_paths = TAPE_LONG_FILES,
  annot_file = ANNOT_FILE,
  min_tapes  = MIN_TAPES
)
cell_list          <- filter_result$cell_list
cell_list_filtered <- filter_result$cell_list_filtered
cell_annot_filtered <- filter_result$cell_annot_filtered

# ==============================================================================
# Step 2: Build lineage tree
# ==============================================================================
message("\n=== Step 2: Building tree ===")
tree <- build_tree(
  dm_file            = DM_FILE,
  cell_list          = cell_list,
  cell_list_filtered = cell_list_filtered,
  method             = HCLUST_METHOD
)

# Save Newick format
nwk_path <- file.path(OUTPUT_DIR, "E8_tree.nwk")
ape::write.tree(tree, file = nwk_path)
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
print(p)

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
