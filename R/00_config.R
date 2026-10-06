# ==============================================================================
# 00_config.R - Configuration for DNA Typewriter lineage tree reconstruction
# ==============================================================================
#
# External dependency:
#   The Choi et al. helper accompanies the TAPE distance analysis.
#   This entry point accepts a precomputed distance matrix for rebuilding.
#   The helper is included at CHOI_SCRIPT_PATH; precomputed tree exports
#   do not require distance-matrix recomputation.
#
# ==============================================================================

library(ape)
library(dplyr)
library(readr)
library(ggtree)
library(ggnewscale)
library(phangorn)
library(tibble)

# ---- Project root (set this to your project directory) ----
BASE_DIR <- here::here()  # or set manually: "/path/to/dna-typewriter-tree"

# ---- External script path ----
CHOI_SCRIPT_PATH <- file.path(BASE_DIR, "external", "extra_script_choi_brief.R")

# ---- Input data paths (relative to BASE_DIR) ----
DATA_DIR <- file.path(BASE_DIR, "data")

DM_FILE    <- file.path(DATA_DIR, "Step7_E8_DM_v2.csv")
TAPE_FILE  <- file.path(DATA_DIR, "Step6_E8_v2_TapeBCpivot.csv")
ANNOT_FILE <- file.path(DATA_DIR, "E8_annotations.csv")

# TAPE v2 uses records from one selected library per merged cell.
# These optional raw inputs are not distributed with this repository.
TAPE_LONG_FILES <- file.path(DATA_DIR, "Step5_E8_v2.csv")

# By default, re-export the included reference tree. Set to TRUE only when
# rebuilding from the matching v2 distance matrix and TAPE tables above.
REBUILD_TREE <- identical(Sys.getenv("LINEAGE_REBUILD_TREE", "false"), "true")
TREE_FILE <- file.path(BASE_DIR, "outputs", "E8_tree.nwk")
EXCLUDED_CELL_TYPES <- "Extraembryonic visceral endoderm"

# ---- Output directory ----
OUTPUT_DIR <- Sys.getenv("LINEAGE_OUTPUT_DIR", file.path(BASE_DIR, "outputs"))
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ---- Parameters ----
MIN_TAPES      <- 8       # Minimum number of edited TAPEs per cell
HCLUST_METHOD  <- "average"  # Hierarchical clustering method
K_CLADES       <- 50      # Number of clades for cutree

# Automatically generate enough colors for all K_CLADES groups.
CLADE_PALETTE <- NULL
