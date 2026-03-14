# ==============================================================================
# 00_config.R - Configuration for DNA Typewriter lineage tree reconstruction
# ==============================================================================
#
# External dependency:
#   This pipeline requires `extra_script_choi_brief.R` from Choi et al.
#   for computing the normalized distance matrix from TAPE barcodes.
#   Please obtain this script from the original authors and place it
#   at the path specified by CHOI_SCRIPT_PATH below.
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

DM_FILE    <- file.path(DATA_DIR, "Step7_E8_DM_v250819_normalized.csv")
TAPE_FILE  <- file.path(DATA_DIR, "Step6_E8_TapeBCpivot_v250819.csv")
ANNOT_FILE <- file.path(DATA_DIR, "E8_annotations_v6.csv")

# Long-format tape selection files (one per library)
TAPE_LONG_FILES <- file.path(DATA_DIR, c(
  "Step5_E8L1_TapeSelect_v250819.csv",
  "Step5_E8L2_TapeSelect_v250819.csv",
  "Step5_E8L3_TapeSelect_v250819.csv"
))

# ---- Output directory ----
OUTPUT_DIR <- file.path(BASE_DIR, "outputs")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ---- Parameters ----
MIN_TAPES      <- 8       # Minimum number of edited TAPEs per cell
HCLUST_METHOD  <- "average"  # Hierarchical clustering method
K_CLADES       <- 50      # Number of clades for cutree

# ---- Color palette for clade visualization ----
# 16-color palette; extend or modify as needed for your K_CLADES setting
CLADE_PALETTE <- c(
  "1"  = "#BB5A3B", "2"  = "#BC3D46", "3"  = "#F4E480", "4"  = "#866159",
  "5"  = "#6785C0", "6"  = "#C15143", "7"  = "#B3B469", "8"  = "#E4B69E",
  "9"  = "#CDCA96", "10" = "#528185", "11" = "#3288bd", "12" = "#94C4A6",
  "13" = "#FAE7C7", "14" = "#A0E0CF", "15" = "#E6A96E", "16" = "#604344"
)
