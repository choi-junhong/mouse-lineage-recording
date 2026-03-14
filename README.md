# DNA Typewriter Lineage Tree Reconstruction & FateVec Analysis

Reconstruct cell lineage trees from DNA Typewriter TAPE barcode data and analyze cell-fate commitment dynamics.

## Overview

### Step 1: Tree Reconstruction (R)
1. **Filters cells** by minimum number of edited TAPEs recovered
2. **Builds a lineage tree** via UPGMA hierarchical clustering on a pairwise distance matrix
3. **Assigns clade labels** using `cutree` and annotates cell metadata
4. **Visualizes** the tree as a circular (fan) layout colored by clade
5. **Exports** tree structure (tips, edges, node-tip mappings) to CSV

### Step 2: FateVec Analysis (Python)
1. **Reconstructs** the tree from Step 1 CSV exports in Python
2. **Propagates fate vectors** bottom-up: each node gets a normalized cell-type proportion vector
3. **Traces** each tip's fate proportion from root to leaf over pseudo-time
4. **Smooths and differentiates** to find the time of maximum fate commitment for each cell type
5. **Exports** peak commitment times to CSV

## Repository Structure

```
dna-typewriter-tree/
├── R/                          # Step 1: Tree reconstruction
│   ├── 00_config.R             # Paths, parameters, color palette
│   ├── 01_filter_cells.R       # Cell filtering by TAPE threshold
│   ├── 02_build_tree.R         # Distance matrix → hclust → phylo
│   ├── 03_assign_clades.R      # cutree clade assignment
│   ├── 04_plot_tree.R          # Circular tree visualization (ggtree)
│   └── 05_export_tree.R        # Export tree structure to CSV
├── run_pipeline.R              # Step 1 entry point
├── fatevec/                    # Step 2: FateVec analysis
│   ├── __init__.py
│   ├── node.py                 # Node class with bottom-up fate vector propagation
│   ├── tree.py                 # Tree construction from CSV
│   ├── analysis.py             # Fate vector tracing & peak derivative finding
│   └── plotting.py             # Visualization functions
├── run_fatevec.py              # Step 2 entry point
├── data/                       # Input data (not tracked)
├── outputs/                    # Generated outputs (not tracked)
└── README.md
```

## Requirements

### R packages (Step 1)

```r
install.packages(c("ape", "phangorn", "dplyr", "readr", "tibble", "here", "scales"))

# ggtree (Bioconductor)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ggtree", "ggnewscale"))
```

### Python packages (Step 2)

```bash
pip install numpy pandas scipy matplotlib
```

### External dependency

This pipeline requires **`extra_script_choi_brief.R`** from [Choi et al.](https://doi.org/10.1038/s41586-022-04922-8) for computing the normalized distance matrix from raw TAPE barcodes. This script is **not included** in this repository. Place it at `external/extra_script_choi_brief.R`.

> If you already have a precomputed distance matrix CSV, this external script is not needed.

## Input Data Format

Place your input files in the `data/` directory.

### 1. Distance Matrix CSV (`Step7_*_DM_*.csv`)

Square, symmetric matrix of pairwise cell distances (no row/column headers).

### 2. TAPE Barcode Pivot Table (`Step6_*_TapeBCpivot_*.csv`)

Wide-format CSV with a `Cell` column. Each row is a cell; columns represent TAPE barcodes.

### 3. Long-format TAPE Selection Files (`Step5_*_TapeSelect_*.csv`)

One file per library. Must contain columns `Cell` and `Site1`. Rows with `Site1 != "None"` count as edited TAPEs.

### 4. Cell Annotation CSV (`*_annotations_*.csv`)

Cell metadata with at minimum `Cell` and `CellType` columns (CellType at column index 2).

| Cell | ... | CellType | ... |
|------|-----|----------|-----|
| cell_001 | ... | Neuron | ... |

## Usage

### Step 1: Build Lineage Tree

1. Edit `R/00_config.R` to set your data file paths and parameters.

2. Run:
```bash
Rscript run_pipeline.R
```

This generates tree CSVs and visualizations in `outputs/`.

### Step 2: FateVec Analysis

1. Edit the configuration section in `run_fatevec.py`:

```python
TREE_NODES_SUMMARY = "outputs/tree_nodes_summary.csv"
CELL_ANNOT_CSV     = "outputs/cell_annot_with_clades.csv"
CELL_TYPE_COL      = 2          # 0-based column index for CellType
WINDOW_LENGTH      = 3501       # Savitzky-Golay window (must be odd)
POLYORDER          = 4          # Polynomial order for smoothing
```

2. Run:
```bash
python run_fatevec.py
```

## Outputs

### Step 1 (R)

| File | Description |
|------|-------------|
| `E8_tree.nwk` | Newick format tree |
| `cell_annot_with_clades.csv` | Cell annotations with `clone_id` column |
| `tree_nodes_summary.csv` | Node summary (depth, clade size, parent) |

### Step 2 (Python)

| File | Description |
|------|-------------|
| `fatevec_curves.pdf` | Average fate proportion curves over pseudo-time |
| `fatevec_derivatives.pdf` | Smoothed dv/dt curves per cell type |
| `fatevec_peak_summary.csv` | Time and magnitude of max fate commitment rate per cell type |

