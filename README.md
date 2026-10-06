# DNA Typewriter Lineage Tree Reconstruction & FateVec Analysis

Reconstruct cell lineage trees from DNA Typewriter TAPE barcode data and analyze fate-bias acquisition dynamics.

## Overview

### Step 1: Tree Reconstruction (R)

1. **Filters cells** by minimum number of edited TAPEs recovered
2. **Builds a lineage tree** via UPGMA hierarchical clustering on a pairwise distance matrix
3. **Assigns clade labels** using `cutree` and annotates cell metadata
4. **Visualizes** the tree as a circular (fan) layout colored by clade
5. **Exports** tree structure (tips, edges, node-tip mappings) to CSV

### Step 2: FateVec Analysis (Python)

1. **Reconstructs** the tree from Step 1 CSV exports in Python
2. **Propagates fate vectors** bottom-up through the tree
3. **Traces** each tip's fate-vector component from root to leaf over normalized tree depth
4. **Smooths and differentiates** the mean curves to estimate **T_on**, the leading half-maximum crossing of each velocity profile
5. **Exports** onset positions (T_on) to CSV

## Repository Structure

```
mouse-lineage-recording/
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
│   ├── analysis.py             # Fate vector tracing & T_on estimation
│   └── plotting.py             # Visualization functions
├── run_fatevec.py              # Step 2 entry point
├── data/                       # Merged E8.5 cell annotations
├── outputs/                    # Reference tree and analysis outputs
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
python -m pip install numpy pandas scipy matplotlib
```

### External dependency

The supporting R helper **`extra_script_choi_brief.R`** from [Choi et al.](https://doi.org/10.1038/s41586-022-04922-8) is included in `external/`. Tree reconstruction reads a precomputed distance matrix.

> The included E8.5 example can be analyzed directly from the supplied tree CSVs.

## Input Data Format

Place your input files in the `data/` directory. The merged E8.5 annotation is included; the distance matrix and TAPE tables below are needed only for rebuilding a tree.

### 1. Distance Matrix CSV (`Step7_E8_DM_v2.csv`)

Square, symmetric matrix of pairwise cell distances (no row/column headers), ordered to match the TAPE pivot table.

### 2. TAPE Barcode Pivot Table (`Step6_E8_v2_TapeBCpivot.csv`)

Wide-format CSV with a `Cell` column. Each row is a cell; columns represent TAPE barcodes.

### 3. Long-format TAPE Selection File (`Step5_E8_v2.csv`)

One selected library per merged cell. Must contain columns `Cell`, `TapeBC` and `Site1`. Rows with `Site1 != "None"` count as edited TAPEs.

### 4. Cell Annotation CSV (`E8_annotations.csv`)

Cell metadata with `Cell` and `subcluster` columns. `Cell` is the unique 16-base barcode.

| Cell | subcluster |
|------|------------|
| AAACCCAAGGGCAATC | Surface ectoderm |

## Usage

### Step 1: Build Lineage Tree

By default, this step re-exports the included merged E8.5 tree. To rebuild from your own TAPE inputs, set `LINEAGE_REBUILD_TREE=true`.

1. Edit `R/00_config.R` to set your data file paths and parameters.

2. Run:

```bash
Rscript run_pipeline.R
```

This generates tree CSVs and visualizations in `outputs/`.

### Step 2: FateVec Analysis

The supplied tree exports let you run this step directly. The merged E8.5 noExVE example contains **10,579 tips and 19 cell types**.

1. Edit the configuration section in `run_fatevec.py`:

```python
TREE_NODES_SUMMARY = BASE_DIR / "outputs/tree_nodes_summary.csv"
CELL_ANNOT_CSV     = BASE_DIR / "outputs/cell_annot_with_clades.csv"
CELL_TYPE_COL     = "subcluster"
T_POINTS          = 1001       # Normalized tree-depth grid, 0 to 1
WINDOW_LENGTH     = 301        # Savitzky-Golay window (30% of the grid)
POLYORDER         = 4          # Polynomial order for smoothing
```

2. Run:

```bash
python run_fatevec.py
```

T_on is reported on normalized tree depth (root = 0, tips = 1).

## Outputs

### Step 1 (R)

| File | Description |
|------|-------------|
| `E8_tree.nwk` | Newick format tree |
| `cell_annot_with_clades.csv` | Cell annotations with `clone_id` column |
| `tree_edges.csv` | Edge list with branch lengths |
| `tree_nodes_summary.csv` | Node summary (depth, clade size, parent) |

### Step 2 (Python)

| File | Description |
|------|-------------|
| `fatevec_curves.pdf` | Mean FateVec component curves over normalized tree depth |
| `fatevec_derivatives.pdf` | Smoothed velocity curves with T_on markers |
| `fatevec_onset_summary.csv` | T_on and maximum velocity for each cell type |
