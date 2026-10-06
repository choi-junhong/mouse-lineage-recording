# DNA Typewriter Lineage Tree Reconstruction & FateVec Analysis

Reconstruct cell lineage trees from DNA Typewriter TAPE records and summarize fate-bias acquisition with **T_on**, the leading half-maximum crossing of each velocity profile.

## Included E8.5 reference

The E8.5 example uses the updated **merged scRNA-seq annotations and TAPE v2 tree**:

- `data/E8_annotations.csv`: final consensus annotations for **13,832 unique cells**, retaining the existing `Cell,subcluster` format.
- `outputs/E8_tree.nwk`: the **10,579-tip noExVE tree**, containing **19 cell types**. Cells have at least eight Site1-edited TAPEs and a final annotation; extraembryonic visceral endoderm (ExVE) was excluded **before clustering**. Amniotic ectoderm and allantois remain included.
- `outputs/cell_annot_with_clades.csv`: matching tree-tip annotations in the existing `Cell,subcluster,clone_id` format. The 50 clades are recomputed on the merged tree.
- `outputs/tree_nodes_summary.csv` and `outputs/tree_edges.csv`: matching node and edge tables with their original column layouts.

`Cell` is the unique 16-base barcode. The old library suffixes are removed. TAPE v2 retains all records from one selected library per cell, based on its GEX capture assignment; it does not combine site calls from different libraries. The default tree corresponds to `E8_tree_TAPE8_v2_noExVE`, with final `subcluster_consensus` labels.

The update replaces the existing annotation, tree, summary-table and PDF outputs. Raw sequencing reads, expression matrices, Seurat objects and additional datasets are not included.

## T_on definition

Let x be **normalized phylogenetic depth**: the root is 0 and the ultrametric tips are 1. For each cell type k, FateVec propagates count-weighted fate vectors through the tree and averages the corresponding component over that type's tip-to-root paths, giving V_k(x).

Normalize the change from its root-depth value:

```text
F_k(x) = [V_k(x) - V_k(0)] / [1 - V_k(0)]
velocity_k(x) = dF_k(x) / dx
T_on,k = first x at which velocity_k(x) reaches 0.5 × max(velocity_k)
```

The reference calculation uses 1,001 equally spaced points in [0, 1], a Savitzky–Golay derivative with a 301-point window (30% of the grid), and polynomial order 4. The first threshold crossing is linearly interpolated between adjacent points. If the velocity already exceeds the threshold at the root, T_on is 0. Flat or already-pure profiles have undefined (NaN) T_on.

T_on replaces the previous peak-location summary; it is computed from the leading half-maximum crossing, not by relabeling peak values. It is dimensionless and is not converted to embryonic days or calibrated against developmental time. All 19 reference onsets agree with the updated E8.5 analysis within its saved rounding precision (maximum absolute difference < 5 × 10⁻⁷).

## Run FateVec on the included example

Install the Python dependencies:

```bash
pip install numpy pandas scipy matplotlib
python run_fatevec.py
```

The included tree CSVs are sufficient; running the R reconstruction first is optional. Configuration in `run_fatevec.py` selects `subcluster` as the cell-type column and defines the grid and smoothing settings.

| Output | Description |
|---|---|
| `outputs/fatevec_curves.pdf` | Mean FateVec components over normalized tree depth |
| `outputs/fatevec_derivatives.pdf` | dF/dx profiles; dots mark T_on |
| `outputs/fatevec_onset_summary.csv` | `Cell type`, `T_on`, and `Max dF/dx` |

The onset CSV replaces `fatevec_peak_summary.csv`. The two existing PDF filenames are retained.

Run the numerical checks with:

```bash
python -m unittest discover -s tests -v
```

## R tree workflow

The default R workflow reads the included merged Newick tree, matches its tips to `data/E8_annotations.csv`, assigns 50 clades with `cutree`, and regenerates the tree exports and fan plot:

```r
install.packages(c("ape", "phangorn", "dplyr", "readr", "tibble", "here", "scales"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ggtree", "ggnewscale"))
```

```bash
Rscript run_pipeline.R
```

Run the R command from the repository root. Set `LINEAGE_OUTPUT_DIR` to write a separate copy of the exports. The local R workflow also generates a fan PDF, tip list and node-to-tip mapping; these auxiliary files are not part of the distributed reference dataset.

### Optional reconstruction from TAPE v2 inputs

To rebuild the tree rather than use the included Newick, supply the matching files below and enable `REBUILD_TREE` in `R/00_config.R`, or run:

```bash
LINEAGE_REBUILD_TREE=true Rscript run_pipeline.R
```

| Input | Format |
|---|---|
| `data/Step7_E8_DM_v2.csv` | Square symmetric distance matrix without row/column headers; order must match the pivot table |
| `data/Step6_E8_v2_TapeBCpivot.csv` | Wide TAPE table with a `Cell` column and one row per merged cell |
| `data/Step5_E8_v2.csv` | One-library-per-cell long table with `Cell`, `TapeBC` and `Site1` |
| `data/E8_annotations.csv` | Included final annotation with `Cell` and `subcluster` |

The first three inputs are not included. The filter requires at least eight edited TAPEs and an annotation, then excludes ExVE before UPGMA clustering. Use the matching v2 inputs together; the old three-library tables use different cell identifiers and are not interchangeable with the merged example. The included helper `external/extra_script_choi_brief.R` accompanies the [Choi et al. method](https://doi.org/10.1038/s41586-022-04922-8); the reconstruction entry point reads a precomputed distance matrix.

## Repository structure

```text
R/                   Filtering, tree reconstruction, clades, plotting and export
run_pipeline.R       R entry point
fatevec/              Fate-vector propagation, tracing, T_on and plotting
run_fatevec.py        Python entry point
external/             Existing R helper
data/                 E8.5 merged cell-type annotation
outputs/              Matching tree, annotation, node/edge tables and FateVec results
tests/                Numerical tests of the T_on estimator
```
