# Mouse lineage recording

Code and example data to reproduce the E8.5 FateVec results.

## Requirements

Python 3 with `numpy`, `pandas`, `scipy` and `matplotlib`.

## Install

```bash
git clone https://github.com/choi-junhong/mouse-lineage-recording.git
cd mouse-lineage-recording
python -m pip install numpy pandas scipy matplotlib
```

## Quickstart

Run the included merged E8.5 example:

```bash
python run_fatevec.py
```

The example contains **10,579 tree tips and 19 cell types**. Its input files are already included:

- [tree_nodes_summary.csv](outputs/tree_nodes_summary.csv): tree structure and node depths.
- [cell_annot_with_clades.csv](outputs/cell_annot_with_clades.csv): cell-type annotations.

## Outputs

Results are saved to `outputs/`:

| File | Content |
|---|---|
| [fatevec_curves.pdf](outputs/fatevec_curves.pdf) | FateVec curves |
| [fatevec_derivatives.pdf](outputs/fatevec_derivatives.pdf) | Velocity curves with T_on markers |
| [fatevec_onset_summary.csv](outputs/fatevec_onset_summary.csv) | T_on for each cell type |

T_on is reported on normalized tree depth (root = 0, tips = 1).

<details>
<summary>Optional: regenerate tree exports with R</summary>

Install the R packages:

```r
install.packages(c("ape", "phangorn", "dplyr", "readr", "tibble", "here", "scales"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ggtree", "ggnewscale"))
```

Run from the repository directory:

```bash
Rscript run_pipeline.R
python run_fatevec.py
```

The R script uses the included Newick tree and merged annotations. Settings are in [R/00_config.R](R/00_config.R).

</details>

Tree workflow reference: [Choi et al. (2022)](https://doi.org/10.1038/s41586-022-04922-8).
