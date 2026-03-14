#!/usr/bin/env python3
"""
run_fatevec.py - FateVec: Cell-Fate Commitment Analysis on Lineage Trees

Usage:
    1. Run Step 1 (run_pipeline.R) first to generate tree CSVs.
    2. Edit the configuration section below to point to your outputs.
    3. Run: python run_fatevec.py
"""

import os

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from fatevec import (
    construct_tree,
    convert_csv_to_nodes,
    avg_vector_component_vs_time,
    find_peak_derivatives,
)
from fatevec.plotting import plot_fate_curves, plot_derivative_curves

# ==============================================================================
# Configuration
# ==============================================================================

# Input files (outputs from Step 1)
TREE_NODES_SUMMARY = os.path.join("outputs", "tree_nodes_summary.csv")
CELL_ANNOT_CSV     = os.path.join("outputs", "cell_annot_with_clades.csv")

# Column name for cell type in the annotation CSV (groups fate vectors by this)
CELL_TYPE_COL = "subcluster"

# Output directory
OUTPUT_DIR = "outputs"

# Pseudo-time range for fate vector sampling
T_MIN = -0.05
T_MAX = 0.50
T_POINTS = 10001

# Savitzky-Golay smoothing parameters
WINDOW_LENGTH = 3501   # Must be odd
POLYORDER = 4

# ==============================================================================
# Step 1: Parse tree and build node structure
# ==============================================================================
print("=== Step 1: Parsing tree and cell annotations ===")
data, cells, n_cells = convert_csv_to_nodes(
    TREE_NODES_SUMMARY, CELL_ANNOT_CSV, cell_type_col=CELL_TYPE_COL
)
print(f"  Cell types found: {list(cells.keys())}")

print("\n=== Step 2: Constructing tree ===")
nodes, head = construct_tree(data, cells, n_cells)

print("\n=== Step 3: Bottom-up fate vector propagation ===")
head.update()
print(f"  Root fate vector: {head.v}")
print(f"  Root cell counts: {head.n}")

# ==============================================================================
# Step 4: Compute average fate curves
# ==============================================================================
print("\n=== Step 4: Computing average fate curves ===")
t_array = np.linspace(T_MIN, T_MAX, T_POINTS)
avg_v = avg_vector_component_vs_time(nodes, t_array)

# Plot fate curves
ax = plot_fate_curves(avg_v, t_array, cells)
fate_plot_path = os.path.join(OUTPUT_DIR, "fatevec_curves.pdf")
ax.figure.savefig(fate_plot_path, dpi=150, bbox_inches="tight")
print(f"  Saved fate curves to: {fate_plot_path}")
plt.close()

# ==============================================================================
# Step 5: Find peak derivatives
# ==============================================================================
print("\n=== Step 5: Finding peak commitment rates ===")
peak_times, peak_values, derivatives = find_peak_derivatives(
    avg_v, t_array, window_length=WINDOW_LENGTH, polyorder=POLYORDER
)

# Plot derivative curves
ax = plot_derivative_curves(
    avg_v, t_array, cells,
    window_length=WINDOW_LENGTH, polyorder=POLYORDER
)
deriv_plot_path = os.path.join(OUTPUT_DIR, "fatevec_derivatives.pdf")
ax.figure.savefig(deriv_plot_path, dpi=150, bbox_inches="tight")
print(f"  Saved derivative curves to: {deriv_plot_path}")
plt.close()

# ==============================================================================
# Step 6: Export results
# ==============================================================================
print("\n=== Step 6: Exporting results ===")
idx_to_label = {v: k for k, v in cells.items()}
records = []
for key in sorted(peak_times.keys()):
    cell_type = idx_to_label.get(key, f"Type_{key}")
    records.append({
        "Cell type": cell_type,
        "Time of max dv/dt": f"{peak_times[key]:.5f}",
        "Max d(vk)/dt": f"{peak_values[key]:.5f}",
    })

df = pd.DataFrame(records)
output_csv = os.path.join(OUTPUT_DIR, "fatevec_peak_summary.csv")
df.to_csv(output_csv, index=False)
print(f"  Saved peak summary to: {output_csv}")

print("\n=== FateVec analysis complete! ===")
