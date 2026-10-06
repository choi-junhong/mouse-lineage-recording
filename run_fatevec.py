#!/usr/bin/env python3
"""Compute FateVec curves and T_on from the included merged E8.5 tree.

Run `python run_fatevec.py` directly; the required tree CSVs are included.
T_on is the leading half-maximum crossing on normalized tree depth.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from fatevec import (
    construct_tree, convert_csv_to_nodes,
    avg_vector_component_vs_time, find_onset_times,
)
from fatevec.plotting import plot_fate_curves, plot_derivative_curves

BASE_DIR = Path(__file__).resolve().parent
TREE_NODES_SUMMARY = BASE_DIR / 'outputs/tree_nodes_summary.csv'
CELL_ANNOT_CSV = BASE_DIR / 'outputs/cell_annot_with_clades.csv'
CELL_TYPE_COL = 'subcluster'
OUTPUT_DIR = BASE_DIR / 'outputs'

# Dimensionless tree depth: root = 0, tips = 1 (not developmental time).
T_POINTS = 1001
WINDOW_LENGTH = 301  # 30% of the grid, rounded to an odd length
POLYORDER = 4


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    data, cells, n_types = convert_csv_to_nodes(
        TREE_NODES_SUMMARY, CELL_ANNOT_CSV, cell_type_col=CELL_TYPE_COL)
    if any(is_tip and label is None for _, is_tip, label, *_ in data):
        raise ValueError('Every tree tip must have a cell-type annotation.')
    nodes, head = construct_tree(data, cells, n_types)
    tip_depths = np.array([n.depth for n in nodes.values() if n.is_tip])
    tip_depth = float(tip_depths.max())
    if (not np.isfinite(tip_depths).all() or tip_depth <= 0 or
            not np.allclose(tip_depths, tip_depth, atol=1e-8, rtol=0)):
        raise ValueError('The reference analysis requires an ultrametric tree.')
    if not np.isfinite(head.depth):
        head.depth = 0.0
    if not np.isclose(head.depth, 0):
        raise ValueError('Root depth must be zero.')
    for node in nodes.values():
        node.depth /= tip_depth
    head.update()
    print(f'Analyzing {len(tip_depths):,} tips and {n_types} cell types.')
    t_array = np.linspace(0.0, 1.0, T_POINTS)
    avg_v = avg_vector_component_vs_time(nodes, t_array)
    onset, maxima, derivatives = find_onset_times(
        avg_v, t_array, window_length=WINDOW_LENGTH, polyorder=POLYORDER)

    ax = plot_fate_curves(avg_v, t_array, cells)
    ax.figure.savefig(OUTPUT_DIR / 'fatevec_curves.pdf', bbox_inches='tight',
                      metadata={'CreationDate': None, 'ModDate': None})
    plt.close(ax.figure)
    ax = plot_derivative_curves(
        avg_v, t_array, cells, window_length=WINDOW_LENGTH,
        polyorder=POLYORDER, onset_times=onset, derivatives=derivatives)
    ax.figure.savefig(OUTPUT_DIR / 'fatevec_derivatives.pdf', bbox_inches='tight',
                      metadata={'CreationDate': None, 'ModDate': None})
    plt.close(ax.figure)

    labels = {index: label for label, index in cells.items()}
    summary = pd.DataFrame([
        {'Cell type': labels[key], 'T_on': onset[key], 'Max dF/dx': maxima[key]}
        for key in sorted(onset)
    ])
    path = OUTPUT_DIR / 'fatevec_onset_summary.csv'
    summary.to_csv(path, index=False, float_format='%.10f')
    print(summary.to_string(index=False))
    print(f'Saved T_on summary to {path}')


if __name__ == '__main__':
    main()
