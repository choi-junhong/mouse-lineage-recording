"""Visualization functions for FateVec analysis."""

import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import savgol_filter


def plot_fate_curves(avg_v, t_array, cells, ax=None, **kwargs):
    """Plot average fate-proportion curves for each cell type.

    Args:
        avg_v: Dict mapping cell-type index -> 1D array of proportions.
        t_array: 1D array of time values.
        cells: Dict mapping cell-type label -> index.
        ax: Optional matplotlib Axes. Creates a new figure if None.
        **kwargs: Additional keyword arguments passed to ax.plot().

    Returns:
        matplotlib Axes.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))

    idx_to_label = {v: k for k, v in cells.items()}

    for key in sorted(avg_v.keys()):
        label = idx_to_label.get(key, f"Type {key}")
        ax.plot(t_array, avg_v[key], label=label, **kwargs)

    ax.set_xlabel("Pseudo-time (depth)")
    ax.set_ylabel("Fate proportion")
    ax.set_title("Average Fate Vector Components vs. Pseudo-time")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize="small")
    ax.figure.tight_layout()
    return ax


def plot_derivative_curves(avg_v, t_array, cells,
                           window_length=3501, polyorder=4,
                           ax=None, **kwargs):
    """Plot smoothed derivative (dv/dt) curves for each cell type.

    Args:
        avg_v: Dict mapping cell-type index -> 1D array of proportions.
        t_array: 1D array of time values.
        cells: Dict mapping cell-type label -> index.
        window_length: Savitzky-Golay window length.
        polyorder: Polynomial order for smoothing.
        ax: Optional matplotlib Axes.
        **kwargs: Additional keyword arguments passed to ax.plot().

    Returns:
        matplotlib Axes.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))

    dt = np.mean(np.diff(t_array))
    idx_to_label = {v: k for k, v in cells.items()}

    for key in sorted(avg_v.keys()):
        dvdt = savgol_filter(avg_v[key], window_length, polyorder, deriv=1, delta=dt)
        label = idx_to_label.get(key, f"Type {key}")
        ax.plot(t_array, dvdt, label=label, **kwargs)

    ax.set_xlabel("Pseudo-time (depth)")
    ax.set_ylabel("d(fate proportion)/dt")
    ax.set_title("Fate Commitment Rate (dv/dt) vs. Pseudo-time")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize="small")
    ax.figure.tight_layout()
    return ax
