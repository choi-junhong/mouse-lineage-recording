"""FateVec curves and half-maximum onset on normalized tree depth."""
import numpy as np
from matplotlib import pyplot as plt
from .analysis import find_onset_times


def plot_fate_curves(avg_v, t_array, cells, ax=None, **kwargs):
    """Plot each cell type's mean FateVec component, V(x)."""
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))
    labels = {index: label for label, index in cells.items()}
    for i, key in enumerate(sorted(avg_v)):
        style = {'color': plt.get_cmap('tab20')(i % 20), **kwargs}
        ax.plot(t_array, avg_v[key], label=labels.get(key, f'Type {key}'), **style)
    ax.set_xlabel('Normalized tree depth (root = 0, tips = 1)')
    ax.set_ylabel('Mean FateVec component')
    ax.set_title('FateVec profiles — merged E8.5')
    ax.set_xlim(0, 1)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize='small')
    ax.figure.tight_layout()
    return ax


def plot_derivative_curves(avg_v, t_array, cells, window_length=301, polyorder=4,
                           ax=None, onset_times=None, derivatives=None, **kwargs):
    """Plot dF/dx and mark its leading half-maximum crossing, T_on.

    F(x) = (V(x) - V(0)) / (1 - V(0)). Optional precomputed arrays ensure
    the plotted velocities and exported onset values are identical.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))
    if onset_times is None or derivatives is None:
        onset_times, _, derivatives = find_onset_times(
            avg_v, t_array, window_length=window_length, polyorder=polyorder)
    labels = {index: label for label, index in cells.items()}
    for i, key in enumerate(sorted(avg_v)):
        style = {'color': plt.get_cmap('tab20')(i % 20), **kwargs}
        line, = ax.plot(t_array, derivatives[key], label=labels.get(key, f'Type {key}'), **style)
        onset = onset_times[key]
        if np.isfinite(onset):
            ax.plot(onset, np.interp(onset, t_array, derivatives[key]), 'o',
                    color=line.get_color(), markersize=4)
    ax.set_xlabel('Normalized tree depth (root = 0, tips = 1)')
    ax.set_ylabel('Fate-bias acquisition velocity, dF/dx')
    ax.set_title(r'Fate-bias acquisition — dots mark $T_{on}$')
    ax.set_xlim(0, 1)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize='small')
    ax.figure.tight_layout()
    return ax
