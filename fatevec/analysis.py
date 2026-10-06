"""Fate vector analysis: tracing, averaging, and half-maximum onset estimation."""

import numpy as np
from scipy.signal import savgol_filter


def trace_back(node):
    """Trace the path from a tip node back to the root.

    Args:
        node: A tip Node.

    Returns:
        List of Nodes from tip to root.
    """
    path = [node]
    while node.parent is not None:
        node = node.parent
        path.append(node)
    return path


def vector_component_vs_time(node, t_array):
    """Get a tip's own cell-type proportion along its lineage over pseudo-time.

    For each time point in t_array, returns the fate vector component
    corresponding to the tip's cell type, sampled from the ancestor at
    that depth.

    Args:
        node: A tip Node (must have is_tip=True).
        t_array: 1D array of time (depth) values to sample.

    Returns:
        1D array of fate proportions at each time point.
    """
    if not node.is_tip:
        raise ValueError("Node is not a tip")

    node_history = trace_back(node)
    t_history = np.asarray([n.depth for n in node_history])
    cell_idx = node.cell_type
    v_history = np.asarray([n.v[cell_idx] for n in node_history])

    t_array = np.asarray(t_array)

    # t_history is descending (tip -> root); reverse for searchsorted
    t_asc = t_history[::-1]
    v_asc = v_history[::-1]

    idx = np.searchsorted(t_asc, t_array, side="right")
    idx_desc = len(t_history) - idx - 1
    idx_desc = np.clip(idx_desc, 0, len(v_history) - 1)

    return v_history[idx_desc]


def avg_vector_component_vs_time(nodes, t_array):
    """Average the fate-proportion trace across all tips of each cell type.

    Args:
        nodes: Dict of node_id -> Node.
        t_array: 1D array of time (depth) values.

    Returns:
        Dict mapping cell-type index -> 1D array of averaged proportions.
    """
    v_by_type = {}
    for node in nodes.values():
        if not node.is_tip:
            continue
        v = vector_component_vs_time(node, t_array)
        ct = node.cell_type
        if ct in v_by_type:
            v_by_type[ct].append(v)
        else:
            v_by_type[ct] = [v]

    avg_v = {}
    for key, traces in v_by_type.items():
        avg_v[key] = np.mean(np.array(traces), axis=0)
    return avg_v


def find_onset_times(avg_v, t_array, window_length=301, polyorder=4):
    """Find the leading half-maximum crossing of each fate-bias velocity.

    For each mean FateVec component V, normalize its dynamic range as
    F(x) = (V(x) - V(0)) / (1 - V(0)). A Savitzky-Golay derivative gives
    dF/dx. T_on is the FIRST crossing of half its global maximum, with
    linear interpolation between grid points. If the curve starts above
    this threshold, T_on is the first grid point. Flat/uninformative
    curves have NaN onset and peak velocity.

    The reference analysis uses normalized tree depth x in [0, 1], 1,001
    points, a 301-point window (30% of the grid), and polynomial order 4.
    No conversion to developmental time is performed.

    Returns:
        onset_times: Dict mapping cell-type index to T_on.
        peak_values: Dict mapping cell-type index to maximum dF/dx.
        derivatives: Dict mapping cell-type index to its full dF/dx array.
    """
    t_array = np.asarray(t_array, dtype=float)
    if (t_array.ndim != 1 or len(t_array) < 3 or
            not np.isfinite(t_array).all()):
        raise ValueError("Depth grid must be a finite one-dimensional array.")
    steps = np.diff(t_array)
    if (steps <= 0).any() or not np.allclose(steps, steps[0]):
        raise ValueError("Depth grid must be increasing and uniformly spaced.")
    if (not isinstance(window_length, (int, np.integer)) or
            window_length % 2 != 1 or not 0 <= polyorder < window_length or
            window_length > len(t_array)):
        raise ValueError("Use an odd window <= grid length and > polyorder.")

    onset_times, peak_values, derivatives = {}, {}, {}
    for key, values in avg_v.items():
        v = np.asarray(values, dtype=float)
        if v.shape != t_array.shape or not np.isfinite(v).all():
            raise ValueError(f"Invalid fate curve for cell type {key}.")
        onset_times[key] = np.nan
        peak_values[key] = np.nan
        if v[0] >= 1 - 1e-12:
            derivatives[key] = np.full_like(v, np.nan)
            continue
        fraction = (v - v[0]) / (1 - v[0])
        velocity = savgol_filter(fraction, window_length, polyorder,
                                deriv=1, delta=float(steps[0]))
        derivatives[key] = velocity
        maximum = float(np.max(velocity))
        if maximum <= 1e-10:
            continue
        peak_values[key] = maximum
        threshold = 0.5 * maximum
        index = int(np.flatnonzero(velocity >= threshold)[0])
        if index == 0:
            onset_times[key] = float(t_array[0])
        else:
            onset_times[key] = float(np.interp(
                threshold, velocity[index - 1:index + 1],
                t_array[index - 1:index + 1]))
    return onset_times, peak_values, derivatives
