"""Fate vector analysis: tracing, averaging, and peak derivative finding."""

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


def find_peak_derivatives(avg_v, t_array, window_length=3501, polyorder=4):
    """Find the time of maximum fate commitment rate for each cell type.

    Applies Savitzky-Golay smoothing to compute the derivative dv/dt,
    then locates the global maximum.

    Args:
        avg_v: Dict from avg_vector_component_vs_time.
        t_array: 1D array of time values.
        window_length: Savitzky-Golay window length (must be odd).
        polyorder: Polynomial order for smoothing.

    Returns:
        peak_times: Dict mapping cell-type index -> time of max dv/dt.
        peak_values: Dict mapping cell-type index -> max dv/dt value.
        derivatives: Dict mapping cell-type index -> full dv/dt array.
    """
    dt = np.mean(np.diff(t_array))
    peak_times = {}
    peak_values = {}
    derivatives = {}

    for key, v_array in avg_v.items():
        dvdt = savgol_filter(v_array, window_length, polyorder, deriv=1, delta=dt)
        max_index = np.argmax(dvdt)
        peak_times[key] = t_array[max_index]
        peak_values[key] = dvdt[max_index]
        derivatives[key] = dvdt

    return peak_times, peak_values, derivatives
