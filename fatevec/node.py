"""Node class for lineage tree with bottom-up fate vector propagation."""

import numpy as np


class Node:
    """A node in the lineage tree carrying a fate vector.

    Each node stores a vector `v` of length `n_cell` (number of unique cell types),
    representing the normalized cell-type proportions in its subtree, and a count
    vector `n` tracking how many tips of each type are descendants.

    Args:
        child: List of child Node objects.
        n_cell: Total number of unique cell types in the tree.
        depth: Edge-length depth from root.
        is_tip: Whether this node is a leaf (tip).
        cell_type: Index into the fate vector for tip nodes (None for internal).
    """

    def __init__(self, child, n_cell, depth, is_tip, cell_type=None):
        self.parent = None
        self.child = child
        self.n_cell = n_cell
        self.depth = depth
        self.is_tip = is_tip
        self.cell_type = cell_type
        self.v = None
        self.n = None

    def set_v(self, v):
        self.v = v

    def set_n(self, n):
        self.n = n

    def update(self):
        """Bottom-up aggregation of fate vectors.

        For each child, accumulates v * n (weighted sum) and n (counts),
        then computes v = weighted_sum / counts and normalizes to sum to 1.
        """
        if len(self.child) == 0:
            self.v = np.zeros((self.n_cell,), dtype=float)
            self.n = np.zeros((self.n_cell,), dtype=float)
            return

        for node in self.child:
            if node.v is None or node.n is None:
                node.update()

        v_agg = np.zeros((self.n_cell,), dtype=float)
        n_agg = np.zeros((self.n_cell,), dtype=float)
        for node in self.child:
            v_agg += node.v * node.n
            n_agg += node.n

        with np.errstate(divide="ignore", invalid="ignore"):
            v = np.zeros_like(v_agg)
            mask = n_agg > 0
            v[mask] = v_agg[mask] / n_agg[mask]

        total = v.sum()
        if total > 0:
            v = v / total

        self.n = n_agg
        self.v = v

    def draw_tree(self, prefix="", is_last=True):
        """Print an ASCII representation of the subtree (for debugging)."""
        connector = "└── " if is_last else "├── "
        label = f"[tip: type={self.cell_type}]" if self.is_tip else "[internal]"
        print(f"{prefix}{connector}{label} depth={self.depth:.4f}")
        new_prefix = prefix + ("    " if is_last else "│   ")
        for i, c in enumerate(self.child):
            c.draw_tree(new_prefix, i == len(self.child) - 1)
