"""Tree construction from CSV exports of the R pipeline."""

import csv

import numpy as np

from .node import Node


def convert_csv_to_nodes(tree_nodes_summary, cell_annot_csv, cell_type_col="CellType"):
    """Parse tree node summary and cell annotation CSVs.

    Args:
        tree_nodes_summary: Path to tree_nodes_summary.csv (from Step 1).
        cell_annot_csv: Path to cell annotation CSV with a CellType column.
        cell_type_col: Column name (str) or 0-based index (int) for cell type
                       in the annotation CSV. Default: "CellType".

    Returns:
        data: List of [node_id, is_tip, label, parent, depth, clade_size].
        cells: Dict mapping cell-type label -> integer index.
        n_cells: Number of unique cell types.
    """
    # Read cell annotation: map cell barcode -> cell type label
    tips_dict = {}
    with open(cell_annot_csv) as f:
        reader = csv.reader(f)
        header = next(reader)

        # Resolve column name to index
        if isinstance(cell_type_col, str):
            if cell_type_col not in header:
                raise ValueError(
                    f"Column '{cell_type_col}' not found in {cell_annot_csv}. "
                    f"Available columns: {header}"
                )
            col_idx = header.index(cell_type_col)
        else:
            col_idx = cell_type_col

        for row in reader:
            if len(row) > col_idx:
                tips_dict[row[0]] = row[col_idx]

    # Read tree nodes summary
    data = []
    with open(tree_nodes_summary) as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            if not row:
                continue

            node_id = int(row[0])
            is_tip = row[1] == "TRUE"

            if is_tip:
                label = tips_dict.get(row[2])
            else:
                label = None

            try:
                parent = int(row[3])
            except (ValueError, IndexError):
                parent = -1

            try:
                depth = float(row[4])
            except (ValueError, IndexError):
                depth = np.nan

            try:
                clade_size = int(row[5])
            except (ValueError, IndexError):
                clade_size = 1

            data.append([node_id, is_tip, label, parent, depth, clade_size])

    # Build cell-type index dictionary
    cells = {}
    n_cells = 0
    for node_id, is_tip, label, parent, depth, clade_size in data:
        if is_tip and label is not None and label not in cells:
            cells[label] = n_cells
            n_cells += 1

    return data, cells, n_cells


def construct_tree(data, cells, n_cells):
    """Build the Node tree from parsed data.

    Args:
        data: List of [node_id, is_tip, label, parent, depth, clade_size].
        cells: Dict mapping cell-type label -> integer index.
        n_cells: Number of unique cell types.

    Returns:
        nodes: Dict mapping node_id -> Node.
        head: The root Node.
    """
    nodes = {}
    head = None

    for line_data in data:
        node_id, is_tip, label, parent, depth, clade_size = line_data[:6]

        if is_tip:
            cell_idx = cells.get(label)
            node = Node([], n_cells, depth, True, cell_idx)
            node.set_v(np.zeros((n_cells,), dtype=float))
            node.set_n(np.zeros((n_cells,), dtype=float))
            if cell_idx is not None:
                node.v[cell_idx] = 1.0
                node.n[cell_idx] = 1.0
            nodes[node_id] = node

            if parent in nodes:
                nodes[parent].child.append(node)
            else:
                nodes[parent] = Node([node], n_cells, 0, False)
            node.parent = nodes[parent]
        else:
            if node_id not in nodes:
                nodes[node_id] = Node([], n_cells, depth, False)
            else:
                nodes[node_id].depth = depth
            node = nodes[node_id]

            if parent == -1:
                head = node
            else:
                if parent in nodes:
                    nodes[parent].child.append(node)
                else:
                    nodes[parent] = Node([node], n_cells, 0, False)
                node.parent = nodes[parent]

    # Fallback: find root if parent=-1 was not encountered
    if head is None:
        child_ids = set()
        for nd in nodes.values():
            for ch in nd.child:
                for k, v in nodes.items():
                    if v is ch:
                        child_ids.add(k)
                        break
        root_candidates = set(nodes.keys()) - child_ids
        head = nodes[next(iter(root_candidates))] if root_candidates else next(iter(nodes.values()))

    return nodes, head
