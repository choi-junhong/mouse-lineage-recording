"""FateVec: Cell-fate commitment analysis on lineage trees."""

from .node import Node
from .tree import construct_tree, convert_csv_to_nodes
from .analysis import (
    trace_back,
    vector_component_vs_time,
    avg_vector_component_vs_time,
    find_peak_derivatives,
)

__all__ = [
    "Node",
    "construct_tree",
    "convert_csv_to_nodes",
    "trace_back",
    "vector_component_vs_time",
    "avg_vector_component_vs_time",
    "find_peak_derivatives",
]
