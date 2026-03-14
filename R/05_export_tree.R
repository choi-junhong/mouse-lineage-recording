# ==============================================================================
# 05_export_tree.R - Export phylogenetic tree structure to CSV files
# ==============================================================================

#' Export a phylo tree to a set of CSV files
#'
#' Generates four files in output_dir:
#'   - tree_tips.csv:          node_id, tip_label
#'   - tree_edges.csv:         parent, child, branch_length
#'   - tree_node_tip_map.csv:  node_id (internal), tip_label (descendant tips)
#'   - tree_nodes_summary.csv: node_id, is_tip, label, parent, depths, clade_size
#'
#' @param tree        A phylo object
#' @param output_dir  Directory to write CSV files
#' @param prefix      Optional filename prefix (default: "tree")
export_tree_csv <- function(tree, output_dir, prefix = "tree") {

  stopifnot(inherits(tree, "phylo"))
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  Ntip  <- length(tree$tip.label)
  Nnode <- tree$Nnode

  # 1) Tips
  tips_df <- tibble(
    node_id   = seq_len(Ntip),
    tip_label = tree$tip.label
  )
  readr::write_csv(tips_df, file.path(output_dir, paste0(prefix, "_tips.csv")))

  # 2) Edges
  edges_df <- as.data.frame(tree$edge)
  names(edges_df) <- c("parent", "child")
  if (!is.null(tree$edge.length)) {
    edges_df$branch_length <- as.numeric(tree$edge.length)
  }
  readr::write_csv(edges_df, file.path(output_dir, paste0(prefix, "_edges.csv")))

  # 3) Parent vector
  parent_vec <- rep(NA_integer_, Ntip + Nnode)
  parent_vec[edges_df$child] <- edges_df$parent

  # 4) Internal nodes: descendant tips & node-tip mapping
  internal_nodes <- (Ntip + 1):(Ntip + Nnode)

  desc_list <- lapply(internal_nodes, function(n) {
    phangorn::Descendants(tree, n, type = "tips")[[1]]
  })

  # Node-to-tip long table
  node_tip_map <- tibble(
    node_id   = rep(internal_nodes, times = vapply(desc_list, length, 0L)),
    tip_index = unlist(desc_list, use.names = FALSE)
  ) %>%
    mutate(tip_label = tree$tip.label[tip_index]) %>%
    select(node_id, tip_label)
  readr::write_csv(node_tip_map, file.path(output_dir, paste0(prefix, "_node_tip_map.csv")))

  # 5) Node summary
  node_summary <- tibble(
    node_id = seq_len(Ntip + Nnode),
    is_tip  = node_id <= Ntip,
    label   = ifelse(node_id <= Ntip, tree$tip.label[node_id], NA_character_),
    parent  = parent_vec,
    depths  = ape::node.depth.edgelength(tree)
  )

  clade_size_tbl <- tibble(
    node_id    = internal_nodes,
    clade_size = vapply(desc_list, length, 0L)
  )

  node_summary <- node_summary %>%
    left_join(clade_size_tbl, by = "node_id") %>%
    mutate(clade_size = ifelse(is_tip, 1L, clade_size))
  readr::write_csv(node_summary, file.path(output_dir, paste0(prefix, "_nodes_summary.csv")))

  message(sprintf(
    "Exported tree (%d tips, %d internal nodes) to %s/",
    Ntip, Nnode, output_dir
  ))
}
