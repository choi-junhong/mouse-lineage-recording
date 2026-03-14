# ==============================================================================
# 02_build_tree.R - Build lineage tree from distance matrix
# ==============================================================================

#' Load distance matrix and build a phylogenetic tree via hierarchical clustering
#'
#' @param dm_file       Path to the distance matrix CSV
#' @param cell_list     Full cell list matching DM row/column order
#' @param cell_list_filtered  Subset of cells to include in the tree
#' @param method        Clustering method for hclust (default: "average")
#' @return A phylo object
build_tree <- function(dm_file, cell_list, cell_list_filtered, method = "average") {

  n_cells <- length(cell_list)
  message("Loading distance matrix (this may take 1-2 minutes)...")

  # First pass: detect format by reading a few rows
  probe <- read.csv(dm_file, stringsAsFactors = FALSE, header = FALSE,
                    na.strings = c("", "NA"), nrows = 3)
  n_col_raw <- ncol(probe)
  message(sprintf("  DM file: detected %d columns, expecting %d cells.", n_col_raw, n_cells))

  # Heuristic: if ncol == n_cells + 1, the first column is likely row names/index
  has_row_index <- (n_col_raw == n_cells + 1)
  if (has_row_index) {
    message("  Detected extra first column (row names/index) — skipping it.")
  }

  # Full read (always as plain data — drop index column manually to avoid
  # "missing values in row.names" when the index column contains NAs)
  distance_matrix <- read.csv(
    dm_file,
    stringsAsFactors = FALSE,
    header = FALSE,
    na.strings = c("", "NA")
  )
  if (has_row_index) {
    distance_matrix <- distance_matrix[, -1, drop = FALSE]
  }

  # If nrow == n_cells + 1, the first row might be a header
  if (nrow(distance_matrix) == n_cells + 1) {
    message("  Detected extra first row (header) — removing it.")
    distance_matrix <- distance_matrix[-1, , drop = FALSE]
    # Re-coerce columns to numeric after dropping the header row
    distance_matrix <- as.data.frame(lapply(distance_matrix, as.numeric))
  }

  # Final dimension check
  if (nrow(distance_matrix) != n_cells || ncol(distance_matrix) != n_cells) {
    stop(sprintf(
      paste0("Dimension mismatch: distance matrix is %d x %d but cell_list has %d entries.\n",
             "  DM file: %s\n",
             "  Please check that the DM and TAPE pivot table refer to the same cell set."),
      nrow(distance_matrix), ncol(distance_matrix), n_cells, dm_file
    ))
  }

  rownames(distance_matrix) <- cell_list
  colnames(distance_matrix) <- cell_list

  # Subset to filtered cells
  dm_filtered <- distance_matrix[cell_list_filtered, cell_list_filtered, drop = FALSE]
  message(sprintf("  Distance matrix filtered: %d x %d cells.", nrow(dm_filtered), ncol(dm_filtered)))

  # Hierarchical clustering -> phylo
  tree <- ape::as.phylo(hclust(as.dist(dm_filtered), method = method))

  message("Tree built successfully.")
  tree
}
