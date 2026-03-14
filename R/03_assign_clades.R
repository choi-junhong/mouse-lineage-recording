# ==============================================================================
# 03_assign_clades.R - Assign clade labels to tree tips via cutree
# ==============================================================================

#' Generate a clade map from cutree grouping
#'
#' @param grp_tip  Named integer vector from cutree (names = tip labels)
#' @return A tibble with columns: Cell, clone_id
clade_map_from_cutree <- function(grp_tip) {
  stopifnot(!is.null(names(grp_tip)))
  tibble(
    Cell     = names(grp_tip),
    clone_id = paste0("CLADE_", as.integer(unname(grp_tip)))
  )
}

#' Generate a clade map from a ggtree groupOTU-annotated tree
#'
#' @param tr_grp  A phylo object annotated by ggtree::groupOTU
#' @return A tibble with columns: Cell, clone_id
clade_map_from_trgrp <- function(tr_grp) {
  pd <- as.data.frame(ggtree::ggtree(tr_grp)$data)
  tibble(
    Cell  = as.character(pd$label),
    clade = pd$clade
  ) %>%
    filter(!is.na(Cell), !is.na(clade)) %>%
    transmute(Cell, clone_id = paste0("CLADE_", as.character(clade)))
}

#' Assign clade IDs to tree tips and annotate cell metadata
#'
#' @param tree    A phylo object
#' @param k       Number of clades for cutree
#' @return A list with:
#'   - grp_tip:   named integer vector from cutree
#'   - groups:    list of tip-label vectors per clade
#'   - tr_grp:    phylo object annotated with groupOTU
#'   - clade_map: tibble mapping Cell -> clone_id
assign_clades <- function(tree, k = 50) {
  hc      <- as.hclust(tree)
  grp_tip <- cutree(hc, k = k)
  groups  <- split(names(grp_tip), grp_tip)
  tr_grp  <- ggtree::groupOTU(tree, groups, group_name = "clade")

  clade_map <- clade_map_from_trgrp(tr_grp)
  message(sprintf("Assigned %d clades to %d tips.", k, length(grp_tip)))

  list(
    grp_tip   = grp_tip,
    groups    = groups,
    tr_grp    = tr_grp,
    clade_map = clade_map
  )
}

#' Merge clade IDs into cell annotation table
#'
#' @param cell_annot  Data frame of cell annotations
#' @param clade_map   Tibble with Cell and clone_id columns
#' @return Updated data frame with clone_id column
write_clade_into_meta <- function(cell_annot, clade_map) {
  clade_map <- clade_map %>%
    filter(!is.na(Cell), !is.na(clone_id)) %>%
    distinct(Cell, .keep_all = TRUE)

  out <- cell_annot %>%
    select(-any_of("clone_id")) %>%
    left_join(clade_map, by = "Cell")

  n_matched <- sum(!is.na(out$clone_id))
  message(sprintf(
    "Wrote clone_id for %d / %d cells (%.1f%%).",
    n_matched, nrow(out), 100 * n_matched / nrow(out)
  ))
  out
}
