# ==============================================================================
# 04_plot_tree.R - Circular (fan) tree visualization with clade coloring
# ==============================================================================

#' Plot a circular phylogenetic tree colored by clade
#'
#' @param tr_grp         A phylo object annotated by groupOTU (from assign_clades)
#' @param clade_palette  Named character vector mapping clade IDs to hex colors.
#'                       If NULL, uses a default qualitative palette.
#' @return A ggplot (ggtree) object
plot_fan_tree <- function(tr_grp, clade_palette = NULL) {

  # Determine which clades are actually present in the tree
  pdat    <- ggtree::ggtree(tr_grp)$data
  present <- sort(unique(na.omit(as.character(pdat$clade))))

  # Build palette
  if (is.null(clade_palette)) {
    n_clades <- length(present)
    clade_palette <- setNames(
      scales::hue_pal()(n_clades),
      present
    )
  } else {
    clade_palette <- clade_palette[present]
  }

  p <- ggtree::ggtree(
    tr_grp,
    layout = "fan",
    aes(color = factor(clade)),
    linewidth = 0.55
  ) +
    ggplot2::scale_color_manual(
      values = clade_palette,
      breaks = present,
      name   = "Clades"
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 5, 5))

  p
}
