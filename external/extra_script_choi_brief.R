
#print("Loading Monocle3 and Seurat")
#suppressMessages(library(monocle3))
#suppressMessages(library(Seurat))

print("Loading packages for regular data analysis, e.g. dplyr")
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(readr))
suppressMessages(library(tibble))

print("Loading packages for plotting, e.g. ggplot2")
suppressMessages(library(ggplot2))
suppressMessages(library(plotly))
suppressMessages(library(htmlwidgets))
suppressMessages(library(gridExtra))
suppressMessages(library(viridis))
suppressMessages(library(gplots))
suppressMessages(library(ggnewscale))

print("Loading packages for phylo tree analysis")
suppressMessages(library(ggtreeExtra))
suppressMessages(library(ggtree))
#suppressMessages(library(PATH))
suppressMessages(library(ape))
suppressMessages(library(qlcMatrix))
suppressMessages(library(castor))
suppressMessages(library(phytools))
#suppressMessages(library(distances))
#suppressMessages(library(DescTools))


rowNorm <- function(W) {
  s <- Matrix::rowSums(W)
  W/ifelse(s==0, 1, s)
}


#' Phylogenetic weight matrix for a node-depth of one. 
#'
#' This function computes the sum normalized phylogenetic weight matrix
#' using the node-depth of one only weighting function. This function
#' returns an N by N dimensional matrix with a value in the ijth position
#' if cell i and j are separated by one node on the phylogeny and a 0 otherwise.

one_node_tree_dist <- function(tree, norm = TRUE) {
  w <- one_node_depth(tree)$W
  if(norm == TRUE) {
    w <- rowNorm(w)
  }
  w <- w/sum(w)
  return(w)
}

#' Get the inverse tree distances between all cells
#'
#' This function gets the inverse tree or phylogenetic distances, in terms of node distances or branch lengths, between all cell pairs. Inverse tree distances are normalized to sum to one.  

inv_tree_dist <- function(tree, node = TRUE, norm = TRUE) {
  w <- castor::get_all_pairwise_distances(tree, only_clades = 1:length(tree$tip.label), as_edge_counts = node)
  if(node == TRUE) {
    w <- w - 1
  }
  w <- 1/w
  diag(w) <- 0
  if(norm == TRUE) {
    w <- rowNorm(w)
  }
  return(w/sum(w))
}

exp_tree_dist <- function(tree, node = TRUE, norm = TRUE) {
  w <- castor::get_all_pairwise_distances(tree, only_clades = 1:length(tree$tip.label), as_edge_counts = node)
  w <- exp(-w)
  diag(w) <- 0
  if(norm == TRUE) {
    w <- rowNorm(w)
  }
  return(w/sum(w))
}





#' Convert categorical cell state vector into categorical matrix.
catMat0 <- function(state_vector, label_order = NULL, n = NULL,
                    return_sparse = TRUE) {
  
  if(is.null(n) & is.numeric(state_vector)==FALSE) {
    n <- max(length(unique(state_vector)), length(unique(label_order)))
  } else if(is.numeric(state_vector)==TRUE){
    n <- max(state_vector, n)
  }
  nI <- n
  if(is.null(label_order)) {
    if(is.numeric(state_vector)==FALSE) {
      x <- factor(state_vector, ordered = TRUE)
      nI <- length(levels(x))
    } else {
      if(min(state_vector)==0) {
        min <- 0
        nI <- nI + 1
        
      } else {
        min <- 1
      }
      x <- factor(state_vector, ordered = T, levels = seq(min,n, 1))
    }
  } else {
    x <- factor(state_vector, ordered = TRUE, levels = label_order)
    nI <- length(levels(x))
  }
  z <- Matrix::.sparseDiagonal(nI)[as.numeric(x),]
  tryCatch({colnames(z)[1:length(levels(x))] <- levels(x)},
           error=function(e) {warning("Column names not assigned", call. = F); return(z)})
  if(is.null(names(state_vector))==FALSE) {
    rownames(z) <- make.names(names(x), unique = TRUE)
  }
  if(return_sparse==TRUE) {
    return(z)
  } else {
    return(as.matrix(z))
  }
}


#' Convert categorical cell state vector into categorical matrix.
catMat <- function(cell_states, 
                   num_states = NULL, state_order = NULL,
                   sparse = TRUE, unformatted = FALSE) {
  # Converts a vector of categorical states into a categorical matrix.
  if(unformatted==FALSE) {
    catMat0(state_vector = cell_states, label_order = state_order, n=num_states, return_sparse = sparse)
  } else {
    if(is.null(num_states)) {
      stop("if unformatted=TRUE, num_states cannot be NULL", call. = F)
    }
    diag(num_states)[factor(cell_states),]
  }
}

