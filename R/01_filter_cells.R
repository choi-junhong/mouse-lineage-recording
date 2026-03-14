# ==============================================================================
# 01_filter_cells.R - Filter cells by TAPE recovery threshold
# ==============================================================================

#' Filter cells based on minimum number of edited TAPEs recovered
#'
#' @param tape_file     Path to the TAPE barcode pivot CSV (wide format)
#' @param long_paths    Character vector of paths to long-format tape selection CSVs
#' @param annot_file    Path to cell-type annotation CSV
#' @param min_tapes     Minimum number of edited TAPEs required per cell (default: 8)
#' @return A list with:
#'   - cell_list_filtered: character vector of filtered cell IDs
#'   - cell_annot_filtered: data.frame of annotations for filtered cells
filter_cells <- function(tape_file, long_paths, annot_file, min_tapes = 8) {

  # Load cell annotations
  cell_annot_all <- read.csv(
    annot_file,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )

  # Load TAPE barcode pivot table to get the master cell list
  tape_table <- read.csv(
    tape_file,
    stringsAsFactors = FALSE,
    header = TRUE,
    na.strings = c("", "NA")
  )
  cell_list <- tape_table$Cell

  # Load and merge long-format tape selection files
  tape_table_long <- lapply(long_paths, function(p) {
    read.csv(p, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  }) %>% bind_rows()

  # Filter cells: require min_tapes edited TAPEs and presence in annotation
  cell_list_filtered <- tape_table_long %>%
    filter(Site1 != "None") %>%
    count(Cell, name = "Tapes") %>%
    filter(Tapes >= min_tapes) %>%
    semi_join(cell_annot_all, by = "Cell") %>%
    pull(Cell)

  # Intersect with the master cell list from the pivot table
  cell_list_filtered <- intersect(cell_list, cell_list_filtered)
  cell_annot_filtered <- filter(cell_annot_all, Cell %in% cell_list_filtered)

  message(sprintf(
    "Filtered to %d cells (from %d) with >= %d edited TAPEs.",
    length(cell_list_filtered), length(cell_list), min_tapes
  ))

  list(
    cell_list         = cell_list,
    cell_list_filtered = cell_list_filtered,
    cell_annot_filtered = cell_annot_filtered
  )
}
