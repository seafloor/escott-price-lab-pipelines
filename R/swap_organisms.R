mouse_to_human <- function(genes) {
  genomes <- set_genomes()
  human <- genomes[[1]]
  mouse <- genomes[[2]]

  cols_to_get <- ensemble_columns_to_extract

  # getting full list of human with chromosome locations too
  human_from_mouse <- getLDS(
    attributes = c("mgi_symbol"),
    filters = "mgi_symbol",
    values = genes, mart = mouse,
    attributesL = cols_to_get,
    martL = human, uniqueRows = T, verbose = FALSE
  )

  human_from_mouse <- tibble::as_tibble(human_from_mouse)
  cols_clean <- c("mgi_symbol", clean_output_name)
  colnames(human_from_mouse) <- cols_clean

  human_from_mouse <- force_canonical_autosomes(human_from_mouse)

  return(human_from_mouse)
}
