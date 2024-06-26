#' Convert Mouse Genes to Human Genes
#'
#' This function maps mouse genes to their human orthologs using BioMart data. It retrieves
#' specified attributes for the mapped human genes, including chromosome locations, and
#' filters the results to include only canonical autosomes.
#'
#' @param genes A character vector of mouse gene symbols (MGI symbols) to be mapped to human.
#' @param attempt An integer specifying the number of attempts to retrieve the data from BioMart.
#' @param genomes A list of two BioMart objects, the first for the human genome and the second for the mouse genome.
#'
#' @return A tibble containing the mapped human genes and their attributes.
#'         The columns include "mgi_symbol" for the mouse gene symbol and
#'         additional columns specified by `ensemble_columns_to_extract` and
#'         `clean_output_name`. Only genes located on canonical autosomes are included.
#'
#' @examples
#' \dontrun{
#' mouse_genes <- c("Apoc1", "Trem2") # Example mouse gene symbols
#' human_orthologs <- mouse_to_human(mouse_genes)
#' print(human_orthologs)
#' }
#'
#' @export
mouse_to_human <- function(genes, attempt = 0, genomes = NULL) {
  # Load the BioMart data
  if(attempt > 0) {
    host <- read_config()$mart_mirrors[attempt]
  } else {
    host <- NULL
  }

  if (is.null(genomes)){
    genomes <- set_genomes(host)
  }
  human <- genomes[[1]]
  mouse <- genomes[[2]]

  cols_to_get <- read_config()$column_labels$ensemble_columns_to_extract
  clean_output_name <- read_config()$column_labels$clean_output_name

  # getting full list of human with chromosome locations too
  message("Querying ensemble...", appendLF = FALSE)
  tryCatch(
    expr = {
      human_from_mouse <- biomaRt::getLDS(
        attributes = c("mgi_symbol"),
        filters = "mgi_symbol",
        values = genes, mart = mouse,
        attributesL = cols_to_get,
        martL = human, uniqueRows = T,
        verbose = FALSE)
    }, error = function(e) {
      attempt <- attempt + 1
      message("failed.")
      message(conditionMessage(e))
      if(attempt < length(read_config()$mart_mirrors)) {
        human_from_mouse <- mouse_to_human(genes, attempt)
      } else {
        stop("Failed to retrieve data from BioMart.")
      }
    }
  )
  message("success.")

  human_from_mouse <- tibble::as_tibble(human_from_mouse)
  cols_clean <- c("mgi_symbol", clean_output_name)
  colnames(human_from_mouse) <- cols_clean

  human_from_mouse <- force_canonical_autosomes(human_from_mouse)

  return(human_from_mouse)
}
