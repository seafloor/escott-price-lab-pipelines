#' Get chromosomal coordinates for a list of genes.
#'
#' @param genes A character vector of gene names.
#' @param gene_format A character string specifying the gene format. Default is "hgnc_symbol".
#' @param build A character string specifying the genome build. Default is "grch38".
#'
#' @return A tibble with columns "gene", "chromosome", "start", "end", "strand".
#'
#' @examples
#' \dontrun{
#' genes <- c("PABPC4;HEYL", "PCSK9")
#' flat_genes <- unlist(stringr::str_split(genes, ";"))
#' coords <- get_regions_from_genes(flat_genes, gene_format = "hgnc_symbol")
#' }
#'
#' @export
get_regions_from_genes <- function(genes, gene_format = "hgnc_symbol",
                                   build = "grch38") {
  mymart <- set_human_genome(build)

  cols_to_get <- ensemble_columns_to_extract
  cols_clean <- clean_output_name

  regions <- biomaRt::getBM(
    attributes = cols_to_get,
    mart = mymart,
    filters = gene_format,
    values = genes
  )

  colnames(regions) <- cols_clean

  regions <- force_canonical_autosomes(regions)

  return(tibble::as_tibble(regions))
}
