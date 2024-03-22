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

  config <- read_config()

  cols_to_get <- config$column_labels$ensemble_columns_to_extract
  cols_clean <- config$column_labels$clean_output_name

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

parse_gff <- function(f) {
  # colnames/types taken from https://www.ensembl.org/info/website/upload/gff.html
  gff_colnames <- c("seqname", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attribute")
  
  # Define column types for all nine columns explicitly
  gff_coltypes <- readr::cols(
    seqname = readr::col_character(),
    source = readr::col_character(),
    feature = readr::col_character(),
    start = readr::col_integer(),
    end = readr::col_integer(),
    score = readr::col_double(),
    strand = readr::col_character(),
    frame = readr::col_character(),
    attribute = readr::col_character()
  )
  
  # Initialize an empty list to store chunks
  results <- list()
  
  # Define the chunk processing function
  process_chunk <- function(df, pos) {
    filtered <- df %>% dplyr::filter(feature == "gene")
    # Append filtered chunk to results list
    results <<- c(results, list(filtered))
  }
  
  # Read the file in chunks
  readr::read_tsv_chunked(file = f, 
                          col_types = gff_coltypes, 
                          col_names = gff_colnames, 
                          comment = "#",
                          chunk_size = 10000, # Adjust chunk size as needed
                          callback = readr::SideEffectChunkCallback$new(process_chunk))
  
  # Combine the chunks
  combined_data <- dplyr::bind_rows(results)
  
  # Now process the attribute column
  final_data <- combined_data %>%
    tidyr::separate_rows(attribute, sep = ";") %>%
    tidyr::separate(attribute, into = c("key", "value"), sep = "=") %>%
    tidyr::pivot_wider(names_from = key, values_from = value, values_fill = list(value = NA)) %>%
    dplyr::select(-score, -strand, -frame)
  
  return(final_data)
}