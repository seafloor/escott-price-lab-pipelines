library(biomaRt)
source(here("R", "utilities", "standard_reference.R"))
source(here("R", "utilities", "utils.R"))

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

  return(as_tibble(regions))
}
