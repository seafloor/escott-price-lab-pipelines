library(biomaRt)
source("standard_reference.R")
source("utils.R")

get_regions <- function(genes, gene_format = "hgnc_symbol") {
  cols_to_get <- ensemble_columns_to_extract
  
  regions <- biomaRt::getBM(attributes = cols_to_get,
                            mart = mymart,
                            filters = gene_format,
                            values = genes)
  
  return(as_tibble(genes))
}
