install_dependencies <- function() {
  # install bioconductor and biomaRt
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  required_bioc_packages <- c("biomaRt")

  missing_packages <- required_bioc_packages[!sapply(required_bioc_packages, requireNamespace, quietly = TRUE)]

  if (length(missing_packages) > 0) {
    message("The following Bioconductor packages are required but not installed: ", paste(missing_packages, collapse = ", "))
    BiocManager::install(missing_packages)
  }
}

install_databases <- function() {
  warning("Function install_databases() is not implemented yet")
  # prompt user for install location if not /data
  # update the toml if location not /data

  # read toml file with up-to-date database information

  # call separate install database functions
}

download_reference_genome_coordinates <- function() {

}

download_sequence_report <- function() {

}

download_reference_genome <- function() {

}


download_dbsnp <- function(file_url, version = "156") {
  print("--> Installing dbSNP")
  warning("Will use approximately 25GB of disk space")

  # download dbSNP

  # annotate human-readable chromosome names


}

download_liftover <- function() {

}

download_chain_files <- function() {

}
