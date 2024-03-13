#' Install Required Bioconductor Dependencies
#'
#' This function checks for and installs required Bioconductor packages for the package.
#' Currently, it ensures that `biomaRt` is installed, using `BiocManager`.
#'
#' @details The function first checks if `BiocManager` is installed and installs it if necessary.
#' Then, it checks for the presence of specified Bioconductor packages (`biomaRt`) and installs any that are missing.
#' This is particularly useful for setting up the package environment or ensuring that dependencies are met for new users.
#'
#' @examples
#' \dontrun{
#' install_dependencies()
#' }
#'
#' @export
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

#' Install Required Databases
#' 
#' This function downloads and installs the required databases for the package.
#' This must be run interactively, as it requires user input to confirm the download
#' due to the size of some databases.
#' 
#' @export
install_databases <- function() {
  if (!interactive()) {
    message("This process requires user input. Please run interactively.")
    return()
  }
  
  # prompt user to confirm they want to download databases
  message("This will download and install the following databases:")
  message("  - liftover (approx. 50MB)")
  message("Do you want to continue? (y/n)")
  response <- readline()
  if (tolower(response) != "y") {
    message("Aborting.")
    return()
  }

  # install required liftover files
  download_liftover()
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

#' Check for Liftover Files
#' 
#' This function checks for the presence of the required liftover files for the package.
#' If the files are not present, it will attempt to download them by calling download_liftover().
#' If there has already been an attempt to download the files, it will stop and remove the tar.gz file.
#' 
#' @param attempt The number of attempts made to download the liftover files. Default is 0.
#' 
#' @export
check_for_liftover <- function(attempt = 0) {
  # check that liftover files are present and non-empty
  if (!file.exists(here::here("data", "liftover", "hg19ToHg38.over.chain.gz")) ||
      !file.exists(here::here("data", "liftover", "hg38ToHg19.over.chain.gz")) ||
      !file.size(here::here("data", "liftover", "liftOver")) ||
      !file.size(here::here("data", "liftover", "liftOver_linux"))) {
    if (attempt > 0) {
      unlink(here::here("data", "liftover.tar.gz"))
      stop("Liftover files not present. Download failed.")
    } else {
      download_liftover(attempt = 1)
    }
  }
}

#' Install Required Liftover Files
#' 
#' This function downloads the required liftover files for the package,
#' including chain files and the max/linux executables.
#' 
#' @export
download_liftover <- function(attempt = 1) {
  # get liftover download info
  db_info <- read_database_toml()
  
  # download liftover
  message("--> Downloading liftover and chain files (approx. 50Mb space after extracting)...")
  download.file(db_info$liftover$url, destfile = here::here("data", "liftover.tar.gz"))
  md5 <- tools::md5sum(here::here("data", "liftover.tar.gz"))
  if (md5 != db_info$liftover$checksum) {
    unlink(here::here("data", "liftover.tar.gz"))
    stop("Downloaded file does not match expected md5sum")
  } else {
    message("done.")
  }
  
  # extract and check liftover files
  message("--> Extracting liftover files...", appendLF = FALSE)
  tryCatch({
    untar(here::here("data", "liftover.tar.gz"), exdir = here::here("data", "liftover"))
  }, error = function(e) {
    unlink(here::here("data", "liftover.tar.gz"))
    stop("Error extracting liftover files: ", conditionMessage(e))
  })
  
  # check liftover files
  check_for_liftover(attempt = attempt)
  
  # remove tar.gz
  unlink(here::here("data", "liftover.tar.gz"))
  
  message("done.")
}
