#' Check if the File is a Zip Archive
#'
#' Determines whether the specified file has a `.zip` extension.
#'
#' @param f The file name to check.
#'
#' @return `TRUE` if the file is a zip archive, otherwise `FALSE`.
#'
#' @examples
#' \dontrun{
#' is_zip("example.zip") # returns TRUE
#' is_zip("example.txt") # returns FALSE
#' }
#'
#' @export
is_zip <- function(f) {
  if (stringr::str_detect(f, "\\.zip$")) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Download and Unzip Supplemental Files
#'
#' Downloads a file from the specified path and unzips it if necessary. If the
#' output path is not specified, a temporary path is used.
#'
#' @param path The URL of the file to download.
#' @param temp_out_path Optional; the output path for the downloaded file.
#'
#' @return The path to the downloaded (and possibly unzipped) files.
#'
#' @export
download_supplement <- function(path, temp_out_path = "") {
  if (temp_out_path == "") {
    temp_out_path <- tempfile()
  }
  download.file(path, temp_out_path)

  files <- list.files(temp_out_path)
  for (i in 1:length(files)) {
    f <- files[i]

    if (is_zip(f)) {
      zip_file <- file.path(temp_out_path, f)
      unzip(zip_file, exdir = temp_out_path)
    }
  }

  return(temp_out_path)
}

#' Generate a List of Files
#'
#' Creates a list of file names based on a specified range and naming pattern.
#'
#' @param range A vector of values to create file names for.
#' @param prepend Text to prepend to each file name.
#' @param append Text to append to each file name.
#'
#' @return A vector of constructed file names.
#'
#' @export
list_files <- function(range, prepend, append) {
  file_list <- sapply(range, paste_names, prepend, append)

  return(file_list)
}

#' Concatenate Names with Prefix and Suffix
#'
#' Concatenates a given index with prefix and suffix strings.
#'
#' @param i Index to be included between prefix and suffix.
#' @param pre Prefix for the name.
#' @param post Suffix for the name.
#'
#' @return A string resulting from the concatenation of prefix, index, and suffix.
#'
#' @export
paste_names <- function(i, pre, post) {
  return(paste(pre, i, post, sep = ""))
}

#' Filter and Adjust Chromosome Data for Canonical Autosomes
#'
#' Filters a dataframe to include only canonical autosomes (1-22) and converts
#' chromosome identifiers to integers. Assumes input chromosomes are strings.
#'
#' @param df A dataframe containing chromosome data as strings.
#'
#' @return A dataframe filtered for canonical autosomes with integer chromosome identifiers.
#'
#' @export
force_canonical_autosomes <- function(df) {
  valid_chromosomes <- as.character(1:22)
  df <- df %>%
    dplyr::filter(chr %in% valid_chromosomes) %>%
    dplyr::mutate(chr = as.integer(chr))

  return(df)
}


#' Set up BioMart connection
#'
#' Initializes a connection to a specified BioMart database and dataset.
#'
#' @param dataset The dataset to connect to within the BioMart database.
#' @param host The URL of the BioMart host to connect to.
#'
#' @return An object of class Mart representing the BioMart connection.
#' 
#' @examples
#' \dontrun{
#'   mart <- set_mart(dataset = "hsapiens_gene_ensembl", host = "http://www.ensembl.org/biomart")
#' }
#' 
#' @export
set_mart <- function(dataset, host){
  biomaRt::useMart(
    biomart = "ensembl", dataset = dataset,
    host = host
  )
}

#' Attempt to establish BioMart connection with error handling
#'
#' Tries to connect to a specified BioMart database and dataset, with
#' error handling for HTTP connection issues. If the initial connection fails
#' due to a server-side error, it attempts to connect to a mirror site.
#'
#' @param dataset The dataset to connect to within the BioMart database.
#' @param host The URL of the BioMart host to connect to initially.
#'
#' @return An object of class Mart representing the BioMart connection, or
#' `NA` if the connection could not be established after error handling.
#' 
#' @examples
#' \dontrun{
#'   mart <- try_mart(dataset = "hsapiens_gene_ensembl", host = "http://www.ensembl.org/biomart")
#' }
#' 
#' @export
try_mart <- function(dataset, host) {
  tryCatch(
    expr = {
      message("Trying ensembl mart")
      mart <- set_mart(dataset, host)
    },
    error = function(e) {
      http_error = as.integer(
        stringr::str_extract(conditionMessage(e),
                             pattern = "\\(HTTP ([0-9]{3})\\)",
                             group = 1)
      )
      print(http_error)
      message("Error in connection to Ensembl:")
      message(conditionMessage(e))
      
      if(is.na(http_error)) {
        stop("Unknown error")
      }
      if(http_error >= 500) {
        message(paste("Ensembl server-side error with host:", host))
        if (stringr::str_detect(host, "grch37")) {
          stop("Client-side error")
        }
        
        host <- stringr::str_replace(host, "mart", "useast")
        
        message(paste("Trying new mirror:", host))
        mart <- set_mart(dataset, host)
      } else if (http_error >= 400 && http_error < 500) {
        message("Client-side error. Check call/connection and retry")
      } else {
        message("Unknown error")
      }
      
      NA
    }
  )
  
  return(mart)
}

#' Set Up Human Genome Annotation
#'
#' Configures the human genome annotation using BioMart based on the specified genome build.
#'
#' @param human_build The human genome build to use, either "grch38" or "grch37".
#'
#' @return A BioMart object configured for the specified human genome build.
#'
#' @export
set_human_genome <- function(human_build = "grch38") {
  if (human_build == "grch38") {
    human <- try_mart(dataset = "hsapiens_gene_ensembl",
                      host = human_genome_38_host_path)
  } else if (human_build == "grch37") {
    human <- try_mart(dataset = "hsapiens_gene_ensembl",
                      host = human_genome_37_host_path)
  } else {
    stop("Human genome build must be grch37 or grch38")
  }

  return(human)
}

#' Set Up Human and Mouse Genome Annotations
#'
#' Configures BioMart objects for both human and mouse genome annotations.
#'
#' @return A list containing BioMart objects for human and mouse genomes.
#'
#' @export
set_genomes <- function() {
  human <- set_human_genome()
  mouse <- try_mart(dataset = "mmusculus_gene_ensembl",
                    host = mouse_genome_host_path)

  return(list(human, mouse))
}

#' Read Reference Genome Coordinates
#'
#' Reads hardcoded reference genome coordinates, specifically for GRCh38.p14.
#'
#' @return A dataframe containing reference genome coordinates.
#'
#' @export
read_ref_genome_coordinates <- function() {
  warning("Reference genome coordinates are hard-coded to GRCh38.p14")

  reference_file <- "grch38p14_sequence_report.tsv"
  reference_path <- system.file("extdata", "annotations",
    reference_file,
    package = "escottpricelabpipelines"
  )
  if (reference_path == "") {
    stop("File not found in package: ", reference_file)
  }

  ref <- readr::read_tsv(reference_path,
    n_max = 22,
    col_types = readr::cols_only(
      `Chromosome name` = readr::col_integer(),
      `Seq length` = readr::col_integer()
    )
  ) %>%
    dplyr::mutate(chr_weights = `Seq length` / `Seq length`[1])

  return(ref)
}

#' Exclude Genes in APOE Region
#'
#' Filters out genes located in the APOE region based on predefined coordinates.
#'
#' @param data A dataframe containing gene annotations.
#'
#' @return A dataframe excluding genes in the APOE region.
#'
#' @export
exclude_apoe_region <- function(data) {
  apoe_coordinates <- apoe_region

  genes_to_drop <- data %>%
    dplyr::filter(chr == apoe_coordinates["chromosome"]) %>%
    dplyr::select(gene_ensemble, gene_start, gene_end) %>%
    tidyr::pivot_longer(-gene_ensemble,
      names_to = "gene_boundary",
      values_to = "coordinates"
    ) %>%
    dplyr::filter(
      coordinates >= apoe_coordinates["start"],
      coordinates <= apoe_coordinates["end"]
    ) %>%
    dplyr::select(gene_ensemble) %>%
    dplyr::pull() %>%
    unique()

  data <- dplyr::filter(data, !gene_ensemble %in% genes_to_drop)

  return(data)
}

#' Read BIM File
#'
#' Reads a .bim file and extracts SNP information.
#'
#' @param f The filepath to the .bim file.
#'
#' @return A dataframe with SNP information from the .bim file.
#'
#' @export
read_bim_file <- function(f) {
  bim <- readr::read_tsv(f,
    col_names = c("chr", "id", "cm", "pos", "alt", "ref"),
    col_types = list(
      chr = readr::col_character(),
      id = readr::col_character(),
      cm = readr::col_integer(),
      pos = readr::col_integer(),
      alt = readr::col_character(),
      ref = readr::col_character()
    )
  )

  return(bim)
}

#' Read Regions to Search for Enrichment
#'
#' Reads in a file specifying regions to search for enrichment analysis, or
#' generates a dummy region for testing. Adds gene length as a column.
#'
#' @param f The filepath to the regions file or "dummy_region" for a test region.
#'
#' @return A dataframe of regions with calculated gene length in base pairs.
#'
#' @export
read_regions_to_search <- function(f) {
  if (f == "dummy_region") {
    regions <- tibble::tribble(
      ~hgnc_symbol, ~gene_ensemble, ~chr, ~gene_start, ~gene_end, ~gene_length_bp,
      HGNCTESTCHR19, ENSIDTESTCHR19, 19, 1, 50000000, 50000000
    )
  } else {
    regions <- readr::read_csv(f) %>%
      dplyr::mutate(gene_length_bp = gene_start - gene_end)
  }

  return(regions)
}