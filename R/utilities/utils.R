library(tidyr)
library(dplyr)
library(stringr)
library(httr)
library(here)
library(DescTools)
source(here("R", "utilities", "standard_reference.R"))

download_supplement <- function(path, temp_out_path = "") {
  # hard coding temp dir for now
  data_dir <- here("data", "temp")
  if (!dir.exists(data_dir)) {
    dir.create(data_dir)
  }

  if (temp_out_path == "") {
    temp_out_path <- tempfile()
  }
  download.file(path, temp_out_path)

  files <- list.files(temp_out_path)
  for (i in 1:length(files)) {
    f <- files[i]
    splitter <- SplitPath(f)

    if (splitter$extension == "zip") {
      zip_file <- file.path(temp_out_path, f)
      unzip(zip_file, exdir = data_dir)
    } else {
      file.copy(file.path(temp_out_path, f), data_dir)
    }
  }

  unlink(temp_out_path, recursive = TRUE)

  return(data_dir)
}

cleanup_files <- function(file_list) {
  for (i in 1:length(file_list)) {
    fn <- file_list[i]

    if (file.exists(fn)) {
      file.remove(fn)
    } else {
      print(fn)
    }
  }
}

list_files <- function(range, prepend, append) {
  file_list <- sapply(range, paste_names, prepend, append)

  return(file_list)
}

paste_names <- function(i, pre, post) {
  return(paste(pre, i, post, sep = ""))
}

force_canonical_autosomes <- function(df) {
  valid_chromosomes <- as.character(1:22)
  df <- df |>
    dplyr::filter(chr %in% valid_chromosomes) |>
    dplyr::mutate(chr = as.integer(chr))

  return(df)
}

set_human_genome <- function(human_build = "grch38") {
  if (human_build == "grch38") {
    human <- useMart(
      biomart = "ensembl", dataset = "hsapiens_gene_ensembl",
      host = human_genome_38_host_path
    )
  } else if (human_build == "grch37") {
    human <- useMart(
      biomart = "ensembl", dataset = "hsapiens_gene_ensembl",
      host = human_genome_37_host_path
    )
  } else {
    stop("Human genome build must be grch37 or grch38")
  }

  return(human)
}

set_genomes <- function() {
  human <- set_human_genome()
  mouse <- useMart(
    biomart = "ensembl", dataset = "mmusculus_gene_ensembl",
    host = mouse_genome_host_path
  )

  return(list(human, mouse))
}

read_ref_genome_coordinates <- function() {
  warning("Reference genome coordinates are hard-coded to GRCh38.p14")

  ref <- read_tsv(here("data", "annotations", "grch38p14_sequence_report.tsv"),
    n_max = 22,
    col_types = cols_only(
      `Chromosome name` = col_integer(),
      `Seq length` = col_integer()
    )
  ) %>%
    mutate(chr_weights = `Seq length` / `Seq length`[1])

  return(ref)
}

exclude_apoe_region <- function(data) {
  apoe_coordinates <- apoe_region

  genes_to_drop <- data |>
    dplyr::filter(chr == apoe_coordinates["chromosome"]) |>
    dplyr::select(gene_ensemble, gene_start, gene_end) |>
    tidyr::pivot_longer(-gene_ensemble,
      names_to = "gene_boundary",
      values_to = "coordinates"
    ) |>
    dplyr::filter(
      coordinates >= apoe_coordinates["start"],
      coordinates <= apoe_coordinates["end"]
    ) |>
    dplyr::select(gene_ensemble) |>
    pull() |>
    unique()

  data <- dplyr::filter(data, !gene_ensemble %in% genes_to_drop)

  return(data)
}


read_bim_file <- function(f) {
  bim <- read_tsv(f,
    col_names = c("chr", "id", "cm", "pos", "alt", "ref"),
    col_types = list(
      chr = col_character(),
      id = col_character(),
      cm = col_integer(),
      pos = col_integer(),
      alt = col_character(),
      ref = col_character()
    )
  )

  return(bim)
}

# read in file for pathways
# add gene length to check for enrichment
read_regions_to_search <- function(f) {
  if (f == "dummy_region") {
    regions <- tribble(
      ~hgnc_symbol, ~gene_ensemble, ~chr, ~gene_start, ~gene_end, ~gene_length_bp,
      HGNCTESTCHR19, ENSIDTESTCHR19, 19, 1, 50000000, 50000000
    )
  } else {
    regions <- read_csv(f) |>
      mutate(gene_length_bp = gene_start - gene_end)
  }

  return(regions)
}
