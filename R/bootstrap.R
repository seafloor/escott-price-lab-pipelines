#' Calculate proportion of regions covered by SNPs
#'
#' @param n_in_regions A vector of the number of SNPs in each region
#' @param region_lengths A vector of the lengths of each region
#'
#' @return A proportion of the regions covered by SNPs
#'
#' @examples
#' calculate_proportion(c(50, 40, 10), c(6000, 3000, 1000))
#'
#' @export
calculate_proportion <- function(n_in_regions, region_lengths) {
  proportion_in_regions <- sum(n_in_regions) / sum(region_lengths)

  return(proportion_in_regions)
}

#' Sample a region from the reference genome
#'
#' @param region_length The length of the region to sample
#' @param chr_reference A data frame containing the reference genome coordinates
#'
#' @return A vector containing the chromosome, start, and end of the region
#'
#' @examples
#' ref_coords <- read_ref_genome_coordinates()
#' sample_region(1000, ref_coords)
#'
#' @export
sample_region <- function(region_length, chr_reference) {
  chr <- sample(1:22, size = 1, prob = chr_reference[["chr_weights"]])

  bp_max_start <- dplyr::filter(chr_reference, `Chromosome name` == chr)[["Seq length"]] - region_length

  bp_start <- sample(1:bp_max_start, size = 1)
  bp_end <- bp_start + region_length

  return(c(chr, bp_start, bp_end))
}

#' Sample multiple regions from the reference genome
#'
#' @param n_regions The number of regions to sample
#' @param min_length The minimum length of a region
#' @param max_length The maximum length of a region
#'
#' @return A data frame containing the sampled regions
#'
#' @examples
#' sample_dummy_regions(n_regions = 10, min_length = 1000, max_length = 2000)
#'
#' @export
sample_dummy_regions <- function(n_regions = 100, min_length = 1000, max_length = 500000) {
  ref_coords <- read_ref_genome_coordinates()
  m <- matrix(NA, nrow = n_regions, ncol = 3)

  for (i in 1:n_regions) {
    l <- sample(min_length:max_length, 1)
    m[i, ] <- sample_region(l, ref_coords)
  }

  m <- tibble::as_tibble(as.data.frame(m))
  colnames(m) <- c("CHR", "BP_START", "BP_END")

  return(m)
}

#' Sample multiple regions from the reference genome given a vector of region lengths
#'
#' @param region_lengths A vector of the lengths of the regions to sample
#'
#' @return A data frame containing the sampled regions
#'
#' @examples
#' sample_regions_from_list(c(1000, 2000, 3000))
#'
#' @export
sample_regions_from_list <- function(region_lengths) {
  ref_coords <- read_ref_genome_coordinates(warn = FALSE)
  n_regions <- length(region_lengths)
  m <- matrix(NA, nrow = n_regions, ncol = 3)

  for (i in 1:n_regions) {
    m[i, ] <- sample_region(region_lengths[i], ref_coords)
  }

  m <- tibble::as_tibble(as.data.frame(m))
  colnames(m) <- c("chr", "gene_start", "gene_end")

  return(m)
}

#' Get the coverage percentage for a chromosome based on a set of regions
#'
#' @param data A data frame containing the regions
#' @param chr_reference A data frame containing the reference genome coordinates
#' @param chrom The chromosome to calculate coverage for
#'
#' @return The coverage percentage for the chromosome
#'
#' @examples
#' ref_coords <- read_ref_genome_coordinates()
#' data <- sample_dummy_regions(n_regions = 100, min_length = 1000, max_length = 2000)
#' get_coverage_for_chromosome(data, ref_coords, chrom = 1)
#'
#' @export
get_coverage_for_chromosome <- function(data, chr_reference, chrom = 1) {
  chr_df <- dplyr::filter(data, CHR == chrom)
  all_coords <- c()
  for (i in 1:nrow(chr_df)) {
    all_coords <- append(all_coords, chr_df[[i, 2]]:chr_df[[i, 3]])
  }

  bp_covered <- length(unique(all_coords))
  coverage <- (bp_covered / dplyr::filter(chr_reference, `Chromosome name` == chrom)[["Seq length"]]) * 100

  return(coverage)
}

#' Check if any SNPs are in a given region
#'
#' @param region_row A row from a data frame containing the region to check
#' @param snps_to_check A data frame containing the SNPs to check
#' @param upstream_bp The number of base pairs upstream of the region to include
#' @param downstream_bp The number of base pairs downstream of the region to include
#'
#' @return The number of SNPs in the region
#'
#' @examples
#' region_row <- tibble::tibble(
#'  chr = 1,
#'  gene_start = 100000,
#'  gene_end = 200000
#' )
#' snps_to_check <- tibble::tibble(
#'  CHR = c(1, 1, 2, 3),
#'  POS = c(100000, 200000, 150000, 250000)
#' )
#'
#' check_if_snps_in_region(region_row, snps_to_check)
#'
#' @export
check_if_snps_in_region <- function(region_row, snps_to_check, upstream_bp = 35000,
                                    downstream_bp = 10000) {
  row_chr_int <- as.integer(region_row["chr"])
  row_bp_start_int <- as.integer(region_row["gene_start"])
  row_bp_end_int <- as.integer(region_row["gene_end"])

  n_top_snps_in_region <- snps_to_check %>%
    dplyr::filter(
      CHR == row_chr_int,
      POS >= (row_bp_start_int - upstream_bp),
      POS <= (row_bp_end_int + downstream_bp)
    ) %>%
    nrow()

  return(n_top_snps_in_region)
}
