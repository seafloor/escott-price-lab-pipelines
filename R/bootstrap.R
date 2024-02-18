calculate_proportion <- function(n_in_regions, region_lengths) {
  proportion_in_regions <- sum(n_in_regions) / sum(region_lengths)

  return(proportion_in_regions)
}

sample_region <- function(region_length, chr_reference) {
  chr <- sample(1:22, size = 1, prob = chr_reference[["chr_weights"]])

  bp_max_start <- dplyr::filter(chr_reference, `Chromosome name` == chr)[["Seq length"]] - region_length

  bp_start <- sample(1:bp_max_start, size = 1)
  bp_end <- bp_start + region_length

  return(c(chr, bp_start, bp_end))
}

# sample multiple regions for running checks on coverage
sample_dummy_regions <- function(n_regions = 100, min_length = 1000, max_length = 500000) {
  m <- matrix(NA, nrow = n_regions, ncol = 3)

  for (i in 1:n_regions) {
    l <- sample(min_length:max_length, 1)
    m[i, ] <- sample_region(l, ref_coords)
  }

  m <- tibble::as_tibble(as.data.frame(m))
  colnames(m) <- c("CHR", "BP_START", "BP_END")

  return(m)
}

# sample multiple regions given a vector of region lengths
sample_regions_from_list <- function(region_lengths) {
  n_regions <- length(region_lengths)
  m <- matrix(NA, nrow = n_regions, ncol = 3)

  for (i in 1:n_regions) {
    m[i, ] <- sample_region(region_lengths[i], ref_coords)
  }

  m <- as_tibble(as.data.frame(m))
  colnames(m) <- c("chr", "gene_start", "gene_end")

  return(m)
}

# check coverage by chromosome from regions
get_coverage_for_chromosome <- function(data, chr_reference, chrom = 1) {
  chr_df <- dplyr::filter(data, chr == chrom)
  all_coords <- c()
  for (i in 1:nrow(chr_df)) {
    all_coords <- append(all_coords, chr_df[[i, 2]]:chr_df[[i, 3]])
  }

  bp_covered <- length(unique(all_coords))
  coverage <- (bp_covered / dplyr::filter(chr_reference, `Chromosome name` == chrom)[["Seq length"]]) * 100

  return(coverage)
}

check_if_snps_in_region <- function(region_row, upstream_bp = 35000,
                                    downstream_bp = 10000) {
  row_chr_int <- as.integer(region_row["chr"])
  row_bp_start_int <- as.integer(region_row["gene_start"])
  row_bp_end_int <- as.integer(region_row["gene_end"])

  n_top_snps_in_region <- top_snp_pos |>
    dplyr::filter(
      CHR == row_chr_int,
      POS >= (row_bp_start_int - upstream_bp),
      POS <= (row_bp_end_int + downstream_bp)
    ) |>
    nrow()

  return(n_top_snps_in_region)
}


#### checks
# checking sampling function has balanced coverage
# shows simila coverage across chromosomes but takes too long to re-run
# d <- sample_dummy_regions(n_regions = 10000, min_length = 10000, max_length = 500000)
# cov <- c()
# for (chr in 1:22) {
#   if (chr %in% unique(d[['CHR']])) {
#     cov <- append(cov, get_coverage_for_chromosome(d, ref_coords, chrom = chr))
#   } else {
#     cov <- append(cov, NA)
#   }
# }
#
# cov <- as_tibble(data.frame('chr' = 1:22, 'coverage' = cov))
# cov
#
# cov |>
#   ggplot(aes(x=chr, y=coverage)) + geom_bar(stat = "identity")
