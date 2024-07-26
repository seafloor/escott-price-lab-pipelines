# Sample data for testing
top_snp_pos <- tibble::tibble(
  CHR = sample(1:22, 100, replace = TRUE),
  POS = sample(1:250000000, 100)
)

test_that("calculate_proportion calculates correctly", {
  n_in_regions <- c(50, 40, 10)
  region_lengths <- c(6000, 3000, 1000)
  expect_equal(calculate_proportion(n_in_regions, region_lengths), 0.01)
})

test_that("sample_region samples correctly", {
  set.seed(42)
  ref_coords <- read_ref_genome_coordinates()
  
  region <- sample_region(1000, ref_coords)
  
  expect_equal(length(region), 3)
  expect_true(region[1] %in% 1:22)
  expect_true(region[2] >= 1)
  expected_end <- ref_coords[ref_coords[['Chromosome name']] == region[1], ][["Seq length"]]
  expect_true(region[3] <= expected_end)
  expect_true(region[3] >= region[2])
  expect_true(region[3] - region[2] == 1000)
})

test_that("sample_dummy_regions samples correctly", {
  set.seed(42)
  regions <- sample_dummy_regions(n_regions = 10, min_length = 1000, max_length = 2000)
  expect_equal(nrow(regions), 10)
  expect_equal(ncol(regions), 3)
  expect_true(all(regions$CHR %in% 1:22))
  expect_true(all(regions$BP_START < regions$BP_END))
})

test_that("sample_regions_from_list samples correctly", {
  set.seed(42)
  region_lengths <- c(1000, 2000, 3000)
  regions <- sample_regions_from_list(region_lengths)
  expect_equal(nrow(regions), length(region_lengths))
  expect_equal(ncol(regions), 3)
  expect_true(all(regions$chr %in% 1:22))
  expect_true(all(regions$gene_start < regions$gene_end))
})

test_that("get_coverage_for_chromosome calculates percentages", {
  set.seed(42)
  chrom <- 1
  ref_coords <- read_ref_genome_coordinates()
  
  data <- sample_dummy_regions(n_regions = 100, min_length = 1000, max_length = 2000)
  coverage <- get_coverage_for_chromosome(data, ref_coords, chrom)
  
  expect_true(coverage >= 0 && coverage <= 100)
})

test_that("check_if_snps_in_region works correctly", {
  region_row <- tibble::tibble(
    chr = 1,
    gene_start = 100000,
    gene_end = 200000
  )
  snps_to_check <- tibble::tibble(
    CHR = c(1, 1, 1, 1, 1, 1, 1),
    POS = c(50000, 90000, 100000, 150000, 200000, 200010, 250000)
  )
  
  snps_in_region <- check_if_snps_in_region(snps_to_check, region_row)
  expect_true(snps_in_region == 5)
})

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
# cov %>%
#   ggplot(aes(x=chr, y=coverage)) + geom_bar(stat = "identity")
