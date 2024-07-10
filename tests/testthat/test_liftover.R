test_that("to_bed creates a tibble with expected output", {
  # Arrange
  input_data <- tibble::tribble(
    ~chr, ~pos, ~id, ~ref, ~alt,
    "1", 101, "rs12", "A", "G",
    "2", 201, "rs567", "T", "C",
    "3", 301, "rs890", "C", "T"
  )
  expected_data <- tibble::tribble(
    ~chr, ~pos0, ~pos1,
    "chr1", 100, 101,
    "chr2", 200, 201,
    "chr3", 300, 301
  )

  # Act
  result <- to_bed(input_data)

  # Assert
  expect_equal(result, expected_data)
  expect_true(inherits(result, "tbl_df"))
})

test_that("to_bed handles alternative column names", {
  # Arrange
  input_data <- tibble::tribble(
    ~chromosome, ~bp, ~id, ~ref, ~alt,
    "1", 101, "rs12", "A", "G",
    "2", 201, "rs567", "T", "C",
    "3", 301, "rs890", "C", "T"
  )
  expected_data <- tibble::tribble(
    ~chr, ~pos0, ~pos1,
    "chr1", 100, 101,
    "chr2", 200, 201,
    "chr3", 300, 301
  )

  # Act
  result <- to_bed(input_data, chrom = "chromosome", position = "bp")

  # Assert
  expect_equal(result, expected_data)
})

test_that("read_unmapped_variants returns a tibble with expected output", {
  # skip_on_ci()
  
  # Arrange
  f <- testthat::test_path("test_data/liftover_example.unmapped")
  result <- read_unmapped_variants(f)

  # Assert
  expect_true(inherits(result, "tbl_df"))
  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), 4)
  expect_equal(colnames(result), c("CHROM", "POS_0", "POS_1", "chrpos"))
})

test_that("call_liftover returns a list of tibbles", {
  # skip_on_ci()
  
  # Arrange tibble of top AD SNPs
  # grch37 and grch38 coordinates taken manually from dbSNP v156
  df_grch37 <- tibble::tribble(
    ~chr, ~pos, ~id, ~ref, ~alt,
    "19", 45411941, "rs429369", "T", "C",
    "19", 45412079, "rs7412", "C", "T",
    "2", 127892810, "rs6733839", "C", "G",
    "11", 85868640, "rs3851179", "T", "C"
  )

  # Act
  result <- call_liftover(df_grch37, check_liftover_installed = FALSE)

  # Asserts
  expect_true(inherits(result, "list"))
  expect_true(length(result) == 2)
  expect_true(inherits(result[[1]], "tbl_df"))

  # Assert that the second element is NULL as all positions should map correctly
  expect_true(is.null(result[[2]]))
})

test_that("call_liftover correctly maps between grch37 and grch38", {
  # skip_on_ci()
  
  # Arrange tibble of top AD SNPs
  # grch37 and grch38 coordinates taken manually from dbSNP v156
  df_grch37 <- tibble::tribble(
    ~chr, ~pos, ~id, ~ref, ~alt,
    "19", 45411941, "rs429369", "T", "C",
    "19", 45412079, "rs7412", "C", "T",
    "2", 127892810, "rs6733839", "C", "G",
    "11", 85868640, "rs3851179", "T", "C"
  )

  df_grch38 <- tibble::tribble(
    ~chr, ~pos, ~id, ~ref, ~alt,
    "19", 44908684, "rs429369", "T", "C",
    "19", 44908822, "rs7412", "C", "T",
    "2", 127135234, "rs6733839", "C", "G",
    "11", 86157598, "rs3851179", "T", "C"
  )

  # Assert lift to 38 is correct
  result <- call_liftover(df_grch37, from_build = 'grch37', to_build = 'grch38',
                          check_liftover_installed = FALSE)
  expect_equal(result[[1]], df_grch38)
  expect_true(is.null(result[[2]]))

  # Assert lift to 37 is correct
  result <- call_liftover(df_grch38, from_build = 'grch38', to_build = 'grch37',
                          check_liftover_installed = FALSE)
  expect_equal(result[[1]], df_grch37)
  expect_true(is.null(result[[2]]))
})

test_that("call_liftover drops unmapped variants", {
  # skip_on_ci()
  
  # Arrange tibble of top AD SNPs
  # grch37 and grch38 coordinates taken manually from dbSNP v156
  df_grch37 <- tibble::tribble(
    ~chr, ~pos, ~id, ~ref, ~alt,
    "19", 45411941, "rs429369", "T", "C",
    "19", 45412079, "rs7412", "C", "T",
    "22", 60000001, "rs1234", "A", "G",
    "21", 55000001, "rs2345", "T", "C"
  )

  expected_mapped <- tibble::tribble(
    ~chr, ~pos, ~id, ~ref, ~alt,
    "19", 44908684, "rs429369", "T", "C",
    "19", 44908822, "rs7412", "C", "T"
  )

  expected_unmapped <- tibble::tribble(
    ~CHROM, ~POS_0, ~POS_1, ~chrpos,
    "22", 60000000, 60000001, "22_60000001",
    "21", 55000000, 55000001, "21_55000001"
  )

  # Act
  result <- call_liftover(df_grch37, check_liftover_installed = FALSE)

  # Assert
  expect_equal(result[[1]], expected_mapped)
  expect_equal(result[[2]], expected_unmapped)
})
