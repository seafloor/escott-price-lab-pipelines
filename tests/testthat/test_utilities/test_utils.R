library(testthat)
library(mockery)
library(here)
source(here("R", "utilities", "utils.R"))

test_that("download_supplement downloads and unzips files", {
  # Arrange
  # Mock the download.file function to not actually perform a download
  test_zip_path <- here("tests", "testthat", "test_data", "test_files.zip")
  temp_data_path <- here("data", "temp")

  # make temp out_dir
  output_dir <- tempfile()
  dir.create(output_dir)

  # mock download.file so it copies zip instead of accessing a URL
  mock_download <- function(url, destfile, ...) {
    file.copy(test_zip_path, destfile)
  }
  stub(download_supplement, "download.file", mock_download)

  # Act
  # Call function with fake URL included since it won't actually download
  download_supplement("http://fakeurl.com/fake.zip", output_dir)

  # Assert
  expect_true(file.exists(file.path(temp_data_path, "file1.xlsx")))
  expect_true(file.exists(file.path(temp_data_path, "file2.xlsx")))
  expect_true(file.exists(file.path(temp_data_path, "file3.xlsx")))

  unlink(output_dir, recursive = TRUE)
  unlink(temp_data_path, recursive = TRUE)
})

test_that("cleanup_files removes files", {
  # Arrange
  file_list <- c("file1.csv", "file2.csv", "file3.csv")
  sapply(file_list, file.create)

  # Assert
  expect_true(file.exists("file1.csv"))
  expect_true(file.exists("file2.csv"))
  expect_true(file.exists("file3.csv"))

  # Act
  cleanup_files(file_list)

  # Assert
  expect_false(file.exists("file1.csv"))
  expect_false(file.exists("file2.csv"))
  expect_false(file.exists("file3.csv"))
})

test_that("list_files returns a list of file names", {
  # Arrange
  range <- 1:3
  prepend <- "file"
  append <- ".txt"

  # Act
  result <- list_files(range, prepend, append)

  # Assert
  expect_equal(result, c("file1.txt", "file2.txt", "file3.txt"))
})

test_that("paste_names returns a string with the given prefix and suffix", {
  # Arrange
  i <- 1
  pre <- "file"
  post <- ".txt"

  # Act
  result <- paste_names(i, pre, post)

  # Assert
  expect_equal(result, "file1.txt")
})

test_that("force_canonical_autosomes filters out non-canonical chromosomes", {
  # Arrange
  df <- data.frame(chr = c(
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
    "11", "12", "13", "14", "15", "16", "17", "18", "19",
    "20", "21", "22", "X", "Y", "MT"
  ))

  # Act
  result <- force_canonical_autosomes(df)

  # Assert
  expect_equal(result$chr, c(1:22))
})

test_that("set_human_genome returns a biomart object", {
  # Arrange
  human_build <- "grch38"

  # Act
  result <- set_human_genome(human_build)

  # Assert
  expect_true(inherits(result, "Mart"))
})

test_that("set_human_genome returns a biomart object for grch37", {
  # Arrange
  human_build <- "grch37"

  # Act
  result <- set_human_genome(human_build)

  # Assert
  expect_true(inherits(result, "Mart"))
})

test_that("set_human_genome uses default host when passed grch38 argument", {
  # Arrange
  human_build <- "grch38"
  human_genome_default_host_path <- "https://mart.ensembl.org:443/biomart/martservice"

  # Act
  result <- set_human_genome(human_build)

  # Assert
  expect_equal(result@host, human_genome_default_host_path)
})

test_that("set_human_genome uses grch37 host when passed grch37 argument", {
  # Arrange
  human_build <- "grch37"
  human_genome_37_host_path <- "https://grch37.ensembl.org:443/biomart/martservice"

  # Act
  result <- set_human_genome(human_build)

  # Assert
  expect_equal(result@host, human_genome_37_host_path)
})

test_that("set_genomes returns a list of biomart objects", {
  # Act
  result <- set_genomes()

  # Assert
  expect_true(inherits(result, "list"))
  expect_true(inherits(result[[1]], "Mart"))
  expect_true(inherits(result[[2]], "Mart"))
})

test_that("set_genomes returns human and mouse genomes", {
  # Act
  result <- set_genomes()

  # Assert
  expect_true(result[[1]]@dataset == "hsapiens_gene_ensembl")
  expect_true(result[[2]]@dataset == "mmusculus_gene_ensembl")
})

test_that("read_ref_genome_coordinates returns a tibble", {
  # Act
  result <- read_ref_genome_coordinates()

  # Assert
  expect_true(inherits(result, "tbl_df"))
})

test_that("exclude_apoe_region returns a tibble", {
  # Arrange
  data <- tribble(
    ~hgnc_symbol, ~gene_ensemble, ~chr, ~gene_start, ~gene_end,
    "MyGene", "ENSG00000000001", 1, 1000, 2000
  )

  # Act
  result <- exclude_apoe_region(data)

  # Assert
  expect_true(inherits(result, "tbl_df"))
})

test_that("read_bim_file returns a tibble", {
  # Arrange
  f <- here("tests", "testthat", "test_data", "test_bim.bim")

  # Act
  result <- read_bim_file(f)

  # Assert
  expect_true(inherits(result, "tbl_df"))
})

test_that("read_regions_to_search returns a tibble", {
  # Arrange
  f <- here("output", "pathways", "gosselin_et_al_2017_microglia_grch38.csv")

  # Act
  result <- read_regions_to_search(f)

  # Assert
  expect_true(inherits(result, "tbl_df"))
})
