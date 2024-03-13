test_that("to_bed creates a tibble with expected output", {
  # Arrange
  input_data <- tibble::tribble(
    ~chr, ~pos, ~id, ~ref, ~alt,
    "1", 100, "rs12", "A", "G",
    "2", 200, "rs567", "T", "C",
    "3", 300, "rs890", "C", "T"
  )
  expected_data <- tibble::tribble(
    ~chr, ~pos0, ~pos1,
    "1", 100, 100,
    "2", 200, 200,
    "3", 300, 300
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
    "1", 100, "rs12", "A", "G",
    "2", 200, "rs567", "T", "C",
    "3", 300, "rs890", "C", "T"
  )
  expected_data <- tibble::tribble(
    ~chr, ~pos0, ~pos1,
    "1", 100, 100,
    "2", 200, 200,
    "3", 300, 300
  )
  
  # Act
  result <- to_bed(input_data, chrom = "chromosome", position = "bp")
  
  # Assert
  expect_equal(result, expected_data)
})