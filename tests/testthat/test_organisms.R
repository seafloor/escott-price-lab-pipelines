
test_that("mouse_to_human correctly processes mocked gene conversion", {
  # Arrange
  # Define mock output of biomaRt::getLDS for TREM2 and APOC1
  # Definitions extracted manually from the links:
  # APOC1 = https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000130208;r=19:44914247-44919349
  # TREM2 = https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000095970;r=6:41158506-41163186
  expected_output <- tibble::tribble(
    ~mgi_symbol, ~hgnc_symbol, ~gene_ensemble, ~chr, ~gene_start, ~gene_end,
    "Apoc1", "APOC1", "ENSG00000130208", 19, 44914247, 44919349,
    "Trem2", "TREM2", "ENSG00000095970", 6, 41158506, 41163186
  )

  # convert to df with wrong column labels to show the output is tibble with correct labels
  mock_output <- as.data.frame(expected_output)
  colnames(mock_output) <- c("mgi_wrongname", "hgnc_wrongname",
                             "ensemble_wrongname", "chr_wrongname",
                             "start_wrongname", "end_wrongname")

  # Mock biomaRt::getLDS to return the predefined mock output
  mockery::stub(mouse_to_human, "biomaRt::getLDS", mock_output)

  # Act
  # Call the function with a mocked input
  mouse_genes <- c("Apoc1", "Trem2")
  test_output <- mouse_to_human(mouse_genes, genomes = list(NULL, NULL))

  # Assert
  # Validate the output
  expect_true(tibble::is_tibble(test_output), info = "Output should be a tibble.")
  expect_equal(nrow(test_output), nrow(expected_output),
               info = "The output should have the same number of rows as the expected output.")
  expect_equal(colnames(test_output), colnames(expected_output),
               info = "Column names of output should match the expected output.")

  # Additional checks can be added here to validate the contents of the output
})

