# contains the defaults for builds, extracted columns from ensemble etc.
# all scripts should read from here to ensure reproducibility

ensemble_columns_to_extract <- c(
  "hgnc_symbol", "ensembl_gene_id",
  "chromosome_name", "start_position",
  "end_position"
)
clean_output_name <- c(
  "hgnc_symbol",
  "gene_ensemble", "chr",
  "gene_start", "gene_end"
)

human_genome_38_host_path <- "https://mart.ensembl.org/"
human_genome_37_host_path <- "https://grch37.ensembl.org/"
mouse_genome_host_path <- "https://mart.ensembl.org/"

apoe_region <- c("chromosome" = 19, "start" = 44400000, "end" = 46500000)
