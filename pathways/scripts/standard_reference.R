# contains the defaults for builds, extracted columns from ensemble etc.
# all scripts should read from here to ensure reproducibility

ensemble_columns_to_extract <- c("hgnc_symbol", "ensembl_gene_id",
                                 "chromosome_name", "start_position",
                                 "end_position")

human_genome_38_host_path = "https://mart.ensembl.org/"
human_genome_37_host_path = "https://grch37.ensembl.org/"
mouse_genome_host_path = "https://mart.ensembl.org/"