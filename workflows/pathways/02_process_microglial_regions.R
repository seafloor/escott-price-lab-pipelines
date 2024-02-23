# load packages with box
box::use(eplp = escottpricelabpipelines[get_regions_from_genes, exclude_apoe_region],
         rdr = readr[read_csv, write_csv],
         here = here[here],
         stringr = stringr[str_detect])

# load data from Galatro et al. 2017 and Gosselin et al. 2017
in_file_galatro_2017 <- here$here("output", "raw_data",
                                  "galatro_et_al_2017_microglia_supp_table_2.csv")
in_file_gosselin_2017 <- here$here("output", "raw_data",
                                   "gosselin_et_al_2017_microglia_supp_table_2.csv")

# max of 1296 required for Galatro et al. 2017 to match
# the highlighted rows in the original file for microglial core genes
raw_df_galatro <- rdr$read_csv(in_file_galatro_2017, n_max = 1296)
raw_df_gosselin <- rdr$read_csv(in_file_gosselin_2017)

# check overlap between the two datasets
matches <- c()
for(i in 1:nrow(raw_df_gosselin)) {
  matches <- c(matches, sum(stringr$str_detect(
    raw_df_galatro[['external_gene_name']],
    paste("^", raw_df_gosselin[['gene_name']][i],
          "$", sep = ""))))
}

print(paste("Number of genes in Gosselin et al. 2017 ",
            "that are also in Galatro et al. 2017: ",
            sum(matches > 0), sep = ""))

# annotate coordinates
# note that this causes a drop in genes from Gosselin et al. 2017
# from 881 to 761 so re-processing of the raw_df_gosselin may
# be needed if you expect to match non-canonical microglial genes
galatro_annotated <- eplp$get_regions_from_genes(raw_df_galatro[['ensembl_gene_id']],
                                                 gene_format = "ensembl_gene_id")

gosselin_annotated <- eplp$get_regions_from_genes(raw_df_gosselin[['gene_name']],
                                                  gene_format = "hgnc_symbol")

combined_annotated <- rbind(galatro_annotated, gosselin_annotated) %>%
combined_annotated <- combined_annotated[!duplicated(combined_annotated), ]

# saved all three as grch38 separately
# However, only Gosselin et al. 2017 will be used for processing in pipelines
rdr$write_csv(galatro_annotated,
              here$here("output", "pathways",
                        "galatro_et_al_2017_microglia_grch38.csv"))
rdr$write_csv(gosselin_annotated,
              here$here("output", "pathways",
                        "gosselin_et_al_2017_microglia_grch38.csv"))
rdr$write_csv(combined_annotated,
              here$here("output", "pathways",
                        "gal_goss_merged_microglia_grch38.csv"))

# saving grch38 without APOE
rdr$write_csv(eplp$exclude_apoe_region(galatro_annotated),
              here$here("output", "pathways",
                        "galatro_et_al_2017_microglia_grch38_no_apoe.csv"))
rdr$write_csv(eplp$exclude_apoe_region(gosselin_annotated),
              here$here("output", "pathways",
                        "gosselin_et_al_2017_microglia_grch38_no_apoe.csv"))
rdr$write_csv(eplp$exclude_apoe_region(combined_annotated),
              here$here("output", "pathways",
                        "gal_goss_merged_microglia_grch38_no_apoe.csv"))

