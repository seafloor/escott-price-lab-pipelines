library(janitor)
library(dplyr)
library(readxl)
library(here)
source(here("R", "utilities", "utils.R"))
source(here("R", "utilities", "genes_to_coordinates.R"))

# load data from Galatro et al. 2017 and Gosselin et al. 2017
in_file_galatro_2017 <- here("output", "raw_data",
                             "galatro_et_al_2017_microglia_supp_table_2.csv")
in_file_gosselin_2017 <- here("output", "raw_data",
                              "gosselin_et_al_2017_microglia_supp_table_2.csv")

# max of 1296 required for Galatro et al. 2017 to match 
# the highlighted rows in the original file for microglial core genes
raw_df_galatro <- read_csv(in_file_galatro_2017, n_max = 1296)
raw_df_gosselin <- read_csv(in_file_gosselin_2017)

# check overlap between the two datasets
matches <- c()
for(i in 1:nrow(raw_df_gosselin)) {
  matches <- c(matches, sum(str_detect(raw_df_galatro[['external_gene_name']],
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
galatro_annotated <- get_regions_from_genes(raw_df_galatro[['ensembl_gene_id']],
                                            gene_format = "ensembl_gene_id")

gosselin_annotated <- get_regions_from_genes(raw_df_gosselin[['gene_name']],
                                             gene_format = "hgnc_symbol")

combined_annotated <- rbind(galatro_annotated, gosselin_annotated) |>
combined_annotated <- combined_annotated[!duplicated(combined_annotated), ]

# saved all three as grch38 separately
# However, only Gosselin et al. 2017 will be used for processing in pipelines
write_csv(galatro_annotated,
          here("output", "pathways",
               "galatro_et_al_2017_microglia_grch38.csv"))
write_csv(gosselin_annotated,
          here("output", "pathways",
               "gosselin_et_al_2017_microglia_grch38.csv"))
write_csv(combined_annotated,
          here("output", "pathways",
               "gal_goss_merged_microglia_grch38.csv"))

# saving grch38 without APOE
write_csv(exclude_apoe_region(galatro_annotated),
          here("output", "pathways",
               "galatro_et_al_2017_microglia_grch38_no_apoe.csv"))
write_csv(exclude_apoe_region(gosselin_annotated),
          here("output", "pathways",
               "gosselin_et_al_2017_microglia_grch38_no_apoe.csv"))
write_csv(exclude_apoe_region(combined_annotated),
          here("output", "pathways",
               "gal_goss_merged_microglia_grch38_no_apoe.csv"))

