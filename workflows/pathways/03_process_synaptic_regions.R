# genes taken from SynGo release on 01/12/2023
# see https://www.syngoportal.org/ and DOI:https://doi.org/10.1016/j.neuron.2019.05.002

# load data from syngo
in_file <- here("output", "raw_data",
                "syngo_database_01_12_2023_release_synaptic_genes.csv")

raw_df <- read_csv(in_file)

# get regions for genes
df_annotated <- get_regions_from_genes(raw_df[['ensembl_id']],
                                       gene_format = "ensembl_gene_id")

# save grch38
write_csv(df_annotated,
          here("output", "pathways",
               "syngo_database_01_12_2023_release_synaptic_genes_grch38.csv"))

# save grch38 without APOE
write_csv(exclude_apoe_region(df_annotated),
          here("output", "pathways",
               "syngo_database_01_12_2023_release_synaptic_genes_grch38_no_apoe.csv"))
