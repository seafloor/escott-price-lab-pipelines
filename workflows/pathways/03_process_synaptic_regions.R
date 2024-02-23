# genes taken from SynGo release on 01/12/2023
# see https://www.syngoportal.org/ and DOI:https://doi.org/10.1016/j.neuron.2019.05.002
# load packages with box
box::use(eplp = escottpricelabpipelines[get_regions_from_genes, exclude_apoe_region],
         rdr = readr[read_csv, write_csv],
         here = here[here])

# load data from syngo
in_file <- here$here("output", "raw_data",
                     "syngo_database_01_12_2023_release_synaptic_genes.csv")

raw_df <- rdr$read_csv(in_file)

# get regions for genes
df_annotated <- eplp$get_regions_from_genes(raw_df[['ensembl_id']],
                                            gene_format = "ensembl_gene_id")

# save grch38
rdr$write_csv(df_annotated,
              here$here("output", "pathways",
                        "syngo_01_12_2023_synaptic_genes_grch38.csv"))

# save grch38 without APOE
rdr$write_csv(eplp$exclude_apoe_region(df_annotated),
              here$here("output", "pathways",
                        "syngo_01_12_2023_synaptic_genes_grch38_no_apoe.csv"))
