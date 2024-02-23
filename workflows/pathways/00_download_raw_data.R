# initiated to avoid spamming journals with requests
# supplementary data is downloaded from the journal website and
# saved to the raw_data folder in the output directory

# load packages with box
box::use(eplp = escottpricelabpipelines[download_supplement],
         rdxl = readxl[read_excel],
         rdr = readr[write_csv],
         jtr = janitor[make_clean_names],
         here = here[here])

# set global output dir for raw data
data_dir <- here$here("output", "raw_data")

################# download astrocyte supplement #################
# DOI:10.1126/science.adc90
in_file <- "https://www.science.org/doi/suppl/10.1126/science.adc9020/suppl_file/science.adc9020_tables_s1_to_s4.zip"
out_file <- "endo_et_al_2022_astrocytes_supp_table_2.csv"

print(paste("--> Checking for file:", out_file))
if (!file.exists(here$here(data_dir, out_file))) {
  temp <- eplp$download_supplement(in_file)

  raw_df <- rdxl$read_excel(here$here(data_dir, "science.adc9020_table_s2.xlsx"),
                            sheet = "825 shared enriched genes",
                            range = "B3:Q828",
                            .name_repair = jtr$make_clean_names)

  unlink(temp, recursive = TRUE)
  rdr$write_csv(raw_df, here$here(data_dir, out_file))
} else {
  print("File already exists")
}

################# download microglial supplement #################
# download data from Galatro et al. 2017
# DOI https://doi.org/10.1038/nn.4597
in_file <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnn.4597/MediaObjects/41593_2017_BFnn4597_MOESM4_ESM.xlsx"
out_file <- "galatro_et_al_2017_microglia_supp_table_2.csv"

print(paste("--> Checking for file:", out_file))
if (!file.exists(here$here(data_dir, out_file))) {
  temp <- eplp$download_supplement(in_file_galatro_2017)
  raw_df <- rdxl$read_excel(temp,
                            sheet = "RPKM microglia_cortex_biopsy",
                            skip = 1,
                            .name_repair = jtr$make_clean_names)

  unlink(temp, recursive = TRUE)
  rdr$write_csv(raw_df, here$here(data_dir, out_file))
} else {
  print("File already exists")
}

# download data from Gosselin et al. 2017
# DOI: 10.1126/science.aal3222
in_file <- "https://www.science.org/doi/suppl/10.1126/science.aal3222/suppl_file/aal3222_gosselin_tables2.xlsx"
out_file <- "gosselin_et_al_2017_microglia_supp_table_2.csv"

print(paste("--> Checking for file:", out_file))
if (!file.exists(here$here(data_dir, out_file))) {
  temp <- eplp$download_supplement(in_file_gosselin_2017)
  raw_df <- rdxl$read_excel(temp,
                            sheet = "Human microglia gene signature",
                            .name_repair = jtr$make_clean_names)

  unlink(temp, recursive = TRUE)
  rdr$write_csv(raw_df, here$here(data_dir, out_file))
} else {
  print("File already exists")
}

################# download synaptic supplement #################
# download from syngo portal
# DOI: https://doi.org/10.1016/j.neuron.2019.05.002
in_file <- "https://syngoportal.org/data/SynGO_bulk_download_release_20231201.zip"
out_file <- "syngo_database_01_12_2023_release_synaptic_genes.csv"

print(paste("--> Checking for file:", out_file))
if (!file.exists(here$here(data_dir, out_file))) {
  temp <- eplp$download_supplement(in_file)
  raw_df <- rdxl$read_excel("syngo_genes.xlsx",
                            sheet = "Sheet 1")

  unlink(temp, recursive = TRUE)
  rdr$write_csv(raw_df, here$here(data_dir, out_file))
} else {
  print("File already exists")
}
