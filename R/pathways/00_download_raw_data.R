# initiated to avoid spamming journals with requests
# supplementary data is downloaded from the journal website and 
# saved to the raw_data folder in the output directory
library(readxl)
library(readr)
library(here)
library(janitor)
source(here("R", "utilities", "utils.R"))

################# download astrocyte supplement #################
# DOI:10.1126/science.adc90
in_file <- "https://www.science.org/doi/suppl/10.1126/science.adc9020/suppl_file/science.adc9020_tables_s1_to_s4.zip"
out_file <- "endo_et_al_2022_astrocytes_supp_table_2.csv"

if (!file.exists(here("output", "raw_data", out_file))) {
  temp <- download_supplement(in_file)
  
  raw_df <- read_excel(here(data_dir, "science.adc9020_table_s2.xlsx"),
                       sheet = "825 shared enriched genes",
                       range = "B3:Q828"
                       .name_repair = make_clean_names)
  
  unlink(temp, recursive = TRUE)
  write_csv(raw_df, here("output", "raw_data", out_file))
}

################# download microglial supplement #################
# download data from Galatro et al. 2017
# DOI https://doi.org/10.1038/nn.4597
in_file <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnn.4597/MediaObjects/41593_2017_BFnn4597_MOESM4_ESM.xlsx"
out_file <- "galatro_et_al_2017_microglia_supp_table_2.csv"

if (!file.exists(here("output", "raw_data", out_file))) {
  temp <- download_supplement(in_file_galatro_2017)
  raw_df <- read_excel(temp,
                       sheet = "RPKM microglia_cortex_biopsy",
                       skip = 1,
                       .name_repair = make_clean_names)
  
  unlink(temp, recursive = TRUE)
  write_csv(raw_df, here("output", "raw_data", out_file))
}

# download data from Gosselin et al. 2017
# DOI: 10.1126/science.aal3222
in_file <- "https://www.science.org/doi/suppl/10.1126/science.aal3222/suppl_file/aal3222_gosselin_tables2.xlsx"
out_file <- "gosselin_et_al_2017_microglia_supp_table_2.csv"

if (!file.exists(here("output", "raw_data", out_file))) {
  temp <- download_supplement(in_file_gosselin_2017)
  raw_df <- read_excel(temp,
                       sheet = "Human microglia gene signature",
                       .name_repair = make_clean_names)
  
  unlink(temp, recursive = TRUE)
  write_csv(raw_df, here("output", "raw_data", out_file)))
}

################# download synaptic supplement #################
# download from syngo portal
# DOI: https://doi.org/10.1016/j.neuron.2019.05.002
in_file <- "https://syngoportal.org/data/SynGO_bulk_download_release_20231201.zip"
out_file <- "syngo_database_01_12_2023_release_synaptic_genes.csv"

if (!file.exists(here("output", "raw_data", out_file))) {
  temp <- download_supplement(in_file)
  raw_df <- read_excel("syngo_genes.xlsx",
                       sheet = "Sheet 1")
  
  unlink(temp, recursive = TRUE)
  write_csv(raw_df, here("output", "raw_data", out_file))
}
