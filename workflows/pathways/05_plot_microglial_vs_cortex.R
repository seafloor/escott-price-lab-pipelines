library(here)
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(forcats)

# load Galatro et al. 2017 data
raw_df_galatro <- read_csv(here("output", "raw_data",
                                "galatro_et_al_2017_microglia_supp_table_2.csv"))

# load gene list from study
xgb_genes <- read_excel("~/Downloads/eadb_ml_supplementary_tables_in_progress_12_02_24_MBS.xlsx",
                  sheet = "Table S1", range = "A1:M109") %>%
  pull('Consensus Locus') %>%
  stringr::str_split(';') %>%
  unlist() %>%
  unique()

dplyr::select(raw_df_galatro, ensembl_gene_id, external_gene_name,
              glia_brain_log_fc, epi_brain_brain_log_fc) %>%
  filter(external_gene_name %in% xgb_genes) %>%
  tidyr::pivot_longer(-c(ensembl_gene_id, external_gene_name),
                      names_to = "region",
                      values_to = "log_fc") %>%
  ggplot(aes(x = fct_reorder(external_gene_name, log_fc), y = log_fc, color = region)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6)) +
  labs(title = "Microglial gene expression in Galatro et al. 2017",
       x = "Region",
       y = "Log fold change")