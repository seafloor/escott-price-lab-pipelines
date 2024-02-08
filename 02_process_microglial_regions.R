library(janitor)
library(readxl)
library(here)
source(here("R", "utilities", "utils.R"))

in_file_galatro_2017 <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnn.4597/MediaObjects/41593_2017_BFnn4597_MOESM4_ESM.xlsx"
in_file_gosselin_2017 <- "https://www.science.org/doi/suppl/10.1126/science.aal3222/suppl_file/aal3222_gosselin_tables2.xlsx"

temp <- download_supplement(in_file_galatro_2017)
raw_df_galatro <- read_excel(temp,
                             sheet = "RPKM microglia_cortex_biopsy",
                             skip = 1,
                             .name_repair = make_clean_names)
unlink(temp)

temp <- download_supplement(in_file_gosselin_2017)
raw_df_gosselin <- read_excel(temp,
                              sheet = "Human microglia gene signature",
                              .name_repair = make_clean_names)
unlink(temp)

matches <- c()
for(i in 1:nrow(raw_df_gosselin)) {
  matches <- c(matches, sum(str_detect(raw_df_galatro[['external_gene_name']], paste("^", raw_df_gosselin[['gene_name']][i], "$", sep = ""))))
}
