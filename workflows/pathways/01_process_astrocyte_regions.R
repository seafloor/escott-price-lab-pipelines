library(biomaRt)
library(readr)
library(dplyr)
library(tibble)
library(here)
source(here("R", "utilities", "swap_organisms.R"))
source(here("R", "utilities", "genes_to_coordinates.R"))
source(here("R", "utilities", "utils.R"))

# load data from endo et al. 2022
in_file <- "endo_et_al_2022_astrocytes_supp_table_2.csv"
raw_df <- read_csv(here("output", "raw_data", in_file))

# convert from mouse to human
human_df <- mouse_to_human(raw_df[['gene']])

# printing duplicated ensemble IDs
print("Duplicated ensemble IDs:")
print(human_df[duplicated(human_df[['gene_ensemble']]),])

# dropping duplicates
human_df <- human_df |>
  select(-mgi_symbol) |>
  distinct() 

# saving grch38
write_csv(human_df,
          here("output", "pathways",
               "endo_et_al_2022_astrocytes_grch38.csv")

# saving grch38 without APOE
write_csv(exclude_apoe_region(human_df),
          here("output", "pathways",
               "endo_et_al_2022_astrocytes_grch38_no_apoe.csv")
