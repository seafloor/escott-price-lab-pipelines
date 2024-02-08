library(biomaRt)
library(readr)
library(readxl)
library(dplyr)
library(tibble)
library(here)
source(here("R", "utilities", "swap_organisms.R"))
source(here("R", "utilities", "genes_to_coordinates.R"))
source(here("R", "utilities", "utils.R"))

# read raw data
in_file <- "https://www.science.org/doi/suppl/10.1126/science.adc9020/suppl_file/science.adc9020_tables_s1_to_s4.zip"

temp <- download_supplement(in_file)

raw_df <- read_excel("science.adc9020_table_s2.xlsx",
                     sheet = "825 shared enriched genes",
                     range = "B3:B828")

# cleanup temp files
unlink(temp)
file_list <- list_files(1:4, prepend = "science.adc9020_table_s",
                        append = ".xlsx")
cleanup_files(file_list)

# convert from mouse to human
human_df <- mouse_to_human(raw_df[['gene']])

# printing duplicated ensemble IDs
print(human_df[duplicated(human_df[['gene_ensemble']]),])

# dropping duplicates
human_df <- human_df |>
  select(-mgi_symbol) |>
  distinct() 

# saving grch38
write_csv(human_df, 'output/pathways/endo_et_al_2022_astrocytes_grch38.csv')

# saving grch38 without APOE
human_df_no_apoe <- exclude_apoe_region(human_df)
write_csv(human_df_no_apoe, '../output/endo_et_al_2022_astrocytes_grch38_no_apoe.csv')
