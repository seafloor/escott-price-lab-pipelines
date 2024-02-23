# load packages with box
box::use(eplp = escottpricelabpipelines[mouse_to_human, exclude_apoe_region],
         rdr = readr[read_csv, write_csv],
         dplyr = dplyr[select, distinct],
         here = here[here])

# load data from endo et al. 2022
in_file <- "endo_et_al_2022_astrocytes_supp_table_2.csv"
raw_df <- rdr$read_csv(here$here("output", "raw_data", in_file))

# convert from mouse to human
human_df <- eplp$mouse_to_human(raw_df[['gene']])

# printing duplicated ensemble IDs
print("Duplicated ensemble IDs:")
print(human_df[duplicated(human_df[['gene_ensemble']]),])

# dropping duplicates
human_df <- human_df %>%
  dplyr$select(-mgi_symbol) %>%
  dplyr$distinct()

# saving grch38
rdr$write_csv(human_df,
              here$here("output", "pathways",
                        "endo_et_al_2022_astrocytes_grch38.csv"))

# saving grch38 without APOE
rdr$write_csv(eplp$exclude_apoe_region(human_df),
              here$here("output", "pathways",
                        "endo_et_al_2022_astrocytes_grch38_no_apoe.csv"))
