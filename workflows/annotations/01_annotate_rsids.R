# load packages with box
box::use(eplp = escottpricelabpipelines[rsids_to_functional_genes],
         rdr = readr[read_csv, write_csv])

rsids <- rdr$read_csv("~/Desktop/All_bellenguez_hits_80SNPs.txt",
                      col_names = c("rsid"))

annotations <- eplp$rsids_to_functional_genes(rsids[['rsid']])

rdr$write_csv(annotations, "~/Desktop/all_bellenguez_hits_annotations.csv")

rsids <- rdr$read_csv("~/Desktop/29_GWS_SNP_not_input_data.txt",
                      col_names = c("rsid"))

annotations <- eplp$rsids_to_functional_genes(rsids[['rsid']])

rdr$write_csv(annotations, "~/Desktop/remaining_bellenguez_hits_annotations.csv")