# load packages with box
box::use(eplp = escottpricelabpipelines[rsids_to_functional_genes],
         rdr = readr[read_csv, write_csv],
         magrittr[`%>%`])

bellenguez_snps <- system.file("extdata", "annotations",
                               "all_bellenguez_hits_common.txt",
                               package = "escottpricelabpipelines")
rsids <- rdr$read_csv(bellenguez_snps,
                      col_names = c("rsid"))

output <- eplp$rsids_to_functional_genes(rsids[['rsid']][1:3])

annotations <- output[[1]]
missing_rsids <- output[[2]]

print(paste("--> Missing", length(missing_rsids),
            "rsids in functional annotations"))
print(missing_rsids)

print("--> Functional annotations for the following rsids:")
print(annotations)