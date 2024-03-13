# imports
box::use(eplp = escottpricelabpipelines[simulate_y],
         jntr = janitor[clean_names, row_to_names],
         dplyr = dplyr[filter, pull, mutate, select],
         stringr = stringr[str_c, str_count],
         here = here[here],
         rdr = readr[read_tsv],
         tibble = tibble[as_tibble])

# set path for your bcftools binary
# will be different depending on where bcftools is installed
# e.g.
Sys.setenv(PATH = paste(Sys.getenv("PATH"), here$here("../bcftools/"), sep = ":"))

# download chr 22
ids_1kg <- here$here("inst", "extdata", "annotations", "igsr_1kg_sample_list_grch38.tsv")
genotypes_1kg <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz"
eur_ids <- rdr$read_tsv(ids_1kg) %>%
  jntr$clean_names() %>%
  dplyr$filter(supserpopulation_code == "EUR") %>%
  dplyr$pull(sample_name) %>%
  paste(collapse = ",")

# call BCFtools, filter to EUR, and MAF > 0.01
system2(command = 'bcftools',
        args = c('view',
                 genotypes_1kg,
                 '--force-samples',
                 '--samples',
                 eur_ids,
                 '--min-af 0.01',
                 '--output-type',
                 'v',
                 '--output',
                 here$here("data", "chr22_1kg_eur_grch38.vcf")))

# potentially delete files afterwards


# read our new vcf file, which is variants * individuals
vcf <- rdr$read_tsv(here$here("data", "chr22_1kg_eur_grch38.vcf"), comment = "##") %>%
                       dplyr$mutate(ID = stringr$str_c(`#CHROM`, "_", POS, "_", REF, "_", ALT)) %>%
                       dplyr$select(-`#CHROM`, -POS, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT)

# swap to individuals * variants and count number of ALT alleles
df <- vcf %>%
  t() %>%
  jntr$row_to_names(row_number = 1) %>%
  apply(2, stringr$str_count, "1") %>%
  tibble$as_tibble()

# drop to 100k SNPs
n_snp <- 100000
df <- df[, as.logical(sample(c(rep(1, n_snp), rep(0, (ncol(df) - n_snp))), ncol(df), replace = FALSE))]

# sample rows with replace = TRUE to get 300k individuals
n_obs <- 500
df_upsampled <- df %>%
  sapply(sample, size=n_obs, replace = TRUE)

y <- eplp$simulate_y(df_upsampled, p_causal = 0.1, h2_l = 0.2, k = 0.05)

# to do:
# test on server
# add apoe simulation
