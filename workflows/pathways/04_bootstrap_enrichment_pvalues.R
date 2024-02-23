# a minimal and non-optimised implementation of boostrapping enrichment p-values
# to be run in scenarios where p-values from summary statistics are not available
# if p-values are available, alternative methods such as MAGMA are available

################# methods description #################
# A – set of all SNPs (N in total)
# B = (b1, b2, b3) – collection of m(=3) regions (genes) of different lengths (N1, N2, N3)
# C – set of top SNPs from ML: c1 – SNPs in B, c2 – SNPs not in B (C=c1+c2)
# AP = (ACTUAL PROPORTION of top SNPs in B) = c1/( N1, N2, N3)

# Simulations:
# From A extract randomly m regions of the same length as in B (i.e. N1, N2, N3) and
# calculate SP_k (SIMULATED PROPORTION at iteration k) as above where c1 is now
# the number of SNPs form C in these randomly selected regions.
# Final p-value is proportion of random simulations where SP_k>=AP.

# Notes:
# - The number of simulations should be large enough to get a stable p-value
#######################################################

################# load libraries #################
library(readr)
library(dplyr)
library(ggplot2)
library(parallel)
library(here)
source(here("R", "utilities", "utils.R"))
source(here("R", "utilities", "boostrap.R"))
# load packages with box
box::use()

################# functions for region matching #################
# read in rsid list of top SNPs as vector
read_top_snps <- function(f){
  top_snp_rsids <- read_csv(f,
                            col_names = FALSE,
                            col_types = cols(.default = col_character()))[['X1']]
  positions <- filter(all_snps, ID %in% top_snp_rsids)

  return(positions)
}

#### run bootstrapping
# set simulation parameters
n_bootstrap <- 1000
pathway_file <- here("output", "pathways",
                     "syngo_database_01_12_2023_release_synaptic_genes_grch38.csv")

# start parallel cluster and load essential libraries on it
cl <- parallel::makeCluster(detectCores())
clusterEvalQ(cl, { library(dplyr); })

# read data into cluster
ref_coords <- read_ref_genome_coordinates()
all_snps <- read_bim_file() %>%
  force_canonical_autosomes()

regions_to_search <- read_regions_to_search("dummy_region")
top_snp_pos <- read_top_snps()

# export data to cluster
clusterExport(cl, 'ref_coords')
clusterExport(cl, 'all_snps')
clusterExport(cl, 'regions_to_search')
clusterExport(cl, 'top_snp_pos')

# get actual proportion of top SNPs in regions
n_in_true_regions <- parApply(cl = cl, X = regions_to_search, MARGIN = 1, FUN = check_if_snps_in_region)
actual_proportion <- calculate_proportion(n_in_true_regions, regions_to_search[['bp_len']])

# empty mat for storing bootstrap results
null_proportions <- rep(NA, times=n_bootstrap)

# get null distribution of top SNPs in regions
for(k in 1:n_bootstrap) {

  # randomly sample regions of equal length and export to cluster
  null_regions_to_search <- sample_regions_from_list(regions_to_search[['bp_len']])
  clusterExport(cl, 'null_regions_to_search')

  # get number of top SNPs present in random region
  n_in_null_regions <- parApply(cl = cl, X = null_regions_to_search, MARGIN = 1, FUN = check_if_snps_in_region)

  # save proportion of SNPs in regions
  null_proportions[k] <- calculate_proportion(n_in_null_regions, regions_to_search[['bp_len']])
}

# get p from null distribution
bootstrap_p <- mean(null_proportions >= actual_proportion)
stopCluster(cl)
