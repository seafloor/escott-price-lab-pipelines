# a minimal and non-optimised implementation of boostrapping enrichment p-values
# to be run in scenarios where p-values from summary statistics are not available
# if p-values are available, other methods such as MAGMA can be used

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
#######################################################

################# load libraries #################
# load packages with box
box::use(rdr = readr[write_csv, cols, col_character],
         dplyr = dplyr[filter],
         magrittr[`%>%`],
         here = here[here],
         prll = parallel[makeCluster, detectCores, stopCluster,
                         clusterEvalQ, clusterExport, parApply],
         eplp = escottpricelabpipelines[
           read_ref_genome_coordinates, force_canonical_autosomes,
           read_bim_file, read_regions_to_search,
           calculate_proportion, sample_regions_from_list])

################# functions for region matching #################
# read in rsid list of top SNPs as vector
read_top_snps <- function(f, all_snps){
  top_snp_rsids <- rdr$read_csv(f,
                                col_names = FALSE,
                                col_types = rdr$cols(.default = rdr$col_character()))[['X1']]
  positions <- dplyr$filter(all_snps, ID %in% top_snp_rsids)

  return(positions)
}

#### run bootstrapping
# set simulation parameters
# example number - increase for more stable p-value
n_bootstrap <- 100
pathway_file <- here$here("output", "pathways",
                          "syngo_01_12_2023_synaptic_genes_grch38.csv")

# start parallel cluster and load essential libraries on it
cl <- prll$makeCluster(prll$detectCores())
prll$clusterEvalQ(cl, { library(dplyr) })

# read the reference genome coordinates (supplied here in a convenience function)
# these are pre-processed from the grch38p14 sequence report 
# available here: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/
ref_coords <- eplp$read_ref_genome_coordinates()

# read in a bim file of all SNPs in the dataset under consideration
# should be everything that passed QC and went into the analysis
# here we use a non-standard file name (no .bim extension) purely because
# otherwise it is automatically tracked by git LFS
bf <- here$here("inst", "extdata", "annotations", "minimal_bim_file.txt")
all_snps <- eplp$force_canonical_autosomes(eplp$read_bim_file(bf))

# read in the regions to search for top SNPs (genes to check for enrichment in)
# set f="dummy_region" to use a dummy region file for testing
regions_to_search <- eplp$read_regions_to_search(f=pathway_file)
top_snp_pos <- eplp$read_top_snps()

# export data to cluster
prll$clusterExport(cl, 'ref_coords')
prll$clusterExport(cl, 'all_snps')
prll$clusterExport(cl, 'regions_to_search')
prll$clusterExport(cl, 'top_snp_pos')

# get actual proportion of top SNPs in regions
n_in_true_regions <- prll$parApply(cl = cl, X = regions_to_search, MARGIN = 1, FUN = check_if_snps_in_region)
actual_proportion <- eplp$calculate_proportion(n_in_true_regions, regions_to_search[['bp_len']])

# empty mat for storing bootstrap results
null_proportions <- rep(NA, times=n_bootstrap)

# get null distribution of top SNPs in regions
for(k in 1:n_bootstrap) {

  # randomly sample regions of equal length and export to cluster
  null_regions_to_search <- eplp$sample_regions_from_list(regions_to_search[['bp_len']])
  prll$clusterExport(cl, 'null_regions_to_search')

  # get number of top SNPs present in random region
  n_in_null_regions <- prll$parApply(cl = cl, X = null_regions_to_search,
                                     MARGIN = 1, FUN = check_if_snps_in_region)

  # save proportion of SNPs in regions
  null_proportions[k] <- eplp$calculate_proportion(n_in_null_regions, regions_to_search[['bp_len']])
}

# get p from null distribution
bootstrap_p <- mean(null_proportions >= actual_proportion)
prll$stopCluster(cl)
