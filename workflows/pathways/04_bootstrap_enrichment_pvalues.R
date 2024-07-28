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
         dplyr = dplyr[filter, select, rename_all],
         stringr = stringr[str_to_upper],
         magrittr[`%>%`],
         here = here[here],
         prll = parallel[makeCluster, detectCores, stopCluster,
                         clusterEvalQ, clusterExport, parApply],
         eplp = escottpricelabpipelines[
           read_ref_genome_coordinates, force_canonical_autosomes,
           read_bim_file, read_regions_to_search,
           calculate_proportion, sample_regions_from_list,
           check_if_snps_in_region])

################# functions for region matching #################s
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

# read in top SNP rsids as a vector
# here using a random selection of 100 snps from the bim file
top_snp_rsids <- sample(all_snps[['id']], 100, replace=FALSE)
top_snp_pos <- all_snps %>%
  dplyr$filter(id %in% top_snp_rsids) %>%
  dplyr$select(chr, pos) %>%
  dplyr$rename_all(stringr$str_to_upper)

# export data to cluster
prll$clusterExport(cl, 'ref_coords')
prll$clusterExport(cl, 'all_snps')
prll$clusterExport(cl, 'regions_to_search')
prll$clusterExport(cl, 'top_snp_pos')

# get actual proportion of top SNPs in regions
n_in_true_regions <- prll$parApply(cl = cl, X = regions_to_search,
                                   MARGIN = 1, FUN = eplp$check_if_snps_in_region,
                                   snps_to_check = top_snp_pos)
actual_proportion <- eplp$calculate_proportion(n_in_true_regions,
                                               abs(regions_to_search[['gene_length_bp']]))

# empty mat for storing bootstrap results
null_proportions <- rep(NA, times=n_bootstrap)

# get null distribution of top SNPs in regions
for(k in 1:n_bootstrap) {

  # randomly sample regions of equal length and export to cluster
  null_regions_to_search <- eplp$sample_regions_from_list(abs(regions_to_search[['gene_length_bp']]))
  prll$clusterExport(cl, 'null_regions_to_search')

  # get number of top SNPs present in random region
  n_in_null_regions <- prll$parApply(cl = cl, X = null_regions_to_search,
                                     MARGIN = 1, FUN = check_if_snps_in_region,
                                     snps_to_check = top_snp_pos)

  # save proportion of SNPs in regions
  null_proportions[k] <- eplp$calculate_proportion(n_in_null_regions, abs(regions_to_search[['gene_length_bp']]))
}

# get p from null distribution
# expect null given randomly sampled example data
bootstrap_p <- mean(null_proportions >= actual_proportion)
prll$stopCluster(cl)
