##### imports #####
devtools::document()
devtools::load_all()
box::use(
  eplp = escottpricelabpipelines[
    simulate_y, read_config,
    binary_from_liability,
    upsample_individuals_independently,
    simulate_e2_e4,
    downsample_individuals
  ],
  jntr = janitor[clean_names, row_to_names],
  dplyr = dplyr[filter, pull, mutate, select],
  stringr = stringr[str_c, str_count],
  here = here[here],
  rdr = readr[read_tsv],
  tibble = tibble[as_tibble],
  bgm = bigmemory[big.matrix],
  prgrss = progress[progress_bar]
)

##### parse 1kg data #####
# download 1kg chr 22
ids_1kg <- here$here("inst", "extdata", "annotations", "igsr_1kg_sample_list_grch38.tsv")
genotypes_1kg <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr22.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz"
eur_ids <- rdr$read_tsv(ids_1kg) %>%
  jntr$clean_names() %>%
  dplyr$filter(superpopulation_code == "EUR") %>%
  dplyr$pull(sample_name) %>%
  paste(collapse = ",")

# call BCFtools, filter to EUR, and MAF > 0.01
system2(
  command = "bcftools",
  args = c(
    "view",
    genotypes_1kg,
    "--force-samples",
    "--samples",
    eur_ids,
    "--min-af 0.01",
    "--output-type",
    "v",
    "--output",
    here$here("data", "chr22_1kg_eur_grch38.vcf")
  )
)

# read our new vcf file, which is variants * individuals
vcf <- rdr$read_tsv(here$here("data", "chr22_1kg_eur_grch38.vcf"), comment = "##") %>%
  dplyr$mutate(ID = stringr$str_c(`#CHROM`, "_", POS, "_", REF, "_", ALT)) %>%
  dplyr$select(-`#CHROM`, -POS, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT)

# swap to individuals * variants and count number of ALT alleles
# requires swapping all row/column counts downstream
# df <- vcf %>%
#   t() %>%
#   jntr$row_to_names(row_number = 1) %>%
#   apply(2, stringr$str_count, "1")

df <- vcf %>%
  dplyr$select(-ID) %>%
  apply(2, stringr$str_count, "1")

snp_ids <- vcf[["ID"]]
rm(vcf)

print(paste("Read in", nrow(df), "SNPs and", ncol(df), "individuals"))

##### up/down sample to desired dimensions for simulation #####
# drop to 100k SNPs, plus 2 for the apoe alleles
n_snp <- 10000 + 2
snps_subsample_bool <- as.logical(sample(
  c(
    rep(1, n_snp),
    rep(0, (nrow(df) - n_snp))
  ),
  nrow(df),
  replace = FALSE
))
df <- df[snps_subsample_bool, ]
snp_ids_subsampled <- snp_ids[snps_subsample_bool]

# note 300k required due to low prevalence and because we first simulate
# the overall populations from which we later draw a case-control sample
ncols <- 5000 # n individuals as matrix transposed
nrows <- 10000 + 2 # n SNPs as matrix transposed
k <- 0.15
p_causal <- 0.1
h2_l <- 0.1
n_cases <- 5000 # used for down-sampling to a case-control study later
n_controls <- 5000

m <- eplp$upsample_individuals_independently(df, nrows, ncols,
  is_transposed = TRUE
)
rm(df)

if (nrows < n_snp) {
  snp_ids_subsampled <- snp_ids_subsampled[1:nrows]
}

##### simulate additive polygenic background #####
out <- eplp$simulate_y(m,
  p_causal = p_causal,
  h2_l = h2_l, k = k,
  is_transposed = TRUE
)

# unpack results
X <- out[[1]]
y <- out[[2]]
polygenic_liability <- out[[3]]
causal_idx <- out[[4]]
simulated_betas <- out[[5]]

##### add APOE effects #####
# read in APOE effect sizes and allele frequencies from config.toml
# note this uses the kunkle et al. (2019) effect sizes and UKBB (2021 release) allele frequencies
# both are Europeans only and specifically the APOE e2/e4 alleles, not the leads APOE SNPs
apoe_sim <- eplp$simulate_e2_e4(ncols)
e2_e4_score_scaled <- scale(apoe_sim[["e2_e4_score"]])

# add additional noise to scale APOE effect sizes
# assume h2 for apoe alone is 0.15
var_e <- (var(e2_e4_score_scaled) / 0.15) - var(e2_e4_score_scaled) # variance in e (noise) based on variance in g (genotypic risk) and heritability on the liability scale
e <- rnorm(ncols, 0, sqrt(var_e)) # simulate error for phenotype for all individuals

# rescale polygenic risk
sim_l_var <- as.vector((var(e2_e4_score_scaled) / (1 - 0.15)) - var(e2_e4_score_scaled))
polygenic_liability_scaled <- scale(polygenic_liability) * sqrt(sim_l_var)
l <- e2_e4_score_scaled + polygenic_liability_scaled + e

hist(e2_e4_score_scaled)
hist(polygenic_liability_scaled)
hist(l)

# create binary phenotype - overwrite y
y <- eplp$binary_from_liability(l, k = k)

# overwrite two random non-causal SNPs with the APOE alleles
apoe_idx <- sample(which(!causal_idx), size = 2, replace = FALSE)
X[apoe_idx, ] <- t(dplyr::select(apoe_sim, e2, e4))

# downsample to make case-control study
case_control_sim <- eplp$downsample_individuals(X, y, n_cases, n_controls,
  is_transposed = TRUE
)

X <- case_control_sim[[1]]
y <- case_control_sim[[2]]

##### optional: check AUC with PRS #####
pb <- prgrss$progress_bar$new(
  format = "  Running GWAS [:bar] :percent eta: :eta",
  total = 100,
  force = TRUE
)

# train/test split to avoid overfitting
set.seed(123)
all_idx <- 1:ncol(X)
train_idx <- sample(all_idx, as.integer(ncol(X) * 0.8))
test_idx <- all_idx[-train_idx]

# run a GWAS
betas <- c()
for (i in 1:nrow(X)) {
  if (i %% (as.integer(nrow(X) / 100)) == 0) {
    pb$tick()
  }
  mod <- glm(y[train_idx] ~ X[i, train_idx], family = binomial)
  betas <- c(betas, coef(mod)[2])
}

# make a PRS with and without APOE on whole data using betas from train split
# note all calculations e.g. dot product are transposed because we're working
# on a genotypes*individuals matrix
betas[is.na(betas)] <- 0
b_mat <- matrix(betas, nrow = 1)
df_prs <- data.frame(prs = scale(as.vector((b_mat %*% X) / nrow(X))), y = y)
df_prs_no_apoe <- data.frame(prs = scale(as.vector((b_mat[, -apoe_idx] %*% X[-apoe_idx, ]) / (nrow(X) - length(apoe_idx)))), y = y)
df_prs_apoe_only <- data.frame(prs = scale(as.vector((b_mat[, apoe_idx] %*% X[apoe_idx, ]) / length(apoe_idx))), y = y)

hist(df_prs$prs)
hist(df_prs_no_apoe$prs)
hist(df_prs_apoe_only$prs)

# check APOE effect sizes in whole datasets
config <- eplp$read_config()
apoe_effects <- unlist(config$apoe_effect_sizes_europeans$kunkle_2019)

e2_mod <- glm(y ~ X[apoe_idx[1], ], family = binomial)
paste(
  "e2 beta expected: ", apoe_effects["e2"],
  ", e2 beta actual: ", coef(e2_mod)[2]
)

e4_mod <- glm(y ~ X[apoe_idx[2], ], family = binomial)
paste(
  "e4 beta expected: ", apoe_effects["e4"],
  ", e4 beta actual: ", coef(e4_mod)[2]
)

# check fit and AUC with APOE only
mod_apoe_only <- glm(y ~ prs,
  data = df_prs_apoe_only[train_idx, ],
  family = binomial
)
summary(mod_apoe_only)
paste(
  "Test AUC for APOE only: ",
  pROC::auc(pROC::roc(
    df_prs_apoe_only[test_idx, "y"],
    predict(mod_apoe_only,
      newdata = df_prs_apoe_only[test_idx, ],
      type = "response"
    )
  ))
)

# check fit and AUC for PRS without APOE
mod_no_apoe <- glm(y ~ prs,
  data = df_prs_no_apoe[train_idx, ],
  family = binomial
)
summary(mod_no_apoe)
paste(
  "Test AUC for PRS without APOE: ",
  pROC::auc(pROC::roc(
    df_prs_no_apoe[test_idx, "y"],
    predict(mod_no_apoe,
      newdata = df_prs_no_apoe[test_idx, ],
      type = "response"
    )
  ))
)

# check fit and AUC with APOE
mod_apoe <- glm(y ~ prs,
  data = df_prs[train_idx, ],
  family = binomial
)
summary(mod_apoe)
paste(
  "Test AUC for PRS with APOE included: ",
  pROC::auc(pROC::roc(
    df_prs[test_idx, "y"],
    predict(mod_apoe,
      newdata = df_prs[test_idx, ],
      type = "response"
    )
  ))
)

##### optional write to plink file #####

# make fake fam files - all fake ids as individuals oversampled
# sex random and unassociated with phenotype
fam <- tibble::tibble(
  fam = stringr$str_c("fid", 1:ncol(X)),
  id = stringr$str_c("iid", 1:ncol(X)),
  pat = 0,
  mat = 0,
  sex = sample(c(1, 2),
    size = ncol(X),
    replace = TRUE
  ),
  pheno = y + 1
)

# make fake bim file from real 1kg marker IDs
bim <- t(data.frame(stringr::str_split(snp_ids_subsampled, "_")))
rownames(bim) <- NULL
colnames(bim) <- c("chr", "pos", "ref", "alt")

bim <- bim %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    chr = stringr::str_replace(chr, "chr", ""),
    pos = as.numeric(pos),
    id = snp_ids_subsampled,
    posg = 0
  ) %>%
  dplyr::select(chr, id, posg, pos, alt, ref)

# correct APOE changes
# positions set to grch38 for rs7412/rs429358 for e2/e4 respectively
bim[apoe_idx, 1] <- c("19", "19")
bim[apoe_idx, 2] <- c("e2", "e4")
bim[apoe_idx, 4] <- c(44908822, 44908684)
bim[apoe_idx, 5] <- c("T", "C")
bim[apoe_idx, 6] <- c("C", "T")

# check APOE SNPs are in the bim file!
dplyr::filter(bim, chr == "19")

# write files
genio::write_plink(here$here("data", "simulated_data_tmp"),
  X = X,
  bim = bim,
  fam = fam
)

# fix split chromosomes
system2(
  command = "plink",
  args = c(
    "--bfile",
    here$here("data", "simulated_data_tmp"),
    "--make-bed",
    "--out",
    here$here("data", "simulated_data")
  )
)

# drop temp files
unlink(here$here("data", "simulated_data_tmp.*"))

##### cleanup files from big.matrix #####
# should drop backing files according to docs for bigmemory
rm(out, X, m, y)
gc()

# additionally remove file backing in case not handled by rm/gc
unlink(here$here("data", "big_matrix.bin"), recursive = TRUE)
unlink(here$here("data", "big_matrix.desc"), recursive = TRUE)
