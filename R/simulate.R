#' Simulate genotypes for genetic association studies
#'
#' This script simulated unnasociated genotypes only. It will draw minor allele frequency for
#' SNPs uniformly from maf_range.
#'
#' @param n number of individuals
#' @param p number of SNPs
#' @param maf_range range of minor allele frequencies
#'
#' @return A n * p matrix of genotypes
#'
#' @export
add_noise <- function(n = 1000, p = 100, maf_range = c(0.05, 0.5)) {
  maf <- runif(p, maf_range[1], maf_range[2])
  genotypes <- sapply(maf, function(freq) rbinom(n, 2, freq))
  return(genotypes)
}

#' Create a new SNP in linkage disequiblibrium (LD) with a given SNP
#'
#' @param causal_snp A vector of genotypes for a single SNP
#' @param r2 The r-squared value for LD between the new SNP and the causal SNP
#'
#' @return A vector of genotypes for the new SNP
#'
#' @export
add_simple_ld <- function(causal_snp, r2) {
  indices_to_shuffle <- sample(1:length(causal_snp),
                               size = (1 - r2) * length(causal_snp),
                               replace = FALSE
  )
  ld_snp <- causal_snp

  if (length(indices_to_shuffle) > 1) {
    shuffledpeople <- sample(causal_snp[indices_to_shuffle])
    ld_snp[indices_to_shuffle] <- shuffledpeople
  }


  return(ld_snp)
}

#' Replace a proportion of causal SNPs with new SNPs in LD
#'
#' @param g A matrix of genotypes
#' @param m The proportion of causal SNPs to keep
#' @param ld_range The range of r-squared values for LD between the new SNPs and the causal SNPs
#'
#' @return A matrix of genotypes with a proportion of causal SNPs replaced by new SNPs in LD
#'
#' @export
replace_with_ld <- function(g, m, ld_range = c(0.1, 0.9)) {
  # g is our genotypes
  # m is the proportion of causal SNPs we want to keep
  # ld_range is the min/max values for a uniform distribution
  # for us to sample from for
  #          the LD for the new SNPs we will create
  #          Note we replace the causal SNPs with the LD ones
  #          (rather than adding them in)
  #          This has the added benefit of stopping dimension size
  #          from blowing-up
  snp_indices_to_replace <- sample(1:dim(g)[2],
                                   as.integer((1 - m) * dim(g)[2]),
                                   replace = FALSE
  )
  ld_for_replacements <- runif(
    length(snp_indices_to_replace),
    ld_range[1], ld_range[2]
  )

  for (i in 1:length(snp_indices_to_replace)) {
    j <- snp_indices_to_replace[i]
    g[, j] <- add_simple_ld(g[, j], ld_for_replacements[i])
  }

  return(g)
}

#' Add a block of SNPs in LD with a given SNP
#'
#' @param snp A vector of genotypes for a single SNP
#' @param block_size The number of SNPs in the block
#'
#' @return A matrix of genotypes for the block of SNPs
#'
#' @export
add_simple_ld_block <- function(snp, block_size = 10) {
  ld_values <- runif(block_size, 0.1, 0.9)
  ld_block <- matrix(0, nrow = length(snp), ncol = block_size)

  for (i in 1:block_size) {
    ld_block[, i] <- add_simple_ld(snp, ld_values[i])
  }

  return(ld_block)
}

#' Simulate a matrix of genotypes with a proportion of causal SNPs in LD
#'
#' Will create n_blocks LD blocks of a size sample uniformly from the
#' range block_size_range.
#'
#' @param p The number of SNPs
#' @param n The number of individuals
#' @param n_blocks The number of LD blocks. Default is 10.
#' @param block_size_range The range of block sizes. Default is c(3, 10).
#'
#' @return A matrix of genotypes with a proportion of causal SNPs in LD
#'
#' @export
simulate_complex_ld_blocks <- function(p = 10, n = 1000, n_blocks = 10, block_size_range = c(3, 10)) {
  g <- add_noise(p = p, n = n)
  block_sizes <- sample(block_size_range[1]:block_size_range[2],
                        n_blocks,
                        replace = TRUE)

  for (i in 1:ncol(g)) {
    block <- add_simple_ld_block(g[, i], block_size = block_sizes[i])

    if (i == 1) {
      ld_blocks <- block
    } else {
      ld_blocks <- cbind(ld_blocks, block)
    }
  }

  return(ld_blocks)
}

#' Simulate a phenotype, given a matrix of genotypes
#'
#' For a given heritability (on the liability scale) and a proportion of causal
#' SNPs, this function will simulate a phenotype based on the genotypes. Note
#' that prevalence, k, and proportion causal SNPs, m, will both influence
#' how easy it is to predict the phenotype from the genotypes, not just
#' heritability. If n_cases and n_controls are not NULL, the function will
#' downsample the cases and controls to the desired numbers. This function
#' is suitable for simulating an additive polygenic background without
#' LD or exceptionally strong main effects from a single locus (i.e.
#' APOE will need to be simulated separately)
#'
#' @param X A matrix or big.matrix of genotypes
#' @param p_causal The proportion of causal SNPs. Default is 0.01.
#' @param h2_l The heritability on the liability scale. Default is 0.2.
#' @param k The prevalence. Default is 0.2.
#' @param n_cases The number of cases. Default is NULL.
#' @param n_controls The number of controls. Default is NULL.
#' @param is_transposed Whether the matrix is transposed (SNPs as rows, individuals as columns).
#' Default is FALSE (individuals as rows, SNPs as columns).
#' @param verbose Whether to print the expected and actual heritability. Default is TRUE.
#'
#' @return A list with five elements:
#' - [[1]] the matrix or big.matrix of individuals*genotypes (X)
#' - [[2]] the binary phenotype (y)
#' - [[3]] the liability-scale phenotype (l)
#' - [[4]] a boolean vector indicating which SNPs are causal
#' - [[5]] a vector of SNP effect sizes on the liability scale
#'
#' @export
simulate_y <- function(X, p_causal = 0.01, h2_l = 0.2, k = 0.2,
                       n_cases = NULL, n_controls = NULL,
                       is_transposed = FALSE, verbose = TRUE) {
  if(any(is.data.frame(X), inherits(X, "tbl_df"))) {
    stop("X must be a matrix")
  }

  if(k == 0 || k == 1) {
    stop("k must be between 0 and 1")
  }

  if(h2_l == 0 || h2_l == 1) {
    stop("h2_l must be between 0 and 1")
  }

  if (is_transposed) {
    n_snps <- nrow(X)
    n_obs <- ncol(X)
    betas_dim <- c(1, n_snps)
  } else {
    n_snps <- ncol(X)
    n_obs <- nrow(X)
    betas_dim <- c(n_snps, 1)
  }

  m <- as.integer(n_snps * p_causal)
  if(m == 0) {
    stop(paste("0 causal SNPs. p_causal must be between",
               1/n_snps,
               "and 1 for", n_snps, "SNPs"))
  }
  betas_causal <- rnorm(m)
  betas_null <- rep(0, n_snps - m)

  # shuffle the betas so causal SNPs are randomly distributed across the dataset
  if (bigmemory::is.big.matrix(X)) {
    betas <- bigmemory::big.matrix(nrow = betas_dim[1], ncol = betas_dim[2], init = 0,
                                   type = "double", separated = FALSE)
  } else {
    betas <- matrix(0, betas_dim[1], betas_dim[2])
  }
  betas_1d <- sample(c(betas_causal, betas_null),
                     size = n_snps, replace = FALSE)

  # derive the genotypic risk from the causal SNPs
  if (is_transposed) {
    betas[1, ] <- betas_1d
    causal_bool <- !dplyr::near(betas_1d, 0) # get index of causal SNPs
    g <- as.vector(betas[, causal_bool] %*% X[causal_bool, ]) # genotypic value
  } else {
    betas[, 1] <- betas_1d
    causal_bool <- !dplyr::near(betas_1d, 0) # get index of causal SNPs
    g <- as.vector(X[, causal_bool] %*% betas[causal_bool, ]) # genotypic value
  }

  # derive the noise (environmental risk, e) to be big enough that we achieve
  # the desired heritability in the final phenotype
  var_g <- var(g)
  var_e <- (var_g / h2_l) - var_g # variance in e (noise) based on variance in g (genotypic risk) and heritability on the liability scale
  e <- rnorm(n_obs, 0, sqrt(var_e)) # simulate error for phenotype for all individuals

  # combine genotypic risk and noise
  l <- g + e
  y <- binary_from_liability(l, k = k)

  # print heritability
  if (verbose) {
    message("h2 liability-scale expected: ", round(h2_l, 2),
            ", h2 liability-scale actual: ", round(var(g) / var(l), 2))
  }

  # handle downsampling of cases/controls to desired numbers
  if(any(!is.null(n_cases), !is.null(n_controls))) {
    ds_ind <- downsample_individuals(X, y, n_cases, n_controls,
                                     is_transposed)
    X <- ds_ind[[1]]
    y <- ds_ind[[2]]
  }

  # return list of outputs
  return(list(X, y, l, causal_bool, betas_1d))
}

#' Convert a liability-scale phenotype to a binary phenotype
#'
#' Requires a liability-scale phenotype and a prevalence (k). The function
#' will use the prevalence as a threshold to convert to a binary phenotype.
#'
#' @param l A liability-scale phenotype
#' @param k The prevalence. Default is 0.2.
#'
#' @return A binary phenotype
#'
#' @export
binary_from_liability <- function(l, k = 0.2) {
  l_z <- (l - mean(l)) / sd(l)
  threshold <- qnorm(1 - k)
  y <- as.integer(l_z > threshold)

  return(y)
}

#' Randomly downsample individuals for a case-control study
#'
#' @param X A matrix or big.matrix of individuals*genotypes
#' @param y A phenotype (binary)
#' @param n_cases The number of cases to downsample to. Default is NULL.
#' @param n_controls The number of controls to downsample to. Default is NULL.
#' @param is_transposed Whether the passed matrix is transposed (SNPs as rows,
#' individuals as columns) or not transposed (individuals as rows, SNPs as columns).
#'
#' @return A list with two elements:
#' - [[1]] the downsampled matrix or big.matrix of individuals*genotypes (X)
#' - [[2]] the downsampled binary phenotype (y)
#'
#' @export
downsample_individuals <- function(X, y, n_cases = NULL, n_controls = NULL,
                                   is_transposed = FALSE) {
  if (!is.null(n_cases)) {
    if (n_cases > sum(y == 1)) {
      stop("n_cases must be less than or equal to the number of cases")
    }

    # downsample cases
    case_idx <- sample(which(y == 1), n_cases)
    y_case <- y[case_idx]
    if (is_transposed) {
      X_case <- X[, case_idx]
    } else {
      X_case <- X[case_idx, ]
    }
  }

  if(!is.null(n_controls)) {
    if (n_controls > sum(y == 0)) {
      stop("n_controls must be less than or equal to the number of controls")
    }

    # downsample controls
    control_idx <- sample(which(y == 0), n_controls)
    y_control <- y[control_idx]
    if (is_transposed) {
      X_control <- X[, control_idx]
    } else {
      X_control <- X[control_idx, ]
    }
  }

  # handle concatenation of transposed matrices
  merge_matrices <- function(X1, X2, is_transposed) {
    if (is_transposed) {
      return(cbind(X1, X2))
    } else {
      return(rbind(X1, X2))
    }
  }

  # combine cases and controls
  if (all(!is.null(n_cases), !is.null(n_controls))) {
    X <- merge_matrices(X_case, X_control, is_transposed)
    y <- c(y_case, y_control)
  } else if (is.null(n_cases)) {
    X <- merge_matrices(X[which(y == 1),], X_control, is_transposed)
    y <- c(y[which(y == 1)], y_control)
  } else if(is.null(n_controls)) {
    X <- merge_matrices(X_case, X[which(y == 0),], is_transposed)
    y <- c(y_case, y[which(y == 0)])
  } else {
    warning("downsample_individuals called with n_cases and n_controls",
            " both NULL - returning original X and y")
  }

  return(list(X, y))
}

#' Upsample individuals in a matrix of genotypes
#'
#' This function will upsample individuals in a matrix of genotypes. It will
#' sample with replacement from the existing individuals to create new
#' individuals. All SNPs are sampled independent instead of as a group,
#' which results in decorrelated SNPs. Note that is_transposed must be
#' set to TRUE if a matrix/data.frame is passed in the classic PLINK
#' format (SNPs as rows, individuals as columns). Note also that nrows
#' and ncols passed have to be swapped if is_transposed is changed too.
#'
#' @param data A matrix, data.frame or big.matrix of individuals*genotypes
#' @param nrows The number of rows to be in the output data
#' @param ncols The number of columns to be in the output data
#' @param n_ticks The number of ticks for the progress bar. Default is 100.
#' @param is_transposed Whether the data is transposed (SNPs as rows,
#' individuals as columns), or not transposed (individuals as rows,
#' SNPs are columns). Default is FALSE.
#'
#' @return A big.matrix of upsampled genotypes of size nrows*ncols.
#'
#' @export
upsample_individuals_independently <- function(data, nrows, ncols, n_ticks = 100,
                                               is_transposed = FALSE) {
  init_matrix <- function(nrows, ncols) {
    # pre-define a file-backed big.matrix to hold upsampled genotypes
    m <- bigmemory::big.matrix(nrow = nrows, ncol = ncols, init = 0,
                               type = "double", separated = FALSE,
                               backingpath = here::here("data"),
                               backingfile = "big_matrix.bin",
                               descriptorfile = "big_matrix.desc")
    options(bigmemory.typecast.warning=FALSE)

    return(m)
  }

  try_init_big_matrix <- function(nrows, ncols) {
    tryCatch({
      m <- init_matrix(nrows, ncols)
      return(m)
    }, error = function(e) {
      message("Error creating big.matrix, attempting to overwite backing files")
      unlink(here::here("data/big_matrix.bin"))
      unlink(here::here("data/big_matrix.desc"))
      m <- init_matrix(nrows, ncols)
      return(m)
    })
  }

  message("Trying to create file-backed big.matrix to hold genotypes...")
  m <- try_init_big_matrix(nrows, ncols)
  if(bigmemory::is.big.matrix(m)) {
    message("Done.\n")
  } else {
    stop("Error creating big.matrix")
  }

  # add progress bar for upsampling
  pb <- progress::progress_bar$new(format = "  Upsampling at each SNP [:bar] :percent eta: :eta",
                                   total = n_ticks,
                                   force = TRUE)

  # upsample rows separately for each column to avoid duplication and drop correlation
  if (is_transposed) {
    max_snps <- nrows
  } else {
    max_snps <- ncols
  }

  for (i in 1:max_snps){
    if (i %% (as.integer(max_snps / n_ticks)) == 0) {
      pb$tick()
    }
    if (is_transposed) {
      m[i, ] <- sample(data[i, ], ncols, replace = TRUE)
    } else {
      m[, i] <- sample(data[, i], nrows, replace = TRUE)
    }
  }

  return(m)
}

#' Simulate APOE alleles and APOE score
#'
#' This function will simulate APOE alleles and APOE score based on the
#' allele frequencies and effect sizes from the config.toml file.
#' The calculated APOE score is the sum of the product of the allele
#' counts and effect sizes for e2 and e4.
#'
#' @param n The number of individuals to simulate
#'
#' @return A tibble with columns for e2, e4, and e2_e4_score
#'
#' @export
simulate_e2_e4 <- function(n) {
  config <- read_config()

  # load empirical effect sizes on the observed scale
  apoe_effects <- unlist(config$apoe_effect_sizes_europeans$kunkle_2019)
  apoe_afreqs <- unlist(config$apoe_allele_freqs_europeans$ukbb_2021_release_rounded)

  # simulate APOE alleles
  e2 <- rbinom(n, 2, apoe_afreqs["e2"])
  e4 <- rbinom(n, 2, apoe_afreqs["e4"])

  e2_e4_score <- (e2 * apoe_effects["e2"]) + (e4 * apoe_effects["e4"])

  results <- tibble::tibble(
    "e2" = e2,
    "e4" = e4,
    "e2_e4_score" = e2_e4_score
  )

  return(results)
}


run_simple_gwas <- function(X, is_transposed = FALSE) {
  NULL
}
