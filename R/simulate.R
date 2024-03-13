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

#' Simulate a phenotype give a matrix of genotypes
#'
#' For a given heritability (on the liability scale) and a proportion of causal
#' SNPs, this function will simulate a phenotype based on the genotypes. Note
#' that prevalence, k, and proportion causal SNPs, m, will both influence
#' how easy it is to predict the phenotype from the genotypes, not just
#' heritability. If n_cases and n_controls are not NULL, the function will
#' downsample the cases and controls to the desired numbers.
#'
#' @param X A matrix of genotypes
#' @param p_causal The proportion of causal SNPs. Default is 0.01.
#' @param h2_l The heritability on the liability scale. Default is 0.2.
#' @param k The prevalence. Default is 0.2.
#' @param n_cases The number of cases. Default is NULL.
#' @param n_controls The number of controls. Default is NULL.
#' @param verbose Whether to print the expected and actual heritability. Default is TRUE.
#'
#' @return A list with three elements: the matrix of genotypes, the phenotype and a logical vector indicating
#' which SNPs are causal.
#'
#' @export
simulate_y <- function(X, p_causal = 0.01, h2_l = 0.2, k = 0.2,
                       n_cases = NULL, n_controls = NULL, verbose = TRUE) {
  if(k == 0 || k == 1) {
    stop("k must be between 0 and 1")
  }

  if(h2_l == 0 || h2_l == 1) {
    stop("h2_l must be between 0 and 1")
  }

  m <- as.integer(ncol(X) * p_causal)
  if(m == 0) {
    stop(paste("0 causal SNPs. p_causal must be between",
               1/ncol(X),
               "and 1 for", ncol(X), "SNPs"))
  }
  betas_causal <- rnorm(m)
  betas_null <- rep(0, ncol(X) - m)

  # shuffle the betas so causal SNPs are randomly distributed across the dataset
  betas <- sample(c(betas_causal, betas_null), ncol(X), replace = FALSE)

  # derive the genotypic risk from the causal SNPs
  causal_bool <- !dplyr::near(betas, 0) # get index of causal SNPs
  g <- as.matrix(X[, causal_bool]) %*% betas[causal_bool] # genotypic value
  var_g <- var(g)

  # derive the noise (environmental risk, e) to be big enough that we achieve
  # the desired heritability in the final phenotype
  var_e <- (var_g / h2_l) - var_g # variance in e (noise) based on variance in g (genotypic risk) and heritability on the liability scale
  e <- rnorm(nrow(X), 0, sqrt(var_e)) # simulate error for phenotype for all individuals

  # combine genotypic risk and noise
  l <- g + e
  l_z <- (l - mean(l)) / sd(l)
  threshold <- qnorm(1 - k)
  p <- as.integer(l_z > threshold)

  # print heritability
  if (verbose) {
    message("h2 liability-scale expected: ", round(h2_l, 2),
            ", h2 liability-scale actual: ", round(var(g) / var(l), 2))
  }

  # handle downsampling of cases/controls to desired numbers
  if (!is.null(n_cases) && !is.null(n_controls)) {
    # handle case where desired n_cases and n_controls are smaller
    # than the observed number of cases/controls
    if (sum(p == 1) < n_cases || sum(p == 0) < n_controls) {
      stop("Desired n_cases/n_controls too small given dim(X). ",
              "Skipping downsampling.")
    }
    # downsample cases
    case_idx <- sample(which(p == 1), n_cases)
    p_case <- p[case_idx]
    X_case <- X[case_idx, ]

    # downsample controls
    control_idx <- sample(which(p == 0), n_controls)
    p_control <- p[control_idx]
    X_control <- X[control_idx, ]

    # combine cases and controls
    X <- rbind(X_case, X_control)
    p <- c(p_case, p_control)

    # handle cases where only one of n_cases or n_controls is not NULL
  } else if (sum(is.null(n_cases), is.null(n_controls)) == 1) {
    stop("Both n_cases and n_controls must be NULL or not NULL")
  }

  # return the phenotype and a logical vector indicating which SNPs are causal
  return(list(X, p, causal_bool))
}
