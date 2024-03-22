test_that("add_noise returns an n x p array", {
  g <- add_noise(p = 10, 1000)

  expect_equal(dim(g), c(1000, 10))
})

test_that("add_noise returns a matrix", {
  g <- add_noise()

  expect_true(is.matrix(g))
})

test_that("add_noise returns genotypes with MAF in maf_range", {
  tol = 0.02
  maf_min = 0.05
  maf_max = 0.5
  g <- add_noise(p = 10, n = 2000, maf_range = c(maf_min, maf_max))


  maf <- colMeans(g) / 2
  expect_true(all(maf >= (maf_min - tol)))
  expect_true(all(maf <= (maf_max + tol)))
})


test_that("add_simple_ld returns a vector of the same length", {
  snp <- rbinom(1000, 2, 0.5)
  ld_snp <- add_simple_ld(snp, 0.5)

  expect_equal(length(ld_snp), length(snp))
})


test_that("add_simple_ld returns a vector with the same mean", {
  snp <- rbinom(1000, 2, 0.5)
  ld_snp <- add_simple_ld(snp, 0.5)

  expect_equal(mean(ld_snp), mean(snp))
})

test_that("add_simple_ld returns with SNP with the correct r2", {
  snp <- rbinom(1000, 2, 0.5)
  tol <- 0.05 # set 5% tolerance for errors due to low number of repeats
  n_checks <- 10

  ld_snp_0 <- c()
  ld_snp_02 <- c()
  ld_snp_04 <- c()
  ld_snp_06 <- c()
  ld_snp_08 <- c()
  ld_snp_1 <- add_simple_ld(snp, 1.0)

  for (i in 1:n_checks) {
    ld_snp_0 <- c(ld_snp_0, cor(add_simple_ld(snp, 0.0), snp))
    ld_snp_02 <- c(ld_snp_02, cor(add_simple_ld(snp, 0.2), snp))
    ld_snp_04 <- c(ld_snp_04, cor(add_simple_ld(snp, 0.4), snp))
    ld_snp_06 <- c(ld_snp_06, cor(add_simple_ld(snp, 0.6), snp))
    ld_snp_08 <- c(ld_snp_08, cor(add_simple_ld(snp, 0.8), snp))
  }

  expect_equal(mean(ld_snp_0), 0, tolerance = tol)
  expect_equal(mean(ld_snp_02), 0.2, tolerance = tol)
  expect_equal(mean(ld_snp_04), 0.4, tolerance = tol)
  expect_equal(mean(ld_snp_06), 0.6, tolerance = tol)
  expect_equal(mean(ld_snp_08), 0.8, tolerance = tol)

  expect_equal(snp, ld_snp_1)
})

test_that("replace_with_ld returns matrix of same shape and MAF as input", {
  n = 1000
  X <- cbind(rbinom(n, 2, 0.2), rbinom(n, 2, 0.4),
             rbinom(n, 2, 0.2), rbinom(n, 2, 0.4))
  X_ld <- replace_with_ld(X, 0.5)

  # check shape is same after replace with LD SNPs
  expect_equal(dim(X_ld), dim(X))

  # check MAF is same after replace with LD SNPs
  expect_equal(colMeans(X_ld) / 2, colMeans(X) / 2)
})

test_that("simulate_y returns list of length 5", {
  X <- cbind(rbinom(1000, 2, 0.2), rbinom(1000, 2, 0.4))
  y <- simulate_y(X, p_causal = 1)

  expect_equal(length(y), 5)
  expect_true(is.list(y))
})

test_that("simulate_y returns vectors of correct shape", {
  n = 1000
  X <- cbind(rbinom(n, 2, 0.2), rbinom(n, 2, 0.4))
  y <- simulate_y(X, p_causal = 1)

  # check phenotype is correct shape and range
  expect_true(is.numeric(y[[2]]))
  expect_true(all(y[[2]] %in% c(0, 1)))
  expect_equal(length(y[[2]]), n)

  # check liability is correct shape
  expect_true(is.numeric(y[[3]]))
  expect_equal(length(y[[3]]), n)

  # check causal SNPs are logical and correct shape
  expect_true(is.logical(y[[4]]))
  expect_equal(length(y[[4]]), 2)

  # check effect sizes are correct shape
  expect_true(is.numeric(y[[5]]))
  expect_equal(length(y[[5]]), 2)
})

test_that("simulate_y returns the correct number of causal SNPs", {
  n = 1000
  X <- cbind(rbinom(n, 2, 0.2), rbinom(n, 2, 0.4),
             rbinom(n, 2, 0.2), rbinom(n, 2, 0.4))

  # check 0 < p_causal < 1
  y <- simulate_y(X, p_causal = 0.5)
  expect_equal(sum(y[[4]]), 2)

  # check p_causal == 0 raises error
  expect_error(simulate_y(X, p_causal = 0))

  # check p_causal == 1 works
  y <- simulate_y(X, p_causal = 1)
  expect_equal(sum(y[[4]]), 4)
})

test_that("simulate_y returns the correct prevalence", {
  tol = 0.02 # set tolerance to a fairly high value (2%)
  n = 5000
  X <- cbind(rbinom(n, 2, 0.2), rbinom(n, 2, 0.4),
             rbinom(n, 2, 0.2), rbinom(n, 2, 0.4))

  # check 0 < k < 1 is within tolerance
  k = 0.5
  y <- simulate_y(X, k = k, p_causal = 0.5)
  expect_true(mean(y[[2]]) >= k - tol)
  expect_true(mean(y[[2]]) <= k + tol)

  # check errors when k == 0 or k == 1
  expect_error(simulate_y(X, k = 0, p_causal = 0.5))
  expect_error(simulate_y(X, k = 1, p_causal = 0.5))
})

test_that("simulate_y returns a matrix of genotypes of expected dim", {
  # check standard run
  n = 1000
  X <- cbind(rbinom(n, 2, 0.2), rbinom(n, 2, 0.4),
             rbinom(n, 2, 0.2), rbinom(n, 2, 0.4))
  y <- simulate_y(X, p_causal = 0.5)

  # check genotypes are correct shape
  expect_true(is.matrix(y[[1]]))
  expect_equal(dim(y[[1]]), dim(X))
})

test_that("simulate_y downsamples cases and controls correctly", {
  # setup
  n = 1000
  X <- cbind(rbinom(n, 2, 0.2), rbinom(n, 2, 0.4),
             rbinom(n, 2, 0.2), rbinom(n, 2, 0.4))

  # no downsampling when both NULL
  y <- simulate_y(X, p_causal = 0.5, k = 0.5,
                  n_cases = NULL, n_controls = NULL)
  expect_equal(dim(X), dim(y[[1]]))

  # set up 500 cases and 500 controls
  n_cases = 100
  n_controls = 200
  y <- simulate_y(X, p_causal = 0.5, k = 0.5,
                  n_cases = n_cases,
                  n_controls = n_controls)

  # downsample cases/controls correctly
  expect_equal(dim(y[[1]]), c(n_cases + n_controls, 4))
  expect_equal(sum(y[[2]] == 1), n_cases)
  expect_equal(sum(y[[2]] == 0), n_controls)

  # error when desired n_cases/n_controls is too high
  expect_error(simulate_y(X, p_causal = 0.5,
                          n_cases = n + 1,
                          n_controls = n_controls))
  expect_error(simulate_y(X, p_causal = 0.5,
                          n_cases = n_cases,
                          n_controls = n + 1))

  # expect correct when only one of n_cases or n_controls is null
  y <- simulate_y(X, p_causal = 0.5, k = 0.5, n_cases = n_cases, n_controls = NULL)
  expect_equal(sum(y[[2]] == 1), n_cases)

  y <- simulate_y(X, p_causal = 0.5, k = 0.5, n_cases = NULL, n_controls = n_controls)
  expect_equal(sum(y[[2]] == 0), n_controls)
})

# test_that("simulate_y returns the correct heritability", {
#   tol = 0.05
#   n = 1000
#   X <- cbind(rbinom(n, 2, 0.2), rbinom(n, 2, 0.4),
#              rbinom(n, 2, 0.2), rbinom(n, 2, 0.4))
#
#   # check 0 < h2_l < 1 is within tolerance
#   h2_l = 0.5
#   y <- simulate_y(X, h2_l = h2_l, p_causal = 0.5)
#
#   # calculate genotypic values as dot product of genotypes and betas for effect of X[, i] on y
#   b <- apply(X, 2, function(x) glm(y[[1]] ~ x, family = binomial(link = "logit"))$coefficients['x'])
#   g <- X %*% l
#
#   # calculate h2 on observed scale and check is within 5% of h2 on liability scale
#   var_g <- var(g)
#   var_p <- var(y[[1]])
#   expect_true(var_g / var_p >= h2_l - tol)
#   expect_true(var_g / var_p <= h2_l + tol)
#
#   # check errors when h2_l == 0 or h2_l == 1
#   expect_error(simulate_y(X, h2_l = 0, p_causal = 0.5))
#   expect_error(simulate_y(X, h2_l = 1, p_causal = 0.5))
# })
