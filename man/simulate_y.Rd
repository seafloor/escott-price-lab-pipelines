% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{simulate_y}
\alias{simulate_y}
\title{Simulate a phenotype, given a matrix of genotypes}
\usage{
simulate_y(
  X,
  p_causal = 0.01,
  h2_l = 0.2,
  k = 0.2,
  n_cases = NULL,
  n_controls = NULL,
  is_transposed = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{X}{A matrix or big.matrix of genotypes}

\item{p_causal}{The proportion of causal SNPs. Default is 0.01.}

\item{h2_l}{The heritability on the liability scale. Default is 0.2.}

\item{k}{The prevalence. Default is 0.2.}

\item{n_cases}{The number of cases. Default is NULL.}

\item{n_controls}{The number of controls. Default is NULL.}

\item{is_transposed}{Whether the matrix is transposed (SNPs as rows, individuals as columns).
Default is FALSE (individuals as rows, SNPs as columns).}

\item{verbose}{Whether to print the expected and actual heritability. Default is TRUE.}
}
\value{
A list with five elements:
\itemize{
\item [\link{1}] the matrix or big.matrix of individuals*genotypes (X)
\item [\link{2}] the binary phenotype (y)
\item [\link{3}] the liability-scale phenotype (l)
\item [\link{4}] a boolean vector indicating which SNPs are causal
\item [\link{5}] a vector of SNP effect sizes on the liability scale
}
}
\description{
For a given heritability (on the liability scale) and a proportion of causal
SNPs, this function will simulate a phenotype based on the genotypes. Note
that prevalence, k, and proportion causal SNPs, m, will both influence
how easy it is to predict the phenotype from the genotypes, not just
heritability. If n_cases and n_controls are not NULL, the function will
downsample the cases and controls to the desired numbers. This function
is suitable for simulating an additive polygenic background without
LD or exceptionally strong main effects from a single locus (i.e.
APOE will need to be simulated separately)
}
