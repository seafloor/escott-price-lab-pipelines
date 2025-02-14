% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{add_noise}
\alias{add_noise}
\title{Simulate genotypes for genetic association studies}
\usage{
add_noise(n = 1000, p = 100, maf_range = c(0.05, 0.5))
}
\arguments{
\item{n}{number of individuals}

\item{p}{number of SNPs}

\item{maf_range}{range of minor allele frequencies}
}
\value{
A n * p matrix of genotypes
}
\description{
This script simulated unnasociated genotypes only. It will draw minor allele frequency for
SNPs uniformly from maf_range.
}
