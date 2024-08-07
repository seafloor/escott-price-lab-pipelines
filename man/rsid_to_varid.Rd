% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsid_to_gene.R
\name{rsid_to_varid}
\alias{rsid_to_varid}
\title{Convert a list of rsids to variant IDs used in Open Targets Genetics.}
\usage{
rsid_to_varid(
  rsids,
  spacing_min = 0,
  spacing_max = 1,
  batch_size = 50,
  batch_spacing = 5
)
}
\arguments{
\item{rsids}{A character vector of rsids}

\item{spacing_min}{A numeric value indicating the minimum number of seconds to wait between requests}

\item{spacing_max}{A numeric value indicating the maximum number of seconds to wait between requests}

\item{batch_size}{A numeric value indicating the number of rsids to be submitted in each batch}

\item{batch_spacing}{A numeric value indicating the number of seconds to wait between batches}
}
\value{
A character vector of variant IDs (formatted as chrom_pos_ref_alt)
}
\description{
Variant IDs are chrom_post_ref_alt, e.g. 1_12345_A_T. Note multiallelic
variants are not handled by this function - it will concatenate all
returned variant IDs into a single string with a semi-colon separator.
Note also that the ref/alt are taken from open targets so may be flipped
compared to other datasets if MAF ~ 0.5.
}
\details{
The function will wait a random number of seconds
(uniformly between spacing_min and spacing_max)
between requests to avoid spamming the server. Note
that the batch size/wait is just to again avoid spamming
servers with a large number of variants. This function
is only suitable for annotating top hits, not whole
summary statistics.
}
