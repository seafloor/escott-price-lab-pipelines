# Escott-Price Lab Pipelines

[![R-CMD-check](https://github.com/seafloor/escott-price-lab-pipelines/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/seafloor/escott-price-lab-pipelines/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/seafloor/escott-price-lab-pipelines/graph/badge.svg?token=IE0NJEL213)](https://codecov.io/gh/seafloor/escott-price-lab-pipelines)
[![lifecycle](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![DOI](https://zenodo.org/badge/751484131.svg)](https://doi.org/10.5281/zenodo.13898716)

## Welcome to the Escott-Price Lab Pipelines repository! :wave: 

This is an r package and associated workflows for standardized pipelines tailored to genetics research in neurodegenerative disorders.

:exclamation: Please note this repository is a work in progress :exclamation: 

We are currently in the process of migrating our scripts and standardising procedures across our group :raised_hands: so everything here should be considered experimental. Anyone stumbling across this should exercise caution and not use the repository as the basis for analysis until the code reaches a more mature and stable phase.

### :cloud: Quick Overview

The aim here is to simplify and standardize our processing of genetic data related to neurodegenerative diseases. Specific workflows are available in the /workflows folder. The status of each is listed below.

#### Planned pipelines:

- **Pathways (experimental):**  Automated downloading and processing of pathway data for annotations.
- **GWAS (planned):** Standardised GWAS procedures. 
- **QC (planned):** Standardised QC procedures. 
- **PRS (planned):** Calculate PRS with standard approaches (C+T, PRScs etc.).
- **Post-GWAS (planned):** Functional interpretation of GWAS hits.
- **Post-PRS (planned):** Interpretation and comparison of models.

### :seedling: Getting Started Locally

Getting started is now as easy as 1, 2, 3.

1. Install the escott-price-lab-pipelines package:

Code can be installed using devtools. Minimum R version is 4.1.0 - see the DESCRIPTION file for details. Paste the code below into the R console to get up and running.

```r
# Ensure devtools is installed
if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools", quiet = TRUE)
}

# Ensure BiocManager is installed for Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", quiet = TRUE)
}

# Install biomaRt using BiocManager if it's not already installed
if (!require("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt", ask = FALSE, quiet = TRUE)
}

# Install escottpricelabpipelines
devtools::install_github("https://github.com/seafloor/escott-price-lab-pipelines", upgrade = "never", quiet = TRUE)
```

2. Install dependencies:

After installing the package, paste the code below to install anything that can't packaged with the normal r package install process (e.g. bioconductor). This can be skipped, but you will be prompted to install these when running a function that needs them. If you're running a script on a server and are not in an interactive session then it will fail, so we recommend running this now.

```r
library(escottpricelabpipelines)
install_dependencies()
```

3. Install required databases:

We use remote servers where possible to query large databases, but sometimes it's faster to run operations locally, or we haven't yet found an optimal server to use. Paste the code below to download all required databases. Again, this can be skipped, and if you don't install now then you will be prompted for consent to download databases as needed. 

```r
install_databases()
```

### :seedling: Getting Started on HPC

Procedure here is similar. First load your chosen R version on the cluster, and run an interactive session.

```r
module load R
R
```

Then install as above. The big difference is in installing dependencies and databases. This is because we need to account for you being able to load module or store data elsewhere in an HPC environment. To do this you can add the `hpc = TRUE` argument and it will check for module availability and prompt the user to give a directory (something like a /scratch or /workdir directory) to install the databases as the home directory on an HPC server usually doesn't allow for much storage.

**Install dependencies for HPC:**

```r
library(escottpricelabpipelines)
install_dependencies(hpc = TRUE)
```

**Install required databases for HPC:**

```r
# change to directory suitable for holding
# large databases on your server. Note it 
# must be accessible from compute nodes.
install_databases("/scratch/user/escottpricelabpipelines", hpc = TRUE)
```

### Genome build

By default we use GRCh38.p14. Build 38 closed gaps and fixed errors in the GRCh37 build, so is recommended to be used as standard, and patch 14 is the latest stable release. However, we use convenience functions to convert to GRCh37 if needed using biomaRt or the comannd line tool from liftover.

### :page_facing_up: License

This project is licensed under the GNU General Public License v3.0.

### :green_heart: How to Contribute

We're thrilled to have you consider contributing to the Escott-Price Lab Pipelines! Whether you're fixing bugs, adding features, or improving documentation, your help is welcome.

Note that the core functions are part of the escottpricelabpipelines package (in /R), and the workflow scripts (in /workflows) use the `box` package [because just importing all functions from a package is mad and doesn't happen in other major languages](https://github.com/klmr/box?tab=readme-ov-file#why-box), so you shouldn't need to ever call `library()` to load functions when contributing.
