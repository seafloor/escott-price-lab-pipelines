![travis-build-status](https://app.travis-ci.com/seafloor/escott-price-lab-pipelines.svg?branch=main)

# Escott-Price Lab Pipelines

Welcome to the Escott-Price Lab Pipelines repository! :wave: 

This is our hub for standardized pipelines tailored to genetics research in neurodegenerative disorders. This includes a variety of tools we use to streamline the analysis of Genome-Wide Association Studies (GWAS), Quality Control (QC), and Polygenic Risk Scores (PRS).

:exclamation: Please note this repository is a work in progress :exclamation: We are currently in the process of migrating our scripts and standardising procedures across our group :raised_hands: :raised_hands:, so everything here should be considered experimental. Anyone stumbling across this should exercise caution and not use the repository as the basis for analysis until the code reaches a more mature and stable phase.

## :cloud: Quick Overview

The aim here is to simplify and standardize our processing of genetic data related to neurodegenerative diseases. We leverage a mix of Bash, R, and Python to this effect.

### Planned Features:

- **Genetic Pathways Pipeline:**  Automate the downloading, processing, and standardization of genetic pathway data for annotating SNPs etc.
- **GWAS Pipeline:** Standardised GWAS procedures. 
- **QC Pipeline:** May be merged into GWAS or PRS pipelines.
- **PRS Pipeline:** Calculate Polygenic Risk Scores with standard approaches (C+T, PRScs)
- **Post GWAS:** Probably a post-GWAS section for functional interpretation of GWAS hits leveraging available tools. May be merged with pathways pipeline.
- **Post PRS:** Probably a post-PRS section to interpret and compare models.

Each pipeline will be developed with the specific needs of neurodegenerative disorder research in mind, such as running separate analysis with and without the *APOE* region.

## :zap: Getting Started

Again, it's a work in progress, so this is a just a "watch this space" plan for now. But eventually, getting up and running with our pipelines should be as easy as:

#### :seedling: Clone this repository:

```bash
git clone https://github.com/seafloor/escott-price-lab-pipelines.git
cd escott-price-lab-pipelines
```

#### :seedling: Install Dependencies:

Our pipelines use a mix of Bash, R, and Python. Please ensure you have these environments set up, along with any necessary libraries or tools specific to each pipeline. Check out the individual pipeline directories for more detailed setup instructions.

#### :seedling: Choose Your Pipeline:

Navigate to the pipeline directory you're interested in (e.g., GWAS, QC, PRS, or Pathways) and follow the README.md instructions there to start your analysis.

## :white_check_mark: Requirements

See the DESCRIPTION file for a full list of requirements. Run install_dependencies.R to install everything. Note that the minimum R version is 4.1.0.

Features like liftover over genome builds require other databases and linke files. Setup for these will be in /data.

## Genome build

By default we use GRCh38.p14. Build 38 closed gaps and fixed errors in the GRCh37 build, so is recommended to be used as standard, and patch 14 is the latest stable release. However, all we use convenience functions to convert to GRCh37 if needed using biomaRt or the comannd line tool from liftover.

## :page_facing_up: License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## :green_heart: How to Contribute

We're thrilled to have you consider contributing to the Escott-Price Lab Pipelines! Whether you're fixing bugs, adding features, or improving documentation, your help is welcome.
