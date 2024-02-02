# Escott-Price Lab Pipelines

Welcome to the Escott-Price Lab Pipelines repository! :wave: 

This is our hub for standardized pipelines tailored to genetics research in neurodegenerative disorders. This includes a variety of tools we use to streamline the analysis of Genome-Wide Association Studies (GWAS), Quality Control (QC), and Polygenic Risk Scores (PRS).

:exclamation: Please note this repository is a work in progress :exclamation: We are currently in the process of migrating our scripts and standardising procedures across our group :raised_hands: :raised_hands:, so everything here should be considered experimental. Anyone stumbling across this should exercise caution and not use the repository as the basis for analysis until the code reaches a more mature and stable phase.

## :cloud: Quick Overview

The aim here is to simplify and standardize our processing of genetic data related to neurodegenerative diseases. We leverage a mix of Bash, R, and Python to this effect.

### Planned Features:

- :seedling: **Genetic Pathways Pipeline:** :seedling: Automate the downloading, processing, and standardization of genetic pathway data for annotating SNPs etc.
- :seedling: **GWAS Pipeline:** :seedling: Standardised GWAS procedures. 
- :seedling: **QC Pipeline:** :seedling: May be merged into GWAS or PRS pipelines.
- :seedling: **PRS Pipeline:** :seedling: Calculate Polygenic Risk Scores with standard approaches (C+T, PRScs)
- :seedling: **Post GWAS:** :seedling: Probably a post-GWAS section for functional interpretation of GWAS hits leveraging available tools. May be merged with pathways pipeline.
- :seedling: **Post PRS:** :seedling: Probably a post-PRS section to interpret and compare models.

Each pipeline will be developed with the specific needs of neurodegenerative disorder research in mind, such as running separate analysis with and without the *APOE* region.

## :zap: Getting Started

Again, it's a work in progress, so this is a just a "watch this space" plan for now. But eventually, getting up and running with our pipelines should be as easy as:

1. **Clone this repository:**

```bash
git clone https://github.com/your-repo/escott-price-lab-pipelines.git
cd escott-price-lab-pipelines
```

2. **Install Dependencies:**

Our pipelines use a mix of Bash, R, and Python. Please ensure you have these environments set up, along with any necessary libraries or tools specific to each pipeline. Check out the individual pipeline directories for more detailed setup instructions.

3. **Choose Your Pipeline:**

Navigate to the pipeline directory you're interested in (e.g., GWAS, QC, PRS, or Pathways) and follow the README.md instructions there to start your analysis.

## :page_facing_up: License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## :green_heart: How to Contribute

We're thrilled to have you consider contributing to the Escott-Price Lab Pipelines! Whether you're fixing bugs, adding features, or improving documentation, your help is welcome.
