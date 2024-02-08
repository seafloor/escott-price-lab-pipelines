# Install CRAN packages
packages <- c("readr", "readxl", "tidyr", "dplyr", "tibble", "janitor", "here")
install.packages(packages)

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt")

# Check for missing dependencies
missing_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(missing_packages)) {
  message("Missing packages: ", paste(missing_packages, collapse = ", "))
  message("Trying to install missing packages...")
  install.packages(missing_packages)
}

message("All dependencies are now installed.")

