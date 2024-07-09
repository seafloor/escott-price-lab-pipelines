#' Install Required Dependencies
#'
#' This function checks for and installs required Bioconductor packages for the package,
#' and any command-line tools needed, e.g. plink, bcftools etc.
#' Currently, it ensures that `biomaRt` is installed, using `BiocManager`.
#'
#' @details The function first checks if `BiocManager` is installed and installs it if necessary.
#' Then, it checks for the presence of specified Bioconductor packages (`biomaRt`) and installs any that are missing.
#' This is particularly useful for setting up the package environment or ensuring that dependencies are met for new users.
#'
#' @param hpc If TRUE, checks modules can be loaded with appropriate tools first.
#'
#' @examples
#' \dontrun{
#' install_dependencies()
#' }
#'
#' @export
install_dependencies <- function(hpc = FALSE) {
  # install bioconductor and biomaRt
  message("Checking for required Bioconductor packages...")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  required_bioc_packages <- c("biomaRt")

  missing_packages <- required_bioc_packages[!sapply(required_bioc_packages, requireNamespace, quietly = TRUE)]

  if (length(missing_packages) > 0) {
    message("The following Bioconductor packages are required but not installed: ", paste(missing_packages, collapse = ", "))
    BiocManager::install(missing_packages)
  }

  # installing genetics packages
  message("Checking for required command-line tools...")
  tools_to_check <- c("plink", "plink2", "bcftools")
  setup_bioinformatics_tools(tools_to_check, hpc = hpc)
}

#' Set Database Directory
#'
#' Creates a directory for holding files and updates the path in the databases.toml file.
#' Useful where the install for R code is in a home directory but large database
#' files need to be stored elsewhere. Note that the location should be accessible
#' from compute nodes in the cluster as it will be accessed during jobs.
#'
#' @param path The full path to the directory to store files in.
#' @param recursive Whether to create any parent directories if they do not exist.
#'
#' @export
set_database_dir <- function(path, recursive = FALSE) {
  # create in database directory for holding files
  if (!dir.exists(path)) {
    message("Creating")
    dir.create(path, recursive = recursive)
  }

  if (dir.exists(path)) {
    message("Done.")
  } else {
    stop("Error creating directory.")
  }

  # read the toml file and overwrite db path
  db_toml <- system.file("inst", "extdata", "databases.toml", package = "escottpricelabpipelines")
  f <- file(db_toml, open = "r+")
  header <- readLines(f, n = 2)
  header[2] <- paste("download_dir = ", "\"", "/my/path/to/stuff", "\"", sep = "")

  # save to the toml file
  t <- tempfile()
  t2 <- tempfile()
  writeLines(header, t)
  writeLines(full[3:length(full)], t)
  close(f)

  # move and cleanup temp files
  unlink(db_toml)
  file.rename(t, db_toml)
  unlink(t)
  unlink(t2)

  # check the update worked
  my_toml <- RcppTOML::parseTOML(db_toml)
  if (my_toml$installation == path) {
    message("Database directory updated in databases.toml to ", path)
  } else {
    stop("Error updating database directory in databases.toml. ",
         "Manually update the download_dir field in ", db_toml)
  }
}

#' Get Install Directory from Database TOML File
#'
#' This function reads the databases.toml file and returns the installation path.
#' Will use the local data directory if the installation path is set to "data".
#'
#' @export
get_database_dir <- function() {
  db_path <- read_database_toml()$installation

  if (db_path == "data") {
    db_path <- system.file("data", package = "escottpricelabpipelines")
  }

  return(db_path)
}

#' Install Required Databases
#'
#' This function downloads and installs the required databases for the package.
#' This must be run interactively, as it requires user input to confirm the download
#' due to the size of some databases.
#'
#' @param path Install location for databases. Will update databases.toml with new location.
#' @param hpc If TRUE, assumes user needs a separate storage location and can
#' load modules from a cluster environment.
#'
#' @export
install_databases <- function(path = NULL, hpc = FALSE) {
  if (!interactive()) {
    message("This process requires user input. Please run interactively.")
    return()
  }

  # set db install location if passed
  if (!is.null(path)) {
    set_database_dir(path)
  }

  # check if manual database directory is set for hpc
  if (hpc) {
    db_path <- get_database_dir()
    if (db_path == "data") {
      message("Database directory not set. ",
              "Please set a directory to store databases.")
      return()
    }
  }

  # prompt user to confirm they want to download databases
  message("This will download and install the following databases:")
  message("  - liftover (approx. 50MB)")
  message("Do you want to continue? (y/n)")
  response <- readline()
  if (tolower(response) != "y") {
    message("Aborting.")
    return()
  }

  # install required liftover files
  download_liftover()
}

download_reference_genome_coordinates <- function() {

}

download_sequence_report <- function() {

}

download_reference_genome <- function() {

}

download_dbsnp <- function(file_url, version = "156") {
  print("--> Installing dbSNP")
  warning("Will use approximately 25GB of disk space")

  # download dbSNP

  # annotate human-readable chromosome names


}

#' Check for Liftover Files
#'
#' This function checks for the presence of the required liftover files for the package.
#' If the files are not present, it will attempt to download them by calling download_liftover().
#' If there has already been an attempt to download the files, it will stop and remove the tar.gz file.
#'
#' @param attempt The number of attempts made to download the liftover files. Default is 0.
#'
#' @export
check_for_liftover <- function(attempt = 0) {
  db_path <- get_database_dir()

  # check that liftover files are present
  if (!file.exists(file.path(db_path, "liftover", "hg19ToHg38.over.chain.gz")) ||
      !file.exists(file.path(db_path, "liftover", "hg38ToHg19.over.chain.gz")) ||
      !file.exists(file.path(db_path, "liftover", "liftOver")) ||
      !file.exists(file.path(db_path, "liftover", "liftOver_linux"))) {
    if (attempt > 0) {
      unlink(file.path(db_path, "liftover.tar.gz"))
      stop("Liftover files not present. Download failed.")
    } else {
      download_liftover(attempt = 1)
    }
  } else {
    message("Liftover files found.")
  }
}

#' Install Required Liftover Files
#'
#' This function downloads the required liftover files for the package,
#' including chain files and the max/linux executables.
#'
#' @export
download_liftover <- function(attempt = 1) {
  # get liftover download info
  db_info <- read_database_toml()
  db_path <- get_database_dir()

  # download liftover
  message("--> Downloading liftover and chain files (approx. 50Mb space after extracting)...")
  download.file(db_info$liftover$url, destfile = file.path(db_path, "liftover.tar.gz"))
  md5 <- tools::md5sum(file.path(db_path, "liftover.tar.gz"))
  if (as.vector(md5) == db_info$liftover$checksum) {
    message("done.")
  } else {
    unlink(file.path(db_path, "liftover.tar.gz"))
    stop("Downloaded file does not match expected md5sum")
  }

  # extract and check liftover files
  message("--> Extracting liftover files...", appendLF = FALSE)
  tryCatch({
    untar(file.path(db_path, "liftover.tar.gz"), exdir = file.path(db_path, "liftover"))
  }, error = function(e) {
    unlink(file.path(db_path, "liftover.tar.gz"))
    stop("Error extracting liftover files: ", conditionMessage(e))
  })

  # check liftover files
  check_for_liftover(attempt = attempt)

  # remove tar.gz
  unlink(file.path(db_path, "liftover.tar.gz"))

  message("done.")
}

#' Check Tool Runs
#'
#' This function checks if a tool is available on the system path
#' and can be executed. It does this by running the tool with the
#' `--version` flag and checking the output.
#'
#' @param tool_command The command to check for on the system path.
#'
#' @export
check_tool_runs <- function(tool_command, load_command = NULL) {
  tool_check_command <- paste(tool_command, "--version")

  if (!is.null(load_command)) {
    tool_check_command <- paste("module purge",
                                paste("module load", load_command),
                                tool_check_command, sep = "; ")
  }

  try_call_tool <- function(tool_check_command) {
    tryCatch({
      tool_check_result <- system(tool_check_command,
                                  intern = TRUE,
                                  ignore.stderr = TRUE)
      return(tool_check_result)
    }, error = function(e) {
      return(c())
    })
  }

  tool_check_result <- try_call_tool(tool_check_command)

  if (length(tool_check_result) > 0 && !grepl("not found|No such file or directory", tool_check_result[1])) {
    message("Running ", tool_command,
            " version: ", tool_check_result[1])
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Check and Load Module
#'
#' This function checks for the availability of a module on the system
#' and attempts to load it. It then checks if the module can be executed
#' successfully.
#'
#' @param module_name The name of the module to check and load.
#'
#' @export
check_and_load_module <- function(module_name) {
  message("\n--> Checking module availability for ", module_name, "...", appendLF = FALSE)

  # check module load
  module_available <- system(paste("module avail",
                                   module_name,
                                   "2>&1 | grep",
                                   module_name),
                             intern = TRUE)

  module_available <- unlist(stringr::str_split(module_available, " +"))

  # check can load, and handle multiple matching modules being available
  if (length(module_available) > 1) {
    message("multiple modules found")

    message("Please select module by entering the number:\n")
    for (i in seq_along(module_available)) {
      cat(i, ": ", module_available[i], "\n", sep = "")
    }

    user_input <- as.integer(readline(prompt = "Option number: "))

    if (user_input %in% seq_along(module_available)) {
      message("\nYou selected:", module_available[user_input], "\n")
      module_available <- module_available[user_input]
    } else {
      message("Invalid selection. Please run the script again and enter a valid number.\n")
      return()
    }

  } else if (length(module_available) == 1) {
    message("found. Loading module...")
  }

  system(paste("module purge; module load", module_available), intern = TRUE)

  # check can execute after load
  if (length(module_available) > 0) {
    if (check_tool_runs(module_name, module_available)) {
      message("loads and executes successfully.")
      return(TRUE)
    } else {
      message("found but cannot be executed.")
      return(FALSE)
    }
  } else {
    message("not found.")
    return(FALSE)
  }
}

#' Setup Bioinformatics Tools
#'
#' This function checks for the availability of required bioinformatics tools
#' on the system path. If the tools are not available, it will attempt to load
#' the tool as a module if the `hpc` argument is set to TRUE. If the module
#' loading fails and a Singularity container is provided, it will attempt to
#' use the container.
#'
#' @param tools A character vector of tool names to check for.
#' @param hpc If TRUE, checks modules can be loaded with appropriate tools first.
#' @param container_paths A named list of Singularity container paths for tools.
#'
#' @export
setup_bioinformatics_tools <- function(tools, hpc = FALSE, container_paths = list()) {
  for (tool in tools) {
    message("--> Checking ", tool, "...")
    tool_available <- check_tool_runs(tool)

    if (!tool_available && hpc) {
      # Attempt to load the tool as a module
      module_loaded <- check_and_load_module(tool)
      if (!module_loaded && !is.null(container_paths[[tool]])) {
        # If module loading fails and a Singularity container is provided, use it
        use_singularity_container(container_paths[[tool]])
      }
    } else if (!tool_available) {
      message("Attempting to install ", tool, "...")
      install_tool(tool)
    }
  }
}

#' Use Singularity Container
#'
#' This function provides guidance on using a Singularity container for a required tool.
#' Function not currently supported
#'
#' @param container_path The path to the Singularity container file.
#'
#' @export
use_singularity_container <- function(container_path) {
  message("Attempting to use Singularity container for the required tool.")
  # Example command to check Singularity availability or pull/load a container
  # Adjust based on your specific needs and container setup
  singularity_check <- system("singularity --version", intern = TRUE, ignore.stderr = TRUE)
  if (length(singularity_check) == 0) {
    message("Singularity is not available. Please contact your system administrator.")
  } else {
    message(paste("Singularity available:", singularity_check))
    # Example of running a command within a Singularity container
    # system("singularity exec /path/to/container.sif some_command")
  }
}

#' Install plink2
#'
#' This function downloads and installs the plink2 binary for the package.
#'
#' @export
install_plink2 <- function() {
  tmp <- tempfile()
  config <- read_config()
  out_dir <- system.file(config$tools[[1]]$location$path, package = "escottpricelabpipelines")

  # check if on mac Darwin
  if (Sys.info()["sysname"] == "Darwin") {
    f <- config$tools[[1]]$mac$plink2
  } else {
    f <- config$tools[[1]]$linux$plink2
  }

  download.file(f, tmp)
  unzip(tmp, exdir = out_dir)
  unlink(tmp)

  update_bash_profile()
  system(paste("chmod u+x", file.path(out_dir, "plink2")))
  runs <- check_tool_runs("plink2")

  if(runs) {
    message(paste("plink2 installed successfully."))
  } else {
    message(paste("Could not install.",
                  "plink2 is not available on this system.",
                  "Please install it manually or check your PATH."))
  }

}

#' Install Tool
#'
#' Wrapper to install specific cli tools.
#'
#' @param tool The name of the tool to install.
#'
#' @export
install_tool <- function(tool) {
  if (tool == "plink2") {
    install_plink2()
  }
}

#' Update .bash_profile
#'
#' This function checks if the CLI install path is in the .bash_profile file.
#' If it is not, it will add the path to the file.
#'
#' @export
update_bash_profile <- function() {
  config <- read_config()
  out_dir <- system.file(config$tools[[1]]$location$path, package = "escottpricelabpipelines")
  cli_tool_line <- paste0("export PATH=$PATH:", out_dir)

  f <- file("~/.bash_profile", open = "r")
  bash_profile <- readLines(f)

  # check if any lines in bash_profile match the cli_tool_line
  path_missing <- !any(grepl(cli_tool_line, bash_profile))
  close(f)

  if(path_missing) {
    message("Adding CLI install path to .bash_profile")
    write(cli_tool_line, file = "~/.bash_profile", append = TRUE)
    system("source ~/.bash_profile")
  }
}