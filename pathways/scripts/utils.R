

cleanup_files <- function(file_list) {
  for (i in 1:length(file_list)) {
    fn <- file_list[i]

    if (file.exists(fn)) {
      file.remove(fn)
    } else {
      print(fn)
    }
  }
}

list_files <- function(range, prepend, append) {
  file_list <- sapply(range, paste_names, prepend, append)
  
  return(file_list)
}

paste_names <- function(i, pre, post) {
  return(paste(pre, i, post, sep = ""))
}

force_canonical_autosomes <- function(df, chromosome="chr") {
  valid_chromosomes <- as.character(1:22)
  df <- filter(df, get(chromosome) %in% valid_chromosomes)
  
  return(df)
}

set_human_genome <- function(human_build = "grch38") {
  if (human_build == "grch38") {
    human <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",
                     host = human_genome_38_host_path)
  } else if (human_build == "grch37") {
    human <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",
                     host = human_genome_37_host_path)
  } else {
    stop("Human genome build must be grch37 or grch38")
  }
  
  return(human)
}

set_genomes <- function() {
  human <- set_human_genome()
  mouse <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl",
                   host = mouse_genome_host_path)
  
  return(list(human, mouse))
}
