#' Submit a query string to the Open Targets Genetics GraphQL API
#'
#' @param query_string A string containing the query to be submitted to the Open Targets Genetics GraphQL API
#' @param variables A list of variables to be passed to the query
#' @param wait A numeric value indicating the number of seconds to wait before submitting the query
#'
#' @return A list containing the results of the query
#'
#' @export
query_opentargets <- function(query_string, variables = NULL, wait = 0) {
  # enforced wait to avoid repetitive spamming of servers
  Sys.sleep(wait)

  # Set base URL of GraphQL API endpoint
  base_url <- "https://api.genetics.opentargets.org/graphql"

  # Set variables object of arguments to be passed to endpoint
  # variables <- list("rsId" = rsId)

  # Construct POST request body object with query string and variables
  if (is.null(variables)) {
    post_body <- list(query = query_string)
  } else {
    post_body <- list(query = query_string, variables = variables)
  }

  # Perform POST request
  r <- httr::POST(url=base_url, body=post_body, encode='json')

  return(httr::content(r))
}

#' Convert a list of rsids to variant IDs used in Open Targets Genetics.
#'
#' Variant IDs are chrom_post_ref_alt, e.g. 1_12345_A_T. Note multiallelic
#' variants are handled by this function - it will take the first returned
#' by default. It can also take the last or concatenate all
#' returned variant IDs into a single string with a semi-colon separator
#' (see multiallelic_strategy). Note that 'first' or 'last' must be passed
#' if the variant is to be later annotated with a gene in OpenTargets - a
#' query for multiple semi-colon separated variant IDs will not return anything.
#' Note also that the ref/alt are taken from open targets so may be flipped
#' compared to other datasets if MAF ~ 0.5.
#'
#' The function will wait a random number of seconds
#' (uniformly between spacing_min and spacing_max)
#' between requests to avoid spamming the server. Note
#' that the batch size/wait is just to again avoid spamming
#' servers with a large number of variants. This function
#' is only suitable for annotating top hits, not whole
#' summary statistics.
#'
#' @param rsids A character vector of rsids
#' @param spacing_min A numeric value indicating the minimum number of seconds to wait between requests
#' @param spacing_max A numeric value indicating the maximum number of seconds to wait between requests
#' @param batch_size A numeric value indicating the number of rsids to be submitted in each batch
#' @param batch_spacing A numeric value indicating the number of seconds to wait between batches
#' @param multiallelic_strategy How to handle return of multiple alleles for a single rsid. Options are: 'first', 'last', 'all'. Default is 'first'.
#'
#' @return A character vector of variant IDs (formatted as chrom_pos_ref_alt)
#'
#' @export
rsid_to_varid <- function(rsids, spacing_min = 0, spacing_max = 0.1, batch_size = 100, batch_spacing = 1,
                          multiallelic_strategy = 'first') {
  # support batches of 100 to avoid spamming servers

  query_string = "
  query rsid2varid($rsId: String!) {
    search(queryString:$rsId) {
      totalVariants
      variants {
        id
      }
    }
  }"

  varids <- c()
  wait_times <- runif(length(rsids), spacing_min, spacing_max)
  pb <- progress::progress_bar$new(format = "  Converting rsids [:bar] :current/:total (:percent)",
                                   total = length(rsids),
                                   force = TRUE)


  for(i in 1:length(rsids)) {
    pb$tick()
    if(i %% batch_size == 0) {
      Sys.sleep(batch_spacing)
    }
    rsId <- rsids[i]
    wait <- wait_times[i]

    varId <- query_opentargets(query_string, variables = c(list("rsId" = rsId)), wait = wait)
    varId <- unlist(varId$data$search$variants)

    # some rsids return e.g. '1_1732685_T_C;1_1732685_T_G' for 'rs2294489'
    # decide which SNP to keep (first, last, or both)
    # note latter will not match any genes in OpenTargets
    if (multiallelic_strategy == 'first') {
      varId <- varId[[1]]
    } else if (multiallelic_strategy == 'last') {
      varId <- varId[[length(varId)]]
    } else if (multiallelic_strategy == 'all') {
      varId <- stringr::str_c(varId, collapse = ";")
    } else {
      stop("multiallelic_strategy must be one of 'first', 'last', or 'all'")
    }

    if (is.null(varId) || length(varId) == 0 || varId == "") {
      varId <- "NA"
    }

    varids <- c(varids, varId)
  }


  return(varids)
}

#' Convert a variant ID to a list of genes and associated functional annotations from Open Targets Genetics
#'
#' The function will return 0 or more rows for each variant.
#' 0 for variants that are not found in the Open Targets Genetics database,
#' 1 for single matches (where it has returned the best prediction from Open Targets)
#' And more than one if there are multiple genes within 1 s.d. of the top gene,
#' or if there are multiple functional annotations or both functional and positional
#' annotations (e.g. QTl data annotates to one genes but the variant is intronic to another gene).
#'
#' @param variantId A character vector of variant IDs (formatted as chrom_pos_ref_alt)
#' @param wait A numeric value indicating the number of seconds to wait before submitting the query
#'
#' @return A tibble containing the results of the query
#'
#' @export
varid_to_genes <- function(variantId, wait = 0) {
  # Build query string
  query_string = "
  query v2g($variantId: String!) {
    genesForVariant(variantId: $variantId) {
      gene {
        id
        symbol
        chromosome
        start
        end
      }
      variant
      overallScore
       qtls {
         sourceId
         tissues {
          tissue {
            name
          }
         }
       }
      functionalPredictions {
        sourceId
        tissues {
          maxEffectLabel
        }
      }
      intervals {
        sourceId
        tissues {
          tissue {
            name
          }
        }
      }
      distances {
        sourceId
        tissues {
          distance
        }
      }
    }
  }"

  # set function for flattening functional annotations
  concat_cols <- function(data, new_col, old_col_prefix) {
    if (ncol(dplyr::select(data, dplyr::starts_with(old_col_prefix))) == 0) {
      out <- dplyr::mutate(data, !!new_col := NA)
    } else {
      out <- tidyr::unite(data,
                          !!new_col,
                          dplyr::starts_with(old_col_prefix),
                          sep = "; ",
                          na.rm = TRUE)
    }

    return(out)
  }

  genes <- query_opentargets(query_string, variables = c(list("variantId" = variantId)), wait = wait)

  # clean output
  if(length(genes$data$genesForVariant) > 0) {
    genes <- lapply(genes$data$genesForVariant, rlist::list.flatten) %>%
      lapply(janitor::clean_names) %>%
      dplyr::bind_rows() %>%
      dplyr::arrange(desc(overall_score)) %>%
      concat_cols("distance_source", "distances_source_id") %>%
      concat_cols("distance_bp", "distances_tissues_distance") %>%
      concat_cols("qtl_source", "qtls_source_id") %>%
      concat_cols("qtl_tissues", "qtls_tissues_tissue_name") %>%
      concat_cols("interval_source", "intervals_source_id") %>%
      concat_cols("interval_tissues", "intervals_tissues_tissue_name") %>%
      concat_cols("function_source", "functional_predictions_source_id") %>%
      concat_cols("function_prediction", "functional_predictions_tissues_max_effect_label") %>%
      dplyr::mutate(interval_source = stringr::str_replace_all(interval_source, "javierre2016", "PCHi-C (Javierre, 2016)"),
                    interval_source = stringr::str_replace_all(interval_source, "jung2019", "PCHi-C (Jung, 2019)"),
                    interval_source = stringr::str_replace_all(interval_source, "andersson2014", "Enhancer-TSS corr (Andersson, 2014)"),
                    interval_source = stringr::str_replace_all(interval_source, "thurman2012", "DHS Promotor corr (Thurman, 2012)"),
                    qtl_source = stringr::str_replace_all(qtl_source, "qtl", "QTL"),
                    function_prediction = stringr::str_replace_all(function_prediction, "(.+)", "Positional (\\1)"))
  } else {
    print(paste("Query for variant", variantId, "in OpenTargets returned empty list"))
    return(tibble::tibble())
  }

  # summarise evidence
  genes <- genes %>%
    dplyr::mutate_if(is.character, list(~dplyr::na_if(.,""))) %>%
    tidyr::unite(functional_evidence_summary,
                 c(qtl_source, interval_source, function_prediction),
                 sep = "; ",
                 na.rm = TRUE,
                 remove = FALSE)

  # get genes within 1 s.d. of top gene and VEP-annotated genes
  score_dist <- as.vector(scale(genes[['overall_score']]))
  genes <- genes %>%
    dplyr::filter(!is.na(function_prediction) | score_dist > score_dist[1] - 1) %>%
    dplyr::rename(variant_id = variant) %>%
    dplyr::select(variant_id, gene_id, gene_symbol, gene_chromosome, gene_start, gene_end, overall_score, functional_evidence_summary, dplyr::everything())

  return(genes)
}

#' Convert a list of rsids to a list of genes and associated functional annotations from Open Targets Genetics
#'
#' This function is a wrapper to call both rsid_to_varid and varid_to_genes.
#' See `?varid_to_genes` for more details on the output.
#'
#' @param rsids A character vector of rsids
#' @param spacing_min A numeric value indicating the minimum number of seconds to wait between requests
#' @param spacing_max A numeric value indicating the maximum number of seconds to wait between requests
#' @param batch_size A numeric value indicating the number of rsids to be submitted in each batch
#' @param batch_spacing A numeric value indicating the number of seconds to wait between batches
#' @param take_top_n_genes A numeric value indicating the max number of genes to return for each variant
#'
#' @return A list containing (index 1), the annotated genes and (index 2), the rsids that were not found in the Open Targets Genetics database
#'
#' @export
rsids_to_functional_genes <- function(rsids, spacing_min = 0, spacing_max = 0.1, batch_size = 100, batch_spacing = 1,
                                      take_top_n_genes = 3) {
  # convert rsids to variants IDs (chrom pos ref alt)
  variant_ids <- rsid_to_varid(rsids, spacing_min, spacing_max)
  assertthat::assert_that(length(variant_ids) == length(rsids))

  # check if any of variant_ids are NA or empty and return list of corresponding rsids
  rsid_missing <- is.na(variant_ids) | variant_ids == ""
  if(any(rsid_missing)) {
    message(paste("Variant ID query for", sum(rsid_missing), "rsids in OpenTargets returned empty list:"))
    message(rsids[rsid_missing])
  }

  # subset to non-missing variant_ids for search
  variant_ids <- variant_ids[!rsid_missing]
  gene_wait_times <- runif(length(variant_ids), spacing_min, spacing_max)
  pb <- progress::progress_bar$new(format = "  Annotating genes [:bar] :current/:total (:percent)",
                                   total = length(variant_ids),
                                   force = TRUE)

  # convert each to a list of genes
  annotations <- tibble::tibble()
  for(i in 1:length(variant_ids)) {
    pb$tick()
    if(i %% batch_size == 0) {
      Sys.sleep(batch_spacing)
    }

    variantId <- variant_ids[i]
    wait <- gene_wait_times[i]
    genes <- varid_to_genes(variantId, wait = wait) %>%
      dplyr::mutate(rsid = rsids[i]) %>%
      dplyr::select(rsid, dplyr::everything())

    if(nrow(genes) > 0) {
      genes <- dplyr::slice_max(genes, n = take_top_n_genes, order_by = overall_score)
      annotations <- rbind(annotations, genes)
    }
  }

  # check if any variant_ids are missing from the annotations
  missing_variants <- setdiff(variant_ids, annotations[['variant_id']])
  if(length(missing_variants) > 0) {
    message(paste("Gene query for", length(missing_variants), "variants in OpenTargets returned empty list:"))
    message(missing_variants)
  }

  # combine the vectors for missing variants for variant and gene annotation
  # make unique to return for positional annotations downstream
  missing_variant_ids <- unique(c(variant_ids[rsid_missing], missing_variants))

  # get the corresponding rsid
  missing_rsids <- rsids[match(missing_variant_ids, variant_ids)]

  return(list(annotations, missing_rsids))
}


rsids_to_positional_genes <- function(rsids) {

}


