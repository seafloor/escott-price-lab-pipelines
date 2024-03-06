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
#' variants are not handled by this function - it will concatenate all 
#' returned variant IDs into a single string with a semi-colon separator. 
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
#' 
#' @return A character vector of variant IDs (formatted as chrom_pos_ref_alt)
#' 
#' @export
rsid_to_varid <- function(rsids, spacing_min = 0, spacing_max = 1, batch_size = 50, batch_spacing = 5) {
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
  
  for(i in 1:length(rsids)) {
    if(i %% batch_size == 0) {
      Sys.sleep(batch_spacing)
    }
    rsId <- rsids[i]
    wait <- wait_times[i]
    
    varId <- query_opentargets(query_string, variables = c(list("rsId" = rsId)), wait = wait)
    varId <- stringr::str_c(unlist(varId$data$search$variants), collapse = ";")
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
      bind_rows() %>% 
      dplyr::arrange(desc(overallScore)) %>%
      concat_cols("distance_source", "distances.sourceId") %>%
      concat_cols("distance_bp", "distances.tissues.distance") %>%
      concat_cols("qtl_source", "qtls.sourceId") %>%
      concat_cols("qtl_tissues", "qtls.tissues.tissue.name") %>%
      concat_cols("interval_source", "intervals.sourceId") %>%
      concat_cols("interval_tissues", "intervals.tissues.tissue.name") %>%
      concat_cols("function_source", "functionalPredictions.sourceId") %>%
      concat_cols("function_prediction", "functionalPredictions.tissues.maxEffectLabel") %>%
      mutate(interval_source = stringr::str_replace_all(interval_source, "javierre2016", "PCHi-C (Javierre, 2016)"),
             interval_source = stringr::str_replace_all(interval_source, "jung2019", "PCHi-C (Jung, 2019)"),
             interval_source = stringr::str_replace_all(interval_source, "andersson2014", "Enhancer-TSS corr (Andersson, 2014)"),
             interval_source = stringr::str_replace_all(interval_source, "thurman2012", "DHS Promotor corr (Thurman, 2012)"),
             qtl_source = stringr::str_replace_all(qtl_source, "qtl", "QTL"),
             function_prediction = stringr::str_replace_all(function_prediction, "(.+)", "Positional (\\1)")) %>%
      janitor::clean_names()
  } else {
    print(paste("Query for variant", variantId, "in OpenTargets returned empty list"))
    return(tibble::tibble())
  }
  
  # summarise evidence
  genes <- genes %>%
    dplyr::mutate_if(is.character, list(~na_if(.,""))) %>%
    tidyr::unite(functional_evidence_summary,
                 c(qtl_source, interval_source, function_prediction),
                 sep = "; ",
                 na.rm = TRUE,
                 remove = FALSE)
  
  # get genes within 1 s.d. of top gene and VEP-annotated genes
  score_dist <- as.vector(scale(genes[['overall_score']]))
  genes <- genes %>%
    filter(!is.na(function_prediction) | score_dist > score_dist[1] - 1) %>%
    rename(variant_id = variant) %>%
    select(variant_id, gene_id, gene_symbol, gene_chromosome, gene_start, gene_end, overall_score, functional_evidence_summary, dplyr::everything())
  
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
#' 
#' @return A tibble containing the results of the query
#' 
#' @export
rsids_to_genes <- function(rsids, spacing_min = 0, spacing_max = 1, batch_size = 50, batch_spacing = 5) {
  # convert rsids to variants IDs (chrom pos ref alt)
  variant_ids <- rsid_to_varid(rsids, spacing_min, spacing_max)
  
  # spacer to avoid spamming
  Sys.sleep(2)
  gene_wait_times <- runif(length(variant_ids), spacing_min, spacing_max)
  
  # convert each to a list of genes
  annotations <- tibble::tibble()
  for(i in 1:length(variant_ids)) {
    if(i %% batch_size == 0) {
      Sys.sleep(batch_spacing)
    }
    variantId <- variant_ids[i]
    wait <- gene_wait_times[i]
    genes <- varid_to_genes(variantId, wait = wait) %>%
      dplyr::mutate(rsid = rsids[i]) %>%
      dplyr::select(rsid, dplyr::everything())
    annotations <- rbind(annotations, genes)
  }
  
  return(annotations)
}

# do mutate(rsid), then rename variant to variant_id, select(rsid, variant_id, everything()), 

##################### QUERY DB VERSION ########################
# query_string = "
# query get_db_version() {
#   meta {
#     name
#     apiVersion {
#       major
#       minor
#       patch
#     }
#     dataVersion
#   }
# }
# "
# 
# db_versions <- query_opentargets(query_string)