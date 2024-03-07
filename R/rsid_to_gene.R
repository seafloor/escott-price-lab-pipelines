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
      mutate(interval_source = stringr::str_replace_all(interval_source, "javierre2016", "PCHi-C (Javierre, 2016)"),
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
#' @return A list containing (index 1), the annotated genes and (index 2), the rsids that were not found in the Open Targets Genetics database
#' 
#' @export
rsids_to_functional_genes <- function(rsids, spacing_min = 0, spacing_max = 1, batch_size = 50, batch_spacing = 5) {
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
  
  # spacer to avoid spamming
  Sys.sleep(2)
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
    annotations <- rbind(annotations, genes)
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

library(biomaRt)
library(escottpricelabpipelines)

#### get pos from rsid
rsids <- c("rs11206378", "rs2479409")

rsid_to_chrpos <- function(rsids) {
  mart <- useEnsembl("snp", dataset = "hsapiens_snp")
  
  # as grch37
  # snp_mart <- useEnsembl(biomart="ENSEMBL_MART_SNP",
  #                        host="grch37.ensembl.org",
  #                        dataset="hsapiens_snp")
  
  # mrcieu_gwas_vcf <- c("#CHROM, POS, ID, REF, ALT, QUAL,
  #                      FILTER, INFO", "FORMAT")
  
  allele_match <- "^(.+)/(.+)$"
  
  annot <- getBM(attributes = c('chrom_start', 'chrom_strand', 'refsnp_id' ,'allele'),
                 filters = c('snp_filter'),
                 values = rsids,
                 mart = mart) %>%
    dplyr::mutate(ref = stringr::str_extract(allele, allele_match, group = 1),
                  alt = stringr::str_extract(allele, allele_match, group = 2)) %>%
    dplyr::rename(rsid = refsnp_id,
                  chr = chrom_strand,
                  pos = chrom_start) %>%
    select(chr, pos, rsid, ref, alt) %>%
    tibble::as_tibble()
  
  return(annot)
}
chrpos <- rsid_to_chrpos(rsids)

# get chrpos from gene list
genes <- c("PABPC4;HEYL", "PCSK9")
flat_genes <- unlist(stringr::str_split(genes, ";"))
coords <- escottpricelabpipelines::get_regions_from_genes(flat_genes, gene_format = "hgnc_symbol")

# get closest mapped gene
myannot <- tribble(
  ~rsid, ~gene,
  "rs11206378", "PABPC4;HEYL",
  "rs2479409", "PCSK9"
)

tmp <- dplyr::inner_join(chrpos, myannot, by = 'rsid')

distance_to_gene <- function(snp_pos, start, end) {
  if (snp_pos >= start && snp_pos <= end) {
    return(0)
  } else if(snp_pos < start) {
    return(start - snp_pos)
  } else {
    return(snp_pos - end)
  }
}

match_closest_gene <- function(row, table_of_gene_coords) {
  rsid_pos <- as.integer(row['pos'])
  gene_options <- row['gene']
  
  # pass gene options as string and convert to vector
  gene_options <- unlist(stringr::str_split(gene_options, ';'))
  
  distances <- c()
  print(length(gene_options))
  for (i in 1:length(gene_options)) {
    current_gene <- gene_options[i]
    current_gene_coords <- filter(table_of_gene_coords, hgnc_symbol == current_gene)
    row_distances <- c()
    if(nrow(current_gene_coords) > 0){
      for (j in 1:nrow(current_gene_coords)) {
        current_row <- current_gene_coords[j, ]
        print(current_row)
        print(rsid_pos)
        row_distances <- c(row_distances,
                           distance_to_gene(rsid_pos, current_row$gene_start,
                                            current_row$gene_end))
      }
      distances <- c(distances, min(row_distances))
    } else {
      distances <- c(1)
    }
  }
  
  # handle the case where genes are equidistant from snp
  gene_idx <- which(distances == min(distances, na.rm = TRUE))
  closest_gene <- stringr::str_c(gene_options[gene_idx], collapse = ";")
  
  return(closest_gene)
}

my_genes <- readr::read_csv("~/Desktop/genes.csv")
unique_snps <- dplyr::distinct(my_genes)

flat_genes <- unlist(stringr::str_split(unique_snps[['Gene (Annovar)']], ";"))

coords <- escottpricelabpipelines::get_regions_from_genes(flat_genes, gene_format = "hgnc_symbol")
unique_genes <- unique(flat_genes)
coords <- escottpricelabpipelines::get_regions_from_genes(unique_genes, gene_format = "hgnc_symbol")

colnames(my_genes) <- c('chr', 'pos', 'rsid', 'gene')

my_closest <- apply(my_genes, 1, match_closest_gene, table_of_gene_coords = coords)
my_genes['closest_biomart'] <- my_closest



# map genes to cytogenic band
get_cytogenic_band <- function(genes) {
  mart <- useEnsembl("hsapiens_gene_ensembl", biomart = "ensembl")
  bands <- biomaRt::getBM(
    attributes = c("hgnc_symbol", "chromosome_name", "band"),
    mart = mart,
    filters = "hgnc_symbol",
    values = genes
  ) %>%
    dplyr::mutate(band_full = stringr::str_c(chromosome_name, band)) %>%
    rename(chr = chromosome_name)
  
  return(bands)
}

flatten_gene_list <- function(genes) {
  return(unlist(stringr::str_split(genes, ";|,")))
}

flat_consensus_genes <- flatten_gene_list(myannot[['Consensus Gene']])
consensus_bands <- get_cytogenic_band(unique(flat_consensus_genes))
consensus_bands <- escottpricelabpipelines::force_canonical_autosomes(consensus_bands)
consensus_bands <- tibble::as_tibble(consensus_bands)

#
all_bands <- c()
for(i in 1:nrow(myannot)){
  curr_genes <- flatten_gene_list(myannot[i, ][['Consensus Gene']])
  curr_bands <- c()
  for(j in 1:length(curr_genes)) {
    b <- filter(consensus_bands, hgnc_symbol == curr_genes[j])[['band_full']]
    curr_bands <- c(curr_bands, b)
  }
  
  if(length(curr_bands) >= 1){
    all_bands <- c(all_bands, stringr::str_c(curr_bands, collapse = ";"))
  } else {
    all_bands <- c(all_bands, "NA")
  }
}

myannot['Consensus Band'] <- all_bands
readr::write_csv(myannot, '~/Desktop/complete_annotations_with_bands.csv')


### merge with classifier list
full_annot <- readr::read_csv('~/Desktop/complete_annotations_with_bands.csv')
supp_xgb <- readxl::read_excel('~/Downloads/eadb final draft/eadb_ml_supplementary_draft1_20_02_24.xlsx',
                               sheet = 'Table S1', range = "A1:N109")
supp_xgb_noapoe <- readxl::read_excel('~/Downloads/eadb final draft/eadb_ml_supplementary_draft1_20_02_24.xlsx',
                                      sheet = 'Table S2', range = "A1:M141")
supp_nn <- readxl::read_excel('~/Downloads/eadb final draft/eadb_ml_supplementary_draft1_20_02_24.xlsx',
                              sheet = 'Table S3', range = "A1:O148")

res <- rbind(select(supp_xgb, rsid, Model),
             select(supp_xgb_noapoe, rsid, Model),
             select(supp_nn, rsid, Model))
res <- res %>%
  group_by(rsid) %>%
  mutate(Model_by_rsid = stringr::str_c(Model, collapse = ";"))

res <- dplyr::distinct(res, rsid, .keep_all = TRUE)

res <- dplyr::left_join(full_annot, select(res, rsid, Model_by_rsid), by = "rsid") %>%
  rename(Model = Model_by_rsid)

readr::write_csv(res, "~/Desktop/gene_annotation_by_model.csv")


# updating mbmdr
library(readr)
library(readxl)

apoe_1d <- readxl::read_excel("withapoe/mbmdr_output_1d_withGENES_mapped.xlsx")
apoe_2d <- readxl::read_excel("withapoe/mbmdr_output_2d_withGENES_mapped.xlsx")
noapoe_1d <- readxl::read_excel("withoutapoe/mbmdr_output_1d_withGENES_mapped.xlsx")
noapoe_2d <- readxl::read_excel("withoutapoe/mbmdr_output_2d_withGENES_mapped.xlsx")

apoe_1d['model'] <- '1d_apoe'
apoe_2d['model'] <- '2d_apoe'
noapoe_2d['model'] <- '2d_noapoe'
noapoe_1d['model'] <- '1d_noapoe'

df_1d <- rbind(apoe_1d, noapoe_1d)
df_2d <- rbind(apoe_2d, noapoe_2d)

dplyr::left_join(apoe_1d, annot, by = join_by(First_Marker == rsid)) |>
  filter(p_value <= 0.999) |>
  group_by(First_Marker) |>
  dplyr::distinct(First_Marker, .keep_all = TRUE) |>
  select(First_Marker, T_value, p_value, chr, pos,
         `Gene (Closest positional; annovar)`,
         `Gene (Functional; OpenTargets)`,
         `Consensus Gene`,
         `Consensus Band`,
         `Annotation Evidence`,
         Model) |>
  rename(rsid = First_Marker,
         t_value = T_value) |>
  readr::write_csv("mbmdr_apoe_1d_annotated.csv")




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

chr_name
chrom_start
chrom_end
ensembl_gene_stable_id
ensembl_gene_name

chr_name
start
band_start
ensembl_gene

tester <- biomaRt::getBM(
  attributes = c("chr_name", "chrom_start", "chrom_end", "ensembl_gene_stable_id", "ensembl_gene_name"),
  mart = ensembl,
  filters = c("snp_filter"),
  values = c("rs3851179")
)
