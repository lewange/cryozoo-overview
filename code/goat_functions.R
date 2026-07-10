# =============================================================================
# GoaT Tetrapod Data Fetcher
# Fetches all species-level taxon + assembly metadata for tetrapod classes
# from the GoaT API (https://goat.genomehubs.org/api/v2)
#
# Parameters derived from GoaT UI URL:
#   &fields=  contig_n50, assembly_date, gene_count, ebp_standard_criteria,
#             assembly_level, assembly_span, scaffold_n50, busco_completeness,
#             c_value, genome_size, genome_size_draft, chromosome_number,
#             haploid_number, sequencing_status, published
#   &names=   assembly_id, assembly_name, common_name
#   &ranks=   subspecies, species, genus, family, order, class, phylum,
#             kingdom, domain
#
# Sharding strategy to stay under the 10k API offset cap:
#   class -> orders -> families (if order > 9500 spp) -> genera (if family > 9500 spp)
#
# Output: tetrapod_data.rds / .csv — one row per species
# =============================================================================

library(httr2)
library(jsonlite)
library(tidyverse)


# Derived from GoaT UI URL — explicit field/name/rank selection
TAXON_FIELDS <- paste(c(
  "contig_n50", "assembly_date", "gene_count", "ebp_standard_criteria",
  "assembly_level", "assembly_span", "scaffold_n50", "busco_completeness",
  "c_value", "genome_size", "genome_size_draft",
  "chromosome_number", "haploid_number",
  "sequencing_status", "published"
), collapse = ",")

TAXON_NAMES <- "assembly_id,assembly_name,common_name"
TAXON_RANKS <- "subspecies,species,genus,family,order,class,phylum,kingdom,domain"

# =============================================================================
# 1.  SHARED UTILITIES
# =============================================================================

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

dir.create(CACHE_DIR, showWarnings = FALSE)

enc       <- function(q) URLencode(URLencode(q, reserved = TRUE), reserved = TRUE)
GOAT_BASE <- "https://goat.genomehubs.org/api/v2/search"

goat_get <- function(url) {
  repeat {
    resp <- tryCatch(request(url) %>% req_perform(), error = function(e) NULL)
    if (!is.null(resp) && resp_status(resp) == 200) return(resp_body_json(resp))
    message("  retrying: ", substr(url, 1, 80), "...")
    Sys.sleep(5 + runif(1, 0, 3))
  }
}

goat_total <- function(query_str, result_type = "taxon") {
  url <- paste0(GOAT_BASE,
                "?taxonomy=ncbi&result=", result_type,
                "&includeEstimates=true",
                "&query=", enc(query_str),
                "&size=0")
  as.integer(goat_get(url)$status$hits %||% 0L)
}

goat_page <- function(query_str, offset, result_type = "taxon",
                      cache_file = NULL) {
  if (!is.null(cache_file) && file.exists(cache_file)) {
    return(jsonlite::fromJSON(cache_file, simplifyVector = FALSE)$results)
  }

  # taxon endpoint: pass explicit fields/names/ranks as in the GoaT UI URL
  extra <- if (result_type == "taxon") {
    paste0(
      "&fields=",             URLencode(TAXON_FIELDS, reserved = TRUE),
      "&names=",              URLencode(TAXON_NAMES,  reserved = TRUE),
      "&ranks=",              URLencode(TAXON_RANKS,  reserved = TRUE),
      "&includeDescendants=true&emptyColumns=false"
    )
  } else {
    ""
  }

  url <- paste0(GOAT_BASE,
                "?taxonomy=ncbi&result=", result_type,
                "&includeEstimates=true",
                "&query=",  enc(query_str),
                "&size=",   PAGE_SIZE,
                "&offset=", offset,
                "&sortBy=scientific_name&sortOrder=asc",
                extra)

  out <- goat_get(url)
  if (!is.null(cache_file)) jsonlite::write_json(out, cache_file, auto_unbox = TRUE)
  Sys.sleep(0.4 + runif(1, 0, 0.3))
  out$results
}

goat_fetch_all <- function(query_str, result_type = "taxon", cache_prefix = NULL) {
  total <- goat_total(query_str, result_type)
  if (total == 0) return(list())
  if (total > 9500) warning("Results capped at 9500 for: ", query_str)
  offsets <- seq(0, min(total - 1, 9499), by = PAGE_SIZE)
  map(offsets, function(offset) {
    cache_file <- if (!is.null(cache_prefix))
      file.path(CACHE_DIR, paste0(cache_prefix, "_off", offset, ".json"))
    else NULL
    goat_page(query_str, offset, result_type = result_type, cache_file = cache_file)
  }) %>%
    unlist(recursive = FALSE)
}

# Get child taxon IDs at a given rank under a parent taxon
get_child_ids <- function(parent_taxid, child_rank, cache_prefix) {
  query <- paste0("tax_tree(", parent_taxid, ") AND tax_rank(", child_rank, ")")
  hits  <- goat_fetch_all(query, result_type = "taxon", cache_prefix = cache_prefix)
  map_chr(hits, ~ .x$result$taxon_id %||% NA_character_) %>%
    na.omit() %>%
    as.character()
}

# =============================================================================
# 2.  SHARDED SPECIES FETCHER
# Recursively shards oversized taxa:
#   order > 9500 spp  -> shard by family
#   family > 9500 spp -> shard by genus
# =============================================================================

fetch_species_for_taxon <- function(taxid, label, class_name, rank = "order") {

  query <- paste0("tax_tree(", taxid, ") AND tax_rank(species)")
  total <- goat_total(query, "taxon")

  if (total == 0) return(list())

  if (total <= 9500) {
    return(goat_fetch_all(query,
                          result_type  = "taxon",
                          cache_prefix = paste0(class_name, "_", rank, taxid)))
  }

  # ── Shard by family ─────────────────────────────────────────────────────────
  message("  ", label, " (", total, " spp) — sharding by family")
  family_ids <- get_child_ids(taxid, "family",
                              cache_prefix = paste0(class_name, "_fam_of_", taxid))

  if (length(family_ids) == 0) {
    warning("No families found for taxid ", taxid, " — truncating at 9500")
    return(goat_fetch_all(query,
                          result_type  = "taxon",
                          cache_prefix = paste0(class_name, "_", rank, taxid)))
  }

  # Fetch all species under each family
  family_results <- map(family_ids, function(fid) {
    fquery <- paste0("tax_tree(", fid, ") AND tax_rank(species)")
    ftotal <- goat_total(fquery, "taxon")

    if (ftotal == 0) return(list())

    if (ftotal <= 9500) {
      return(goat_fetch_all(fquery,
                            result_type  = "taxon",
                            cache_prefix = paste0(class_name, "_fam", fid)))
    }

    # ── Shard by genus ───────────────────────────────────────────────────────
    message("    family ", fid, " (", ftotal, " spp) — sharding by genus")
    genus_ids <- get_child_ids(fid, "genus",
                               cache_prefix = paste0(class_name, "_gen_of_", fid))

    if (length(genus_ids) == 0) {
      warning("No genera found for family taxid ", fid, " — truncating at 9500")
      return(goat_fetch_all(fquery,
                            result_type  = "taxon",
                            cache_prefix = paste0(class_name, "_fam", fid)))
    }

    map(genus_ids, function(gid) {
      gquery <- paste0("tax_tree(", gid, ") AND tax_rank(species)")
      goat_fetch_all(gquery,
                     result_type  = "taxon",
                     cache_prefix = paste0(class_name, "_genus", gid))
    }) %>%
      unlist(recursive = FALSE)

  }) %>%
    unlist(recursive = FALSE)

  # ── Catch species with no family (incertae sedis / unclassified) ─────────────
  # These sit directly under the order and are missed by family sharding
  family_taxon_ids <- map_chr(family_results, ~ .x$result$taxon_id %||% NA_character_) %>%
    na.omit()

  direct_query <- paste0("tax_name(", taxid, ") AND tax_rank(species)")
  direct_total <- goat_total(direct_query, "taxon")

  direct_results <- if (direct_total > 0) {
    message("  ", label, ": fetching ", direct_total,
            " species directly under order (no family assigned)")
    goat_fetch_all(direct_query,
                   result_type  = "taxon",
                   cache_prefix = paste0(class_name, "_", rank, taxid, "_direct"))
  } else {
    list()
  }

  c(family_results, direct_results)
}
# =============================================================================
# 3.  TAXON PARSER  (result=taxon)
# With fields/names/ranks params the response contains:
#   r$fields : all requested fields
#   r$names  : keyed list — common_name, assembly_id, assembly_name
#   r$lineage: full ranks including domain
# =============================================================================

parse_taxon_result <- function(x) {
  r <- x$result

  lineage_rank <- function(rank) {
    hit <- keep(r$lineage %||% list(), ~ .x$taxon_rank == rank)
    if (length(hit) == 0) return(NA_character_)
    as.character(hit[[1]]$scientific_name %||% NA_character_)
  }

  # Always take first value — some fields (e.g. ebp_standard_criteria) return
  # a vector; tibble() requires length-1 scalars
  fval <- function(field, type = "real") {
    v <- r$fields[[field]][["value"]] %||% NA
    v <- if (length(v) > 1) v[[1]] else v
    switch(type, real = as.numeric(v), int = as.integer(v), as.character(v))
  }

  # r$names is a keyed list when &names= param is passed
  nval <- function(field) as.character(r$names[[field]] %||% NA_character_)

  # taxon_names still holds synonym and tolid prefix
  name_of_class <- function(name_class) {
    hit <- keep(r$taxon_names %||% list(), ~ .x[["class"]] == name_class)
    if (length(hit) == 0) return(NA_character_)
    map_chr(hit, ~ as.character(.x[["name"]] %||% NA_character_)) %>%
      na.omit() %>% paste(collapse = "; ")
  }

  tibble(
    # ── Identity ──────────────────────────────────────────────────────────────
    taxon_id        = as.character(r$taxon_id        %||% NA_character_),
    taxon_rank      = as.character(r$taxon_rank      %||% NA_character_),
    scientific_name = as.character(r$scientific_name %||% NA_character_),
    common_name     = nval("common_name"),
    synonym         = name_of_class("synonym"),
    tolid_prefix    = name_of_class("tolid prefix"),
    assembly_id     = nval("assembly_id"),
    assembly_name   = nval("assembly_name"),
    # ── Lineage ───────────────────────────────────────────────────────────────
    species      = lineage_rank("species"),
    genus        = lineage_rank("genus"),
    family       = lineage_rank("family"),
    order        = lineage_rank("order"),
    class        = lineage_rank("class"),
    phylum       = lineage_rank("phylum"),
    kingdom      = lineage_rank("kingdom"),
    superkingdom = lineage_rank("domain"),
    # ── Fields (all confirmed present via GoaT UI URL) ────────────────────────
    genome_size           = fval("genome_size",           "real"),
    genome_size_draft     = fval("genome_size_draft",     "real"),
    c_value               = fval("c_value",               "real"),
    chromosome_number     = fval("chromosome_number",     "int"),
    haploid_number        = fval("haploid_number",        "int"),
    assembly_level        = fval("assembly_level",        "chr"),
    assembly_span         = fval("assembly_span",         "real"),
    assembly_date         = fval("assembly_date",         "chr"),
    contig_n50            = fval("contig_n50",            "real"),
    scaffold_n50          = fval("scaffold_n50",          "real"),
    busco_completeness    = fval("busco_completeness",    "real"),
    gene_count            = fval("gene_count",            "int"),
    ebp_standard_criteria = fval("ebp_standard_criteria", "chr"),
    sequencing_status     = fval("sequencing_status",     "chr"),
    published             = fval("published",             "chr")
  )
}

fetch_taxon_class <- function(class_name, tax_id) {
  message("\n── Taxon: ", class_name, " (tax_id=", tax_id, ") ──────────────────")

  order_hits <- goat_fetch_all(
    paste0("tax_tree(", tax_id, ") AND tax_rank(order)"),
    result_type  = "taxon",
    cache_prefix = paste0(class_name, "_orders")
  )
  order_ids <- map_chr(order_hits, ~ .x$result$taxon_id %||% NA_character_) %>%
    na.omit() %>% as.character()

  message(class_name, ": ", length(order_ids), " orders — fetching species")

  all_results <- map(order_ids, function(oid) {
    fetch_species_for_taxon(oid,
                            label      = paste0("order ", oid),
                            class_name = class_name,
                            rank       = "order")
  }) %>%
    unlist(recursive = FALSE)

  message(class_name, ": ", length(all_results), " species records fetched")

  map_dfr(all_results, parse_taxon_result) %>%
    distinct(taxon_id, .keep_all = TRUE)
}

# =============================================================================
# 4.  ASSEMBLY PARSER  (result=assembly)
# Confirmed fields: assembly_level, assembly_span, contig_n50, scaffold_n50,
#                   chromosome_count, gc_percent, busco_completeness, gene_count
# Confirmed identifiers: genbank_accession, refseq_accession, wgs_accession
# =============================================================================

parse_assembly_result <- function(x) {
  r <- x$result

  fval <- function(field, type = "real") {
    v <- r$fields[[field]][["value"]] %||% NA
    v <- if (length(v) > 1) v[[1]] else v
    switch(type, real = as.numeric(v), int = as.integer(v), as.character(v))
  }

  ival <- function(id_class) {
    hit <- keep(r$identifiers %||% list(), ~ .x[["class"]] == id_class)
    if (length(hit) == 0) return(NA_character_)
    as.character(hit[[1]]$identifier %||% NA_character_)
  }

  tibble(
    taxon_id           = as.character(r$taxon_id    %||% NA_character_),
    assembly_id        = as.character(r$assembly_id %||% NA_character_),
    genbank_accession  = ival("genbank_accession"),
    refseq_accession   = ival("refseq_accession"),
    wgs_accession      = ival("wgs_accession"),
    assembly_level     = fval("assembly_level",     "chr"),
    assembly_span      = fval("assembly_span",      "real"),
    contig_n50         = fval("contig_n50",         "real"),
    scaffold_n50       = fval("scaffold_n50",       "real"),
    chromosome_count   = fval("chromosome_count",   "int"),
    gc_percent         = fval("gc_percent",         "real"),
    busco_completeness = fval("busco_completeness", "real"),
    gene_count         = fval("gene_count",         "int")
  )
}

fetch_assembly_class <- function(class_name, tax_id) {
  message("\n── Assembly: ", class_name, " (tax_id=", tax_id, ") ──────────────")

  order_hits <- goat_fetch_all(
    paste0("tax_tree(", tax_id, ") AND tax_rank(order)"),
    result_type  = "taxon",
    cache_prefix = paste0(class_name, "_orders")   # reuses taxon cache
  )
  order_ids <- map_chr(order_hits, ~ .x$result$taxon_id %||% NA_character_) %>%
    na.omit() %>% as.character()

  all_results <- map(order_ids, function(oid) {
    goat_fetch_all(
      paste0("tax_tree(", oid, ") AND tax_rank(species)"),
      result_type  = "assembly",
      cache_prefix = paste0(class_name, "_asm_order", oid)
    )
  }) %>%
    unlist(recursive = FALSE)

  message(class_name, ": ", length(all_results), " assembly records fetched")
  if (length(all_results) == 0) return(tibble())

  map_dfr(all_results, parse_assembly_result)
}

# =============================================================================
# 5.  ASSEMBLY DEDUPLICATION
# Keep best assembly per taxon_id:
#   chromosome > scaffold > contig, then best busco_completeness
# =============================================================================

LEVEL_RANK <- c(
  "chromosome" = 1L, "complete genome" = 1L,
  "scaffold"   = 2L, "contig"          = 3L
)

best_assembly_per_species <- function(asm_df) {
  asm_df %>%
    mutate(level_rank = LEVEL_RANK[tolower(assembly_level)] %>% replace_na(4L)) %>%
    arrange(taxon_id, level_rank, desc(busco_completeness)) %>%
    distinct(taxon_id, .keep_all = TRUE) %>%
    select(-level_rank)
}

