# Simple phylogenetic tree plots for CryoZoo species
# One circular tree per common_group, with IUCN status as outer ring
# Run interactively in RStudio — not via Rscript

library(tidyverse)
library(rotl)
library(ape)
library(ggtree)
library(ggtreeExtra)

# ── Paths ──────────────────────────────────────────────────────────────────────
path     <- here::here()
out_dir  <- file.path(path, "output", "figures", "trees")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── IUCN colours ──────────────────────────────────────────────────────────────
IUCN_COLOURS <- c(
  LC  = "#60C659",
  NT  = "#CCE226",
  VU  = "#F9E814",
  EN  = "#FC7F3F",
  CR  = "#D81E05",
  EW  = "#542344",
  EX  = "#000000",
  DD  = "#D1D1C6",
  NE  = "#AAAAAA"
)

# Map verbose IUCN categories to short codes
iucn_recode <- c(
  "LEAST CONCERN"         = "LC",
  "NEAR THREATENED"       = "NT",
  "VULNERABLE"            = "VU",
  "ENDANGERED"            = "EN",
  "CRITICALLY ENDANGERED" = "CR",
  "EXTINCT IN THE WILD"   = "EW",
  "EXTINCT"               = "EX",
  "DATA DEFICIENT"        = "DD"
)

# ── 1. Load & tidy species table ──────────────────────────────────────────────
species_list <- readRDS(file.path(path, "output", "Species_List_with_genome_and_cl_info.rds")) %>%
  mutate(
    common_group = case_when(
      class %in% c("Actinopteri", "Chondrichthyes") ~ "Fish",
      class %in% c("Lepidosauria")                  ~ "Reptilia",
      is.na(class)                                  ~ "Reptilia",
      TRUE                                           ~ class
    ),
    iucn_code = recode(iucn_category, !!!iucn_recode, .default = "NE"),
    iucn_code = factor(iucn_code, levels = names(IUCN_COLOURS)),
    # Fix known name mismatches for OTL
    species_TOL = case_match(
      species,
      "Astur gentilis"    ~ "Accipiter gentilis",
      "Caribicus warreni" ~ "Celestus warreni",
      .default = species
    )
  )

cat("Species per common_group:\n")
print(count(species_list, common_group))

# ── 2. Fetch OTL tree for a set of species ────────────────────────────────────
fetch_tree <- function(species_vec, group_name) {
  message("Fetching OTL tree for ", group_name, " (", length(species_vec), " spp)...")

  matched <- tryCatch(
    tnrs_match_names(species_vec, context_name = NULL),
    error = function(e) { message("  OTL match failed: ", e$message); NULL }
  )
  if (is.null(matched)) return(NULL)

  ott_ids <- matched %>%
    filter(!is.na(ott_id), !flags %in% c("barren", "extinct")) %>%
    pull(ott_id) %>%
    unique()

  if (length(ott_ids) < 3) {
    message("  Too few resolved names (", length(ott_ids), "), skipping")
    return(NULL)
  }

  tree <- tryCatch(
    tol_induced_subtree(ott_ids = ott_ids, label_format = "name"),
    error = function(e) { message("  Subtree fetch failed: ", e$message); NULL }
  )
  if (is.null(tree)) return(NULL)

  tree$tip.label <- tree$tip.label %>%
    str_replace_all("_", " ") %>%
    str_remove("\\s+ott\\d+$")

  tree <- compute.brlen(tree, method = "Grafen")
  tree
}

# ── 3. Build one circular tree plot with IUCN ring ────────────────────────────
build_tree_plot <- function(tree, meta_df, group_name) {
  message("Building plot for ", group_name, "...")

  tip_df <- tibble(species_TOL = tree$tip.label) %>%
    left_join(meta_df, by = "species_TOL") %>%
    mutate(iucn_code = factor(iucn_code, levels = names(IUCN_COLOURS)))

  p <- ggtree(tree, layout = "circular", size = 0.4, colour = "grey40") %<+% tip_df +
    geom_tiplab(aes(label = label), size = 2, offset = 0.05,
                fontface = "italic", colour = "grey20") +
    labs(title = group_name) +
    theme(
      plot.title      = element_text(face = "bold", size = 13, hjust = 0.5),
      legend.position = "right",
      legend.title    = element_text(size = 9, face = "bold"),
      legend.text     = element_text(size = 8)
    )

  # IUCN ring
  iucn_ring <- tip_df %>% select(species_TOL, iucn_code)

  p <- p +
    geom_fruit(
      data    = iucn_ring,
      geom    = geom_tile,
      mapping = aes(y = species_TOL, fill = iucn_code),
      width   = 0.06, offset = 0.12,
      colour  = "white"
    ) +
    scale_fill_manual(
      values   = IUCN_COLOURS,
      na.value = "grey90",
      name     = "IUCN Status",
      drop     = FALSE
    )

  p
}

# ── 4. Run for each common_group ──────────────────────────────────────────────
groups <- unique(species_list$common_group)

for (grp in groups) {
  sp_vec <- species_list %>%
    filter(common_group == grp) %>%
    pull(species_TOL)

  tree <- fetch_tree(sp_vec, grp)
  if (is.null(tree)) {
    message("Skipping ", grp, " — no tree returned")
    next
  }

  p <- build_tree_plot(tree, species_list, grp)

  # Size plot by number of species
  n   <- length(tree$tip.label)
  dim <- max(8, round(sqrt(n) * 1.5))

  out_file <- file.path(out_dir, paste0("tree_", grp, ".pdf"))
  ggsave(out_file, p, width = dim, height = dim, units = "in", device = cairo_pdf)
  message("Saved: ", out_file)
}

message("\nDone! All trees saved to: ", out_dir)
