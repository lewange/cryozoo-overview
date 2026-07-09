# Analysis Scripts — Cryozoo Species Selection

## Overview

The scripts follow a linear pipeline: data ingestion and cleaning → inventory updates → summary plots → phylogenetic trees. Run them in the order shown below.

---

## Script descriptions

### 1. `Clean_CryoZoo_Lists.Rmd`
**Starting point of the pipeline.**

Reads the master *List of Species* Excel file and the Google Sheet inventory, tidies taxonomy, resolves misspelled species names, fetches IUCN status and genome/annotation metadata (GOAT/NCBI), and joins everything into two clean summary tables:

- `output/Species_List_with_genome_and_cl_info.rds` — one row per species; includes taxonomy, IUCN category, genome assembly info, TOGA coverage, cell-line status, and number of individuals.
- `output/Sample_List_all_species_clean.rds` — one row per sample/individual; includes tissue of origin, sex, age, and provenance.

Both tables are also written back to the *Cryozoo_clean* Google Sheet.

**Key inputs:**
- `data/List of species.xlsx`
- `data/Wrong_Species_names.xlsx`
- Google Sheets (species inventory)

**R packages:** `tidyverse`, `readxl`, `googlesheets4`, `taxize`, `janitor`, `runner`, `udpipe`, `cowplot`

---

### 2. `update_Cryostock.Rmd`
**Updates the physical cryotube inventory.**

Reads all box sheets from the N₂/−80 °C storage Google Sheet, classifies each box and tube (fibroblast, biopsy, PBMC, iPSC, extracts, other), parses animal IDs and passage numbers, and produces:

- `data/Cryo_stocks_long.rds` — one row per cryotube position.
- `data/Cryo_stocks_summary.rds` — per-animal, per-cell-type summary (total tubes, max passage, number of passages).
- `data/Extractions_long.rds` — inventory of nucleic-acid extraction boxes.

Results are also pushed back to the *Cryo_stocks* Google Sheet.

**Key inputs:** Google Sheet (N₂/−80 inventory)

**R packages:** `tidyverse`, `googlesheets4`

---

### 3. `update_Cell_culture.Rmd`
**Summarises cell culture / thawing outcomes.**

Pulls thawing records from two Google Sheets (Illumina and Silencer project tracking), cleans them, and saves a combined harvest summary:

- `output/cell_culture_summary.rds` — one row per thawing event, with thawing date, passage, cell number, and harvest outcome.

Results are also written back to a consolidated Google Sheet.

**Key inputs:** Google Sheets (cell culture tracking)

**R packages:** `tidyverse`, `googlesheets4`, `purrr`, `janitor`

---

### 4. `update_Extractions.Rmd`
**Summarises nucleic-acid extraction records.**

Reads RNA and DNA extraction tables from a Google Sheet and the extraction location inventory (`data/Extractions_long.rds` produced by `update_Cryostock.Rmd`). Currently used interactively to inspect extraction yields per animal.

**Key inputs:**
- Google Sheet (RNA/DNA extractions)
- `data/Extractions_long.rds` ← `update_Cryostock.Rmd`

**R packages:** `tidyverse`, `googlesheets4`

---

### 5. `update_Sequencing.Rmd`
**Aggregates sequencing status across all library types.**

Pulls per-library status from individual Google Sheets for RNA-seq, ATAC-seq, WGS (Illumina), EM-seq, and WG-Nanopore, merges them into a single table, and pushes it to the *All_Sequenced* tab of the sequencing tracker. Also creates per-species and per-sample wide-format summaries.

Outputs:
- `output/sequencing_summary.rds`

**Key inputs:**
- `output/Sample_List_all_species_clean.rds` ← `Clean_CryoZoo_Lists.Rmd`
- `output/Species_List_with_genome_and_cl_info.rds` ← `Clean_CryoZoo_Lists.Rmd`
- Multiple sequencing Google Sheets

**R packages:** `tidyverse`, `googlesheets4`

---

### 6. `Update_Sample_table.Rmd`
**Joins all inventory data into a master sample table.**

Loads the clean species and sample tables, then left-joins cryotube counts, sequencing status, and cell culture outcomes into a single comprehensive view (displayed with `view()`; not saved to disk yet).

**Key inputs (all must be produced first):**
- `output/Species_List_with_genome_and_cl_info.rds` ← `Clean_CryoZoo_Lists.Rmd`
- `output/Sample_List_all_species_clean.rds` ← `Clean_CryoZoo_Lists.Rmd`
- `data/Cryo_stocks_long.rds` ← `update_Cryostock.Rmd`
- `data/Cryo_stocks_summary.rds` ← `update_Cryostock.Rmd`
- `output/cell_culture_summary.rds` ← `update_Cell_culture.Rmd`
- Google Sheet (sequencing tracker, *All_Sequenced* tab) ← `update_Sequencing.Rmd`

**R packages:** `tidyverse`, `readxl`, `googlesheets4`, `cowplot`

---

### 7. `plot_Cryozoo_summary.Rmd`
**Main summary figures for the collection.**

Reads `Species_List_with_genome_and_cl_info.rds` (and the sample table) and produces bar charts and beeswarm plots covering: species counts by taxonomic class, IUCN status breakdown, fibroblast line establishment, genome availability and assembly quality, TOGA annotation overlap, sample provenance by institution, and a selection of best-genome species written to CSVs.

Also saves an updated species table as:
- `output/Species_List_with_genome_and_cl_info2.rds` (adds a `cat` column used downstream by the tree scripts)

**Key inputs:**
- `output/Species_List_with_genome_and_cl_info.rds` ← `Clean_CryoZoo_Lists.Rmd`
- `output/Sample_List_all_species_clean.rds` ← `Clean_CryoZoo_Lists.Rmd`
- `data/mammal_data_orders.rds`

**Key outputs:** figures in `output/figures/`; `output/Species_List_with_genome_and_cl_info2.rds`; `selection_best_genomes.csv`; `selected_species.csv`

**R packages:** `tidyverse`, `cowplot`, `ggsci`, `ggbeeswarm`, `scales`

---

### 8. `plot_Tree_mammal.Rmd`
**Phylogenetic tree for mammals.**

Uses the MamPhy supertree (vertlife) to place Cryozoo mammal species on a time-calibrated phylogeny. Corrects species name mismatches, subsets the tree to one representative per family (filling in Cryozoo species), and produces circular and rectangular tree layouts with a heatmap annotation (fibroblast status, genome availability, replicates, sequencing status).

Saves the pruned tree object and metadata for reuse:
- `data/mamaltree_object_allfamilies.Rds`
- `data/mamaltree_object_allfamilies_metadata.Rds`
- `data/mamaltree_object_allfamilies_cols.Rds`

**Key inputs:**
- `output/Species_List_with_genome_and_cl_info2.rds` ← `plot_Cryozoo_summary.Rmd`
- `data/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre`
- `data/Wrong_Species_names2.xlsx`
- Google Sheet (species selection: silencer / illumina flags)

**R packages:** `tidyverse`, `cowplot`, `treedataverse` (ggtree), `magrittr`, `RColorBrewer`, `readxl`, `googlesheets4`, `ggsci`

---

### 9. `plot_Tree_birds.Rmd`
**Phylogenetic tree for birds.**

Mirrors `plot_Tree_mammal.Rmd` for birds. Reads a NEXUS tree from vertlife.org (selects one of 100 posterior trees), maps GOAT taxonomic metadata onto the tips, corrects species names, and generates circular and rectangular tree visualisations with the same heatmap annotation approach.

**Key inputs:**
- `output/Species_List_with_genome_and_cl_info2.rds` ← `plot_Cryozoo_summary.Rmd`
- `data/all_birds_tree.nex`
- `data/Wrong_Species_names2.xlsx`
- `data/GenomeOnATree_Vertebrata1.tsv`

**R packages:** `tidyverse`, `cowplot`, `treedataverse` (ggtree), `magrittr`, `RColorBrewer`, `readxl`, `ggsci`

---

## Dependency diagram

```
Clean_CryoZoo_Lists.Rmd
    │
    ├──► update_Cryostock.Rmd ──────────────────────────────────────┐
    │                                                                 │
    ├──► update_Cell_culture.Rmd ───────────────────────────────────►│
    │                                                                 │
    ├──► update_Extractions.Rmd (also uses Cryo_stocks output)       │
    │                                                                 │
    ├──► update_Sequencing.Rmd ─────────────────────────────────────►│
    │                                                                 │
    │                               Update_Sample_table.Rmd ◄────────┘
    │
    └──► plot_Cryozoo_summary.Rmd
              │
              ├──► plot_Tree_mammal.Rmd
              └──► plot_Tree_birds.Rmd
```

## Archived scripts

Older or exploratory versions of the analysis are kept in `analysis/archive/` and are not part of the active pipeline.
