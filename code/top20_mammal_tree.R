##
library(readxl)
library(tidyverse)
library(cowplot)
library(treedataverse)
library(magrittr)
library(RColorBrewer)

## plotting options
theme_set(theme_bw())
here::here()
path <- here::here()


summary_species <-
  readRDS(paste0(path, "/output/Species_List_with_genome_and_cl_info.rds")) %>%
  filter(class == "Mammalia") %>%
  mutate(selected = if_else(cat == "Both", "selected", "FALSE")) %>%
  mutate(
    orders_num = as.integer(as.factor(order)),
    alternative_name_ = gsub(x = alternative_name, replacement = "_", " "),
    annotation = !(is.na(Release.Date) & is.na(Ensembl.Assembly))
  )

best20 <- summary_species %>%
  filter(selected == "selected") %>%
  group_by(order, genus) %>%
  slice_max(scaffold_n50, n = 1) %>%
  group_by(order) %>%
  slice_max(fibro, n = 5)

summary_species <- summary_species %>%
  mutate(top20 = if_else(
    scientific_name_ %in% best20$scientific_name_,
    "selected",
    "FALSE"
  ))

table(summary_species$top20, summary_species$order)

table(summary_species$selected, summary_species$order)

big_tree <-
  drop.tip(read.nexus(
    paste0(
      path,
      "/data/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre"
    )
  ), "_Anolis_carolinensis")

tree_meta <- data.frame(label = big_tree$tip.label) %>%
  separate(
    col = label,
    sep = "_",
    into = c("s1", "s2", "family", "order"),
    remove = F
  ) %>%
  unite(c(s1, s2), col = "scientific_name_", sep = "_") %>%
  mutate(family = str_to_title(family),
         order = str_to_title(order))

table(summary_species$scientific_name_ %in% tree_meta$scientific_name_)

summary_species <- summary_species %>%
  mutate(
    scientific_name_ = if_else(
      alternative_name_ %in% tree_meta$scientific_name_,
      alternative_name_,
      scientific_name_
    )
  )

table(summary_species$scientific_name_ %in% tree_meta$scientific_name_)

summary_species[!(summary_species$scientific_name_ %in% tree_meta$scientific_name_), ]

# replace incorrect species names in tree
Wrong_Species_names2 <- read_excel("data/Wrong_Species_names2.xlsx")

tree_meta[tree_meta$scientific_name_ %in% Wrong_Species_names2$mammal_tree_name, ]$scientific_name_ <-
  Wrong_Species_names2$species_name_

# finlly all names match
table(summary_species$scientific_name_ %in% tree_meta$scientific_name_)
#big_tree$tip.label<-paste(tree_meta$scientific_name_,toupper(tree_meta$family),toupper(tree_meta$order),sep="_")

remove <-
  !(tree_meta$scientific_name_ %in% summary_species$scientific_name_)

table(remove)

big_tree <- drop.tip(big_tree, tree_meta$label[remove])

ggtree(big_tree)

keep <-
  (tree_meta$scientific_name_ %in% summary_species$scientific_name_)

tree_meta <- tree_meta[keep, ]

# tree_meta <- data.frame(label=big_tree$tip.label) %>%
#   separate(col = label,sep = "_",into = c("s1","s2","family","order"),remove = F) %>%
#   unite(c(s1,s2),col = "scientific_name_",sep = "_") %>%
#   mutate(family=str_to_title(family),
#          order=str_to_title(order))
table(big_tree$tip.label == tree_meta$label)

big_tree$tip.label <- tree_meta$scientific_name_

big_tree$tip.label[!(big_tree$tip.label %in% summary_species$scientific_name_)]

summary_species <-
  summary_species[na.omit(match(big_tree$tip.label, summary_species$scientific_name_)), ]

pos <- which(colnames(summary_species) == "scientific_name_")

summary_species <-
  summary_species[, c(pos, 1:(pos - 1), (pos + 1):ncol(summary_species))]

p.tree <- ggtree(big_tree, color = "grey60", lwd = 0.5)

cols <-
  colorRampPalette(brewer.pal("RdPu", n = 9)[-c(1:3)])(length(unique(summary_species$order)))

p.tree2 <- p.tree %<+% summary_species +
  geom_tippoint(aes(color = order), size = 2) +
  #guides(color="none")+
  ggsci::scale_color_futurama()
p.tree2

mat_sel2 <- summary_species %>%
  ungroup() %>%
  filter(!(duplicated(scientific_name_))) %>%
  mutate(fibro = fibro_l,
         replicates = samples > 1) %>%
  dplyr::select(scientific_name_, fibro, genome, replicates, top20) %>%
  column_to_rownames("scientific_name_") %>%
  as.matrix

heatmap.colours <- c("purple", "darkblue", "grey75")

names(heatmap.colours) <- c("selected", "TRUE", "FALSE")

tree_map2 <-
  gheatmap(
    p.tree2,
    mat_sel2,
    width = .2,
    colnames = T,
    legend_title = "",
    colnames_position = "top",
    font.size = 4,
    colnames_angle = 90,
    colnames_offset_y = 7,
    offset = 2
  )   +
  scale_fill_manual(values = heatmap.colours,
                    #breaks=c("Chromosome","Contig","Scaffold",2,4,6,8,10),
                    na.value = "grey75",
                    limits = force) +
  ylim(NA, 110) +
  labs(colour = "",
       fill = "",
       shape = "") +
  theme(legend.position = c(0.4, 0.9)) +
  theme(legend.direction = "horizontal") +
  guides(fill = "none")
tree_map2

