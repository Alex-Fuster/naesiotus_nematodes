library(tidyverse)
library(readxl)
library(ggpubr)
library(phytools)
library(ape)
library(ggtree)

tree <- read.newick(here::here("data/phylogeny_2024/tree2024.nwk"))
df_rads <- read_excel(here::here("data/phylogeny_2024/New_tree_taxa_28Aug24.xlsx"))
species_names_df_raw <- read.csv(here::here("data/species_names_checks.csv"))

species_names_df <- species_names_df_raw %>% 
  rename(species = annotated_name)  %>% 
  filter(species_name != "")

tree_data <- fortify(tree)
ggtree(tree) +
  geom_tiplab() +  # Add tip labels for species
  theme_tree() +  # Remove legend
  xlim(0, max(tree_data$x) + 0.3) 


## Eliminate accents

tree$tip.label <- gsub("'", "", tree$tip.label)
tree$tip.label[grepl("sp._nov._Volcán_Wolf_1", tree$tip.label)] <- "sp. nov. Volcan Wolf 1"
tree$tip.label[grepl("sp._nov._Volcán_Wolf_2", tree$tip.label)] <- "sp. nov. Volcan Wolf 2"
tree$tip.label[grepl("sp._nov._Volcán_Wolf_3", tree$tip.label)] <- "sp. nov. Volcan Wolf 3"
tree$tip.label[which(tree$tip.label == "sp._nov._pinzón_1")]<- "sp. nov. pinzon 1"

tree_data <- fortify(tree)
ggtree(tree) +
  geom_tiplab() +  # Add tip labels for species
  theme_tree() +  # Remove legend
  xlim(0, max(tree_data$x) + 0.3) 



# Tree formatting: delete duplicate names and homogenize names format

# Identify the positions of the duplicates in the tip labels
amastroides_tips <- which(tree$tip.label == "amastroides")
perspectivus_tips <- which(tree$tip.label == "perspectivus")

# Append "delete" to one instance of each duplicated name
tree$tip.label[amastroides_tips[1]] <- "amastroidesdelete"
tree$tip.label[perspectivus_tips[1]] <- "perspectivusdelete"

# Drop the tips with "delete" in their names
tree <- drop.tip(tree, c("amastroidesdelete", "perspectivusdelete"))

# Replace underscores with hyphens (or other suitable character)
tree$tip.label <- gsub("_", "-", tree$tip.label)

# Replace periods with empty space or other suitable character
tree$tip.label <- gsub("\\.", "", tree$tip.label)

# Replace hyphens with an empty string or another suitable character
tree$tip.label <- gsub("-", "", tree$tip.label)

plot(tree)

write.tree(tree, file = here::here("BAMM/cleannames_tree_doublenames_041024.txt"))
write.tree(tree, file = here::here("BAMM/cleannames_tree_doublenames_041024.newick"))