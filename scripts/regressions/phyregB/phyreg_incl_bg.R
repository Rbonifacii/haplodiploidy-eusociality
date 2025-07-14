# ============================
# Phylogenetic Regression Script
# ============================

# ==== Package Setup ====
# Install and load required packages
# install.packages("renv")  # Uncomment if renv is not installed
renv::restore()  # Restore package environment for reproducibility

# Install phyregB manually (if needed)
# install.packages("/Users/user/Downloads/phyregB_0.1.0.tar.gz", repos = NULL, type = "source")

library(phyregB)
library(here)
library(ape)
library(phytools)
library(dplyr)

here::i_am("scripts/regressions/phyregB/phyreg_incl_mw.R")

# ==== Load Data ====
# Load cleaned dataset
insect_data <- read.csv(here("data", "datasetS1.csv")) %>%
  select(species_name, superorganism_bg, haplodiploidy_inclusive) %>%
  filter(!is.na(haplodiploidy_inclusive), !is.na(superorganism_bg))  # Remove rows with missing values

# Load and prune tree to match species in data
tree <- read.tree(here("tree", "chestersklado_pruned.tree"))
species_to_keep <- intersect(tree$tip.label, insect_data$species_name)
tree <- drop.tip(tree, setdiff(tree$tip.label, species_to_keep))

# Convert tree to internal format for phyregB
tree_vector <- phyfromphylo(tree)

# ==== Ensure Data & Tree Alignment ====
# Reorder data to match the tree tip labels
insect_data <- insect_data[match(tree$tip.label, insect_data$species_name), ]

# Confirm alignment
stopifnot(all(tree$tip.label == insect_data$species_name))

# ==== Run Phylogenetic Regression ====
phylo_model <- phyreg(
  control = superorganism_bg ~ 1,                  # Null model (intercept only)
  test = ~ haplodiploidy_inclusive,           # Test predictor
  data = insect_data,                         # Trait data
  phydata = tree_vector$phy                   # Phylogeny (internal format)
)

# ==== Save Model Output ====
saveRDS(phylo_model, file = here("output", "phyreg_incl_bg.rds"))
