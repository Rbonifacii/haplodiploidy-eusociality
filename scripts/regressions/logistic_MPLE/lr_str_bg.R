# ===============================
# Setup

library(ape)
library(phytools)
library(phylolm)
library(dplyr)
library(here)

# ===============================
# Load tree and data

tree <- read.tree(here("tree", "chestersklado_pruned.tree"))

insect_data <- read.csv(here("data", "datasetS1.csv")) %>%
  select(species_name, haplodiploidy_strict, superorganism_bg) %>%
  filter(!is.na(haplodiploidy_strict), !is.na(superorganism_bg))

# ===============================
# Clean and align data with tree

# Keep only species present in both tree and data
species_to_keep <- intersect(tree$tip.label, insect_data$species_name)
tree <- drop.tip(tree, setdiff(tree$tip.label, species_to_keep))

# Order data to match tree tips
insect_data <- insect_data[match(tree$tip.label, insect_data$species_name), ]
stopifnot(all(tree$tip.label == insect_data$species_name))

# Convert columns to numeric
insect_data <- insect_data %>%
  mutate(
    haplodiploidy_strict = as.numeric(as.character(haplodiploidy_strict)),
    superorganism_bg = as.numeric(as.character(superorganism_bg))
  )

# Set rownames to species and drop species_name column
rownames(insect_data) <- insect_data$species_name
insect_data$species_name <- NULL

# ===============================
# Phylogenetic Logistic Regression

# Model 1: with haplodiploidy
model1 <- phyloglm(
  superorganism_bg ~ haplodiploidy_strict,
  data = insect_data,
  phy = tree,
  method = "logistic_MPLE",
  boot = 1000
)

summary(model1)

# Fit a null model (no haplodiploidy)
model_null <- phyloglm(
  superorganism_bg ~ 1,
  data = insect_data,
  phy = tree,
  method = "logistic_MPLE",
  boot = 1000
)

summary(model_null)


# ===============================
# Save model
saveRDS(model1, file = "output/logistic_str_bg.rds")

# ============================
# End of Script
# ============================