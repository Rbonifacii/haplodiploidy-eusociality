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
  select(species_name, haplodiploidy_inclusive, eusocial_cy) %>%
  filter(!is.na(haplodiploidy_inclusive), !is.na(eusocial_cy))

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
    haplodiploidy_inclusive = as.numeric(as.character(haplodiploidy_inclusive)),
    eusocial_cy = as.numeric(as.character(eusocial_cy))
  )

# Set rownames to species and drop species_name column
rownames(insect_data) <- insect_data$species_name
insect_data$species_name <- NULL

# ===============================
# Phylogenetic Logistic Regression

# Model 1: with haplodiploidy
model1 <- phyloglm(
  eusocial_cy ~ haplodiploidy_inclusive,
  data = insect_data,
  phy = tree,
  method = "logistic_MPLE",
  boot = 1000
)

summary(model1)

# Fit a null model (no haplodiploidy)
model_null <- phyloglm(
  eusocial_cy ~ 1,
  data = insect_data,
  phy = tree,
  method = "logistic_MPLE",
  boot = 1000
)

summary(model_null)


# ===============================
# Save model
saveRDS(model1, file = "output/logistic_incl_cy.rds")

# ============================
# End of Script
# ============================
