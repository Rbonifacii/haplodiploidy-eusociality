# ============================
# HMM Correlated Evolution Script
# ============================

# ==== Package Setup ====
# Load required libraries
library(corHMM)   
library(ape)      
library(dplyr)   
library(here)     
library(phytools)  

set.seed(2025)
here::i_am("scripts/HMM/HMIM_str_mw.R")

# ==== Load and Prepare Data ====

# Load phylogeny and apply branch lengths
tree <- read.tree(here("tree","chestersklado_pruned.tree"))

# Load and clean trait data
trait_data <- read.csv(here("data", "datasetS1.csv"), stringsAsFactors = FALSE) %>%
  select(sp = species_name, X = haplodiploidy_strict, S = eusocial_mw) %>%
  mutate(across(c(X, S), as.character)) %>%
  filter(!is.na(as.numeric(X)), !is.na(as.numeric(S))) %>%
  mutate(across(c(X, S), as.numeric))

# Match species between tree and data
common_species <- intersect(tree$tip.label, trait_data$sp)
tree <- drop.tip(tree, setdiff(tree$tip.label, common_species))
trait_data <- trait_data[trait_data$sp %in% common_species, ]
trait_data <- trait_data[match(tree$tip.label, trait_data$sp), ]

# ==== Define Model Structures ====

# Dependent (Pagel's) model
dep_mat <- getStateMat4Dat(trait_data, collapse = FALSE, indep = FALSE)$rate.mat

# Independent model
indep_mat <- getStateMat4Dat(trait_data, collapse = FALSE, indep = TRUE)$rate.mat

# Hidden Markov Independent Model (HMIM)
hmim_state_mats <- list(indep_mat, indep_mat)
rate_cat_mat <- equateStateMatPars(getRateCatMat(2), 1:2)
HMIM_mat <- getFullMat(hmim_state_mats, rate_cat_mat)

# Hidden Markov Dependent Model (HMIM Dependent)
hmim_dep_state_mats <- list(dep_mat, dep_mat)
HMIM_dep_mat <- getFullMat(hmim_dep_state_mats, rate_cat_mat)

# ==== Fit Models ====

fit_dep <- corHMM(
  phy = tree,
  data = trait_data,
  rate.cat = 1,
  rate.mat = dep_mat,
  model = "ARD",
  collapse = TRUE,
  node.states = "none",
  nstarts = 0,
  root.p = c(1, 0, 0, 0)
)

fit_ind <- corHMM(
  phy = tree,
  data = trait_data,
  rate.cat = 1,
  rate.mat = indep_mat,
  model = "ARD",
  collapse = TRUE,
  node.states = "none",
  nstarts = 0,
  root.p = c(1, 0, 0, 0)
)

fit_HMIM <- corHMM(
  phy = tree,
  data = trait_data,
  rate.cat = 2,
  rate.mat = HMIM_mat,
  model = "ARD",
  collapse = TRUE,
  node.states = "none",
  nstarts = 0,
  root.p = c(1, 0, 0, 0)
)

fit_HMIM_dep <- corHMM(
  phy = tree,
  data = trait_data,
  rate.cat = 2,
  rate.mat = HMIM_dep_mat,
  model = "ARD",
  collapse = TRUE,
  node.states = "none",
  nstarts = 0,
  root.p = c(1, 0, 0, 0)
)

# ==== Compare Model Fits ====

aic_table <- tibble::tibble(
  model = c("Dependent", "Independent", "HMIM", "HMIM Dependent"),
  logLik = c(fit_dep$loglik, fit_ind$loglik, fit_HMIM$loglik, fit_HMIM_dep$loglik),
  AIC = c(fit_dep$AIC, fit_ind$AIC, fit_HMIM$AIC, fit_HMIM_dep$AIC),
  delta_AIC = AIC - min(AIC),
  Akaike_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))
)

print(aic_table)

# Save table
readr::write_csv(aic_table, here("scripts", "analyseschesters17", "output", "HMIM_AIC_table_str_mw.csv"))

# ==== Visualize Model Structures ====

plotMKmodel(HMIM_mat, rate.cat = 2, display = "row", text.scale = 0.7)
plotMKmodel(HMIM_dep_mat, rate.cat = 2, display = "row", text.scale = 0.7)

# ==== Stochastic Character Mapping (HMIM Dependent) ====

model_mat <- fit_HMIM_dep$solution
model_mat[is.na(model_mat)] <- 0
diag(model_mat) <- -rowSums(model_mat)

simmap <- makeSimmap(
  tree = fit_HMIM_dep$phy,
  data = fit_HMIM_dep$data,
  model = model_mat,
  rate.cat = 2,
  nSim = 1,
  nCores = 1
)

phytools::plotSimmap(simmap[[1]], fsize = 0.1)

# ==== Save Model Fits ====

saveRDS(
  list(
    fit_dep = fit_dep,
    fit_ind = fit_ind,
    fit_HMIM = fit_HMIM,
    fit_MIM_dep = fit_HMIM_dep
  ),
  file = here("scripts", "analyseschesters17", "output", "HMIM_model_fits_str_mw.rds")
)

# ============================
# End of Script
# ============================
