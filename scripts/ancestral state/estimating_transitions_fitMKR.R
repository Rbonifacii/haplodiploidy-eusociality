# ------------------------------------------------------------------------------
# This script performs ancestral state reconstruction for all trait combinations usinf phytools:: fitMK
# The aim is to estimate transition counts between character states across the tree.
# ------------------------------------------------------------------------------

library(ape)
library(phytools)
library(dplyr)
library(here)

# --- Load phylogenetic tree and trait data ---
tree_master <- read.tree(here("tree", "chestersklado_pruned.tree"))
df <- read.csv(here("data", "datasetS1.csv"), stringsAsFactors = FALSE)

# --- Define state combinations and interaction traits ---
possible_states <- c("0.0", "1.0", "0.1", "1.1")

interaction_defs <- list(
  ploidyIncl_eusocial_mw      = c("haplodiploidy_inclusive", "eusocial_mw"),
  ploidyIncl_eusocial_cy      = c("haplodiploidy_inclusive", "eusocial_cy"),
  ploidyIncl_superorganism_bg = c("haplodiploidy_inclusive", "superorganism_bg"),
  ploidyStr_eusocial_mw       = c("haplodiploidy_strict", "eusocial_mw"),
  ploidyStr_eusocial_cy       = c("haplodiploidy_strict", "eusocial_cy"),
  ploidyStr_superorganism_bg  = c("haplodiploidy_strict", "superorganism_bg")
)

# --- Construct interaction state variables ---
for (iname in names(interaction_defs)) {
  vars <- interaction_defs[[iname]]
  df[[iname]] <- factor(
    interaction(df[[vars[1]]], df[[vars[2]]], drop = TRUE),
    levels = possible_states
  )
}

# --- Initialise results table ---
results <- data.frame(
  trait_col              = character(),
  haplo_eusocial_gains   = integer(),
  diploid_eusocial_gains = integer(),
  haplodiploidy_gains    = integer(),
  total_transitions      = integer(),
  stringsAsFactors       = FALSE
)

model_outputs <- list()

# --- Main loop over all trait combinations ---
for (trait_col in names(interaction_defs)) {
  cat("\n\n=== Trait:", trait_col, "===\n")
  
  trait <- df[[trait_col]]
  names(trait) <- df$species_name
  
  valid_tips <- intersect(tree_master$tip.label, names(trait)[!is.na(trait)])
  tree2 <- drop.tip(tree_master, setdiff(tree_master$tip.label, valid_tips))
  
  trait <- trait[tree2$tip.label]
  cat("Tip state distribution:\n")
  print(table(trait))
  
  # --- Ancestral state reconstruction (ARD model) ---
  stopifnot(all(levels(trait) == possible_states))
  states_vec <- trait
  names(states_vec) <- tree2$tip.label
  
  fit <- fitMk(tree2, states_vec, model = "ARD", root.prior = c(1, 0, 0, 0))
  lik_anc <- ancr(fit)$ace
  
  node_states <- apply(lik_anc, 1, function(x) colnames(lik_anc)[which.max(x)])
  names(node_states) <- rownames(lik_anc)
  
  all_states <- c(as.character(trait), node_states)
  names(all_states) <- c(tree2$tip.label, rownames(lik_anc))
  
  model_outputs[[trait_col]] <- list(
    fit      = fit,
    lik_anc  = lik_anc,
    states   = trait,          # observed states at tips
    node_states = node_states  # most probable states at internal nodes
  )
  
  # --- Map states to tree edges ---
  Ntip <- length(tree2$tip.label)
  Nnode <- tree2$Nnode
  all_labels <- c(tree2$tip.label, as.character((Ntip + 1):(Ntip + Nnode)))
  names(all_labels) <- as.character(1:(Ntip + Nnode))
  
  parent_states <- all_states[all_labels[as.character(tree2$edge[, 1])]]
  child_states  <- all_states[all_labels[as.character(tree2$edge[, 2])]]
  
  # --- Count transitions between states ---
  transitions <- table(
    factor(parent_states, levels = possible_states),
    factor(child_states, levels = possible_states)
  )
  diag(transitions) <- 0
  
  haplo_eusocial_gains   <- transitions["1.0", "1.1"]
  diploid_eusocial_gains <- transitions["0.0", "0.1"]
  haplodiploidy_gains    <- sum(transitions[c("0.0", "0.1"), c("1.0", "1.1")])
  total_transitions      <- sum(transitions)
  
  cat("Haplodiploid gains of eusociality: ", haplo_eusocial_gains, "\n")
  cat("Diploid gains of eusociality:      ", diploid_eusocial_gains, "\n")
  cat("Gains of haplodiploidy:            ", haplodiploidy_gains, "\n")
  cat("Total transitions:                 ", total_transitions, "\n")
  
  results <- rbind(
    results,
    data.frame(
      trait_col              = trait_col,
      haplo_eusocial_gains   = haplo_eusocial_gains,
      diploid_eusocial_gains = diploid_eusocial_gains,
      haplodiploidy_gains    = haplodiploidy_gains,
      total_transitions      = total_transitions,
      stringsAsFactors       = FALSE
    )
  )
}

# --- Output summary results ---
cat("\n\n=== Summary Table ===\n")
print(results)

# --- Save results as .rds file in the output folder ---
saveRDS(
  list(
    summary_table = results,
    models        = model_outputs
  ),
  file = here("output", "transition_models_and_summary.rds")
)

# ============================
# End of Script
# ============================