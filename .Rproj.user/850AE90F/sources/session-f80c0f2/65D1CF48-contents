# ==============================================================================
# 📜 SCRIPT PURPOSE
# This script performs Read & Nee-style phylogenetic contrasts to test whether
# haplodiploid clades show higher prevalence of eusociality than their non-haplodiploid
# sister clades. The analysis includes:
# - Identifying pure haplodiploid clades
# - Finding appropriate sister clades
# - Statistical testing (binomial, Wilcoxon)
# - Robustness analysis using random sister selection
# - This script uses the "inclusive" definition of Haplodiploidy, and Michener/Wilsons definition of eusociality
# ==============================================================================

# --- Libraries ---
library(ape)
library(phangorn)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(here)
library(tidyr)
library(ggplot2)
library(coin)
source(here("functions", "readnee_functions.R"))

# --- Load data ---
insect_data <- read.csv(here("data", "datasetS1"))
tree <- read.tree(here("tree", "chestersklado_pruned_to_data.tree"))
haplo_vec <- setNames(insect_data$haplodiploidy_inclusive, insect_data$species_name)[tree$tip.label]

# --- Identify pure haplodiploid clades ---
root_node <- Ntip(tree) + 1
pure_clades <- find_pure_clades(tree, root_node, haplo_vec)

# --- Add clade depth and sort ---
node_depths <- node.depth(tree)
names(node_depths) <- as.character(1:(Ntip(tree) + tree$Nnode))
pure_clades$node_depth <- node_depths[as.character(pure_clades$clade_node)]
pure_clades <- arrange(pure_clades, node_depth)

# --- Find best sister clades for each pure clade ---
used_nodes <- integer()
used_tips <- character()
final_contrasts <- list()

for (i in seq_len(nrow(pure_clades))) {
  clade_row <- pure_clades[i, ]
  focal_node <- clade_row$node
  focal_n_species <- clade_row$n_species
  
  sister <- find_best_sister_clade(tree, focal_node, haplo_vec, used_nodes, used_tips)
  
  if (!is.null(sister)) {
    best_sister <- arrange(sister, levels_up, abs(sister_n_tips - focal_n_species)) %>% slice(1)
    
    final_contrasts[[length(final_contrasts) + 1]] <- tibble(
      node = clade_row$node,
      n_species = clade_row$n_species,
      tips = clade_row$tips,
      parent_node = best_sister$parent_node,
      sister_nodes = best_sister$sister_nodes,
      sister_n_tips = best_sister$sister_n_tips,
      sister_tips = best_sister$sister_tips,
      levels_up = best_sister$levels_up
    )
    
    used_nodes <- union(used_nodes, Descendants(tree, as.integer(best_sister$sister_nodes), type = "all")[[1]])
    used_tips <- union(used_tips, str_split(best_sister$sister_tips, ",\\s*")[[1]])
  }
}

contrasts_df <- bind_rows(final_contrasts)
write.csv(contrasts_df, here("output", "read_nee_contrasts_incl.csv"), row.names = FALSE)

# --- Statistical Tests ---

# Prepare data for analysis
contrast_results <- contrasts_df %>%
  rowwise() %>%
  mutate(
    haplo_tips = str_split(tips, ",\\s*"),
    sister_tips = str_split(sister_tips, ",\\s*"),
    haplo_eusocial = sum(insect_data$species_name %in% unlist(haplo_tips) & insect_data$eusocial_mw == "1"),
    haplo_total = length(unlist(haplo_tips)),
    sister_eusocial = sum(insect_data$species_name %in% unlist(sister_tips) & insect_data$eusocial_mw == "1"),
    sister_total = length(unlist(sister_tips)),
    haplo_prop = haplo_eusocial / haplo_total,
    sister_prop = sister_eusocial / sister_total,
    informative = haplo_prop != sister_prop,
    outcome = case_when(
      haplo_prop > sister_prop ~ 1,
      haplo_prop < sister_prop ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup()

# Binomial test
informative_contrasts <- filter(contrast_results, informative)
binom_test <- binom.test(sum(informative_contrasts$outcome == 1),
                         nrow(informative_contrasts),
                         p = 0.5, alternative = "greater")
print(binom_test)

# Wilcoxon signed-rank test (permutation-based)
wilcox_data <- contrast_results %>%
  filter(!is.na(haplo_prop), !is.na(sister_prop)) %>%
  mutate(pair_id = row_number()) %>%
  pivot_longer(cols = c(haplo_prop, sister_prop), names_to = "group", values_to = "prop") %>%
  mutate(group = factor(group), pair_id = factor(pair_id))

balanced_ids <- wilcox_data %>% count(pair_id) %>% filter(n == 2) %>% pull(pair_id)
wilcox_clean <- filter(wilcox_data, pair_id %in% balanced_ids)

wilcox_test <- wilcoxsign_test(prop ~ group | pair_id, data = wilcox_clean,
                               distribution = approximate(nresample = 10000),
                               alternative = "greater")
print(wilcox_test)

write.csv(contrast_results, here("output", "read_nee_results_incl_mw.csv"), row.names = FALSE)

# ==============================================================================

# --- Robustness Analysis ---

set.seed(1365)

n_reps <- 100
replicate_results <- vector("list", n_reps)

for (rep in 1:n_reps) {
  used_nodes <- integer()
  used_tips <- character()
  final_contrasts <- list()
  
  for (i in seq_len(nrow(pure_clades))) {
    clade_row <- pure_clades[i, ]
    focal_node <- clade_row$node
    focal_n_species <- clade_row$n_species
    
    sisters <- find_best_sister_clade_expanded(tree, focal_node, haplo_vec, used_nodes, used_tips)
    if (!is.null(sisters)) {
      top_candidates <- arrange(sisters, levels_up, abs(sister_n_tips - focal_n_species)) %>% slice_head(n = 10)
      chosen <- slice_sample(top_candidates, n = 1)
      
      final_contrasts[[length(final_contrasts) + 1]] <- tibble(
        rep = rep,
        focal_node = clade_row$node,
        n_species = clade_row$n_species,
        tips = clade_row$tips,
        parent_node = chosen$parent_node,
        sister_nodes = chosen$sister_nodes,
        sister_n_tips = chosen$sister_n_tips,
        sister_tips = chosen$sister_tips,
        levels_up = chosen$levels_up
      )
      
      used_nodes <- union(used_nodes, Descendants(tree, as.integer(chosen$sister_nodes), type = "all")[[1]])
      used_tips <- union(used_tips, str_split(chosen$sister_tips, ",\\s*")[[1]])
    }
  }
  replicate_results[[rep]] <- bind_rows(final_contrasts)
}

# Combine all replicate results
eusocial_results <- bind_rows(replicate_results) %>%
  rowwise() %>%
  mutate(
    haplo_tips = str_split(tips, ",\\s*"),
    sister_tips = str_split(sister_tips, ",\\s*"),
    haplo_eusocial = sum(insect_data$species_name %in% unlist(haplo_tips) & insect_data$eusocial_mw == 1),
    haplo_total = length(unlist(haplo_tips)),
    sister_eusocial = sum(insect_data$species_name %in% unlist(sister_tips) & insect_data$eusocial_mw == 1),
    sister_total = length(unlist(sister_tips)),
    haplo_prop = ifelse(haplo_total > 0, haplo_eusocial / haplo_total, NA_real_),
    sister_prop = ifelse(sister_total > 0, sister_eusocial / sister_total, NA_real_),
    informative = haplo_prop != sister_prop,
    outcome = case_when(
      haplo_prop > sister_prop ~ 1,
      haplo_prop < sister_prop ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup()

write.csv(eusocial_results, here("output", "read_nee_robustness_incl_mw.csv"), row.names = FALSE)

# Summarise replicate outcomes
replicate_summary <- eusocial_results %>%
  filter(informative) %>%
  group_by(rep) %>%
  summarise(
    n = n(),
    haplo_better = sum(outcome == 1, na.rm = TRUE),
    sister_better = sum(outcome == 0, na.rm = TRUE),
    ties = sum(is.na(outcome)),
    prop_haplo_wins = mean(outcome == 1, na.rm = TRUE),
    majority = case_when(
      haplo_better > sister_better ~ 1,
      haplo_better < sister_better ~ 0,
      TRUE ~ NA_real_
    ), .groups = "drop"
  )

summary_across_reps <- summarise(replicate_summary,
                                 total_reps = n(),
                                 majority_haplo = sum(majority == 1, na.rm = TRUE),
                                 majority_sister = sum(majority == 0, na.rm = TRUE),
                                 ties = sum(is.na(majority)),
                                 prop_reps_haplo_wins = mean(majority == 1, na.rm = TRUE),
                                 avg_prop_haplo_wins_within_rep = mean(prop_haplo_wins, na.rm = TRUE)
)

print(summary_across_reps)
