library(phytools)   # For stochastic mapping and transition counting
library(ape)        # For phylogenetic tree manipulation (e.g. drop.tip)
library(stringr)    # For splitting tip label strings
library(dplyr)      # For tidy manipulation and rowwise evaluation


# ============================
# 📁 FUNCTIONS: Read & Nee Contrast Identification
# ============================

# This script defines two core functions used to identify phylogenetically 
# independent contrasts following the Read & Nee (1995) method. These contrasts 
# allow comparative tests of trait evolution that avoid pseudoreplication caused 
# by shared ancestry.

# FUNCTIONS:

# 1. find_pure_clades(tree, node, haplo_vec)
#    ---------------------------------------
#    Recursively traverses a phylogenetic tree to identify all "pure" 
#    haplodiploid clades — that is, monophyletic groups where all descendant 
#    species have haplodiploid01 == 1. Returns a dataframe of such clades, 
#    including single-species clades (singletons).

# 2. find_best_sister_clade(tree, focal_node, haplo_vec, used_nodes = ..., used_tips = ...)
#    -------------------------------------------------------------------------------------
#    For each pure haplodiploid clade, finds the best-matching sister clade 
#    that contains only non-haplodiploid species (haplodiploid01 == 0), 
#    has not been previously used, and is phylogenetically proximate. 
#    Returns a ranked list of eligible sister clades for each focal clade.

# These functions form the foundation for building independent evolutionary 
# comparisons in downstream Read & Nee-style analyses.

# Required packages: ape, phangorn, dplyr, purrr, tibble, stringr

# === Recursive function to find all pure haplodiploid clades ===
find_pure_clades <- function(tree, node, haplo_vec) {
  desc_tips <- tree$tip.label[Descendants(tree, node, type = "tips")[[1]]]
  trait_vals <- haplo_vec[desc_tips]
  clean_vals <- trait_vals[!is.na(trait_vals)]
  
  if (length(clean_vals) == 0) return(NULL)
  
  # === Singleton
  if (length(desc_tips) == 1 && trait_vals == 1) {
    tip_index <- which(tree$tip.label == desc_tips)
    parent_node <- Ancestors(tree, tip_index, type = "parent")
    
    return(tibble(
      node = tip_index,        # used for contrast logic
      clade_node = parent_node, # used for node depth & plotting
      n_species = 1,
      tips = desc_tips
    ))
  }
  
  # === All haplodiploid
  if (length(desc_tips) >= 2 && all(clean_vals == 1)) {
    return(tibble(
      node = node,
      clade_node = node,
      n_species = length(desc_tips),
      tips = paste(desc_tips, collapse = ", ")
    ))
  }
  
  # === Two-taxon, one haplodiploid
  if (length(desc_tips) == 2 && sum(trait_vals == 1, na.rm = TRUE) == 1) {
    # Identify which tip is haplodiploid
    haplo_tip <- desc_tips[which(trait_vals == 1)]
    haplo_tip_index <- which(tree$tip.label == haplo_tip)
    
    # Get MRCA of both tips
    clade_mrca <- getMRCA(tree, desc_tips)
    
    # Find the two children of the MRCA
    children <- Children(tree, clade_mrca)
    
    # Find which child contains the haplodiploid tip
    matching_child <- NULL
    for (child in children) {
      desc <- Descendants(tree, child, type = "tips")[[1]]
      if (haplo_tip_index %in% desc) {
        matching_child <- child
        break
      }
    }
    
    # Return that child as the 'node' (for contrast finding),
    # and the MRCA as the 'clade_node' (for depth sorting)
    return(tibble(
      node = matching_child,
      clade_node = clade_mrca,
      n_species = 1,
      tips = haplo_tip
    ))
  }
  
  
  # === Recurse
  children <- Children(tree, node)
  if (length(children) == 0) return(NULL)
  
  results <- lapply(children, function(child) find_pure_clades(tree, child, haplo_vec))
  return(bind_rows(results))
}


# === Step 2: Sister-finding function ===
find_best_sister_clade <- function(tree, focal_node, haplo_vec, used_nodes = integer(), used_tips = character()) {
  # Initialize climbing variables
  level_up <- 0
  current_node <- focal_node
  
  # Helper to gather all internal and tip nodes under a given node
  get_all_internal_nodes <- function(tree, node) {
    children <- Children(tree, node)
    if (length(children) == 0) {
      return(node)
    }
    return(c(node, unlist(lapply(children, get_all_internal_nodes, tree = tree))))
  }
  
  repeat {
    # Move up to the parent
    parent <- Ancestors(tree, current_node, type = "parent")
    if (length(parent) == 0) {
      # Reached root without finding a suitable sister
      return(NULL)
    }
    
    # Identify all sibling branches
    siblings <- setdiff(Children(tree, parent), current_node)
    if (length(siblings) == 0) {
      current_node <- parent
      level_up <- level_up + 1
      next
    }
    
    # Collect every candidate node under each sibling (including tips)
    candidate_nodes <- unique(unlist(lapply(siblings, function(s) get_all_internal_nodes(tree, s))))
    
    # Evaluate each candidate for eligibility
    all_candidates <- list()
    for (node in candidate_nodes) {
      subtree_nodes <- Descendants(tree, node, type = "all")[[1]]
      if (any(subtree_nodes %in% used_nodes)) next
      
      # Extract tip labels for this candidate
      node_tips <- tree$tip.label[Descendants(tree, node, type = "tips")[[1]]]
      node_tips <- trimws(node_tips)
      clean_tips <- setdiff(node_tips, used_tips)
      if (length(clean_tips) == 0) next
      
      # Check that all are non-haplodiploid
      trait_vals <- haplo_vec[clean_tips]
      if (all(trait_vals == 0, na.rm = TRUE)) {
        all_candidates[[length(all_candidates) + 1]] <- tibble(
          focal_node         = focal_node,
          parent_node        = parent,
          sister_nodes       = as.character(node),
          sister_n_tips      = length(clean_tips),
          sister_has_non_haplo = TRUE,
          sister_tips        = paste(clean_tips, collapse = ", "),
          levels_up          = level_up
        )
      }
    }
    
    # If we found any candidates, rank and return the best
    if (length(all_candidates) > 0) {
      # Determine focal clade size
      focal_size <- length(Descendants(tree, focal_node, type = "tips")[[1]])
      
      result_df <- bind_rows(all_candidates) %>%
        mutate(
          size_difference = abs(sister_n_tips - focal_size)
        ) %>%
        arrange(levels_up, size_difference) %>%
        slice(1)
      
      return(result_df)
    }
    
    # Otherwise, climb up and try again
    current_node <- parent
    level_up <- level_up + 1
  }
}


# Expanded find_best_sister_clade_expanded (randomized robustness: return top N candidates)
find_best_sister_clade_expanded <- function(tree, focal_node, haplo_vec,
                                            used_nodes = integer(), used_tips = character(),
                                            min_required_candidates = 10, max_levels_up = 10) {
  level_up <- 0
  current_node <- focal_node
  all_candidates <- list()
  focal_size <- length(Descendants(tree, focal_node, type = "tips")[[1]])
  
  get_all_internal_nodes <- function(tree, node) {
    children <- Children(tree, node)
    if (length(children) == 0) return(node)
    return(c(node, unlist(lapply(children, get_all_internal_nodes, tree = tree))))
  }
  
  repeat {
    if (level_up > max_levels_up) break
    parent <- Ancestors(tree, current_node, type = "parent")
    if (length(parent) == 0) break
    
    siblings <- setdiff(Children(tree, parent), current_node)
    if (length(siblings) == 0) {
      current_node <- parent
      level_up <- level_up + 1
      next
    }
    
    # collect all internal and tips under siblings
    candidate_nodes <- unique(unlist(lapply(siblings, function(s) get_all_internal_nodes(tree, s))))
    for (node in candidate_nodes) {
      subtree_all <- Descendants(tree, node, type = "all")[[1]]
      if (any(subtree_all %in% used_nodes)) next
      
      node_tips <- tree$tip.label[Descendants(tree, node, type = "tips")[[1]]]
      clean_tips <- setdiff(trimws(node_tips), used_tips)
      if (length(clean_tips) == 0) next
      
      trait_vals <- haplo_vec[clean_tips]
      if (all(trait_vals == 0, na.rm = TRUE)) {
        all_candidates[[length(all_candidates) + 1]] <- tibble(
          focal_node          = focal_node,
          parent_node         = parent,
          sister_nodes        = as.character(node),
          sister_n_tips       = length(clean_tips),
          sister_has_non_haplo = TRUE,
          sister_tips         = paste(clean_tips, collapse = ", "),
          levels_up           = level_up,
          size_diff           = abs(length(clean_tips) - focal_size)
        )
      }
    }
    
    # stop when we have enough candidates
    if (length(all_candidates) >= min_required_candidates) break
    current_node <- parent
    level_up <- level_up + 1
  }
  
  if (length(all_candidates) == 0) return(NULL)
  
  # return all sorted by levels_up and then size diff
  return(bind_rows(all_candidates) %>% arrange(levels_up, size_diff))
}
