# -----------------------------------------
# Script for Robust Phylogenetic Regression Using ROBRT
# -----------------------------------------

# ==== Load Libraries ====

library(ROBRT)
library(geiger)
library(phytools)
library(ape)
library(dplyr)
library(here)

# ==== Load and Prepare Tree & Data ====
here::i_am("scripts/regressions/robust_regression/ROBRT_str_mw.R")

tree <- read.tree(here("tree", "chestersklado_pruned.tree"))

insect_data <- read.csv(here("data", "datasetS1.csv")) %>%
  select(species_name, haplodiploidy_strict, eusocial_mw) %>%
  filter( !is.na(haplodiploidy_strict),
          !is.na(eusocial_mw)
 )

# ==== Match Tree and Data ====
species_to_keep <- intersect(tree$tip.label, insect_data$species_name)
tree <- drop.tip(tree, setdiff(tree$tip.label, species_to_keep))

insect_data <- insect_data[match(tree$tip.label, insect_data$species_name), ]
rownames(insect_data) <- insect_data$species_name
insect_data$species_name <- NULL

# ==== Run Regression ====
estimators <- "L1"

robust_results <- Conduct.Robust_PhylogeneticRegression_PGLS(
  tree,
  insect_data$eusocial_mw,
  insect_data$haplodiploidy_strict,
  estimators
)

summary(robust_results$L1)

# ==== Save Model Output ====
saveRDS(robust_results, here("output", "robust_regression_str_mw.rds"))

# ==== Diagnostic Plots ====
fit_L1 <- robust_results$L1
res <- fit_L1$residuals
fitted_vals <- fit_L1$fitted.values

if (all(is.finite(res)) && all(is.finite(fitted_vals))) {
  par(mfrow = c(2, 2))  # 2x2 layout for plots
  
  plot(fitted_vals, res, main = "Residuals vs Fitted", xlab = "Fitted", ylab = "Residuals")
  abline(h = 0, col = "red")
  
  qqnorm(res)
  qqline(res, col = "red")
  
  hist(res, main = "Histogram of Residuals", xlab = "Residuals")
  
  plot(fitted_vals, abs(res), main = "Scale-Location", xlab = "Fitted", ylab = "|Residuals|")
} else {
  cat("Residuals or fitted values contain NA/NaN/Inf — skipping diagnostic plots.\n")
}
