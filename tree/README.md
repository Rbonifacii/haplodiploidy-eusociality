# Phylogenetic trees

This folder contains the phylogenetic trees associated with Bonifacii et al. (2026), *No evidence that haplodiploidy favors the evolution of eusociality*.

Two trees are provided:

- **`chestersklado_pruned.tree`** — the final uncalibrated phylogeny, including Allodapini, used as the starting point for divergence time estimation.
- **`chesters_pruned_dated.tre`** — the final time-calibrated phylogeny used in the main analyses, with branch lengths expressed in millions of years.

### Divergence time estimation

Divergence times were estimated using secondary calibrations from the dated family-level insect phylogeny of Rainford et al. (2014). Shared nodes between the phylogenies were identified and calibrations transferred using `congruify` in the R package `geiger`. These calibrations were then applied using penalized likelihood dating in `treePL`, with cross-validation used to select the optimal smoothing parameter.

The resulting time-calibrated phylogeny was used for the main analyses. Additional analyses using Grafen-scaled branch lengths were conducted to assess the robustness of the results and are reported in SI Appendix, Tables S3 and S4(b).

### Repository update

The time-calibrated tree and updated uncalibrated tree (including Allodapini) were added to this repository in October 2026. These are the original trees used in preparing the published analyses, rather than newly generated phylogenies.

The `archive/` folder contains an earlier version of the uncalibrated phylogeny, prior to the addition of Allodapini and divergence-time calibration during peer review.
