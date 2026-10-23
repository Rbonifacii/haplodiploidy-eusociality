# Phylogenetic trees

This folder contains the phylogenetic trees used in Bonifacii et al. (2026), *No evidence that haplodiploidy favors the evolution of eusociality*.

Two trees are provided:

- chestersklado_pruned.tree — final uncalibrated phylogeny, including Allodapini, used for divergence time estimation.
- chesters_pruned_dated.tre — final time-calibrated phylogeny used in the main analyses, with branch lengths in millions of years.

### Divergence time estimation

Divergence times were estimated using secondary calibrations from the dated family-level insect phylogeny of Rainford et al. (2014). Shared nodes were identified and calibrations transferred using `congruify` in the R package `geiger`, followed by penalized likelihood dating in `treePL`. Cross-validation was used to select the smoothing parameter.

Analyses using Grafen-scaled branch lengths are reported in SI Appendix, Tables S3 and S4(b).

### Repository update

These trees were used in preparing the published analyses and subsequently added to this repository in October 2026.

The `archive/` folder contains an earlier version of the uncalibrated phylogeny, before Allodapini was added and divergence-time calibration was performed during peer review.
