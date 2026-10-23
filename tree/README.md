# Phylogenetic trees

This folder contains the phylogenetic trees associated with Bonifacii et al. (2026), *No evidence that haplodiploidy favors the evolution of eusociality*.

Two trees are provided:

- **`chestersklado_pruned.tree`** — original, uncalibrated phylogeny used as the starting point for divergence time estimation.
- **`chesters_pruned_dated.tre`** — time-calibrated phylogeny used in the main analyses, with branch lengths expressed in millions of years.

Divergence times were estimated using secondary calibrations transferred from the dated family-level insect phylogeny of Rainford et al. (2014), using congruify in the R package geiger, followed by penalized likelihood dating using treePL. using `congruify` in the R package `geiger`, followed by penalized likelihood dating using `treePL`.

The time-calibrated tree was added to this repository in October 2026 to make the phylogeny used in the published analyses directly available.

Analyses using alternative Grafen-scaled branch lengths are reported in SI Appendix, Tables S3 and S4(b).
