# Phylogenetic turnover thresholds in Atlantic Forest birds and anurans

**Associated manuscript:**  
Anunciação, P.R., Barros, F.M., Ribeiro, M.C. *et al.* Phylogenetic Shifts in Anuran and Avian Communities Driven by Land-Use Change in the Atlantic Forest. *Journal of Applied Ecology* (under review).

\---

## Overview

This repository contains all R code used to analyse phylogenetic community turnover thresholds in birds and anurans across land-use gradients in the Brazilian Atlantic Forest. The analysis uses GradientForest (Ellis et al. 2012) as the primary framework, with Generalised Dissimilarity Modelling (Ferrier et al. 2007) as a supplementary robustness check at the species-level.

\---

## Repository structure

```
README.md
01\_vif\_analysis.R
02\_gf\_analysis\_combined.R
03\_anuran\_tree\_matching.R
04\_bird\_tree\_matching.R
05\_gdm\_phylogenetic\_analysis.R
data/
  site\_PG\_birds\_new.txt
  site\_PG\_Anura.txt
  sites\_envir\_birds\_coord.txt
  sites\_envir\_anura\_coord.txt
  sites\_species\_anura\_names.txt
  sites\_species\_birds\_names\_new.txt
```

### Scripts

|Script|Description|
|-|-|
|`01\_vif\_analysis.R`|Collinearity diagnostics — reduces predictors from 9 to 7 per taxon using `vifcor` (Pearson \|r\| > 0.7) and `vifstep` (VIF > 5)|
|`02\_gf\_analysis\_combined.R`|Main GradientForest analysis for both taxa — ntree stability, full model, 100-bootstrap threshold uncertainty and variable importance, Figure 2 and Figure 3|
|`03\_anuran\_tree\_matching.R`|Anuran taxonomy matching and tree pruning — synonymises species names, prunes Jetz \& Pyron (2018) consensus tree to 36 study species, makes ultrametric via chronos|
|`04\_bird\_tree\_matching.R`|Bird taxonomy matching and tree pruning — synonymises species names, excludes 15 species absent from BirdTree, merges 7 congeneric pairs, prunes Jetz et al. (2012) Hackett backbone to 145 tips|
|`05\_gdm\_phylogenetic\_analysis.R`|Supplementary GDM analysis — computes PhyloSor and betaMNTD, fits GDM models with geo = FALSE, permutation importance tests, Mantel tests, I-spline figures|

### Data files

|File|Description|
|-|-|
|`site\_PG\_birds\_new.txt`|Bird phylogenetic group presence/absence matrix (32 sites × 130 PGs)|
|`site\_PG\_Anura.txt`|Anuran phylogenetic group presence/absence matrix (15 sites × 19 PGs)|
|`sites\_envir\_birds\_coord.txt`|Environmental predictors for bird sampling units (32 sites × 7 predictors + coordinates)|
|`sites\_envir\_anura\_coord.txt`|Environmental predictors for anuran sampling units (15 sites × 7 predictors + coordinates)|
|`sites\_species\_anura\_names.txt`|Anuran species-level presence/absence matrix (15 sites × 37 species; used for GDM only)|
|`sites\_species\_birds\_names\_new.txt`|Bird species-level presence/absence matrix (32 sites × 166 species; used for GDM only)|

**Land-use classification data** (input to `01\_vif\_analysis.R`) are available from the corresponding author upon reasonable request (p.ribeiroanunciacao@lancaster.ac.uk).

**Phylogenetic trees** required by scripts 03 and 04 must be downloaded separately:

* Anurans: Jetz \& Pyron (2018) consensus tree (`amph\_shl\_new\_Consensus\_7238.tre`) from Dryad: https://doi.org/10.5061/dryad.db265
* Birds: BirdTree Hackett backbone posterior trees (`hacket\_vertlife\_seq.nex`) from https://birdtree.org (download 1000 trees using the Hackett backbone)

\---

## Reproducibility

All stochastic analyses use `set.seed(42)`. To reproduce the analysis from scratch:

1. Set `rerun\_analysis <- TRUE` in `02\_gf\_analysis\_combined.R`
2. Set `rerun\_distances <- TRUE` in `05\_gdm\_phylogenetic\_analysis.R`

To regenerate figures only from saved outputs (much faster):

1. Set `rerun\_analysis <- FALSE` in `02\_gf\_analysis\_combined.R` (requires `gf\_model\_birds` and `gf\_model\_anurans` in the R environment from a saved session)
2. Set `rerun\_distances <- FALSE` in `05\_gdm\_phylogenetic\_analysis.R`

**Expected runtime** (full analysis from scratch):

* Script 01: < 1 minute
* Script 02: approximately 2–3 hours (bootstrap replicates)
* Scripts 03–04: < 5 minutes each
* Script 05: approximately 30–60 minutes (permutation tests)

\---

## Dependencies

All analyses were conducted in **R version 4.3.1** (R Core Team 2023).

|Package|Version|Use|
|-|-|-|
|`gradientForest`|0.1-18|GradientForest analysis|
|`gdm`|1.5|Generalised Dissimilarity Modelling|
|`picante`|1.8.2|PhyloSor, betaMNTD, Faith's PD|
|`ape`|5.7|Phylogenetic tree manipulation, chronos|
|`vegan`|2.6-4|Bray-Curtis dissimilarity, Mantel test|
|`usdm`|2.1-7|VIF collinearity diagnostics|
|`readr`|2.1.4|CSV reading|
|`corrplot`|0.92|Correlation matrix visualisation|

Install all dependencies:

```r
install.packages(c("gdm", "picante", "ape", "vegan",
                   "usdm", "readr", "corrplot"))

# gradientForest is not on CRAN — install from R-Forge:
install.packages("gradientForest",
                 repos = "http://R-Forge.R-project.org")
```

\---

## Data sources

**Species occurrence data:**

* Anuran occurrence data: Anunciação et al. (2021) *Biological Conservation* https://doi.org/10.1016/j.biocon.2021.109137
* Bird occurrence data: Barros et al. (2019) *Landscape Ecology* https://doi.org/10.1007/s10980-019-00812-z

**Phylogenies:**

* Birds: Jetz, W., Thomas, G.H., Joy, J.B., Hartmann, K. \& Mooers, A.O. (2012) The global diversity of birds in space and time. *Nature* 491, 444–448. https://doi.org/10.1038/nature11631
* Anurans: Jetz, W. \& Pyron, R.A. (2018) The interplay of past diversification and evolutionary isolation with present imperilment across the amphibian tree of life. *Nature Ecology \& Evolution* 2, 850–858. https://doi.org/10.1038/s41559-018-0515-5

\---

## Contact

Paula Ribeiro Anunciação  
Marie Skłodowska-Curie Postdoctoral Fellow  
Lancaster Environment Centre, Lancaster University, UK  
p.ribeiroanunciacao@lancaster.ac.uk

\---

## Licence

Code is released under the MIT Licence. See `LICENSE` for details.

