# Benthic communities in Southwestern Atlantic oceanic islands
https://doi.org/10.5281/zenodo.21890794

Data and R scripts for analyzing benthic communities across Southwestern Atlantic oceanic islands


[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


## 📖 About

Benthic communities are essential components of reef ecosystems, contributing to their biodiversity, functioning, and structural complexity. This repository contains the data, analytical scripts, and outputs used to investigate benthic community structure across four Brazilian oceanic islands:

- **St. Peter and St. Paul Archipelago** (SP)
- **Rocas Atoll** (RA)
- **Fernando de Noronha Archipelago** (FN)
- **Trindade Island** (TR)

## 🎯 Objective
This study aims to:
- Quantify spatial and temporal variation in benthic community structure across four Southwestern Atlantic oceanic islands (SP, RA, FN, TR).
- Identify environmental variables associated with community composition (wave power, SST, PAR, POC) using multivariate analyses (PERMANOVA, PCoA, dbRDA).
- Assess temporal trends (2013–2019) in dominant benthic functional groups using Mann–Kendall tests.
- Identify benthic categories associated with each island or combinations of islands.

## 📊 Methods

Shallow reefs were sampled using photo-quadrats collected along semi-fixed transects. Benthic cover was identified at the image level and aggregated at the transect level, with transects used as the sampling units for community-level analyses.

### Community composition and ordination

PERMANOVA, PCoA, envfit, and indicator analyses were performed using **18 benthic categories**, comprising **17 biological categories plus abiotic substrate (ABI)**.

The analyses included:

- **PERMANOVA** (`vegan::adonis2`) to test differences among islands, years, and their interaction.
- **PCoA** for visualization of benthic community differences among islands.
- **envfit** to evaluate associations between individual benthic categories and the ordination.
- **Indicator analysis** (`indicspecies`) to identify benthic categories associated with individual islands or combinations of islands.

### Environmental drivers

Distance-based Redundancy Analysis (**dbRDA**) was used to evaluate relationships between benthic community composition and environmental variables (wave power, SST, PAR, and POC).

For the dbRDA, benthic cover was summarized into **10 broader benthic groups**:

- Epilithic Algal Matrix (EAM)
- Macroalgae (MAL)
- Crustose Coralline Algae (CCA)
- Articulated Coralline Algae (ACA)
- Scleractinian corals (SCL)
- Abiotic substrate (ABI)
- Invertebrates (INV)
- Zoanthids (ZOA)
- Cyanobacteria (CYA)
- Suspension feeders (SUS)

The dbRDA was based on mean annual benthic cover for each island × year combination.

### Temporal trends
Temporal trends were assessed using:
- LOESS smoothing for visualization.
- Mann–Kendall trend tests applied to observed (non-standardized) annual mean cover values derived from transect-level aggregation, ensuring independence among observations.
- Sen’s slope estimator to quantify trend magnitude.

This approach ensures that trends reflect real ecological changes rather than artifacts of data standardization.


## Main Results

The analyses revealed strong spatial differentiation in benthic community structure among the four oceanic islands, with island identity explaining substantially more variation than temporal differences.

The environmental analysis indicated distinct relationships between benthic community composition and environmental conditions among islands. Wave exposure was particularly associated with the benthic structure of Trindade, whereas other environmental variables, including POC, SST, and PAR, contributed to differences among islands.

Temporal trajectories were island-specific. Significant temporal changes were detected for some dominant groups, whereas other apparent trajectories were not statistically significant over the monitoring period.

These results emphasize that Brazilian oceanic islands should not be considered a single homogeneous reef system, but rather as distinct benthic communities shaped by local environmental conditions and temporal dynamics.

###  Final Remarks
As sentinels of the South Atlantic, these islands record not only the history of their impacts, but also the potential for their recovery. The unfolding of this story will depend on the continuity of monitoring and the effectiveness of recently implemented protection measures.

## 📁 Repository Structure
```text
Benthic_communities_in_Oceanic_islands/
├── data/
│   ├── raw/
│   │   ├── benthic_complete_data.csv
│   │   └── environment.csv
│   └── processed/
├── code/
│   ├── 01_benthic_composition.R
│   ├── 02_ordination.R
│   ├── 03_dbRDA.R
│   ├── 04_temporal_trends_SP_transect.R
│   ├── 05_temporal_trends_RA_transect.R
│   ├── 06_temporal_trends_FN_transect.R
│   ├── 07_temporal_trends_TR_transect.R
│   └── functions_transect.R
├── results/
│   ├── figures/
│   └── tables/
└── README.md
```
## 🚀 How to Reproduce
Set the working directory to the project root before running the scripts.

Example:
```r
setwd("path/to/Benthic_communities_in_Oceanic_islands")
```
The scripts should be run in numerical order:

1. `01_benthic_composition.R` — Data processing and benthic composition.
2. `02_ordination.R` — PERMANOVA, PCoA, envfit, and indicator analysis using 18 benthic categories.
3. `03_dbRDA.R` — Distance-based Redundancy Analysis using 10 broader benthic groups and environmental variables.
4. `04_temporal_trends_SP_transect.R` — Temporal trends for SP.
5. `05_temporal_trends_RA_transect.R` — Temporal trends for RA.
6. `06_temporal_trends_FN_transect.R` — Temporal trends for FN.
7. `07_temporal_trends_TR_transect.R` — Temporal trends for TR.

Custom functions used by the temporal analyses are available in `functions_transect.R`.

## 📦 Dependencies

Analyses were performed in R. Key packages include:

- `vegan`
- `ggplot2`
- `dplyr`
- `tidyr`
- `readr`
- `ggrepel`
- `indicspecies`
- `patchwork`
- `zyp`
- `zoo`

## 📝 Citation

If you use these data or scripts, please cite the associated Zenodo repository:

DOI: 10.5281/zenodo.21890794

## 📄 License

This project is licensed under the MIT License. See the `LICENSE` file for details.

## 📬 Contact

**Ana Monteiro-Leonel**  
anamonteiroleonel@alumni.usp.br

## 🙏 Acknowledgments

This study was financed in part by the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior - Brasil (CAPES) - Finance Code 001. Ana Monteiro-Leonel acknowledges funding from CAPES Processo PROEX: 88887358027/2019-00. Tito Lotufo acknowledges CNPq grant number 443318/2019-0.
We thank the members of the PELD ILOC program for providing the benthic data collected over the years, as well as the many researchers involved in the sampling efforts. We are grateful to ICMBio for granting the sampling authorization (SISBio #41327-54, CELF) and to the Brazilian Navy for logistical support at the St. Peter and St. Paul Archipelago and Trindade Island.
Financial support was also provided by CNPq for sampling (grant numbers #441750/2024-9 and #446005/2024-0, CELF). We thank the anonymous reviewers for their valuable comments.
