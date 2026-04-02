# Benthic_communities_in_Oceanic_islands
Data and R scripts for analyzing benthic communities across Southwestern Atlantic oceanic islands
# Benthic communities in Southwestern Atlantic oceanic islands

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

## 📖 About

The benthic community is essential to reef ecosystems, contributing to their biodiversity and complexity. In the Southwest Atlantic, the Epilithic Algae Matrix (EAM) and macroalgae are the most abundant components, including in better-conserved areas such as the Brazilian oceanic islands: 

- **São Pedro and São Paulo Archipelago** (SP)
- **Rocas Atoll** (RA)  
- **Fernando de Noronha Archipelago** (FN)
- **Trindade Island** (TR)

## 🎯 Objective
This study aims to:
•	Quantify spatial and temporal variation in benthic community structure
•	Identify environmental drivers shaping community composition
•	Assess temporal trends in dominant benthic functional groups

## 📊 Methods

Shallow reefs were sampled annually between 2013-2019 using photo-quadrats collected along semi-fixed transects. 
Benthic cover was:
•	Identified at the image level
•	Aggregated to the transect level (used as the sampling unit to avoid pseudoreplication)
•	Grouped into morphofunctional categories
Temporal trends were assessed using:
•	LOESS smoothing for visualization
•	Mann–Kendall trend test applied to observed (non-standardized) annual mean cover values derived from transect-level aggegation, ensuring independence among observations
•	Sen’s slope estimator to quantify trend magnitude
This approach ensures that trends reflect real ecological changes rather than artifacts of data standardization.


## 🔬 Key Results
📊 Benthic Richness and Composition
Island	                      Richness	            Dominant Groups	                        Distinctive Feature
SP (St. Peter and St. Paul)	  23 taxa	              EAM, Macroalgae, Suspension feeders	    Uniquely structured by Palythoa caribaeorum (10.2% cover)
RA (Rocas Atoll)	            62 taxa	              EAM (29.6%), Macroalgae (18.6%)	        High richness, similar to FN
FN (Fernando de Noronha)	    48 taxa	              EAM (43.7%), Macroalgae (26.3%)	        Coral resilience despite regional stressors
TR (Trindade)	                35 taxa	              Macroalgae (~27%), EAM (16%)	          Most distinct assemblage, high abiotic cover

🌊 Environmental Drivers
* Wave power emerged as the primary factor structuring benthic communities at Trindade, overriding temperature effects
* Temperature, photosynthetic active radiation (PAR), and water transparency explained compositional differences among islands
* Geographic proximity and shared environmental forcing link FN and RA assemblages

📈 Temporal Trajectories (2013-2019)
Island	Trend	Key Observation
SP	  Significant increase in macroalgae; Decline in zoanthids
RA	  Consistent but non-significant trends
FN	  Benthic decline	Episodic cyanobacterial blooms; corals remain resilient
TR	  Consistent but non-significant trends

🎯 Main Takeaways
✅ Island-specific trajectories: Each island follows its own benthic path
✅ MAL expands in SP
✅ Algal reorganization (RA, TR)
✅ Corals resilient in Noronha despite environmental stress
✅ Wave power, not temperature, drives Trindade's benthic composition

🎯 Final Remarks
As sentinels of the South Atlantic, these islands record not only the history of their impacts, but also the potential for their recovery. The unfolding of this story will depend on the continuity of monitoring and the effectiveness of recently implemented protection measures.

## 📁 Repository Structure
Benthic_communities_in_Southwestern_Atlantic_oceanic_islands/

├── data/

│   ├── raw/

│   │   ├── benthic_complete_data.csv

│   │   └── environment.csv

│   └── processed/

│       ├── benthic_cover_photos.csv

│       ├── benthic_cover_transects.csv

│       └── benthic_cover_summary.csv

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

## 🚀 How to Reproduce
⚠️ Important: Set your working directory to the project root before running scripts.
Example:
setwd("path/to/Benthic_communities_in_Southwestern_Atlantic_oceanic_islands")

All analyses were performed in R version 4.2.3. Scripts should be run in numerical order:

1. `01_benthic_composition.R` — Data processing and stacked bar plot
2. `02_ordination.R` — PCoA, PERMANOVA, envfit, and indicator species
3. `03_dbRDA.R` — Distance-based Redundancy Analysis
4. `04_07_temporal_trends_*_transect.R` — Temporal trends per island (LOESS + Mann-Kendall)

Custom functions are available in `functions_transect.R`.

## 📦 Dependencies

- R (version ≥ 4.0.0)
- Key packages: "vegan", "ggplot2", "dplyr", "tidyr", "readr", "ggrepel", "indicspecies", "patchwork", "zyp"

## 📝 Citation

If you use this data or code, please cite:

[To be added upon publication]

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📬 Contact

[Ana Monteiro-Leonel] - [anamonteiroleonel@alumni.usp.br]

## 🙏 Acknowledgments

This study was financed in part by the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior - Brasil (CAPES) - Finance Code 001. Ana Monteiro-Leonel acknowledges funding from CAPES Processo PROEX: 88887358027/2019-00. Tito Lotufo CNPq grant number 443318/2019-0. We also thank the members of the PELD ILOC program for providing the benthic data collected over the years, as well as the many researchers involved in the sampling efforts. We are grateful to ICMBio for granting the sampling authorization (SISBio #41327-54, CELF). We also thank the Brazilian Navy for logistical support at the São Pedro and São Paulo Archipelago and Trindade Island. Financial support was also provided by CNPq for sampling (grant numbers #441750/2024-9 and #446005/2024-0, CELF). We thank the anonymous referees for their valuable comments.
