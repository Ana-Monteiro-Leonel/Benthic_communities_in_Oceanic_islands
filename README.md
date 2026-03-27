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

This study analyzes the composition and spatio-temporal structure of benthic reef communities on these islands and assesses which environmental variables better explain similarities or dissimilarities between sites.

## 📊 Methods

Shallow reefs were sampled annually from 2013-2019 using photo-quadrats on semi-fixed transects. We estimated the cover of benthic organisms and classified them into morphofunctional groups.

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
SP	  EAM expansion	Progressive consolidation of epilithic algal matrix
RA	  EAM expansion	Similar trajectory to SP
FN	  Benthic decline	Episodic cyanobacterial blooms; corals remain resilient
TR	  CCA increase	Shift toward calcareous algae dominance

🎯 Main Takeaways
✅ Island-specific trajectories: Each island follows its own benthic path
✅ EAM expands in equatorial islands (SP and RA)
✅ CCA increases in Trindade over time
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
│   └── functions.R
├── results/
│   ├── figures/
│   │   ├── Figure_2_benthic_composition.png
│   │   ├── Figure_2B_PCoA_ordination.png
│   │   ├── Figure_3_dbRDA_ordination.png
│   │   ├── Figure_4_temporal_trends_SP.png
│   │   ├── Figure_5_temporal_trends_RA.png
│   │   ├── Figure_6_temporal_trends_FN.png
│   │   └── Figure_7_temporal_trends_TR.png
│   └── tables/
│       ├── Table_1_PERMANOVA_results.csv
│       ├── Table_2_envfit_results.csv
│       ├── Table_3_indicator_species.csv
│       ├── Table_4_dbRDA_summary.csv
│       ├── Table_5_dbRDA_axes_summary.csv
│       ├── Table_6_env_variables_significance_sequential.csv
│       ├── Table_6_env_variables_significance_marginal.csv
│       ├── Table_VIF_results.csv
│       ├── Table_rankindex_results.csv
│       └── Table_S1_relative_abundance.csv
└── README.md

## 🚀 How to Reproduce
All analyses were performed in R version 4.2.3. Scripts should be run in numerical order:

1. `01_benthic_composition.R` — Data processing and stacked bar plot
2. `02_ordination.R` — PCoA, PERMANOVA, envfit, and indicator species
3. `03_dbRDA.R` — Distance-based Redundancy Analysis
4. `04_temporal_trends_*.R` — Temporal trends per island (LOESS + Mann-Kendall)

Custom functions are available in `functions.R`.

## 📦 Dependencies

- R (version ≥ 4.0.0)
- Key packages: "vegan", "ggplot2", "dplyr", "tidyr", "readr", "ggrepel", "indicspecies", "patchwork"

## 📝 Citation

If you use this data or code, please cite:

[To be added upon publication]

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📬 Contact

[Ana Monteiro-Leonel] - [anamonteiroleonel@alumni.usp.br]

## 🙏 Acknowledgments

This study was financed in part by the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior - Brasil (CAPES) - Finance Code 001. Ana Monteiro-Leonel acknowledges funding from CAPES Processo PROEX: 88887358027/2019-00. Tito Lotufo CNPq grant number 443318/2019-0. We also thank the members of the PELD ILOC program for providing the benthic data collected over the years, as well as the many researchers involved in the sampling efforts. We are grateful to ICMBio for granting the sampling authorization (SISBio #41327-54, CELF). We also thank the Brazilian Navy for logistical support at the São Pedro and São Paulo Archipelago and Trindade Island. Financial support was also provided by CNPq for sampling (grant numbers #441750/2024-9 and #446005/2024-0, CELF). We thank the anonymous referees for their valuable comments.
