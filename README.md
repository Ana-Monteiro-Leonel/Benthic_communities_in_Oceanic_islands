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
   
   │   ├── raw/                 # Original field data (immutable)
   
   │   ├── processed/           # Cleaned datasets in CSV format
   
   │   └── metadata/            # Data dictionary and sampling protocols
   
   ├── code/
   
   │   ├── 01_clean-data.Rmd    # Data cleaning workflow
   
   │   ├── 02_analyze.Rmd       # Statistical analyses and figures
   
   │   └── functions.R          # Custom R functions
   
   ├── results/
   
   │   ├── figures/             # Generated plots
   
   │   └── tables/              # Summary tables
   
   └── docs/                    # Additional documentation
   


## 🚀 How to Reproduce

1. Clone this repository
2. Open the R project file
3. Run the scripts in numerical order:
   - `code/01_clean-data.Rmd`
   - `code/02_analyze.Rmd`

## 📦 Dependencies

- R (version ≥ 4.0.0)
- Key packages: vegan, tidyverse, ggplot2, knitr

## 📝 Citation

If you use this data or code, please cite:

[To be added upon publication]

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📬 Contact

[Ana Monteiro-Leonel] - [anamonteiroleonel@alumni.usp.br]

## 🙏 Acknowledgments

This study was financed in part by the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior - Brasil (CAPES) - Finance Code 001. Ana Monteiro-Leonel acknowledges funding from CAPES Processo PROEX: 88887358027/2019-00. Tito Lotufo CNPq grant number 443318/2019-0. We also thank the members of the PELD ILOC program for providing the benthic data collected over the years, as well as the many researchers involved in the sampling efforts. We are grateful to ICMBio for granting the sampling authorization (SISBio #41327-54, CELF). We also thank the Brazilian Navy for logistical support at the São Pedro and São Paulo Archipelago and Trindade Island. Financial support was also provided by CNPq for sampling (grant numbers #441750/2024-9 and #446005/2024-0, CELF). We thank the anonymous referees for their valuable comments.
