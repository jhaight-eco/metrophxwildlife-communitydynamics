[![DOI](https://zenodo.org/badge/631009107.svg)](https://doi.org/10.5281/zenodo.15043780)

# A repository for:
Haight, J. D.,; S. J. Hall; J. S. Lewis. Landscape modification and species traits shape seasonal wildlife community dynamics within an arid metropolitan region

This repository contains code and data used for fitting a dynamic community occupancy model to mammal occurrence data across multiple seasons, conducting a meta-analysis of species richness estimates, and producing visualizations based on these results. Occurrence data were collected via a network of 50 wildlife camera traps sited across an urbanization gradient within the urbanized Sonoran Desert region of central Arizona, centered around Phoenix, AZ, USA.


This `README` file includes information on the various scripts and datasets used for this analysis, though not every data source is saved in this repository due to data size limitations (the dynamic occupancy model output; GIS data). The manuscript includes citation of all locations from which environmental covariates were derived. This `README` file is organized into three sections corresponding to the three main folders into which the repository is organized
1. [code](#code)  
2. [data](#data)
3. [figures](#figures)

Please direct all questions to Jeffrey Haight jdhaight.eco(at)gmail.com


---
---
<div align="center"> <h3>code</h3> </div>

This folder contains all .R script and RMarkdown (.rmd) files used for conducting the analyses that produce the `model_outputs` and `model_summaries` data and the `figures`. There are a total of eight R scripts, numbered by the order in which they are to be run. Scripts 1_0 and 1_1 are for preparing data for modeling and examining environmental relationships, 2_0 through 2_4 are for conducting the analyses, and 3_0 through 3_2 are for examining model outputs and producing figures components.

| File  | Description  |
|---|---|
|**./code/1_0_ModelDataPrep.R**   | Scripts for cleaning and organizing raw datasets, in preparation for the fitting of occupancy models |
|**./code/1_1_EnvironmentalRelationships.R**   | Scripts for examining spatial and temporal relationships among environmental conditions, including urbanization, vegetation greenness, and landscape metrics | 
|**./code/2_0_model_DCM.R**   |  The Bayesian dynamic community occupancy model that we fit to multi-season species occurrence data   |
|**./code/2_1_FittingUnivariateDCMs_ScaleOptimization.R**   | Script for using JAGS to fit the univariate Bayesian dynamic community occupancy models within R, using the R package `jagsUI`, for the purpose of optimizing the spatial scale of habitat use covariates    |
|**./code/2_2_FittingGlobalDCM.R**   | Script for using JAGS to fit the global Bayesian dynamic community occupancy models containing scale-optimized covariates within R, using the R package `jagsUI` |
|**./code/2_3_model_alphadiversity_multiseason**  | The log-normal model used for the Bayesian meta-analysis of species richness estimates   |
|**./code/2_4_FittingSpeciesRichnessMetaanalysis.R**  | Script for using JAGS to fit the Bayesian log-normal meta-analysis model of species richness within R, using the R package `jagsUI`  |
|**./code/3_0_Prediction_HabitatUse.R**   |  Script for extracting posterior parameter (intercept and slopes) estimates of species occupancy (i.e., habitat use), environmental relationships, and trait effect and using them to predict species occupancy across gradients of environmental variables |
|**./code/3_1_Visualization_HabitatUse.R**  | Script for creating figures visualizing relationships between community-level and species-level occupancy and environmental variables (urbanization, patch diversity, and vegetation greenness), as well as among-species variation in relationships associated with functional traits |
|**./code/3_2_Visualization_SpeciesRichness.R**| Script for creating figures visualizing relationships between seasonal species richness and environmental variables, particularly urbanization |  


Code was run using the following R packages:  
| Package  | Version  |
|---|---|
| jagsUI | 1.6.2 |
|tidyverse| 2.0.0 |
| dplyr | 1.1.4 |
| beepr | 2.0 |
| ggplot2 | 3.5.1 |
| ggcorrplot | 0.1.4.1 |
| gghighlight | 0.4.1 |
| GGally | 2.2.1 |
|ggrepel| 0.9.5 |
|ggpubr| 0.6.0 |
| png | 0.1-8 |
| scales | 1.3.0 |
|lterpalettefinder|1.1.0.900|
|RColorBrewer| 1.1-3 |
|viridis| 0.6.5 |


---
<div align="center"> <h3>data</h3> </div>
This contains three subfolders with all the data files used for conducting the analyses and producing the `figures`:  
**`model_inputs`**  
**`model_outputs`**   
**`model_summaries`** 

 

Within **`model_inputs`**, there are five files:  

**./data/model_inputs/arizonapalettes.Rdata**  
Contains color palettes used for creating species richness figures. Colors were manually selected by Jeffrey Haight from personal photos of Arizona flora, including wildflowers, lichens, and mosses.

**./data/model_inputs/covariates_50camsite.csv**
Contains data of environmental variables for 50 wildlife camera sites across the Phoenix metropolitan area, Arizona, USA


**./data/model_inputs/EltonTraits_22mammalbird**    
Contains data of species-level attributes, including each species' names, taxonomic groupings, and species traits (body mass, log-transformed as 'logmass'; carnivory). Species traits were derived from the EltonTraits 1.0 database (https://doi.org/10.1890/13-1917.1).

 

**./data/model_inputs/ModelInputData_DCM.Rdata**  
Contains all the cleaned datasets necessary for fitting the Bayesian multi-city community occupancy model, including the following R objects:  

| Object Name	| Description   |
|---------------------------|--------|
| data_site	| A dataframe containing environmental variables at 50 wildlife camera sites in a wide tabular format (with seasonal variables as separate columns), including the modeled covariates of urbanization (e.g., `imp_1km`), patch diversity (e.g., `sdhi_5km`), and seasonal vegetation greenness (e.g., `ndvi_wd1km`, `ndvi_ww1km`, `ndvi_cw1km`) |
| data_site	| A dataframe containing environmental variables at 50 wildlife camera sites in a long tabular format (with `season` and `scale` as variables), including the modeled covariates of urbanization (`imp`), patch diversity (`sdhi`), and vegetation greenness (`ndvi`) |
| data.spp		| A dataframe of traits for the 20 focal mammal species (plus two ground-dwelling bird species), derived from the EltonTraits database (Wilman et al. 2016).	|
| K 	        	| A matrix including the number of days in which each site was surveyed by camera traps during each season (max = 96 days), as well as total by-site sampling days across all seasons (`effort`) and standardized effort (`effort_std`)|
| K.survey         	| A matrix including the number of repeat survey occassions at each site during each season (max = 6 surveys )	|
| K.std         	| A matrix including the standardized by-site variable of sampling `effort`, used as a detection covariate in the fitting of occupancy models	|
| n.season		| The number of sampling seasons|
| n.site		| The number of sites surveyed in each season|
| n.spp		| The number of species observed across all seasons	|
| n.survey		| The number of survey occasions per season	|
| names.common		| A vector of common names for each modeled species	|
| spp.model		| A vector of shortened common names for each modeled species, the names by which species occurrence data was alphabetized	|
| y	|A four-dimensional (site X survey X season X species) array of species detections (1 = detected, 0 = not detected), including only modeled species|
| ysum	|A three-dimensional (site X season X species) array of species detections summed across surveys (max value = n.survey = 6). This is the occurrence dataset formatted for fitting the dynamic community occupancy model model (**./code/model_DCM.R**) using JAGS|

Additional explanation of each R object is provided as comments within the script preparing the data (**./code/1_0_ModelDataPrep.R**).  
 

**./data/model_inputs/y_s22j50_occ16day_3season.rds**  
An R object formatted as a four-dimensional (site X survey X season X species) array of species detections (1 = detected, 0 = not detected). This precursory dataset contains detections for additional species not included in our finalized occupancy analyses (two ground-dwelling bird species), whereas the `y` object of **./data/model_inputs/ModelInputData_DCM.Rdata** contains only the 20 modeled species.  



---
<div align="center"> <h3>figures</h3> </div>

This folder contains all images utilized in the production of the manuscript figures. This folder includes visualizations produced using the R scripts (see **code** above). These visualizations were then combined within with one another within Inkscape (https://inkscape.org/) to create the figures in the published manuscript. 

The `figures` folder also includes the subfolder **./figures/mammalgraphics** which contains image files used to represent mammals species in Figure 1 of the manuscript. All mammal graphics were sourced from PhyloPic (https://www.phylopic.org/) and were utilized as part of the public domain (https://creativecommons.org/publicdomain/zero/1.0/) or under the Creative Commons Attribution 3.0 license (https://creativecommons.org/licenses/by/3.0/). Mammal graphics are accompanied by the text file **./figures/mammalgraphics/imageattributions.txt**, which specifies each image's source and provides attribution for each image, as outlined below.
| File Name  | Creative Commons License | Attribution  |  Source  |
|---|---|---|---|
| Canis latrans  | CC0 1.0 Universal Public Domain Dedication   |  Margot Michaud  |  https://www.phylopic.org/images/5a0398e3-a455-4ca6-ba86-cf3f1b25977a/canis-latrans  |
| Lepus americanus  | CC0 1.0 Universal Public Domain Dedication   | Maha Ghazal   | https://www.phylopic.org/images/8e61e166-11f4-4377-a923-9b5b597b6eba/lepus-americanus   |
| Lepus europaeus  |  CC0 1.0 Universal Public Domain Dedication  |  Anthony Caravaggi  |  https://www.phylopic.org/images/630a20e2-2d73-49ec-8cff-80f212b638c6/lepus-europaeus  |
| Lynx rufus  | CC0 1.0 Universal Public Domain Dedication   |  Margot Michaud  |  http://phylopic.org/image/20442a12-596d-4987-a668-509c19a155da/  |
| Odocoileus virginianus  | CC0 1.0 Universal Public Domain Dedication   |  Tracy A. Heath  |  http://phylopic.org/image/56f6fdb2-15d0-43b5-b13f-714f2cb0f5d0/  |
| Pecari tajacu  | CC0 1.0 Universal Public Domain Dedication   |  Steven Traver  | http://phylopic.org/image/44fb7d4f-6d59-432b-9583-a87490259789/   |
| Procyon lotor  | CC0 1.0 Universal Public Domain Dedication   |  Mathieu Basille  |  http://phylopic.org/image/e805d164-21e7-4657-979a-226f6ccc7f15/  |
| Rodent, unknown  |  N/A  |  Kinley Ragan  |  Personal communication  |
|  Urocyon cinereoargenteus | CC0 Attribution 3.0 Unported | Brian Gratwicke (photo) and T. Michael Keesey (vectorization)   |  http://phylopic.org/image/20da6c7c-2584-4cee-921b-ebd09384567b/  |

---
 
