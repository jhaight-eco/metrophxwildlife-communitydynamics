# metrophxwildlife-communitydynamics

# A repository for:
Haight, J. D.,; S. J. Hall; J. S. Lewis. Landscape modification and species traits shape seasonal wildlife community dynamics within an arid metropolitan region

This repository contains code and data used for fitting a dynamic community occupancy model to mammal occurrence data across multiple seasons, conducting a meta-analysis of species richness estimates, and producing visualizations based on these results. Occurrence data were collected via a network of 50 wildlife camera traps sited across an urbanization gradient within the urbanized Sonoran Desert region of central Arizona, centered around Phoenix, AZ, USA.


This `README` file includes information on the various scripts and datasets used for this analysis, though not every data source is saved in this repository due to data size limitations (the dynamic occupancy model output; GIS data). The manuscript includes citation of all locations from which environmental covariates were derived. This `README` file is organized into three sections corresponding to the three main folders into which the repository is organized
3. [code](#code)  
1. [data](#data)
2. [figures](#figures)

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
| jagsUI | 1.5.2 |
| dplyr | 1.0.8 |
| beepr | 1.3 |
| reshape2 | 1.4.4 |
| ggplot2 | 3.3.5 |
| [peRReo](https://github.com/jbgb13/peRReo) | 0.1.0 |
| png | 0.1-7 |
| scales | 1.2.0 |
| gghighlight | 0.3.2 |
| GGally | 2.1.2 |

 

