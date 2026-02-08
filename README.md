# Balancing-COVID-19-deaths-and-suicides-in-Japan-A-critical-trade-off-of-the-pandemic
Code and data for the analysis.
## Overview
This repository contains the data and R code used in the analysis of the trade-off between COVID-19-related deaths and excess suicides in Japan during the pandemic period.

The code provided here is intended to reproduce the main empirical analyses and figures reported in the associated manuscript.
## Data
All datasets used in this study are publicly available (open data) and are provided in the `data/` directory.

- `data/raw/`  
  Contains the original data sources used in the analysis (e.g., COVID-19 cases and deaths, suicide counts).

- `data/processed/`  
  Contains pre-fitted Bayesian model objects saved as `.rds` files (posterior draws).
  
## Code
All analysis scripts are written in R and are located in the `code/` directory.  
The scripts should be executed in the order listed below. Note that Scripts `05_Random_effects_model.R` and `06_Correlation_analysis.R` need to be run after Scripts `01`–`04` within the same R session, as they rely on objects created in earlier scripts.

1. `01_analysis_excess_suicide.R`  
   Estimates baseline and excess suicides by sex and age group using time-series models.

2. `02_analysis_excess_unemployment.R`  
   Estimates baseline and excess unemployment by sex and age group using time-series models.

3. `03_analysis_excess_depression.R`  
   Estimates baseline and excess Social Mood Index (depression) for the overall population using time-series models.

4. `04_SEIR_model.R`  
   Implements the SEIR model and performs epidemic-related simulations and analyses.

5. `05_Random_effects_model.R`  
   Estimates the random effects model to account for heterogeneity across groups or time periods and performs post-processing of the corresponding results.

6. `06_Correlation_analysis.R`  
   Performs correlation analyses among key variables to inform the selection of time lags.

## Reproducibility
### Software environment
- R (version 4.5.0).
- Key R packages: deSolve (version 1.40), nloptr (version 2.2.1), forecast (version 8.24.0), brms (version 2.22.0), cmdstanr (version 0.9.0).

### Notes
- Some components of the analysis pipeline are computationally intensive and may require substantial computing time (several hours) to run from scratch.
- To avoid numerical discrepancies arising from MCMC sampling across environments, pre-fitted `.rds` files are provided to ensure exact reproduction of the reported results and figures.

## Contact
For questions regarding the code or data, please contact the corresponding author.
