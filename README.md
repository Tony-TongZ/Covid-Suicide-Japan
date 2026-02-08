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
The scripts are intended to be executed in the following order:

1. `01_prepare_data.R`  
   Data cleaning, preprocessing, and construction of analysis variables.

2. `02_fit_model.R`  
   Model specification and estimation (e.g., Bayesian models and SEIR-related analyses).

3. `03_postprocess.R`  
   Post-estimation processing, including uncertainty quantification and summary statistics.

4. `04_make_figures.R`  
   Generation of tables and figures reported in the manuscript.

## Reproducibility
### Software environment
- R (version 4.5.0)
- Key R packages: deSolve (version 1.40), nloptr (version 2.2.1), forecast (version 8.24.0), brms (version 2.22.0), cmdstanr (version 0.9.0)

### Notes
- Some components of the analysis pipeline are computationally intensive and may require substantial computing time (several hours) to run from scratch.
- To avoid numerical discrepancies arising from MCMC sampling across environments, pre-fitted `.rds` files are provided to ensure exact reproduction of the reported results and figures.

## Contact
For questions regarding the code or data, please contact the corresponding author.
