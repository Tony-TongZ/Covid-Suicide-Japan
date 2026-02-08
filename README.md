# Balancing-COVID-19-deaths-and-suicides-in-Japan-A-critical-trade-off-of-the-pandemic
Code and data for the analysis
## Overview
This repository contains the data and R code used in the analysis of the trade-off between COVID-19-related deaths and excess suicides in Japan during the pandemic period.

The code provided here is intended to reproduce the main empirical analyses and figures reported in the associated manuscript.
## Data
All data files are stored in the `data/` directory.

- `data/raw/`  
  Contains the original data sources used in the analysis (e.g., COVID-19 cases and deaths, suicide counts, population statistics).

- `data/processed/`  
  Contains processed datasets and intermediate objects generated from the raw data and used for model estimation and visualization.

If certain data cannot be publicly shared, instructions for access or data sources are documented within the corresponding directory.

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
- R (version 4.x)
- Key R packages: brms, cmdstanr, deSolve, tidyverse

### Notes
- Some model estimation procedures are computationally intensive and may require several hours to run.
- Where applicable, precomputed intermediate results are provided to facilitate replication without re-running all models from scratch.

## Contact
For questions regarding the code or data, please contact the corresponding author.
