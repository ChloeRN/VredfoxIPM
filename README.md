# VredfoxIPM: an integrated population modelling workflow for supporting mesopredator management    

This repository contains code for a reproducible, modular workflow for Bayesian integrated population analysis (using Nimble, <https://github.com/nimble-dev/nimble>) for hunted red foxes on the Varanger penninsula, Norway. 

The workflow is available as both a manual implementation ("RedFox_IPM_Analysis.R") and as a semi-automated  "targets" pipeline (see <https://books.ropensci.org/targets/>) and contains the following parts: 
- Download of raw data from COAT data portal (<https://data.coat.no/>)
- Wrangling, formatting, and preparation of observational data and prior information
- Definition and setup of integrated population model (IPM)
- Complete simulation of initial values
- Running of IPM using MCMC
- Extraction and visualization of primary results
- Follow-up analyses (incl. vizualizations)

A schematic overview of the workflow and its parts is shown below: 
![](Diagrams/RedFox_targets_light.png)

The workflow is fully function-based, and each function is documented in its respective script within the folder "R".
The folder "vignettes" contains a richly documented "walkthrough" of the manual workflow ("workflow_walkthrough.Rmd"), including instructions for installing necessary dependencies and accessing all requisite data. 

Auxiliary data, computed results, and finished plots and other visualizations can be found on OSF: 
<https://osf.io/756re/>

The scientific article associated with this repository and presenting the workflow and results on Varanger red foxes is available as a pre-print here: <https://ecoevorxiv.org/repository/view/7489/>

We extensively used GitHub issues to document model development, sensitivity analyses, etc.  
Some topics potentially of interest include: 

- <https://github.com/ChloeRN/VredfoxIPM/issues/78> (incl. sub-issues): Summary of model revisions v1 to v2
- <https://github.com/ChloeRN/VredfoxIPM/issues/81> & <https://github.com/ChloeRN/VredfoxIPM/issues/83>: Stepwise addition of density-dependence and compensatory mechanisms
- <https://github.com/ChloeRN/VredfoxIPM/issues/82>: Comparing models with/without secondary summer census and relaxing assumption that there is no natural mortality during summer
- <https://github.com/ChloeRN/VredfoxIPM/issues/87>: Optimizing density-dependent IPMs for PVA projections
- <https://github.com/ChloeRN/VredfoxIPM/issues/86>: Modelling immigration by rate vs. numbers
- <https://github.com/ChloeRN/VredfoxIPM/issues/80>: Exploring different options for rodent covariate (species-aggregated vs. lemmings and voles separately)
- <https://github.com/ChloeRN/VredfoxIPM/issues/79>: Testing the effects of including/excluding additional environmental covariate for reindeer carcass availability
- <https://github.com/ChloeRN/VredfoxIPM/issues/90>: Relaxing assumption that all immigrants are foxes in their first year of life

