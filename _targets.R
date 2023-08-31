## Load packages required to define the pipeline
library(targets)
library(nimble)


## Source all functions in "R" folder
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir('R')


#-----------------------------#
# Workflow and model settings #
#-----------------------------#

## Seed
mySeed <- 0

## Test run vs. full run
#testRun <- TRUE # Runs a test with only 10 MCMC iterations for model fitting
testRun <- FALSE # Runs the full MCMC with 30000 iterations unless otherwise specified

## General parameters
Amax <- 5 # Number of age classes
Tmax <- 18  # Number of years
minYear <- 2004 # First year to consider
maxAge_yrs <- 10 # Age of the oldest female recorded
summer_removal <- c(6, 7, 8, 9) # Months to be removed from age-at-harvest data (summer harvest is not currently included)
area_selection <- c("Inner", "BB", "Tana") # Varanger sub-area to include in analyses (BB = Batsfjord and Berlevag areas)

# start and end of placental scars and embryo sample periods (julian day)
plac_start <- 180 #including
plac_end   <- 80  #until, not including
embr_start <- 100 #including
embr_end   <- 140 #until, not including

## Dataset names, versions, and directories
carcass.dataset.name <- "v_redfox_carcass_examination_v3"
carcass.dataset.version <- 3

rodent.dataset.name <-"v_rodents_snaptrapping_abundance_regional_v5"
rodent.dataset.version <- 5

genetics.datapath <- "Data/RedFox_genetics_immigrant_probabilities.txt"

meta.datafile <- "Data/RedFox_LiteratureData.csv"

hoening.datafile <- "Data/HoenigMod_Posteriors_fromTomPorteus.txt"


## Credentials for accessing the COAT database

# Stijn
#shapefile.dir <- "C:\\Users\\sho189\\OneDrive - UiT Office 365\\PhD\\RedfoxIPM\\Fox areas shapefile\\tana rest"
#COAT_key <- Sys.getenv("API_COAT_Stijn") # Stijn's API key for the COAT dataportal is saved as an environmental variable on the computer 

# Chloe
shapefile.dir <- "C:/Users/chloe.nater/OneDrive - NINA/Documents/Projects/RedFox_IPM/Data/shapefiles"
COAT_key <- Sys.getenv("COAT_API")


## "Switches" for running different model versions

# Covariate toggles
fitCov.mH <- FALSE # Fit covariates on mH (harvest effort)
fitCov.mO <- TRUE # Fit covariates on mO (rodent abundance x reindeer carcasses)
fitCov.Psi <- TRUE # Fit covariates on Psi (rodent abundance)
fitCov.rho <- TRUE # Fit covariates on rho (rodent abundance)
fitCov.immR <- TRUE # Fit covariates on immigration rate (rodent abundance) - only if immigration is estimated as a rate
rCov.idx <- FALSE # Use discrete vs. continuous rodent covariate
nLevels.rCov <- 2 # 2-level discrete rodent covariate
#nLevels.rCov <- 3 # 3-level discrete rodent covariate (data not currently prepared)
standSpec.rCov <- TRUE # standardize different rodent species before summing (offset catchability) v.s. simply sum all numbers

# Random year effect toggles
mO.varT <- TRUE

# Annual survival prior type toggles
HoeningPrior <- FALSE # Use prior on natural mortality derived from Hoening model
#sPriorSource <- "Bristol" # Base survival prior on data from Bristol (not hunted)
sPriorSource <- "NSweden" # Base survival prior on data from North Sweden (lightly hunted)
#sPriorSource <- "metaAll" # Base survival prior on meta-analysis including all populations
#sPriorSource <- "metaSub" # Base survival prior on meta-analysis including only not/lightly hunted populations

# Immigration parameters toggle
imm.asRate <- TRUE # Estimating immigration as a rate as opposed to numbers

# Genetic immigration data toggles (details in documentation of wrangleData_gen function
poolYrs.genData <- TRUE # Pool data across all years
useData.gen <- TRUE # Use genetic data for estimation of immigration rate
indLikelihood.genData <- FALSE # Apply an individual-level likelihood for genetic data
threshold <- 0.2
#pImm.type <- "original"
pImm.type <- "rescaled"
#pImm.type <- "LL-based"


## Toggles for LTRE analyses
HazardRates <- TRUE
PopStructure <- TRUE


#-------------------#
# Workflow pipeline #
#-------------------#

## Set target-specific options such as packages.
tar_option_set(packages = c("tidyverse", "sf", "reshape2", "remotes", "ckanr", "purrr", "dplyr", "metafor", "ggplot2", "patchwork"),
               format = "qs",
               memory = "transient", 
               garbage_collection = TRUE)


## Define Targets List
list(
  # Download raw redfox carcass data from COAT database
  tar_target(
    carcass.data.raw,
    downloadData_COAT(COAT_key = COAT_key, 
                      COATdataset.name = carcass.dataset.name,
                      COATdataset.version = carcass.dataset.version)
  ),
  
  # Reformat redfox carcass data
  tar_target(
    carcass.data,
    reformatData_carcass(Amax = Amax,   
                         summer_removal = summer_removal ,
                         area_selection = area_selection,
                         plac_start = plac_start,
                         plac_end = plac_end ,
                         embr_start = embr_start ,
                         embr_end = embr_end,
                         carcass.dataset = carcass.data.raw,
                         shapefile.dir = shapefile.dir)
  ),
  
  # Extract winter age-at-harvest matrix
  tar_target(
    wAaH.data,
    wrangleData_winterAaH(wAaH.datafile = carcass.data$AaH.matrix, 
                          Amax = Amax)
  ),

  # Extract reproduction data
  tar_target(
    rep.data,
    wrangleData_rep(P1.datafile = carcass.data$P1var, 
                    P2.datafile = carcass.data$P2var,
                    Amax = Amax, 
                    minYear = minYear)
  ),
  
  # Load and format genetic data
  tar_target(
    gen.data,
    wrangleData_gen(datapath = genetics.datapath,
                    minYear, 
                    onlyFemales = FALSE, 
                    poolYrs.genData = poolYrs.genData, 
                    threshold = threshold)
  ),
  
  # Download and reformat data on hunting effort
  tar_target(
    hunter.data,
    reformatData_hunters(area_selection = area_selection,
                         carcass.dataset = carcass.data.raw,
                         shapefile.dir = shapefile.dir)
  ),
  
  # Download raw rodent abundance data from COAT database
  tar_target(
    rodent.data.raw,
    downloadData_COAT(COAT_key = COAT_key, 
                      COATdataset.name = rodent.dataset.name,
                      COATdataset.version = rodent.dataset.version)
  ),
  
  # Reformat rodent data
  tar_target(
    rodent.data,
    reformatData_rodent(rodent.dataset = rodent.data.raw,
                        minYear = minYear)
  ),
  
  # Load and reformat reindeer data
  tar_target(
    reindeer.data,
    reformatData_reindeer(minYear = minYear,
                          Tmax = Tmax)
  ),
  
  # Collate conceptual year information
  tar_target(
    YearInfo,
    collate_yearInfo(minYear = minYear,
                     Tmax = Tmax)
  ),
  
  # Simulate and collate prior information on survival
  tar_target(
    surv.priors,
    collate_priorInfo(meta.datafile = meta.datafile,
                      simulateSD = TRUE,
                      hoening.datafile = hoening.datafile, 
                      nsim = 30, 
                      mu.t.max = 22.61062, 
                      maxAge = maxAge_yrs)
  ),
  
  # Define type of prior to use in analyses
  tar_target(
    survPriorType,
    definePriorType_AnnSurv(HoeningPrior = HoeningPrior, 
                            sPriorSource = sPriorSource)
  ),
  
  # Write model code
  tar_target(
    redfox.code,
    writeCode_redfoxIPM(indLikelihood.genData = indLikelihood.genData)
  ),
  
  # Assemble model input data
  tar_target(
    input.data,
    assemble_inputData(Amax = Amax, 
                       Tmax = Tmax, 
                       minYear = minYear,
                       maxPups = 14,
                       uLim.N = 800,
                       uLim.Imm = 3000,
                       nLevels.rCov = nLevels.rCov,
                       standSpec.rCov = standSpec.rCov,
                       poolYrs.genData = poolYrs.genData,
                       pImm.type = pImm.type,
                       wAaH.data = wAaH.data, 
                       rep.data = rep.data, 
                       gen.data = gen.data,
                       rodent.data = rodent.data, 
                       reindeer.data = reindeer.data,
                       hunter.data = hunter.data, 
                       surv.priors = surv.priors,
                       survPriorType = survPriorType)
  ),
  
  # Set up model
  tar_target(
    model.setup,
    setupModel(modelCode = redfox.code, 
               nim.data = input.data$nim.data, 
               nim.constants = input.data$nim.constants, 
               minN1 = c(600, 50, 50, 50, 50), 
               maxN1 = c(800, 400, 400, 400, 400), 
               minImm = 50, 
               maxImm = 600,
               fitCov.mH = fitCov.mH, 
               fitCov.mO = fitCov.mO,
               fitCov.Psi = fitCov.Psi, 
               fitCov.rho = fitCov.rho,
               fitCov.immR = fitCov.immR,
               rCov.idx = rCov.idx,
               mO.varT = mO.varT,
               HoeningPrior = HoeningPrior,
               testRun = testRun,
               initVals.seed = mySeed
    )
    
  ),
  
  # Run model
  tar_target(
    IPM.out,
    nimbleMCMC(code = model.setup$modelCode,
               data = input.data$nim.data, 
               constants = input.data$nim.constants,
               inits = model.setup$initVals, 
               monitors = model.setup$modelParams,
               nchains = model.setup$mcmcParams$nchains, 
               niter = model.setup$mcmcParams$niter, 
               nburnin = model.setup$mcmcParams$nburn, 
               thin = model.setup$mcmcParams$nthin, 
               samplesAsCodaMCMC = TRUE, 
               setSeed = mySeed)
  ),
  
  # Save model output as .rds
  tar_target(
    IPM.out.saveRDS,
    saveRDS(IPM.out, file = "RedfoxIPM_ModelRun.rds"),
    format = "file"
  ),
  
  # Plot basic IPM results
  tar_target(
    basePlots,
    plotIPM_basicOutputs(MCMC.samples = IPM.out,
                         nim.data = input.data$nim.data,
                         Amax = Amax, Tmax = Tmax, minYear = minYear),
    format = "file"
    
  ),
  
  # Plot vital rate - covariate relationships
  tar_target(
    covariatePlots,
    plotIPM_covariateEffects(MCMC.samples = IPM.out,
                             rCov.idx = rCov.idx,
                             rodentMIN = -1.75, rodentMAX = 4,
                             reindeerMIN = -1.5, reindeerMAX = 1.5,
                             AgeClass = 1),
    format = "file"
    
  ),
  
  # Extract parameter samples for further analyses
  tar_target(
    paramSamples,
    extractParamSamples(MCMC.samples = IPM.out,
                        Amax = Amax, Tmax = Tmax)
  ),
  
  # Calculate sensitivities and elasticities
  tar_target(
    sensitivities,
    calculateSensitivities(paramSamples = paramSamples,
                           Amax = Amax)
  ),
  
  # Plot sensitivities and elasticities
  tar_target(
    sensitivityPlots,
    plotSensitivities(sensitivities = sensitivities,
                      Amax = Amax),
    format = "file"
  ),
  
  # Run random design transient LTRE analysis
  tar_target(
    randomLTRE,
    runLTRE_randomDesign(paramSamples = paramSamples, 
                         sensitivities = sensitivities, 
                         Amax = Amax, Tmax = Tmax, 
                         HazardRates = HazardRates, 
                         PopStructure = PopStructure)
  ),
  
  # Plot results from random design transient LTRE
  tar_target(
    LTREPlots.random,
    plotLTRE_randomDesign(LTRE_results = randomLTRE,
                          Amax = Amax,
                          HazardRates = HazardRates,
                          PopStructure = PopStructure),
    format = "file"
  ),
  
  # Run fixed design transient LTRE analysis
  tar_target(
    fixedLTRE,
    runLTRE_fixedDesign_allYears(paramSamples = paramSamples, 
                                 Amax = Amax, Tmax = Tmax, 
                                 HazardRates = HazardRates, 
                                 PopStructure = PopStructure)
  ),
  
  # Plot results from fixed design transient LTRE
  tar_target(
    LTREPlots.fixed,
    plotLTRE_fixedDesign(LTRE_results = fixedLTRE, 
                         Amax = Amax, Tmax = Tmax, minYear = minYear, 
                         HazardRates = HazardRates, 
                         PopStructure = PopStructure)
  )
)

# Test by running tar_manifest(fields = all_of("command")) and tar_visnetwork() in the console

# Run workflow using tar_make() in the console
# Check with tar_network() in the console

# We can then use tar_read() and tar_load() to inspect and work with results
