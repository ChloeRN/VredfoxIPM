library(ggplot2)
library(nimble)
library(sf)
library(reshape2)
library(remotes)
library(ckanr)
library(purrr)
library(dplyr)
library(metafor)
library(patchwork)
library(coda)

#**********#
# 0) SETUP #
#**********#

## Set seed
mySeed <- 57

## Set general parameters
Amax <- 5 # Number of age classes
Tmax <- 20  # Number of years
minYear <- 2004 # First year to consider
maxAge_yrs <- 10 # Age of the oldest female recorded
summer_removal <- c(6,7,8,9)    #removal of summer months: numerical months to be removed from winter age at harvest data
winter_removal <- c(1:6, 10:12) #removal of winter months: numerical months to be removed from summer age at harvest data
area_selection <- c("Inner", "BB",  "Tana")# choosing varanger sub area ("Inner" / "BB" / "Tana)     ((BB = Batsfjord and Berlevag areas))
# start and end of placental scars and embryo sample periods (julian day)
plac_start <- 180 #including
plac_end   <- 80  #until, not including
embr_start <- 100 #including
embr_end   <- 140 #until, not including

# Normalizing value for population size when modelling density-dependence
normN <- 500 # Based on mean/median of estimated N.tot-Imm 

## set dataset names, versions, and directories, and access
hunting.dataset.name <- "v_redfox_hunting_v3"
hunting.dataset.version <- 3

carcass.dataset.name <- "v_redfox_carcass_examination_v4"
carcass.dataset.version <- 4

rodent.dataset.name <-"v_rodents_snaptrapping_abundance_regional_v7"
rodent.dataset.version <- 7

# Stijn
shapefile.dir <- "C:\\Users\\sho189\\OneDrive - UiT Office 365\\PhD\\RedfoxIPM\\Fox areas shapefile\\tana rest"
COAT_key <- Sys.getenv("API_COAT_Stijn") # Stijn's API key for the COAT dataportal is saved as an environmental variable on the computer 

# Chloe
#shapefile.dir <- "C:/Users/chloe.nater/OneDrive - NINA/Documents/Projects/RedFox_IPM/Data/shapefiles"
shapefile.dir <- "Data/shapefiles"
COAT_key <- Sys.getenv("COAT_API")

## Source all functions in "R" folder
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir('R')

## Set "switches" for running different model versions

# Covariate toggles
fitCov.mH <- FALSE # Fit covariates on mH (harvest effort)
fitCov.mO <- TRUE # Fit covariates on mO (rodent abundance)
fitCov.Psi <- TRUE # Fit covariates on Psi (rodent abundance)
fitCov.rho <- TRUE # Fit covariates on rho (rodent abundance)
fitCov.immR <- TRUE # Fit covariates on immigration rate (rodent abundance) - only if immigration is estimated as a rate
rCov.idx <- FALSE # Use discrete vs. continuous rodent covariate
nLevels.rCov <- 2 # 2-level discrete rodent covariate
#nLevels.rCov <- 3 # 3-level discrete rodent covariate (data not currently prepared)
standSpec.rCov <- TRUE # standardize different rodent species before summing (offset catchability) v.s. simply sum all numbers
reinCov.VarTana <- TRUE # Calculate the reindeer carcass data count covariate using Varanger (+Tana) municipalities as geographical area. FALSE is for whole of Eastern Finnmark

# Random year effect toggles
mO.varT <- TRUE

# Age-at-harvest data toggles
add.sumr.unaged <- TRUE # Add summer harvested individuals as un-aged individuals to the total harvested individuals in winter
saAH.years <- c(2005:2012) # Years for which the summer age at harvest matrix should be constructed (e.g. years in which summer harvest was aged consistently)

# Annual survival prior type toggles
HoenigPrior <- FALSE # Use prior on natural mortality derived from Hoenig model
#sPriorSource <- "Bristol" # Base survival prior on data from Bristol (not hunted)
#sPriorSource <- "NSweden" # Base survival prior on data from North Sweden (lightly hunted)
sPriorSource <- "metaAll" # Base survival prior on meta-analysis including all populations
#sPriorSource <- "metaSub" # Base survival prior on meta-analysis including only not/lightly hunted populations

# Immigration parameters toggle
imm.asRate <- TRUE # Estimating immigration as a rate as opposed to numbers

# Genetic immigration data toggles (details in documentation of wrangleData_gen function
poolYrs.genData <- TRUE # Pool data across all years
useData.gen <- TRUE # Use genetic data for estimation of immigration rate
indLikelihood.genData <- FALSE # Apply an individual-level likelihood for genetic data
threshold <- 0.05
#threshold <- 0.2
#pImm.type <- "original"
pImm.type <- "rescaled"
#pImm.type <- "LL-based"

# Den survey prior and data toggles
useData.pup <- TRUE
useInfPrior.S0 <- FALSE

## Changes to denning survival prior
S0.mean.offset <- 0
S0.sd.factor <- 1

## Density effects toggles
DD.mO <- TRUE
DD.immR <- TRUE
DDxRodent <- FALSE

## Compensation toggles
comp.mO <- TRUE
comp.immR <- FALSE
comp.RE <- FALSE


#*********************#
# 1) DATA PREPARATION #
#*********************#

# 1a) Download and reformat carcass data #
#----------------------------------------#

## Download hunting data (this is the record of foxes hunted, before they end up in the carcass examination lab)
hunting.data.raw <- downloadData_COAT(COAT_key = COAT_key, 
                                      COATdataset.name = hunting.dataset.name,
                                      COATdataset.version = hunting.dataset.version)

## Reformat hunting data
hunting.data  <- reformatData_hunting(summer_removal = summer_removal,
                                      hunting.dataset = hunting.data.raw)


## Download carcass data
carcass.data.raw <- downloadData_COAT(COAT_key = COAT_key, 
                                     COATdataset.name = carcass.dataset.name,
                                     COATdataset.version = carcass.dataset.version)

## Reformat carcass data
carcass.data <- reformatData_carcass(Amax = Amax,   
                                     summer_removal = summer_removal ,
                                     winter_removal = winter_removal ,
                                     area_selection = area_selection,
                                     plac_start = plac_start,
                                     plac_end = plac_end ,
                                     embr_start = embr_start ,
                                     embr_end = embr_end,
                                     carcass.dataset = carcass.data.raw,
                                     shapefile.dir = shapefile.dir,
                                     add.sumr.unaged = add.sumr.unaged, 
                                     saAH.years = saAH.years,
                                     hunting.data = hunting.data)



# 1b) Age-at-Harvest data #
#-------------------------#

## Winter AaH data
wAaH.data <- wrangleData_AaH(AaH.datafile = carcass.data$WAaH.matrix, 
                             Amax = Amax)
## Summer AaH data
sAaH.data <- wrangleData_AaH(AaH.datafile = carcass.data$SAaH.matrix, 
                             Amax = Amax)


# 1c) Reproduction data #
#-----------------------#

## Set data paths/filenames
P1.datafile <- carcass.data$P1var # Placental scar/embryo count
P2.datafile <- carcass.data$P2var # Presence of placental scars/embryos/pregnancy signs

## Prepare reproduction data
rep.data <- wrangleData_rep(P1.datafile = P1.datafile, 
                            P2.datafile = P2.datafile,
                            Amax = Amax, 
                            minYear = minYear)


# 1d) Genetic data #
#------------------#

## Set data paths
genetics.datapath <- "Data/RedFox_genetics_immigrant_probabilities.txt"
#genetics.datapath <- "Data/RedFox_genetics_immigrant_probabilities_LvarLother.txt"

## Prepare genetic data
gen.data <- wrangleData_gen(datapath = genetics.datapath,
                            minYear, 
                            onlyFemales = FALSE, 
                            poolYrs.genData = poolYrs.genData, 
                            threshold = threshold)


# 1e) Opportunistic pup observation data #
#----------------------------------------#

## Set data path
pups.datapath <- "Data/Rfox_early_litter_sizes.csv"

## Prepare pup observation data
pup.data <- wrangleData_pup(datapath = pups.datapath,
                            minYear = minYear)


# 1f) Harvest effort data #
#-------------------------#

## Prepare harvest effort data
hunter.data <- reformatData_hunters(area_selection = area_selection,
                                    carcass.dataset = carcass.data.raw,
                                    shapefile.dir = shapefile.dir)


# 1g) Environmental data #
#------------------------#

## Download rodent data
rodent.data.raw <- downloadData_COAT(COAT_key = COAT_key, 
                                     COATdataset.name = rodent.dataset.name,
                                     COATdataset.version = rodent.dataset.version)

## Reformat rodent data
rodent.data <- reformatData_rodent(rodent.dataset = rodent.data.raw,
                                          minYear = minYear)

## Reformat reindeer data
reindeer.data <- reformatData_reindeer(minYear = minYear,
                                       Tmax = Tmax,
                                       reinCov.VarTana = reinCov.VarTana)


# 1h) Conceptual year information #
#---------------------------------#

YearInfo <- collate_yearInfo(minYear = minYear,
                             Tmax = Tmax)


#**********************#
# 2) PRIOR INFORMATION #
#**********************#

## Parameters/paths for making informative priors for survival based on meta-analysis of literature data
meta.datafile <- "Data/RedFox_LiteratureData.csv"
simulateSD <- TRUE

## Parameters/paths for making informative priors for natural mortality using Tom Porteus' Hoenig model approach
mu.t.max <- 22.61062
hoenig.datafile <- "Data/HoenigMod_Posteriors_fromTomPorteus.txt"
nsim <- 30

## Collate all prior information
surv.priors <- collate_priorInfo(meta.datafile = meta.datafile,
                                 simulateSD = simulateSD,
                                 hoenig.datafile = hoenig.datafile, 
                                 nsim = nsim, 
                                 mu.t.max = mu.t.max, 
                                 maxAge = maxAge_yrs,
                                 S0.mean.offset = S0.mean.offset,
                                 S0.sd.factor = S0.sd.factor)

## Define type of prior to use for annual survival
survPriorType <- definePriorType_AnnSurv(HoenigPrior = HoenigPrior, 
                                         sPriorSource = sPriorSource)

#****************#
# 3) MODEL SETUP #
#****************#

# 3a) Write model code #
#----------------------#

redfox.code <- writeCode_redfoxIPM(indLikelihood.genData = indLikelihood.genData)
#redfox.code <- writeCode_redfoxIndepMod(indLikelihood.genData = indLikelihood.genData)


# 3b) Assemble IPM input data #
#-----------------------------#

input.data <- assemble_inputData(Amax = Amax, 
                                 Tmax = Tmax, 
                                 minYear = minYear,
                                 maxPups = 14,
                                 uLim.N = 800,
                                 uLim.Imm = 3000,
                                 normN = normN,
                                 nLevels.rCov = nLevels.rCov,
                                 standSpec.rCov = standSpec.rCov,
                                 poolYrs.genData = poolYrs.genData,
                                 pImm.type = pImm.type,
                                 wAaH.data = wAaH.data, 
                                 sAaH.data = sAaH.data,
                                 rep.data = rep.data, 
                                 gen.data = gen.data,
                                 pup.data = pup.data,
                                 rodent.data = rodent.data, 
                                 reindeer.data = reindeer.data,
                                 hunter.data = hunter.data, 
                                 surv.priors = surv.priors,
                                 survPriorType = survPriorType)



# Adjustments for running independent models
#input.data$nim.constants$uLim.N_1 <- 3500
#input.data$nim.constants$uLim.N_2to5 <- 1500


# 3c) Set up for model run (incl. simulating initial values) #
#------------------------------------------------------------#

model.setup <- setupModel(modelCode = redfox.code, 
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
                          HoenigPrior = HoenigPrior,
                          imm.asRate = imm.asRate,
                          DD.mO = DD.mO, 
                          DD.immR = DD.immR,
                          DDxRodent = DDxRodent,
                          comp.mO = comp.mO,
                          comp.immR = comp.immR,
                          comp.RE = comp.RE,
                          testRun = FALSE,
                          initVals.seed = mySeed
                          )


# Adjustments for running independent models
#for(i in 1:model.setup$mcmcParams$nchains){
#  model.setup$initVals[[i]]$meanLS <- c(0, runif(Tmax - 1, 2, 10))
#}
#
#model.setup$modelParams <- model.setup$modelParams[which(!(model.setup$modelParams %in% c("initN", "N.tot", "B.tot", "R.tot", "B", "L", "R", "Imm")))]

####################
# 4) MODEL FITTING #
####################

t1 <- Sys.time()
IPM.out <- nimbleMCMC(code = model.setup$modelCode,
                      data = input.data$nim.data, 
                      constants = input.data$nim.constants,
                      inits = model.setup$initVals, 
                      monitors = model.setup$modelParams,
                      nchains = model.setup$mcmcParams$nchains, 
                      niter = model.setup$mcmcParams$niter, 
                      nburnin = model.setup$mcmcParams$nburn, 
                      thin = model.setup$mcmcParams$nthin, 
                      samplesAsCodaMCMC = TRUE, 
                      setSeed = rep(mySeed, model.setup$mcmcParams$nchains))
Sys.time() - t1


saveRDS(IPM.out, file = "RedFoxIPM_main.rds") 

#saveRDS(IPM.out, file = "RedFoxIPM_genData1.rds")
#saveRDS(IPM.out, file = "RedFoxIPM_genData2.rds")
#saveRDS(IPM.out, file = "RedFoxIPM_survPrior1.rds")
#saveRDS(IPM.out, file = "RedFoxIPM_survPrior2.rds") 
#saveRDS(IPM.out, file = "RedFoxIPM_survPrior3.rds")
#saveRDS(IPM.out, file = "RedFoxIndepModels.rds")


#MCMCvis::MCMCtrace(IPM.out)


########################
# 5) MODEL COMPARISONS #
########################

## Updated data models
compareModels(Amax = Amax, 
              Tmax = Tmax, 
              minYear = minYear, 
              post.filepaths = c("RedFoxIPM_main.rds",
                                 "RedFoxIPM_main_huntingData.rds",
                                 "RedFoxIPM_main_2004-2021.rds"), 
              model.names = c("Updated data (-2024)", 
                              "Updated data (-2024) + corr.",
                              "Original data (-2022)"), 
              plotFolder = "Plots/Comp_DataUpdate2024")

## Genetic data likelihoods
compareModels(Amax = Amax, 
              Tmax = Tmax, 
              minYear = minYear, 
              post.filepaths = c("RedFoxIPM_main.rds",
                                 "RedFoxIPM_genData1.rds",
                                 "RedFoxIPM_genData2.rds"), 
              model.names = c("sum. likelihood (th = 0.05)", 
                              "sum. likelihood (th = 0.2)",
                              "ind. likelihood (rescaled p)"), 
              plotFolder = "Plots/CompFinal_GenData")


## Survival priors
compareModels(Amax = Amax, 
              Tmax = Tmax, 
              minYear = minYear, 
              post.filepaths = c("RedFoxIPM_main.rds", 
                                 "RedFoxIPM_survPrior1.rds",
                                 "RedFoxIPM_survPrior2.rds",
                                 "RedFoxIPM_survPrior3.rds"), 
              model.names = c("Meta-analysis", 
                              "Hoenig model",
                              "North Sweden",
                              "Bristol"), 
              plotFolder = "Plots/CompFinal_SurvPriors")


## Integrated vs. independent model
compareModels(Amax = Amax, 
              Tmax = Tmax, 
              minYear = minYear, 
              post.filepaths = c("RedFoxIPM_main_noDD&comp.rds",
                                 "RedFoxIndepModels.rds"), 
              model.names = c("Integrated", 
                              "Independent"), 
              plotFolder = "Plots/CompFinal_IndepModels")


###########################################
# 6) IPM RESULTS - STUDY PERIOD ESTIMATES #
###########################################

IPM.out <- readRDS("RedfoxIPM_main.rds")
#IPM.out <- readRDS("RedfoxIPM_ModelRun.rds")

## Write posterior summaries to files
writePostSummaries(MCMC.samples = IPM.out,
                   Amax = Amax,
                   minYear = minYear)

## Plot basic IPM outputs (vital rate & population size estimates)
plotIPM_basicOutputs(MCMC.samples = IPM.out,
                     nim.data = input.data$nim.data,
                     Amax = Amax, Tmax = Tmax, minYear = minYear)

## Plot covariate relationships
plotIPM_covariateEffects(MCMC.samples = IPM.out,
                         rCov.idx = rCov.idx,
                         rodentMIN = -1.75, rodentMAX = 4,
                         mHdevMIN = -0.6, mHdevMAX = 0.6,
                         densityMIN = -0.5, densityMAX = 1,
                         normN = normN,
                         AgeClass = 1) 

#########################
# 7) FOLLOW-UP ANALYSES #
#########################

## Extract parameter samples
paramSamples <- extractParamSamples(MCMC.samples = IPM.out,
                                    Amax = Amax, Tmax = Tmax)

## Calculate sensitivities and elasticities
sensitivities <- calculateSensitivities(paramSamples = paramSamples,
                                        Amax = Amax)

## Plot sensitivities
plotSensitivities(sensitivities = sensitivities,
                  Amax = Amax)

## Set LTRE options
HazardRates <- TRUE
PopStructure <- TRUE

## Run random design LTRE
randomLTRE <- runLTRE_randomDesign(paramSamples = paramSamples, 
                                   sensitivities = sensitivities, 
                                   Amax = Amax, Tmax = Tmax, 
                                   HazardRates = HazardRates, 
                                   PopStructure = PopStructure)

## Plot results from random design LTRE
plotLTRE_randomDesign(LTRE_results = randomLTRE,
                      Amax = Amax,
                      HazardRates = HazardRates,
                      PopStructure = PopStructure)

## Run fixed design LTRE
fixedLTRE <- runLTRE_fixedDesign_allYears(paramSamples = paramSamples, 
                                          Amax = Amax, Tmax = Tmax, 
                                          HazardRates = HazardRates, 
                                          PopStructure = PopStructure)

## Plot results from fixed design LTRE
plotLTRE_fixedDesign(LTRE_results = fixedLTRE, 
                     Amax = Amax, Tmax = Tmax, minYear = minYear, 
                     HazardRates = HazardRates, 
                     PopStructure = PopStructure)
  
## Plot decomposition of mO into covariates and random effect
plotVariance_comp_mO(MCMC.samples = IPM.out, 
                     Tmax = Tmax,
                     minYear = minYear,
                     normN = normN)

## Calculate post-hoc parameter correlations to check for signs of density dependence
calculate_p.hoc_param.corr(MCMC.samples = IPM.out, 
                           Tmax = Tmax)
