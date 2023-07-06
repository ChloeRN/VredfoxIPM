library(ggplot2)
library(nimble)
library(sf)
library(reshape2)
library(remotes)
library(ckanr)
library(purrr)
library(dplyr)
library(metafor)

#**********#
# 0) SETUP #
#**********#

## Set seed
mySeed <- 0

## Set general parameters
Amax <- 5 # Number of age classes
Tmax <- 18  # Number of years
minYear <- 2004 # First year to consider
maxAge_yrs <- 10 # Age of the oldest female recorded
summer_removal <- c(6,7,8,9) #removal of summer months: numerical months to be removed from age at harvest data
area_selection<- c("Inner", "BB",  "Tana")# choosing varanger sub area ("Inner" / "BB" / "Tana)     ((BB = Batsfjord and Berlevag areas))
# start and end of placental scars and embryo sample periods (julian day)
plac_start <- 180 #including
plac_end   <- 80  #until, not including
embr_start <- 100 #including
embr_end   <- 140 #until, not including

## set dataset names, versions, and directories, and access
carcass.dataset.name <- "v_redfox_carcass_examination_v3"
carcass.dataset.version <- 3

rodent.dataset.name <-"v_rodents_snaptrapping_abundance_regional_v5"
rodent.dataset.version <- 5

# Stijn
shapefile.dir <- "C:\\Users\\sho189\\OneDrive - UiT Office 365\\PhD\\RedfoxIPM\\Fox areas shapefile\\tana rest"
COAT_key <- Sys.getenv("API_COAT_Stijn") # Stijn's API key for the COAT dataportal is saved as an environmental variable on the computer 

# Chloe
shapefile.dir <- "C:/Users/chloe.nater/OneDrive - NINA/Documents/Projects/RedFox_IPM/Data/shapefiles"
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
fitCov.Psi <- TRUE # Fit covariates on Psi (rodent abundance)
fitCov.rho <- TRUE # Fit covariates on rho (rodent abundance)
rCov.idx <- FALSE # Use discrete vs. continuous rodent covariate
nLevels.rCov <- 2 # 2-level discrete rodent covariate
#nLevels.rCov <- 3 # 3-level discrete rodent covariate (data not currently prepared)
standSpec.rCov <- TRUE # standardize different rodent species before summing (offset catchability) v.s. simply sum all numbers

# Annual survival prior type toggles
HoeningPrior <- FALSE # Use prior on natural mortality derived from Hoening model
#sPriorSource <- "Bristol" # Base survival prior on data from Bristol (not hunted)
sPriorSource <- "NSweden" # Base survival prior on data from North Sweden (lightly hunted)
#sPriorSource <- "metaAll" # Base survival prior on meta-analysis including all populations
#sPriorSource <- "metaSub" # Base survival prior on meta-analysis including only not/lightly hunted populations

# Genetic immigration data toggles (details in documentation of wrangleData_gen function
#GeneClass.approach <- 1 # Using first approach for GeneClass analysis 
GeneClass.approach <- 2 # Using second approach for GeneClass analysis
poolYrs.genData <- TRUE # Pool data across all years
imm.asRate <- TRUE # Estimating immigration as a rate as opposed to numbers
useData.gen <- TRUE # Use genetic data for estimation of immigration rate


#*********************#
# 1) DATA PREPARATION #
#*********************#

# 1a) Download and reformat carcass data
#-------------------------------#

## Download carcass data
carcass.data.raw <- downloadData_COAT(COAT_key = COAT_key, 
                                     COATdataset.name = carcass.dataset.name,
                                     COATdataset.version = carcass.dataset.version)

## Reformat carcass data
carcass.data <- reformatData_carcass(Amax = Amax,   
                                     summer_removal = summer_removal ,
                                     area_selection = area_selection,
                                     plac_start = plac_start,
                                     plac_end = plac_end ,
                                     embr_start = embr_start ,
                                     embr_end = embr_end,
                                     carcass.dataset = carcass.data.raw,
                                     shapefile.dir = shapefile.dir)

# 1b) Winter Age-at-Harvest data #
#--------------------------------#

## Set data path/filename
wAaH.datafile <- carcass.data$AaH.matrix

## Prepare winter AaH data
wAaH.data <- wrangleData_winterAaH(wAaH.datafile = wAaH.datafile, 
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

## Prepare genetic data
gen.data <- wrangleData_gen(datapath = genetics.datapath,
                            minYear, 
                            onlyFemales = FALSE, 
                            GeneClass.approach = GeneClass.approach, 
                            poolYrs.genData = poolYrs.genData)


# 1e) Harvest effort data #
#-------------------------#

## Prepare harvest effort data
hunter.data <- reformatData_hunters(area_selection = area_selection,
                                    carcass.dataset = carcass.data.raw,
                                    shapefile.dir = shapefile.dir)


# 1f) Environmental data #
#------------------------#

## Download rodent data
rodent.data.raw <- downloadData_COAT(COAT_key = COAT_key, 
                                     COATdataset.name = rodent.dataset.name,
                                     COATdataset.version = rodent.dataset.version)

## Reformat rodent data
rodent.data <- reformatData_rodent(rodent.dataset = rodent.data.raw,
                                          minYear = minYear)


# 1g) Conceptual year information #
#---------------------------------#

YearInfo <- collate_yearInfo(minYear = minYear,
                             Tmax = Tmax)


#**********************#
# 2) PRIOR INFORMATION #
#**********************#

## Parameters/paths for making informative priors for survival based on meta-analysis of literature data
meta.datafile <- "Data/RedFox_LiteratureData.csv"
simulateSD <- TRUE

## Parameters/paths for making informative priors for natural mortality using Tom Porteus' Hoening model approach
mu.t.max <- 22.61062
hoening.datafile <- "Data/HoenigMod_Posteriors_fromTomPorteus.txt"
nsim <- 30

## Collate all prior information
surv.priors <- collate_priorInfo(meta.datafile = meta.datafile,
                                 simulateSD = simulateSD,
                                 hoening.datafile = hoening.datafile, 
                                 nsim = nsim, 
                                 mu.t.max = mu.t.max, 
                                 maxAge = maxAge_yrs)

## Define type of prior to use for annual survival
survPriorType <- definePriorType_AnnSurv(HoeningPrior = HoeningPrior, 
                                         sPriorSource = sPriorSource)

#****************#
# 3) MODEL SETUP #
#****************#

# 3a) Write model code #
#----------------------#

redfox.code <- writeCode_redfoxIPM()


# 3b) Assemble IPM input data #
#-----------------------------#

input.data <- assemble_inputData(Amax = Amax, 
                                 Tmax = Tmax, 
                                 minYear = minYear,
                                 maxPups = 14,
                                 uLim.N = 800,
                                 uLim.Imm = 3000,
                                 nLevels.rCov = nLevels.rCov,
                                 standSpec.rCov = standSpec.rCov,
                                 poolYrs.genData = poolYrs.genData,
                                 wAaH.data = wAaH.data, 
                                 rep.data = rep.data, 
                                 gen.data = gen.data,
                                 rodent.data = rodent.data, 
                                 hunter.data = hunter.data, 
                                 surv.priors = surv.priors,
                                 survPriorType = survPriorType)


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
                          fitCov.Psi = fitCov.Psi, 
                          fitCov.rho = fitCov.rho,
                          rCov.idx = rCov.idx, 
                          HoeningPrior = HoeningPrior,
                          testRun = FALSE,
                          initVals.seed = mySeed)


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
                      setSeed = 0)
Sys.time() - t1

saveRDS(IPM.out, file = "immR_rodentsEff&mO_varT.rds")
MCMCvis::MCMCtrace(IPM.out)


#######################
# 5) MODEL COMPARISON #
#######################

compareModels(Amax = Amax, 
              Tmax = Tmax, 
              minYear = minYear, 
              post.filepaths = c("ImmNum_naive.rds", "ImmRate_naive.rds", "ImmRate_genData.rds"), 
              model.names = c("Number, naive", "Rate, naive", "Rate, pooled gen data"), 
              plotFolder = "Plots/Comp_ImmModels")

