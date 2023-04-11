library(ggplot2)
library(nimble)
library(sf)
library(reshape2)

#**********#
# 0) SETUP #
#**********#

## Set seed
mySeed <- 0

## Set general parameters
Amax <- 5 # Number of age classes
Tmax <- 15  # Number of years
minYear <- 2004 # First year to consider

maxAge_yrs <- 10 # Age of the oldest female recorded

#removal of summer months:
summer_removal <- c(6,7,8,9) #numerical months to be removed from age at harvest data
# choosing varanger sub area ("Inner" / "BB" / "Tana)
area_selection<- c("Inner", "BB",  "Tana")
# start and end of placental scars and embryo sample periods (julian day)
plac_start <- 140 #including
plac_end   <- 80  #until, not including
embr_start <- 100 #including
embr_end   <- 140 #until, not including
## set directories
carcass.dir        <-"C:\\Users\\sho189\\OneDrive - UiT Office 365\\PhD\\RedfoxIPM\\Data from google disk\\carcass_examination"
shapefile.dir      <-"C:\\Users\\sho189\\OneDrive - UiT Office 365\\PhD\\RedfoxIPM\\Fox areas shapefile\\tana rest"

  
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
fitCov.Psi <- FALSE # Fit covariates on Psi (rodent abundance)
rCov.idx <- TRUE # Use discrete vs. continuous rodent covariate

# Annual survival prior type toggles
HoeningPrior <- FALSE # Use prior on natural mortality derived from Hoening model
sPriorSource <- "Bristol" # Base survival prior on data from Bristol (not hunted)
#sPriorSource <- "NSweden" # Base survival prior on data from North Sweden (lightly hunted)
#sPriorSource <- "metaAll" # Base survival prior on meta-analysis including all populations
#sPriorSource <- "metaSub" # Base survival prior on meta-analysis including only not/lightly hunted populations

#*********************#
# 1) DATA PREPARATION #
#*********************#

# 1aa) reformatting carcass data
#-------------------------------#
carcassData<- data_reformatting_carcass (
  Amax               = Amax,   
  summer_removal     = c(6,7,8,9) ,
  area_selection     = c("Inner", "BB",  "Tana"),
  plac_start         = 140,
  plac_end           = 80 ,
  embr_start         = 100 ,
  embr_end           = 140,
  carcass.dir        ="C:\\Users\\sho189\\OneDrive - UiT Office 365\\PhD\\RedfoxIPM\\Data from google disk\\carcass_examination",
  shapefile.dir      ="C:\\Users\\sho189\\OneDrive - UiT Office 365\\PhD\\RedfoxIPM\\Fox areas shapefile\\tana rest"
)

# 1a) Winter Age-at-Harvest data #
#--------------------------------#

## Set data path/filename
wAaH.datafile <- carcassData$AaH.matrix

## Prepare winter AaH data
wAaH.data <- wrangleData_winterAaH(wAaH.datafile = wAaH.datafile, 
                                   Amax = Amax)

# 1b) Reproduction data #
#-----------------------#

## Set data paths/filenames
P1.datafile <- carcassData$P1var # Placental scar/embryo count
P2.datafile <- carcassData$P2var # Presence of placental scars/embryos/pregnancy signs

## Prepare reproduction data
rep.data <- wrangleData_rep(P1.datafile  = P1.datafile, 
                            P2.datafile  = P2.datafile,
                            Amax = Amax, 
                            minYear = minYear)


# 1c) Environmental data #
#------------------------#

## Set data paths/filenames
rodent.datafile <- "Data/stor_intensiv_04_20-year-var.txt"

## Prepare rodent abundance data
rodent.data <- wrangleData_rodent(datafile = rodent.datafile,
                                  minYear = minYear,
                                  adjust = TRUE)

## Prepare harvest effort data
hunter.data <- data_reformatting_hunters(area_selection = area_selection,
                          carcass.dir = carcass.dir,
                          shapefile.dir = shapefile.dir)

#TODO: Add option where hunting effort not taken into account in analysis (assumed equal)

# 1d) Conceptual year information #
#---------------------------------#

YearInfo <- collate_yearInfo(minYear = minYear,
                             Tmax = Tmax)


#**********************#
# 2) PRIOR INFORMATION #
#**********************#

## Make informative priors for natural mortality using Tom Porteus' Hoening model approach
mu.t.max <- 22.61062
hoening.datafile <- "Data/HoenigMod_Posteriors_fromTomPorteus.txt"
nsim <- 30


## Collate all prior information
surv.priors <- collate_priorInfo(datafile = hoening.datafile, 
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
                                 uLim.Imm = 800,
                                 #nLevels.rCov = 2,
                                 wAaH.data = wAaH.data, 
                                 rep.data = rep.data, 
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
                          rCov.idx = rCov.idx, 
                          HoeningPrior = HoeningPrior,
                          testRun = TRUE,
                          initVals.seed = mySeed)

####################
# 4) MODEL FITTING #
####################

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

