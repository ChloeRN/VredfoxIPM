library(ggplot2)
library(nimble)
library(tidyverse)

## Set toggle combos for models to run
cov_toggles <- list(
  fitCov.mH = rep(FALSE, 4),
  fitCov.Psi = c(FALSE, rep(TRUE, 3)),
  fitCov.rho = c(FALSE, rep(TRUE, 3)),
  rCov.idx = rep(c(FALSE, TRUE), each = 2),
  nLevels.rCov = c(rep(2, 3), 3)
)

## List model names
model_names <- c("Base",
                 "rodCont",
                 "rodIdx2",
                 "rodIdx3")


## Altered function for initial value simulation
simulateInitVals_repModels <- function(nim.data, nim.constants, fitCov.Psi, fitCov.rho, rCov.idx){
  
  Amax <- nim.constants$Amax
  Tmax <- nim.constants$Tmax
  
  #---------------------------------------------------#
  # Set initial values for vital rate base parameters #
  #---------------------------------------------------#
  
  ## Pregnancy rate
  Mu.Psi <- c(0, runif(Amax-1, 0.2, 0.8))
  
  ## Placental scars
  Mu.rho <- c(0, runif(Amax-1, 3, 6))
  
  ## Random effect standard deviations
  sigma.Psi <- runif(1, 0.05, 0.5)
  sigma.rho <- runif(1, 0.05, 0.5)
  
  ## Random effects (initialize at to 0)
  epsilon.Psi <- rep(0, Tmax+1)
  epsilon.rho <- rep(0, Tmax+1)
  
  ## Covariate effects
  
  # Rodent abundance on Psi
  if(fitCov.Psi){
    if(rCov.idx){
      betaR.Psi <- rep(0, nLevels.rCov)
      for(x in 2:nim.constants$nLevels.rCov){
        betaR.Psi[x] <- runif(1, -5, 5)
      }
    }else{
      betaR.Psi <- runif(1, 0, 2)
    }
  }else{
    betaR.Psi <- 0
  }
  
  
  # Rodent abundance on rho
  if(fitCov.rho){
    if(rCov.idx){
      betaR.rho <- rep(0, nLevels.rCov)
      for(x in 2:nim.constants$nLevels.rCov){
        betaR.rho[x] <- runif(1, -5, 5)
      }
    }else{
      betaR.rho <- runif(1, 0, 2)
    }
  }else{
    betaR.rho <- 0
  }
  
  #-------------------------------------#
  # Calculate year-specific vital rates #
  #-------------------------------------#
  
  Psi <- rho <- matrix(NA, nrow = Amax, ncol = Tmax + 1)
  
  for(t in 1:(Tmax+1)){
    
    
    ## Pregnancy rate
    Psi[1, t] <- 0
    
    if(fitCov.Psi & rCov.idx){
      Psi[2:Amax, t] <- plogis(qlogis(Mu.Psi[2:Amax]) + betaR.Psi[nim.constants$RodentIndex[t]+1] + epsilon.Psi[t])
    }else{
      Psi[2:Amax, t] <- plogis(qlogis(Mu.Psi[2:Amax]) + betaR.Psi*nim.constants$RodentAbundance[t] + epsilon.Psi[t])
    }
    
    ## Placental scars
    rho[1, t] <- 0
    
    if(fitCov.rho & rCov.idx){
      rho[2:Amax, t] <- exp(log(Mu.rho[2:Amax]) + betaR.rho[nim.constants$RodentIndex[t]+1] + epsilon.rho[t])
    }else{
      rho[2:Amax, t] <- exp(log(Mu.rho[2:Amax]) + betaR.rho*nim.constants$RodentAbundance[t] + epsilon.rho[t])
    }
  }


## List all initial values
InitVals <- list(
  Mu.Psi = Mu.Psi,
  Mu.rho = Mu.Psi,
  
  sigma.Psi = sigma.Psi,
  sigma.rho = sigma.rho,
  
  epsilon.Psi = epsilon.Psi,
  epsilon.rho = epsilon.Psi,
  
  betaR.Psi = betaR.Psi,
  betaR.rho = betaR.rho,
  
  Psi = Psi, 
  rho = rho
)


## Return initial values
return(InitVals)
}

## Set general parameters
Amax <- 5 # Number of age classes
Tmax <- 15  # Number of years
minYear <- 2004 # First year to consider

maxAge_yrs <- 10 # Age of the oldest female recorded

## Source all functions in "R" folder
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir('R')


## Run model variants
for(i in 1:4){
  
  #**********#
  # 0) SETUP #
  #**********#
  
  ## Set seed
  mySeed <- 0
  
  ## Set "switches" for running different model versions
  
  # Covariate toggles
  fitCov.Psi <- cov_toggles$fitCov.Psi[i]
  fitCov.rho <- cov_toggles$fitCov.rho[i]
  rCov.idx <- cov_toggles$rCov.idx[i]
  nLevels.rCov <- cov_toggles$nLevels.rCov[i]
  
  
  #*********************#
  # 1) DATA PREPARATION #
  #*********************#
  #*
  
  # 1b) Reproduction data #
  #-----------------------#
  
  ## Set data paths/filenames
  rep.datafiles <- c("Data/P1var_tot.txt", # Placental scar/embryo count
                     "Data/P2var_tot.txt") # Presence of placental scars/embryos/pregnancy signs
  
  ## Prepare reproduction data
  rep.data <- wrangleData_rep(datafiles = rep.datafiles, 
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
  
 
  # 1d) Conceptual year information #
  #---------------------------------#
  
  YearInfo <- collate_yearInfo(minYear = minYear,
                               Tmax = Tmax)
  
  
  #**********************#
  # 2) PRIOR INFORMATION #
  #**********************#

  
  #****************#
  # 3) MODEL SETUP #
  #****************#
  
  # 3a) Write model code #
  #----------------------#
  
  redfox.code <- writeCode_repModels()
  
  
  # 3b) Assemble IPM input data #
  #-----------------------------#

  # Reproduction data
  P1 <- subset(rep.data$P1, repryear %in% c(minYear + 0:Tmax))
  P2 <- subset(rep.data$P2, repryear %in% c(minYear + 0:Tmax))
  
  
  ## Select relevant categorical rodent covariate
  if(is.na(nLevels.rCov)){
    RodentIndex <- NA
  }else{
    
    if(nLevels.rCov == 2){
      RodentIndex <- rodent.data$cat2
    }else{
      RodentIndex <- rodent.data$cat3
    }
  }
  
  input.data <- list(
    
    nim.data = list(
      P1 = P1$P1,
      P2 = P2$P2
    ),
    
    nim.constants = list(
      Amax = Amax,
      Tmax = Tmax,
      minYear = minYear,
      
      maxPups = 14,
      
      P1_age = P1$age_adj,
      P1_year = P1$RepYearIndex,
      X1 = length(P1$P1),
      
      P2_age = P2$age_adj,
      P2_year = P2$RepYearIndex,
      X2 = length(P2$P2),
      
      RodentAbundance = rodent.data$cont,
      RodentIndex = RodentIndex,
      nLevels.rCov = nLevels.rCov
    )
  )
  
  # 3c) Set up for model run (incl. simulating initial values) #
  #------------------------------------------------------------#
  
  ## Simulate initial values
  set.seed(mySeed)

  initVals <- list()
  for(x in 1:3){
    initVals[[x]] <- simulateInitVals_repModels(nim.data = input.data$nim.data, 
                                                nim.constants = input.data$nim.constants, 
                                                fitCov.Psi = fitCov.Psi,
                                                fitCov.rho = fitCov.rho,
                                                rCov.idx = rCov.idx)

  }
  
  ## Set parameters to monitor
  params <- c("Mu.Psi", "Mu.rho", "sigma.Psi", "sigma.rho", "Psi", "rho")
  
  if(fitCov.Psi){
    params <- c(params, "betaR.Psi")
  }
  
  if(fitCov.rho){
    params <- c(params, "betaR.rho")
  }
  
  ####################
  # 4) MODEL FITTING #
  ####################
  
  IPM.out <- nimbleMCMC(code = redfox.code,
                        data = input.data$nim.data, 
                        constants = input.data$nim.constants,
                        inits = initVals, 
                        monitors = params,
                        nchains = 3, 
                        niter = 10000, 
                        nburnin = 2000, 
                        thin = 2, 
                        samplesAsCodaMCMC = TRUE, 
                        setSeed = 0)
  
  saveRDS(IPM.out, file = paste0("RepModels_", model_names[i], ".rds"))
  
}


#######################
# 5) MODEL COMPARISON #
#######################


compareModels(Amax = Amax, 
              Tmax = Tmax, 
              minYear = minYear,
              post.filepaths = paste0("RepModels_", model_names, ".rds"), 
              model.names = model_names, 
              plotFolder = "Plots/Comp_RepModels")


#######################################
# EXTRA) FREQUENTIST MODEL COMPARISON #
#######################################

library(lme4)

## Rodent data
rod.data <- data.frame(
  Year = 1:length(rodent.data$cont),
  cont = rodent.data$cont,
  idx2 = rodent.data$cat2,
  idx3 = rodent.data$cat3 
)

## Reproduction proportions
propRep.data <- data.frame(
  P2 = input.data$nim.data$P2,
  Age = input.data$nim.constants$P2_age,
  Year = input.data$nim.constants$P2_year
) %>%
  dplyr::left_join(rod.data, by = "Year")

## Number of placental scars
placScar.data <- data.frame(
  P1 = input.data$nim.data$P1,
  Age = input.data$nim.constants$P1_age,
  Year = input.data$nim.constants$P1_year
) %>%
  dplyr::left_join(rod.data, by = "Year")

## Fit alternative models and compare AIC/BIC/residual variation

# Reproduction proportions
modP2.cont <- glmer(P2 ~ Age + cont + (1|Year), family = binomial(link = "logit"), data = propRep.data)
modP2.idx2 <- glmer(P2 ~ Age + idx2 + (1|Year), family = binomial(link = "logit"), data = propRep.data)
modP2.idx3 <- glmer(P2 ~ Age + idx3 + (1|Year), family = binomial(link = "logit"), data = propRep.data)

AIC(modP2.cont)
AIC(modP2.idx2)
AIC(modP2.idx3)

BIC(modP2.cont)
BIC(modP2.idx2)
BIC(modP2.idx3)

sd(resid(modP2.cont))
sd(resid(modP2.idx2))
sd(resid(modP2.idx3))

# Placental scars
modP1.cont <- glmer(P1 ~ Age + cont + (1|Year), family = poisson(link = "log"), data = placScar.data)
modP1.idx2 <- glmer(P1 ~ Age + idx2 + (1|Year), family = poisson(link = "log"), data = placScar.data)
modP1.idx3 <- glmer(P1 ~ Age + idx3 + (1|Year), family = poisson(link = "log"), data = placScar.data)

AIC(modP1.cont)
AIC(modP1.idx2)
AIC(modP1.idx3)

BIC(modP1.cont)
BIC(modP1.idx2)
BIC(modP1.idx3)

sd(resid(modP1.cont))
sd(resid(modP1.idx2))
sd(resid(modP1.idx3))

