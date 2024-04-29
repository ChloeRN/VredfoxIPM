library(targets)
library(tidyverse)
library(nimble)

tar_load(rodent.data)

Tmax <- length(rodent.data$cont.wintvar)

#-------------------------------------#
# EXPLORATORY FREQUENTIST LM ANALYSIS #
#-------------------------------------#

## Assemble rodent data (including lags) in a data frame
covData <- data.frame(Index = 1:Tmax,
                      wintvar = rodent.data$cont.wintvar.stsp,
                      wintvar_lag1 = c(NA, rodent.data$cont.wintvar.stsp[1:(Tmax-1)]),
                      wintvar_lag2 = c(NA, NA, rodent.data$cont.wintvar.stsp[1:(Tmax-2)]),
                      fallstor = c(NA, rodent.data$cont.fallstor.stsp))

## Test simple linear models including time-lags
fullMod_wintvar <- lm(wintvar ~ wintvar_lag1 + wintvar_lag2 + wintvar_lag1:wintvar_lag2, data = covData)
autoregMod_wintvar <- lm(wintvar ~ wintvar_lag1, data = covData)

summary(fullMod_wintvar)
summary(autoregMod_wintvar)

AIC(fullMod_wintvar)
AIC(autoregMod_wintvar)

BIC(fullMod_wintvar)
BIC(autoregMod_wintvar)
# --> The model using both 1- and 2-year time-lags appears to be substantially better

## Test correlation between fall and winter variables
cor.test(covData$wintvar, covData$fallstor)
# --> Correlation is fairly decent


#-------------------#
# BAYESIAN ANALYSIS #
#-------------------#

## Assemble data and constants (add 2 "dummy years")
nim.data <- list(wintvar = c(rep(NA, 2), rodent.data$cont.wintvar.stsp),
                 fallstor = c(rep(NA, 3), rodent.data$cont.fallstor.stsp))

nim.constants <- list(Tmax = Tmax + 2)

## Write model code
rodMod.code <- nimbleCode({
  
  # Impute missing values
  for(t in 1:2){
    wintvar[t] ~ dnorm(mean = 0, sd = 1)
  }
  
  # Linear model for rodent dynamics in Varanger
  for(t in 3:Tmax){
    wintvar_pred[t] <- beta.RodMod[1]*wintvar[t-1] + beta.RodMod[2]*wintvar[t-2] + beta.RodMod[3]*wintvar[t-1]*wintvar[t-2]
    wintvar[t] ~ dnorm(mean = wintvar_pred[t], sd = sigmaT.wintvar)
  }
  
  # Correlative model for rodent dynamics in the larger area
  for(t in 1:Tmax){
    fallstor[t] ~ dnorm(mean = beta.RodCorr*wintvar[t], sd = sigmaT.fallstor)
  }

  # Priors
  for(i in 1:3){
    beta.RodMod[i] ~ dunif(-5, 5)
  }
  beta.RodCorr ~ dunif(0, 1)
  sigmaT.wintvar ~ dunif(0, 5)
  sigmaT.fallstor ~ dunif(0, 5)
})

## Function for simulating initial values
sim_Inits <- function(Tmax){
  
  wintvar <- fallstor <- rep(NA, Tmax)
  wintvar[c(1:2, Tmax)] <- rnorm(3, 0, 1)
  fallstor[1:3] <- rnorm(3, 0, 1)
  
  beta.RodMod <- runif(3, -1, 1)
  beta.RodCorr <- runif(1, 0, 1)
  sigmaT.wintvar <- runif(1, 0, 0.2)
  sigmaT.fallstor <- runif(1, 0, 0.2)
  
  return(list(wintvar = wintvar,
              fallstor = fallstor,
              beta.RodMod = beta.RodMod,
              beta.RodCorr = beta.RodCorr,
              sigmaT.wintvar = sigmaT.wintvar,
              sigmaT.fallstor = sigmaT.fallstor))
}

## Set up for running model
params <- c("wintvar", "fallstor", 
            "beta.RodMod", "beta.RodCorr", 
            "sigmaT.wintvar", "sigmaT.fallstor")

Inits <- list(sim_Inits(Tmax = nim.constants$Tmax),
              sim_Inits(Tmax = nim.constants$Tmax),
              sim_Inits(Tmax = nim.constants$Tmax))

## Run model
modelTest <- nimbleMCMC(code = rodMod.code,
                        data = nim.data, 
                        constants = nim.constants,
                        inits = Inits, 
                        monitors = params,
                        nchains = 3, 
                        niter = 10000, 
                        nburnin = 3000, 
                        thin = 1, 
                        samplesAsCodaMCMC = TRUE)

## Check results
plot(modelTest, ask = TRUE)
# --> Estimates look good

