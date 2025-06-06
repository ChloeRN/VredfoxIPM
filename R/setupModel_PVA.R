#' Set up model, data, and initial values for running MCMC (PVA-IPM)
#'
#' @param modelCode string. Relative path to the model file to be used
#' @param nim.data list of input objects representing data
#' @param nim.constants list of input objects representing constants.
#' @param minN1 integer vector. Lower bound for initial population size 
#' (per age class) to use in initial value simulation.
#' @param maxN1 integer vector. Upper bound for initial population size 
#' (per age class) to use in initial value simulation. 
#' @param minImm integer. Lower bound for the annual number of immigrants to 
#' use in initial value simulation.  
#' @param maxImm integer. Upper bound for the annual number of immigrants to 
#' use in initial value simulation. 
#' @param fitCov.mH logical. If TRUE, sets up model including covariate
#' effects on harvest mortality.
#' @param fitCov.mO logical. If TRUE, sets up model including covariate
#' effects on natural mortality.
#' @param fitCov.Psi logical. If TRUE, sets up model including covariate
#' effects on pregnancy rates. 
#' @param fitCov.rho logical. If TRUE, sets up model including covariate
#' effects on litter size. 
#' @param fitCov.immR logical. If TRUE, sets up model including covariate
#' effects on limmigration rate. 
#' @param rCov.idx logical. Only required if fitCov.Psi = TRUE. If TRUE, assumes
#' a categorical rodent abundance covariate. If FALSE, assumes a continuous rodent
#' abundance covariate.
#' @param mO.varT logical. If TRUE, sets up model with time variation in 
#' natural mortality. 
#' @param HoenigPrior logical. If TRUE, sets up a model using informative natural 
#' mortality priors based on the Hoenig model. If FALSE, sets up a model using 
#' informative survival priors based on literature. 
#' @param imm.asRate logical. If TRUE, sets up a model that estimates immigration
#' as a rate. 
#' @param Mu.mO_fixInits logical. If TRUE (default), sets initial values for
#' age-specific average natural mortality hazard rates to pre-defined values
#' taken from the North Sweden red fox population as presented in Devenish-Nelson
#' et al. 2017. Using these values seems to produce good sets of initial values
#' for the entire model. If set to FALSE, initial values for Mu.mO parameters
#' are instead simulated fron uniform distributions. 
#' @param DD.mO logical. If TRUE, models natural mortality (of age class 1) as
#' density-dependent. 
#' @param DD.immR logical. If TRUE, models immigration rate as density-dependent.
#' @param DDxRodent logical. If TRUE and if any density-dependence is included,
#' models both a direct density effect and a rodent x density interaction on the
#' relevant vital rate(s).
#' @param comp.mO logical. If TRUE, models natural mortality as potentially
#' compensatory for harvest (correlation with harvest).
#' @param comp.immR logical. If TRUE, models immigration rate as potentially
#' compensatory for harvest (correlation with harvest).
#' @param comp.RE logical. If TRUE and if any compensation is included, models
#' compensatory mechanisms via correlated random effects (after the approximation
#' approach from Nater et al. 2020). If FALSE and if any compensation is included,
#' models compensatory mechanisms as an effect of log-deviance of annual harvest
#' mortality from average harvest mortality (= harvest mortality RE plus any 
#' potential covariate effects on harvest mortality).
#' @param niter integer. Number of MCMC iterations (default = 100000)
#' @param nthin integer. Thinning factor (default = 8)
#' @param nburn integer. Number of iterations to discard as burn-in (default = 37500)
#' @param nchains integer. Number of chains to run (default = 3).
#' @param testRun logical. If TRUE, sets up for a test run with 10 iterations,
#' no thinning, and no burn-in (default = FALSE)
#' @param initVals.seed integer. Seed to use for inital value simulation.
#'
#' @return list of list containing all components necessary for running model 
#' with `nimble::nimbleMCMC()`
#' @export
#'
#' @examples

setupModel_PVA <- function(modelCode,
                           nim.data, nim.constants,
                           minN1, maxN1, minImm, maxImm,
                           fitCov.mH, fitCov.mO, fitCov.Psi, fitCov.rho, fitCov.immR, rCov.idx, 
                           mO.varT, HoenigPrior, imm.asRate, Mu.mO_fixInits = TRUE,
                           DD.mO, DD.immR, DDxRodent,
                           comp.mO, comp.immR, comp.RE,
                           niter = 200000, nthin = 20, nburn = 75000, nchains = 4,
                           testRun = FALSE, initVals.seed){
  
  
  ## Set parameters to monitor in all model versions
  params <- c("Mu.mH", "Mu.mO", "Mu.Psi", "Mu.rho", "Mu.S0",
              "sigma.mH", "sigma.mO", "sigma.Psi", "sigma.rho",
              "Psi", "rho", "mH", "mO", "S",
              "initN",
              "N.tot", "B.tot", "R.tot", 
              "N", "B", "L", "R", "Imm", "immR",
              "beta.RodMod", "beta.RodCorr", "sigmaT.RodAbun", "sigmaT.RodAbun2",
              "RodentAbundance", "RodentAbundance2", "pertFac.mH.flex")
  
  ## Add additional parameters to monitor depending on model version
  if(HoenigPrior){
    params <- c(params, "JuvAdRatio", "Mu.mO.ad")
  }else{
    params <- c(params, "Mu.Snat")
  }
  
  if(fitCov.mH){
    params <- c(params, "betaHE.mH", "HarvestEffort")
  }
  
  if(fitCov.mO){
    params <- c(params, "betaR.mO")
    
    if(DD.mO){
      params <- c(params, "betaD.mO")
      if(DDxRodent){
        params <- c(params, "betaRxD.mO")
      }
    }
    
    if(comp.mO & !comp.RE){
      params <- c(params, "gamma.mO")
    }
  }
  
  if(comp.mO & comp.RE){
    params <- c(params, "tau.mO", "C.mO")
  }
  
  if(fitCov.Psi){
    params <- c(params, "betaR.Psi")
  }
  
  if(fitCov.rho){
    params <- c(params, "betaR.rho")
  }
  
  if(imm.asRate){
    params <- c(params, "Mu.immR", "sigma.immR", "tau.immR", "C.immR")
    
    if(!poolYrs.genData){
      params <- c(params, "immR_pre")
    }
  }else{
    params <- c(params, "Mu.Imm", "Mu.immR", "sigma.immR", "tau.immR", "C.immR")
  } 
  
  if(fitCov.immR){
    params <- c(params, "betaR.immR", "RodentAbundance2")
    
    if(DD.immR){
      params <- c(params, "betaD.immR")
      if(DDxRodent){
        params <- c(params, "betaRxD.immR")
      }
    }
    
    if(comp.immR & !comp.RE){
      params <- c(params, "gamma.immR")
    }
  }
  
  if(comp.immR & comp.RE){
    params <- c(params, "tau.immR", "C.immR")
  }
  
  ## Simulate initial values
  set.seed(initVals.seed)
  initVals <- list()
  for(c in 1:nchains){
    initVals[[c]] <- simulateInitVals_PVA(nim.data = nim.data, 
                                          nim.constants = nim.constants, 
                                          minN1 = minN1, maxN1 = maxN1, 
                                          minImm = minImm, maxImm = maxImm, 
                                          fitCov.mH = fitCov.mH, 
                                          fitCov.mO = fitCov.mO,
                                          fitCov.Psi = fitCov.Psi, 
                                          fitCov.rho = fitCov.rho, 
                                          fitCov.immR = fitCov.immR,
                                          rCov.idx = rCov.idx, 
                                          mO.varT = mO.varT,
                                          HoenigPrior = HoenigPrior,
                                          imm.asRate = imm.asRate,
                                          Mu.mO_fixInits = Mu.mO_fixInits)
  }
  
  ## Adjust MCMC parameters if doing a test run
  if(testRun){
    niter <- 50
    nthin <- 1
    nburn <- 0
  }
  
  ## Collate model setup variables in a list
  setup <- list(
    modelCode = modelCode,
    modelParams = params,
    initVals = initVals,
    mcmcParams = list(niter = niter, nthin = nthin, 
                      nburn = nburn, nchains = nchains)
  )
  
  return(setup) 
}
