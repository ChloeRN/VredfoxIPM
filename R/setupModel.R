#' Set up model, data, and initial values for running MCMC
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
#' @param fitCov.Psi logical. If TRUE, sets up model including covariate
#' effects on pregnancy rates. 
#' @param rCov.idx logical. Only required if fitCov.Psi = TRUE. If TRUE, assumes
#' a categorical rodent abundance covariate. If FALSE, assumes a continuous rodent
#' abundance covariate.
#' @param HoeningPrior logical. If TRUE, sets up a model using informative natural 
#' mortality priors based on the Hoening model. If FALSE, sets up a model using 
#' informative survival priors based on literature. 
#' @param niter integer. Number of MCMC iterations (default = 30000)
#' @param nthin integer. Thinning factor (default = 4)
#' @param nburn integer. Number of iterations to discard as burn-in (default = 5000)
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

setupModel <- function(modelCode,
                       nim.data, nim.constants,
                       minN1, maxN1, minImm, maxImm,
                       fitCov.mH, fitCov.Psi, rCov.idx, HoeningPrior,
                       niter = 30000, nthin = 4, nburn = 5000, nchains = 3,
                       testRun = FALSE, initVals.seed){
  
  
  ## Set parameters to monitor in all model versions
  params <- c("Mu.mH", "Mu.mO", "Mu.Psi", "Mu.rho", "Mu.S0",
              "sigma.mH", "sigma.Psi", "sigma.rho",
              "Psi", "rho", "mH",
              "initN",
              "N.tot", "B.tot", "R.tot", 
              "N", "B", "L", "Imm")
  
  ## Add additional parameters to monitor depending on model version
  if(HoeningPrior){
    params <- c(params, "JuvAdRatio", "Mu.mO.ad")
  }else{
   params <- c(params, "Mu.Snat")
  }
  
  if(fitCov.mH){
    params <- c(params, "betaHE.mH")
  }
  
  if(fitCov.Psi){
    params <- c(params, "betaR.Psi")
  }
  
  ## Simulate initial values
  set.seed(initVals.seed)
  initVals <- list()
  for(c in 1:nchains){
    initVals[[c]] <- simulateInitVals(nim.data = nim.data, 
                                      nim.constants = nim.constants, 
                                      minN1 = minN1, maxN1 = maxN1, 
                                      minImm = minImm, maxImm = maxImm, 
                                      fitCov.mH = fitCov.mH, 
                                      fitCov.Psi = fitCov.Psi, 
                                      rCov.idx = rCov.idx, 
                                      HoeningPrior = HoeningPrior)
  }
  
  ## Adjust MCMC parameters if doing a test run
  if(testRun){
    niter <- 20
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