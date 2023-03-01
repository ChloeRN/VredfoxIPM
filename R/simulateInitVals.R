#' Simulate a complete set of initial values
#'
#' @param nim.data list of input objects representing data
#' @param nim.constants list of input objects representing constants
#' @param minN1 integer vector. Lower bound for initial population size (per age class).
#' @param maxN1 integer vector. Upper bound for initial population size (per age class).
#' @param minImm integer. Lower bound for the annual number of immigrants. 
#' @param maxImm integer. Upper bound for the annual number of immigrants.
#' @param fitCov.mH logical. If TRUE, simulates initial values including covariate
#' effects on harvest mortality.
#' @param fitCov.Psi logical. If TRUE, simulates initial values including covariate
#' effects on pregnancy rates. 
#' @param rCov.idx logical. Only required if fitCov.Psi = TRUE. If TRUE, assumes
#' a categorical rodent abundance covariate. If FALSE, assumes a continuous rodent
#' abundance covariate.
#' @param HoeningPrior logical. If TRUE, simulates initial values for a model 
#' using informative natural mortality priors based on the Hoening model. If 
#' FALSE, simulates initial values for a model using informative survival priors
#' based on literature. 
#'
#' @return a list containing a complete set of initial values for all parameters
#' in the IPM. 
#' @export
#'
#' @examples

simulateInitVals <- function(nim.data, nim.constants, minN1, maxN1, minImm, maxImm, fitCov.mH, fitCov.Psi, rCov.idx, HoeningPrior){
  
  Amax <- nim.constants$Amax
  Tmax <- nim.constants$Tmax
  
  #-------------------------------------------------#
  # Set initial values for missing covariate values #
  #-------------------------------------------------#
  
  ## Number of successful hunters
  NHunters <- nim.data$HarvestEffort
  if(NA %in% NHunters){
    NHunters[which(is.na(NHunters))] <- mean(NHunters, na.rm = TRUE)
  }
  
  #---------------------------------------------------#
  # Set initial values for vital rate base parameters #
  #---------------------------------------------------#
  
  ## Harvest and natural mortality
  Mu.mH <- runif(Amax, 0.05, 0.2)
  
  if(HoeningPrior){
    
    JuvAdRatio <- exp(nim.constants$ratioJA.logmean)
    Mu.mO.ad <- exp(nim.constants$mnat.logmean)
    
    Mu.mO <- c(Mu.mO.ad*JuvAdRatio, rep(Mu.mO.ad, Amax-1))
    
  }else{
    Mu.Snat <- nim.constants$Snat.mean
    Mu.mO <- -log(Mu.Snat)
  }
  
  ## Annual survival
  Mu.S <- exp(-(Mu.mH + Mu.mO))
  
  ## Pregnancy rate
  Mu.Psi <- c(0, runif(Amax-1, 0.2, 0.8))
  
  ## Placental scars
  Mu.rho <- c(0, runif(Amax-1, 3, 6))
  
  ## Early survival
  Mu.S0 <- nim.constants$S0.mean
  
  ## Immigration
  Mu.Imm <- runif(1, minImm, maxImm)
  sigma.Imm <- runif(1, 10, 40)
  
  ## Random effect standard deviations
  sigma.mH <- runif(1, 0.05, 0.5)
  sigma.Psi <- runif(1, 0.05, 0.5)
  sigma.rho <- runif(1, 0.05, 0.5)
  
  ## Random effects (initialize at to 0)
  epsilon.mH <- rep(0, Tmax)
  epsilon.Psi <- rep(0, Tmax+1)
  epsilon.rho <- rep(0, Tmax+1)
  
  ## Covariate effects
  
  # Harvest effort on mH
  if(fitCov.mH){
    betaHE.mH <- runif(1, 0, 0.2)
  }else{
    betaHE.mH <- 0
  }
  
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

  
  #-------------------------------------#
  # Calculate year-specific vital rates #
  #-------------------------------------#
  
  mH <- mO <- matrix(NA, nrow = Amax, ncol = Tmax)
  Psi <- rho <- matrix(NA, nrow = Amax, ncol = Tmax + 1)
  S0 <- rep(NA, Tmax + 1)
  
  for(t in 1:(Tmax+1)){
    
    if(t <= Tmax){
      ## Harvest and natural mortality
      mH[1:Amax, t] <- exp(log(Mu.mH[1:Amax]) + betaHE.mH*NHunters[t] + epsilon.mH[t])
      
      ## Other (natural) mortality hazard rate
      mO[1:Amax, t] <- exp(log(Mu.mO[1:Amax]))
    }
    
    ## Pregnancy rate
    Psi[1, t] <- 0
    
    if(fitCov.Psi & rCov.idx){
      Psi[2:Amax, t] <- plogis(qlogis(Mu.Psi[2:Amax]) + betaR.Psi[nim.constants$RodentIndex[t]+1] + epsilon.Psi[t])
    }else{
      Psi[2:Amax, t] <- plogis(qlogis(Mu.Psi[2:Amax]) + betaR.Psi*nim.constants$RodentAbundance[t] + epsilon.Psi[t])
    }

    ## Placental scars
    rho[1, t] <- 0
    rho[2:Amax, t] <- exp(log(Mu.rho[2:Amax]) + epsilon.rho[t])
    
    ## Early survival
    S0[t] <- Mu.S0
  }

  ## Survival probability
  S <- exp(-(mH + mO))
  
  ## Proportion harvest mortality
  alpha <- mH/(mH + mO)
  
  ## Harvest rate
  h <- (1 - S)*alpha

  ## Immigration
  Imm <- round(truncnorm::rtruncnorm(Tmax+1, a = 0, b = maxImm, mean = Mu.Imm, sd = sigma.Imm))
  Imm[1] <- 0
  

  #---------------------------------------------------------#
  # Simulate initial values for population-level quantities #
  #---------------------------------------------------------#
  
  ## Prepare empty vectors and matrices
  H <- N <- B <- L <- R <- matrix(NA, nrow = Amax, ncol = Tmax)
  
  ## Set initial population sizes
  for(a in 1:Amax){
    N[a, 1] <- round(runif(1, minN1[a], maxN1[a]))
  }
  
  ## Set age class 0 reproductive contributions to 0
  B[1, 1:Tmax] <- L[1, 1:Tmax] <- R[1, 1:Tmax] <- 0
  
  
  for (t in 1:(Tmax-1)){

    # a) Project local survivors to the next year
    #---------------------------------------------
      
    ## Age classes 1 to 3 (indeces = 2, 3, 4): age classes 0, 1, and 2 survivors    
    for(a in 1:(Amax-2)){
      N[a+1, t+1] <- rbinom(1, size = N[a, t], prob = S[a, t])
    }			
      
    # Age class 4+ (index = 5): age class 3 and 4+ survivors
    N[Amax, t+1] <- rbinom(1, size = N[Amax-1, t] + N[Amax, t], prob = S[Amax, t])
      
      
    # b) Sample through reproductive season 
    #--------------------------------------
      
    for (a in 2:Amax){
        
      ## Breeding Population Size: Number of females that reproduce
      B[a, t+1] <- rbinom(1, size = N[a, t+1], prob = Psi[a, t+1])
        
      ## Litter Size (in utero): Number of pups produced by females of age class a
      L[a, t+1] <- rpois(1, lambda = B[a, t+1] * rho[a, t+1] * 0.5)
        
      ## Number Recruits: Number of pups surviving to emerge from the den
      R[a, t+1] <- rbinom(1, size = L[a, t+1], prob = S0[t+1])
    }
      
      
    # c) Add new recruits and immigrants
    #------------------------------------
      
    N[1, t+1] <- sum(R[1:Amax, t+1]) + Imm[t+1]
 }
  
  
  # d) Check for years with more harvests than alive individuals
  #-------------------------------------------------------------
  
  if(any(nim.data$C > N)){
    stop('Simulation resulted in less alive than harvested. Retry.')
  }
  
  # e) Assemble and return list of initial values
  #----------------------------------------------  
  
  ## Fill out NA values in B, L, and R
  for (a in 2:Amax){
    #B[a,1] <- rbinom(1, size = N[a, 1], prob = Psi[a, 1])
    #L[a,1] <- rpois(1, lambda = B[a, 1] * rho[a, 1] * 0.5)
    #R[a,1] <- rbinom(1, size = L[a, 1], prob = S0[t])
    B[a,1] <- 0
    L[a,1] <- 0
    R[a,1] <- 0
  }
  # NOTE: These nodes do not appear in the model and it therefore does not 
  #       matter what numbers they contain. Filling them in prevents a warning
  #       about NA nodes when building the model. 
  
  
  ## List all initial values
  InitVals <- list(
    N = N, 
    initN = N[, 1],
    B = B, 
    L = L,
    R = R,
    Imm = Imm,
    
    Mu.mH = Mu.mH,
    Mu.mO = Mu.mO,
    Mu.Psi = Mu.Psi,
    Mu.rho = Mu.Psi,
    Mu.S0 = Mu.S0,
    
    sigma.mH = sigma.mH,
    sigma.Psi = sigma.Psi,
    sigma.rho = sigma.rho,
    
    epsilon.mH = epsilon.mH,
    epsilon.Psi = epsilon.Psi,
    epsilon.rho = epsilon.Psi,
    
    betaHE.mH = betaHE.mH,
    betaR.Psi = betaR.Psi,
    
    mH = mH,
    mO = mO, 
    S = S,
    alpha = alpha, 
    h = h, 
    Psi = Psi, 
    rho = rho,
    S0 = S0
    
    #Mu.Imm = Mu.Imm,
    #sigma.Imm = sigma.Imm
  )
  
  ## Add initial values for parameters specific to survival prior model versions
  if(HoeningPrior){
    InitVals$JuvAdRatio <- JuvAdRatio
    InitVals$Mu.mO.ad <- Mu.mO.ad
  }else{
    InitVals$Mu.Snat <- Mu.Snat
  }
  
  ## Add initial values for missing covariate values (if applicable)
  if(fitCov.mH & (NA %in% nim.data$HarvestEffort)){
    Inits_NHunters <- rep(NA, length(NHunters))
    Inits_NHunters[which(is.na(nim.data$HarvestEffort))] <- NHunters[which(is.na(nim.data$HarvestEffort))]
    InitVals$HarvestEffort <- Inits_NHunters
  }
  
  ## Return initial values
  return(InitVals)
}
