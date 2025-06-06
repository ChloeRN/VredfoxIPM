#' Simulate a complete set of initial values for IPM-PVA
#'
#' @param nim.data list of input objects representing data
#' @param nim.constants list of input objects representing constants
#' @param minN1 integer vector. Lower bound for initial population size (per age class).
#' @param maxN1 integer vector. Upper bound for initial population size (per age class).
#' @param minImm integer. Lower bound for the annual number of immigrants. 
#' @param maxImm integer. Upper bound for the annual number of immigrants.
#' @param fitCov.mH logical. If TRUE, simulates initial values including covariate
#' effects on harvest mortality.
#' @param fitCov.mO logical. If TRUE, simulates initial values including covariate
#' effects on natural mortality.
#' @param fitCov.Psi logical. If TRUE, simulates initial values including covariate
#' effects on pregnancy rates. 
#' @param fitCov.rho logical. If TRUE, simulates initial values including covariate
#' effects on litter size. 
#' @param fitCov.immR logical. If TRUE, adds 0 inital values for covariate effects
#' on immigration rates (these effects are not currently formally incorporated into
#' the initial value simulation).
#' @param rCov.idx logical. Only required if fitCov.Psi = TRUE. If TRUE, assumes
#' a categorical rodent abundance covariate. If FALSE, assumes a continuous rodent
#' abundance covariate.
#' @param mO.varT logical. If TRUE, simulates initial values for a model with 
#' time variation in natural mortality.
#' @param HoenigPrior logical. If TRUE, simulates initial values for a model 
#' using informative natural mortality priors based on the Hoenig model. If 
#' FALSE, simulates initial values for a model using informative survival priors
#' based on literature. 
#' @param imm.asRate logical. If TRUE, returns initial values associated with 
#' immigration rate.
#' @param Mu.mO_fixInits logical. If TRUE (default), sets initial values for
#' age-specific average natural mortality hazard rates to pre-defined values
#' taken from the North Sweden red fox population as presented in Devenish-Nelson
#' et al. 2017. Using these values seems to produce good sets of initial values
#' for the entire model. If set to FALSE, initial values for Mu.mO parameters
#' are instead simulated from uniform distributions. 
#' 
#' @return a list containing a complete set of initial values for all parameters
#' in the IPM. 
#' @export
#'
#' @examples

simulateInitVals_PVA <- function(nim.data, nim.constants, minN1, maxN1, minImm, maxImm, 
                                 fitCov.mH, fitCov.mO, fitCov.Psi, fitCov.rho, fitCov.immR, rCov.idx, 
                                 mO.varT, HoenigPrior, imm.asRate, Mu.mO_fixInits = TRUE){
  
  Amax <- nim.constants$Amax
  Tmax <- nim.constants$Tmax + nim.constants$Tmax_sim
  
  pertFac.mH <- nim.data$pertFac.mH
  pertFac.mO <- nim.data$pertFac.mO
  pertFac.S0 <- nim.data$pertFac.S0
  pertFac.immR <- nim.data$pertFac.immR
  pertFac.rodent <- nim.data$pertFac.rodent
  
  if(any(c(pertFac.rodent) != 1)){
    warning("Initial value simulation for scenarios with perturbations to covariates (rodents) have not been tested throroughly yet and may not be implemented correctly.")
  }
  
  #-------------------------------------------------#
  # Set initial values for missing covariate values #
  #-------------------------------------------------#
  
  ## Number of successful hunters
  NHunters <- nim.data$HarvestEffort
  if(NA %in% NHunters){
    NHunters[which(is.na(NHunters))] <- mean(NHunters, na.rm = TRUE)
  }
  
  ## Rodent abundance (categorical)
  RodentIndex <- nim.data$RodentIndex
  if(NA %in% RodentIndex){
    RodentIndex[which(is.na(RodentIndex))] <- sample(1:nLevels.rCov, length(which(is.na(RodentIndex))), replace = TRUE)
  }
  
  RodentIndex2 <- nim.data$RodentIndex2
  if(NA %in% RodentIndex2){
    RodentIndex2[which(is.na(RodentIndex2))] <- sample(1:nLevels.rCov, length(which(is.na(RodentIndex2))), replace = TRUE)
  }
  
  ## Rodent abundance (continuous)
  RodentAbundance <- nim.data$RodentAbundance
  if(NA %in% RodentAbundance){
    #RodentAbundance[which(is.na(RodentAbundance))] <- mean(RodentAbundance, na.rm = TRUE)
    RodentAbundance[which(is.na(RodentAbundance))] <- 0
  }
  RodentAbundance_pert <- RodentAbundance + (1-pertFac.rodent)
  
  RodentAbundance2 <- nim.data$RodentAbundance2
  if(NA %in% RodentAbundance2){
    #RodentAbundance2[which(is.na(RodentAbundance2))] <- mean(RodentAbundance2, na.rm = TRUE)
    RodentAbundance2[which(is.na(RodentAbundance2))] <- 0
  }
  RodentAbundance2_pert <- RodentAbundance2 + (1-pertFac.rodent)
  
  
  #--------------------------------------------------------#
  # Set initial values for conditional perturbation factor #
  #--------------------------------------------------------#
  
  ## Harvest perturbation based on rodent abundance
  if(nim.data$thresholdAbove){
    pertFac.mH.flex <- ifelse(RodentAbundance_pert > nim.data$threshold.rodent.mH, nim.data$factor.mH.rodent, 1)
  }else{
    pertFac.mH.flex <- ifelse(RodentAbundance_pert < nim.data$threshold.rodent.mH, nim.data$factor.mH.rodent, 1)
  }
  
  pertFac.mH.flex[1:nim.constants$Tmax] <- 1
  pertFac.mH.flex <- pertFac.mH.flex[1:(Tmax+1)]
  
  #---------------------------------------------------#
  # Set initial values for vital rate base parameters #
  #---------------------------------------------------#
  
  ## Harvest and natural mortality
  Mu.mH.juv <- runif(1, 0.05, 0.2)
  Mu.mH.ad <- runif(1, 0.05, 0.2)
  
  Mu.mH <- c(Mu.mH.juv, rep(Mu.mH.ad, length(2:Amax)))
  
  if(Mu.mO_fixInits){
    
    if(HoenigPrior){
      
      Mu.mO.ad <- -log(mean(c(0.71, 0.5, 0.59, 0.59)))
      Mu.mO <- c(-log(0.33), rep(Mu.mO.ad, Amax-1))
      JuvAdRatio <- Mu.mO[1]/Mu.mO.ad
      
    }else{
      
      Mu.Snat <- c(0.33, 0.71, 0.5, 0.59, 0.59)
      Mu.mO <- -log(Mu.Snat)
    }
    
  }else{
    
    if(HoenigPrior){
      
      JuvAdRatio <- exp(nim.constants$ratioJA.logmean)
      Mu.mO.ad <- exp(nim.constants$mnat.logmean)
      Mu.mO <- c(Mu.mO.ad*JuvAdRatio, rep(Mu.mO.ad, Amax-1))
      
    }else{
      #Mu.Snat <- nim.constants$Snat.mean
      Mu.Snat <- rep(NA, Amax)
      Mu.Snat[1] <- runif(1, 0.4, 0.6)
      Mu.Snat[2:Amax] <- runif(Amax-1, 0.6, 1)
      Mu.mO <- -log(Mu.Snat)
    }
    
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
  
  logsigma.Imm <- sd(log(truncnorm::rtruncnorm(10000, a = 0, b = maxImm, mean = Mu.Imm, sd = sigma.Imm)))
  
  
  ## Random effect standard deviations
  sigma.mH <- runif(1, 0.05, 0.5)
  sigma.Psi <- runif(1, 0.05, 0.5)
  sigma.rho <- runif(1, 0.05, 0.5)
  
  if(mO.varT){
    sigma.mO <- runif(1, 0.05, 0.5)
  }else{
    sigma.mO <- 0
  }
  
  ## Random effects (initialize at to 0)
  epsilon.mH <- rep(0, Tmax+1)
  epsilon.mO <- rep(0, Tmax+1)
  epsilon.Psi <- rep(0, Tmax+1)
  epsilon.rho <- rep(0, Tmax+1)
  
  ## Covariate effects
  
  # Harvest effort on mH
  if(fitCov.mH){
    betaHE.mH <- runif(1, 0, 0.2)
  }else{
    betaHE.mH <- 0
  }
  
  # Rodent abundance on mO
  if(fitCov.mO){
    betaR.mO <- runif(1, -0.2, 0)
  }else{
    betaR.mO <- 0
  }
  
  # Rodent abundance on Psi
  if(fitCov.Psi){
    if(rCov.idx){
      betaR.Psi <- rep(0, nLevels.rCov)
      for(x in 2:nim.constants$nLevels.rCov){
        betaR.Psi[x] <- runif(1, 0, 2)
      }
    }else{
      betaR.Psi <- runif(1, 0, 0.2)
    }
  }else{
    
    if(rCov.idx){
      betaR.Psi <- rep(0, nLevels.rCov)
    }else{
      betaR.Psi <- 0
    }
  }
  
  # Rodent abundance on rho
  if(fitCov.rho){
    
    if(rCov.idx){
      betaR.rho <- rep(0, nLevels.rCov)
      for(x in 2:nim.constants$nLevels.rCov){
        betaR.rho[x] <- runif(1, 0, 2)
      }
    }else{
      betaR.rho <- runif(1, 0, 0.2)
    }
    
  }else{
    
    if(rCov.idx){
      betaR.rho <- rep(0, nLevels.rCov)
    }else{
      betaR.rho <- 0
    }
  }
  
  # Rodent abundance on immR
  if(fitCov.immR){
    if(rCov.idx){
      betaR.immR <- rep(0, nLevels.rCov)
    }else{
      betaR.immR <- 0
    }
  }
  
  #-------------------------------------#
  # Calculate year-specific vital rates #
  #-------------------------------------#
  
  mH <- mO <- Psi <- rho <- matrix(NA, nrow = Amax, ncol = Tmax + 1)
  S0 <- rep(NA, Tmax + 1)
  
  for(t in 1:(Tmax+1)){
      
    ## Winter harvest mortality hazard rate
    mH[1:Amax, t] <- exp(log(Mu.mH[1:Amax]) + betaHE.mH*NHunters[t] + epsilon.mH[t])*pertFac.mH[t]*pertFac.mH.flex[t]
      
    ## Other (natural) mortality hazard rate
    mO[1:Amax, t] <- exp(log(Mu.mO[1:Amax]) + betaR.mO*RodentAbundance_pert[t+1] + epsilon.mO[t])*pertFac.mO[t]
    
    ## Pregnancy rate
    Psi[1, t] <- 0
    
    if(fitCov.Psi & rCov.idx){
      Psi[2:Amax, t] <- plogis(qlogis(Mu.Psi[2:Amax]) + betaR.Psi[RodentIndex[t]] + epsilon.Psi[t])
    }else{
      Psi[2:Amax, t] <- plogis(qlogis(Mu.Psi[2:Amax]) + betaR.Psi*RodentAbundance_pert[t] + epsilon.Psi[t])
    }
    
    ## Placental scars
    rho[1, t] <- 0
    if(fitCov.Psi & rCov.idx){
      rho[2:Amax, t] <- exp(log(Mu.rho[2:Amax]) + betaR.rho[RodentIndex[t]] + epsilon.rho[t])
    }else{
      rho[2:Amax, t] <- exp(log(Mu.rho[2:Amax]) + betaR.rho*RodentAbundance_pert[t] + epsilon.rho[t])
    }
    
    ## Early survival
    S0[t] <- Mu.S0*pertFac.S0[t]
  }
  
  ## Survival probability
  S <- exp(-(mH + mO))
  
  ## Proportion harvest mortality
  alpha <- mH/(mH + mO)
  
  ## Harvest rate
  h <- (1 - S)*alpha
  
  ## Immigrant numbers
  Imm <- round(truncnorm::rtruncnorm(Tmax+1, a = 0, b = maxImm, mean = Mu.Imm, sd = sigma.Imm)*pertFac.immR)
  Imm[1] <- 0
  
  
  #---------------------------------------------------------#
  # Simulate initial values for population-level quantities #
  #---------------------------------------------------------#
  
  ## Prepare empty vectors and matrices
  H <- N <- B <- L <- R <- matrix(NA, nrow = Amax, ncol = Tmax+1)
  
  ## Set initial population sizes
  for(a in 1:Amax){
    N[a, 1] <- round(runif(1, minN1[a], maxN1[a]))
  }
  
  ## Set age class 0 reproductive contributions to 0
  B[1, 1:(Tmax+1)] <- L[1, 1:(Tmax+1)] <- R[1, 1:(Tmax+1)] <- 0
  
  ## Set first year reproductive contributions to 0
  B[2:Amax, 1] <- L[2:Amax, 1] <- R[2:Amax, 1] <- 0
  
  
  for (t in 1:Tmax){
    
    # a) Project local survivors to the next year
    #---------------------------------------------
    
    ## Age classes 1 to 3 (indeces = 2, 3, 4): age classes 0, 1, and 2 survivors    
    for(a in 1:(Amax-2)){
      N[a+1, t+1] <- rbinom(1, size = N[a, t], prob = S[a, t])
    }			
    
    ## Age class 4+ (index = 5): age class 3 and 4+ survivors
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
  
  if(any(nim.data$C_w[,1:nim.constants$Tmax] > N[, 1:nim.constants$Tmax])){
    stop('Simulation resulted in less alive than harvested (winter). Retry.')
  }
  
  
  # e) Assemble and return list of initial values
  #----------------------------------------------  
  
  ## Fill out NA values in N, survN1, B, L, and R
  for (a in 2:Amax){
    #B[a,1] <- rbinom(1, size = N[a, 1], prob = Psi[a, 1])
    #L[a,1] <- rpois(1, lambda = B[a, 1] * rho[a, 1] * 0.5)
    #R[a,1] <- rbinom(1, size = L[a, 1], prob = S0[t])
    B[a, 1] <- 0
    L[a, 1] <- 0
    R[a, 1] <- 0
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
    localN.tot = colSums(R[2:Amax, 1:(Tmax+1)]) + colSums(N[2:Amax, 1:(Tmax+1)]) + 1,
    
    Mu.mH = Mu.mH,
    Mu.mO = Mu.mO,
    Mu.Psi = Mu.Psi,
    Mu.rho = Mu.rho,
    Mu.S0 = Mu.S0,
    
    sigma.mH = sigma.mH,
    sigma.mO = sigma.mO,
    sigma.Psi = sigma.Psi,
    sigma.rho = sigma.rho,
    
    epsilon.mH = epsilon.mH,
    epsilon.mO = epsilon.mO,
    epsilon.Psi = epsilon.Psi,
    epsilon.rho = epsilon.Psi,
    
    eta.mH = epsilon.mH,
    eta.mO = epsilon.mO,
    tau.mO = 0,
    C.mO = 0,
    
    mH = mH,
    mO = mO, 
    S = S,
    alpha = alpha, 
    h = h, 
    Psi = Psi, 
    rho = rho,
    S0 = S0,
    
    pertFac.mH.flex = pertFac.mH.flex,
    
    logDev.mH = log(mH[1, ]) - log(Mu.mH[1])
  )
  
  ## Add initial values for parameters specific to survival prior model versions
  if(HoenigPrior){
    InitVals$JuvAdRatio <- JuvAdRatio
    InitVals$Mu.mO.ad <- Mu.mO.ad
  }else{
    InitVals$Mu.Snat <- Mu.Snat
  }
  
  ## Add initial values for covariate effects
  if(fitCov.mH){
    InitVals$betaHE.mH <- betaHE.mH
  }
  
  if(fitCov.mO){
    InitVals$betaR.mO <- betaR.mO
    InitVals$betaD.mO <- 0
    InitVals$betaRxD.mO <- 0
    InitVals$gamma.mO <- 0
  }
  
  if(fitCov.Psi){
    InitVals$betaR.Psi <- betaR.Psi
  }
  
  if(fitCov.rho){
    InitVals$betaR.rho <- betaR.rho
  }
  
  if(fitCov.immR){
    InitVals$betaR.immR <- betaR.immR
    InitVals$betaR.immR <- betaR.immR
    InitVals$betaD.immR <- 0
    InitVals$betaRxD.immR <- 0
    InitVals$gamma.immR <- 0
  }
  
  ## Add initial values specific to immigration model versions
  if(imm.asRate){
    
    if(useData.gen){
      ImmData <- rbinom(n = nim.constants$Xgen, size = 1, p = mean(nim.data$pImm))
      ImmData[which(nim.data$pImm == 0)] <- 0
      ImmData[which(nim.data$pImm == 1)] <- 1
      
      InitVals$ImmData <- ImmData
      InitVals$Mu.immR <- sum(InitVals$ImmData)/(nim.constants$Xgen-sum(InitVals$ImmData))
      
      if(!poolYrs.genData){
        ImmData_pre <- rbinom(n = nim.constants$Xgen_pre, size = 1, p = mean(nim.data$pImm_pre))
        ImmData_pre[which(nim.data$pImm_pre == 0)] <- 0
        
        InitVals$ImmData_pre <- ImmData_pre
        InitVals$immR_pre <- rep(InitVals$Mu.immR, nim.constants$Tmax_Gen_pre)
      }
      
    }else{
      InitVals$Mu.immR <- mean(Imm[2:(Tmax+1)]/colSums(R)[2:(Tmax+1)])
    }
    
    InitVals$immR <- c(0, rep(InitVals$Mu.immR, Tmax))
    InitVals$sigma.immR <- runif(1, 0, 0.5)
    InitVals$epsilon.immR <- rep(0, Tmax+1)
    InitVals$eta.immR <- rep(0, Tmax+1)
    InitVals$tau.immR <- 0
    InitVals$C.immR <- 0
    
  }else{
    
    InitVals$Mu.Imm <- Mu.Imm
    InitVals$sigma.immR <- logsigma.Imm
    InitVals$epsilon.immR <- rep(0, Tmax+1)
    InitVals$eta.immR <- rep(0, Tmax+1)
    InitVals$tau.immR <- 0
    InitVals$C.immR <- 0
    InitVals$ImmExp <- Imm
  }
  
  ## Add initial values for missing covariate values (if applicable)
  if(fitCov.mH & (NA %in% nim.data$HarvestEffort)){
    Inits_NHunters <- rep(NA, length(NHunters))
    Inits_NHunters[which(is.na(nim.data$HarvestEffort))] <- NHunters[which(is.na(nim.data$HarvestEffort))]
    InitVals$HarvestEffort <- Inits_NHunters
  }
  
  if(rCov.idx & (NA %in% nim.data$RodentIndex)){
    Inits_RodentIndex <- rep(NA, length(RodentIndex))
    Inits_RodentIndex[which(is.na(nim.data$RodentIndex))] <- RodentIndex[which(is.na(nim.data$RodentIndex))]
    InitVals$RodentIndex <- Inits_RodentIndex
  }
  if(!rCov.idx & (NA %in% nim.data$RodentAbundance)){
    Inits_RodentAbundance <- rep(NA, length(RodentAbundance))
    Inits_RodentAbundance[which(is.na(nim.data$RodentAbundance))] <- RodentAbundance[which(is.na(nim.data$RodentAbundance))]
    InitVals$RodentAbundance <- Inits_RodentAbundance
  }
  
  if(rCov.idx & (NA %in% nim.data$RodentIndex2)){
    Inits_RodentIndex2 <- rep(NA, length(RodentIndex2))
    Inits_RodentIndex2[which(is.na(nim.data$RodentIndex2))] <- RodentIndex2[which(is.na(nim.data$RodentIndex2))]
    InitVals$RodentIndex2 <- Inits_RodentIndex2
  }
  if(!rCov.idx & (NA %in% nim.data$RodentAbundance2)){
    Inits_RodentAbundance2 <- rep(NA, length(RodentAbundance2))
    Inits_RodentAbundance2[which(is.na(nim.data$RodentAbundance2))] <- RodentAbundance2[which(is.na(nim.data$RodentAbundance2))]
    InitVals$RodentAbundance2 <- Inits_RodentAbundance2
  }
  
  if(indLikelihood.genData & rCov.idx & !poolYrs.genData){
    InitVals$RodentIndex2_pre <- sample(0:nim.constants$nLevels.rCov, size = nim.constants$Tmax_Gen_pre, replace = TRUE)
  }
  if(!rCov.idx & !poolYrs.genData){
    InitVals$RodentAbundance2_pre <- rnorm(nim.constants$Tmax_Gen_pre, mean = 0, sd = 1)
  }
  
  ## Add perturbed covariate values
  InitVals$RodentAbundance_pert <- RodentAbundance_pert
  InitVals$RodentAbundance2_pert <- RodentAbundance2_pert
  
  ## Add initial values for rodent model parameters
  InitVals$beta.RodMod <- runif(3, -1, 1)
  InitVals$beta.RodCorr <- runif(1, 0, 1)
  InitVals$sigmaT.RodAbun <- runif(1, 0, 0.2)
  InitVals$sigmaT.RodAbun2 <- runif(1, 0, 0.2)
  
  ## Return initial values
  return(InitVals)
}
