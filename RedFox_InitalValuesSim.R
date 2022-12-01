#*******************************************************************************#
#* INITIAL VALUE FUNCTIONS
#*******************************************************************************#

## Function to simulate and assemble initial values for IPM
RF.IPM.inits <- function(IPM.data, IPM.constants, minN1, maxN1, minImm, maxImm){
  
  A <- IPM.constants$A
  Tmax <- IPM.constants$Tmax
  
  #------------------------------------#
  # Set initial values for vital rates #
  #------------------------------------#
  
  ## Harvest and natural mortality
  Mu.mH <- runif(A, 0.05, 0.2)
  
  # Literature prior
  Mu.Snat <- Snat.mean
  Mu.mO <- -log(Mu.Snat)
  #Mu.mnat <- Mu.mO <- mnat.mean
  
  # Hoening model prior
  #Mu.mO <- c(exp(mnat.logmean)*exp(ratioJA.logmean),
  #           rep(exp(mnat.logmean), 4))
  
  ## Annual survival
  Mu.S <- exp(-(Mu.mH + Mu.mO))
  
  ## Pregnancy rate
  Mu.Psi <- c(0, runif(A-1, 0.2, 0.8))
  
  ## Placental scars
  Mu.rho <- c(0, runif(A-1, 3, 6))
  
  ## Early survival
  mean.S0 <- Mu.S0 <- S0.mean
  
  ## Random effect standard deviations
  sigma.mH <- runif(1,0.05,0.5)
  #sigma.mO <- runif(1,0.05,0.5)
  sigma.Psi <- runif(1,0.05,0.5)
  sigma.rho <- runif(1,0.05,0.5)
  
  ## Random effects (set to 0)
  epsilon.mH <- rep(0, RF.constants$Tmax)
  #epsilon.mO <- rep(0, RF.constants$Tmax),
  epsilon.Psi <- rep(0, RF.constants$Tmax+1)
  epsilon.rho <- rep(0, RF.constants$Tmax+1)
  
  ## Covariate effects
  #betaR.Psi <- runif(1, 0, 2),
  #betaRI2.Psi <- runif(1, 0, 2),
  #betaRI3.Psi <- c(0, runif(2, 0, 2),
  
  Mu.Imm <- runif(1, minImm, maxImm)
  #sigma.Imm <- (runif(1, 10, 40)),
  

  #---------------------------------------------------------#
  # Simulate initial values for population-level quantities #
  #---------------------------------------------------------#
  
  ## Prepare empty vectors and matrices
  H <- N <- B <- L <- R <- matrix(NA, nrow = A, ncol = Tmax)
  
  ## Sample preliminary immigrant numbers
  Imm <- c(0, rpois(Tmax-1, Mu.Imm))
  
  ## Set initial population sizes
  for(a in 1:A){
    N[a, 1] <- round(runif(1, minN1[a], maxN1[a]))
  }
  
  ## Set age class 0 reproductive contributions to 0
  B[1,1:Tmax] <- L[1,1:Tmax] <- R[1,1:Tmax] <- 0
  
  
  for (t in 1:(Tmax-1)){

    # a) Project local survivors to the next year
    #---------------------------------------------
      
    ## Age classes 1 to 3 (indeces = 2, 3, 4): age classes 0, 1, and 2 survivors    
    for(a in 1:(A-2)){
      N[a+1,t+1] <- rbinom(1, size = N[a,t], prob = Mu.S[a])
    }			
      
    # Age class 4+ (index = 5): age class 3 and 4+ survivors
    N[A,t+1] <- rbinom(1, size = N[A-1,t] + N[A,t], prob = Mu.S[A])
      
      
    # b) Sample through reproductive season 
    #--------------------------------------
      
    for (a in 2:A){
        
      ## Breeding Population Size: Number of females that reproduce
      B[a,t+1] <- rbinom(1, size = N[a,t+1], prob = Mu.Psi[a])
        
      ## Litter Size (in utero): Number of pups produced by females of age class a
      L[a,t+1] <- rpois(1, lambda = B[a,t+1]*Mu.rho[a]*0.5)
        
      ## Number Recruits: Number of pups surviving to emerge from the den
      R[a,t+1] <- rbinom(1, size = L[a,t+1], prob = Mu.S0)
    }
      
      
    # c) Add new recruits and immigrants
    #------------------------------------
      
    N[1,t+1] <- sum(R[1:A,t+1]) + Imm[t+1]
 }
  
  
  # d) Check for years with more harvests than alive individuals
  #-------------------------------------------------------------
  
  if(any(IPM.data$C > N)){
    stop('Simulation resulted in less alive than harvested. Retry.')
  }
  
  # e) Assemble and return list of initial values
  #----------------------------------------------  
  
  ## Fill out NA values in B, L, and R
  for (a in 2:A){
    #B[a,1] <- rbinom(1, size = N[a,1], prob = Mu.Psi[a])
    #L[a,1] <- rpois(1, lambda = B[a,1]*Mu.rho[a]*0.5)
    #R[a,1] <- rbinom(1, size = L[a,1], prob = Mu.S0)
    B[a,1] <- 0
    L[a,1] <- 0
    R[a,1] <- 0
  }
  # NOTE: These nodes do not appear in the model and it therefore does not 
  #       matter what numbers they contain. Filling them in prevents a warning
  #       about NA nodes when building the model. 
  
  
  ## Return list
  return(list(
    N = N, 
    initN = N[,1],
    B = B, 
    L = L,
    R = R,
    Imm = Imm,
    
    Mu.mH = Mu.mH,
    
    # Literature prior
    Mu.Snat = Snat.mean,
    #Mu.mnat = mnat.mean,
    
    # Hoening prior
    #Mu.mO.ad = exp(mnat.logmean),
    #JuvAdRatio = exp(ratioJA.logmean),
    
    Mu.Psi = Mu.Psi,
    Mu.rho = Mu.Psi,
    
    mean.S0 = mean.S0,
    
    sigma.mH = sigma.mH,
    #sigma.mO = sigma.mO,
    sigma.Psi = sigma.Psi,
    sigma.rho = sigma.rho,
    
    epsilon.mH = epsilon.mH,
    #epsilon.mO = epsilon.mO,
    epsilon.Psi = epsilon.Psi,
    epsilon.rho = epsilon.Psi
    
    #betaR.Psi = betaR.Psi,
    #betaRI2.Psi = betaRI2.Psi,
    #betaRI3.Psi = betaRI3.Psi,
    
    #Mu.Imm = Mu.Imm,
    #sigma.Imm = sigma.Imm,
    
  ))
  
}
