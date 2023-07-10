#' Write NIMBLE code for red fox IPM
#'
#' @return an R call object specifying the model structure for the red fox IPM. 
#' @export
#'
#' @examples

writeCode_redfoxIPM <- function(){
  
  ## Check for incompatible toggles
  if(!imm.asRate & fitCov.immR){
    stop("Incompatible model settings. Rodent covariate effect on immigration can only be fit (fitCov.immR = TRUE) if immigration is estimated as a rate (imm.asRate = TRUE). ")
  }
  
  ## Write model code
  redfox.code <- nimbleCode({
    
    
    ##########################  
    #### POPULATION MODEL ####
    ##########################
    
    ### Likelihood (age classes: 1, 2, 3+)
    
    ## Survival
    
    for(t in 1:Tmax){ 
      
      # Age class 0 (index = 1): sum of local reproduction & immigrants
      N[1, t+1] <- sum(R[2:Amax, t+1]) + Imm[t+1]     
      
      # Age classes 1 to 3 (indeces = 2, 3, 4): age classes 0, 1, and 2 survivors    
      for(a in 1:(Amax-2)){
        N[a+1, t+1] ~ dbin(S[a, t], N[a, t])
      }			
      
      # Age class 5+ (index = Amax = 5): age class 4 and 5+ survivors
      N[Amax, t+1] ~ dbin(S[Amax, t], N[Amax-1, t] + N[Amax, t])
    }
    
    ## Reproduction
    
    # Age class 0 (young of the year --> do not reproduce in year of birth)
    B[1, 1:(Tmax+1)] <- 0
    L[1, 1:(Tmax+1)] <- 0
    R[1, 1:(Tmax+1)] <- 0
    
    # Age classes 1 to 3+    	    
    for(t in 1:(Tmax+1)){        				
      for(a in 2:Amax){
        
        # Breeding Population Size: Number of females that reproduce
        B[a, t] ~ dbin(Psi[a, t], N[a, t])
        
        # Litter Size (in utero): Number of pups produced by females of age class a
        L[a, t] ~ dpois(B[a, t]*rho[a, t]*0.5)
        
        # Number Recruits: Number of pups surviving to emerge from the den
        R[a, t] ~ dbin(S0[t], L[a, t])
      } 
    }
    
    #===============================================================================================
    
    
    
    ############################
    #### DERIVED QUANTITIES ####
    ############################
    
    for(t in 1:(Tmax+1)){
      N.tot[t] <- sum(N[1:Amax, t])
      R.tot[t] <- sum(R[1:Amax, t])		
      B.tot[t] <- sum(B[1:Amax, t])
    }
    
    #===============================================================================================
    
    
    
    ######################################
    #### WINTER AGE-AT-HARVEST MODULE ####
    ######################################
    
    ### Parameters:
    # N = number of individuals in a given age class at a given time
    # h = time-dependent probability of dying from hunting ([1] = adults, [2] = juveniles)
    
    ### Data:
    # C = age-at-harvest matrix
    # pData = annual proportion of harvests with (complete) carcass data
    
    ### Likelihood
    
    for(t in 1:Tmax){
      for(a in 1:Amax){
        C[a, t] ~ dbin(h[a, t]*pData[t], N[a, t])
      }
    }
    
    #===============================================================================================
    
    
    
    ###############################
    #### PLACENTAL SCAR MODULE ####
    ###############################
    
    ### Parameters:
    # rho = expected number of placental scars (fetuses)
    # Psi = pregnancy rate
    
    ## Data:
    # P1 = individual placental scar counts
    # P1_age = individual ages associated with P1
    # P1_year = year associated with P1
    
    # P2 = individual presence/absence of placental scars
    # P2_age = individual ages associated with P2
    # P2_year = year associated with P2
    
    
    ### Likelihood (litter size)
    
    for(x in 1:X1){
      P1[x] ~ dpois(rho[P1_age[x], P1_year[x]])
    }
    
    ### Likelihood (pregnancy rate)
    
    for(x in 1:X2){
      P2[x] ~ dbern(Psi[P2_age[x], P2_year[x]])
    }
    
    #===============================================================================================
    
    
    ###########################################
    #### GENETIC IMMIGRATION STATUS MODULE ####
    ###########################################
    
    ### Parameters:
    # Mu.immR = average immigration rate
    # immR = annual immigration rates
    
    ## Data:
    # pImm = individual probability of being an immigrant

    
    ### Likelihood (immigration status of sampled individuals)
    if(imm.asRate){
      if(useData.gen){
        
        ## Likelihood for individuals (within study period) to be immigrants
        for(x in 1:Xgen){
          ImmData[x] ~ dbern(pImm[x])
        }
        
        if(poolYrs.genData){
          
          ## Derivation of average immigration rate
          Mu.immR <- sum(ImmData[1:Xgen]) / (Xgen - sum(ImmData[1:Xgen]))
        
        }else{
          
          ## Likelihood for individuals outside the study period to be immigrants
          for(x in 1:Xgen_pre){
            ImmData_pre[x] ~ dbern(pImm_pre[x])
          }
          
          ## Derivation of year-specific immigration rates
          # Within study period
          immR[1:Tmax_Gen] <- calculateImmR(ImmData = ImmData[1:Xgen], 
                                            yearIdx = pImm_yrs[1:Xgen],
                                            Tmax = Tmax_Gen, skip_t1 = FALSE)
          
          # Outside study period
          immR_pre[1:Tmax_Gen_pre] <- calculateImmR(ImmData = ImmData_pre[1:Xgen_pre], 
                                                    yearIdx = pImm_yrs_pre[1:Xgen_pre],
                                                    Tmax = Tmax_Gen_pre, skip_t1 = FALSE)
          
        }
        
      }else{
        
        Mu.immR ~ dunif(0, 10)
        
      }
    }

    
    ################################
    #### PRIORS AND CONSTRAINTS ####
    ################################
    
    ## Survival and mortality
    
    for(t in 1:Tmax){ 
      
      # Harvest mortality hazard rate
      if(fitCov.mH){
        log(mH[1:Amax, t]) <- log(Mu.mH[1:Amax]) + betaHE.mH*HarvestEffort[t] + epsilon.mH[t]
      }else{
        log(mH[1:Amax, t]) <- log(Mu.mH[1:Amax]) + epsilon.mH[t]
      }
      
      # Other (natural) mortality hazard rate
      log(mO[1:Amax, t]) <- log(Mu.mO[1:Amax]) + betaRd.mO*Reindeer[t] + betaR.mO*RodentAbundance[t] + betaRxRd.mO*Reindeer[t]*RodentAbundance[t] + epsilon.mO[t]
      #log(mO[1:Amax, t]) <- log(Mu.mO[1:Amax]) + epsilon.mO[t]
      
      # Survival probability
      S[1:Amax, t] <- exp(-(mH[1:Amax, t] + mO[1:Amax,t]))
      
      # Proportion harvest mortality
      alpha[1:Amax, t] <- mH[1:Amax, t]/(mH[1:Amax, t] + mO[1:Amax, t])
      
      # Harvest rate
      h[1:Amax, t] <- (1-S[1:Amax, t])*alpha[1:Amax, t]
      
    }
    
    # Median harvest mortality hazard rates
    
    # Age-dependent
    for(a in 1:Amax){
      Mu.mH[a] ~ dunif(0, 5) 
    }
    
    # Age-independent   
    #Mu.mH.all ~ dunif(0, 5) 
    #Mu.mH[1:Amax] <- Mu.mH.all
    
    # Median other (natural) cause mortality hazard rates
    #* INFORMATIVE PRIOR REQUIRED: LITERATURE VALUE / HOENING MODEL CALCULATION
    
    if(HoeningPrior){
      # Using prior distributions calculated with Hoening model
      Mu.mO.ad ~ dlnorm(mnat.logmean, sdlog = mnat.logsd)
      Mu.mO[2:5] <- Mu.mO.ad
      Mu.mO[1] <- Mu.mO.ad*JuvAdRatio #* NOTE: Can be provided as constant or distribution
      
      #JuvAdRatio <- exp(JAratio.logmean)
      JuvAdRatio ~ dlnorm(ratioJA.logmean, sdlog = ratioJA.logsd)
      
    }else{
      # Using literature values on age-specific survival
      for(a in 1:Amax){
        Mu.mO[a] <- -log(Mu.Snat[a])
        Mu.Snat[a] ~ T(dnorm(Snat.mean[a], sd = Snat.sd[a]), 0, 1)   
      }
    }
    
   
    ## Covariate effects
    if(fitCov.mH){
      betaHE.mH ~ dunif(0, 5) # Effect of harvest effort on mH
    }
    
    betaRd.mO ~ dunif(-5, 0)
    betaR.mO ~ dunif(-5, 0)
    betaRxRd.mO ~ dunif(-5, 5)
    
    #---------------------------------------------------------------------------------------------
    
    
    ## Pregnancy rate
    
    for(t in 1:(Tmax+1)){
      Psi[1, t] <- 0
      
      if(fitCov.Psi){
        if(rCov.idx){
          logit(Psi[2:Amax,t]) <- logit(Mu.Psi[2:Amax]) + betaR.Psi[RodentIndex[t]] + epsilon.Psi[t] # Reindeer.rodent interaction not (yet) written in
        }else{
          logit(Psi[2:Amax,t]) <- logit(Mu.Psi[2:Amax]) + betaR.Psi*RodentAbundance[t] + epsilon.Psi[t]
        }
      }else{
        logit(Psi[2:Amax, t]) <- logit(Mu.Psi[2:Amax]) + epsilon.Psi[t]
      }
    }
    
    Mu.Psi[1] <- 0
    for(a in 2:Amax){	
      Mu.Psi[a] ~ dunif(0, 1)
    }
    
    if(fitCov.Psi){
     if(rCov.idx){
       betaR.Psi[1] <- 0 # --> Lowest level corresponds to intercept
       for(x in 2:nLevels.rCov){
         betaR.Psi[x] ~ dunif(-5, 5)
       }
     }else{
       betaR.Psi ~ dunif(-5, 5)
     }
    }
    
    
    #---------------------------------------------------------------------------------------------
    
    
    ## Litter size
    
 
    for(t in 1:(Tmax+1)){
      rho[1, t] <- 0
      
      if(fitCov.rho){
        if(rCov.idx){
          log(rho[2:Amax, t]) <- log(Mu.rho[2:Amax]) + betaR.rho[RodentIndex[t]] + epsilon.rho[t]
        }else{
          log(rho[2:Amax, t]) <- log(Mu.rho[2:Amax]) + betaR.rho*RodentAbundance[t] + epsilon.rho[t]
        }
      }else{
        log(rho[2:Amax, t]) <- log(Mu.rho[2:Amax]) + epsilon.rho[t]
      }
    }
    
    Mu.rho[1] <- 0
    for(a in 2:Amax){
      Mu.rho[a] ~ dunif(0, maxPups) # Baseline number of pups 
    }
    
    if(fitCov.rho){
      if(rCov.idx){
        betaR.rho[1] <- 0 # --> Lowest level corresponds to intercept
        for(x in 2:nLevels.rCov){
          betaR.rho[x] ~ dunif(-5, 5)
        }
      }else{
        betaR.rho ~ dunif(-5, 5)
      }
    }
    
    #---------------------------------------------------------------------------------------------  
    
    
    ## Denning survival
    #* INFORMATIVE PRIOR REQUIRED: LITERATURE VALUE
    
    for(t in 1:(Tmax+1)){ 
      S0[t] <- Mu.S0
      #S0[t] <- exp(-m0[t])
      #log(m0[t]) <- log(-log(Mu.S0)) + epsilon.m0[t]
    }
    
    Mu.S0 ~ T(dnorm(S0.mean, sd = S0.sd), 0, 1)  
    #---------------------------------------------------------------------------------------------
    
    
    ## Immigration
    
    if(imm.asRate){
      
      if(useData.gen & !poolYrs.genData){
        
        ## Extraction of log mean and sd for immigration rate
        log(Mu.immR) <- log(mean(c(immR[1:Tmax_Gen], immR_pre[1:Tmax_Gen_pre]))) 
        sigma.immR <- calculateLogSD(immR = c(immR[1:Tmax_Gen], immR_pre[1:Tmax_Gen_pre]),
                                     replace0 = 0.01)
        
        ## Projection of immigration rates beyond genetic data coverage
        for(t in (Tmax_Gen+1):(Tmax+1)){
          immR[t] ~ dlnorm(meanlog = log(Mu.immR), sdlog = sigma.immR)
        }
        
      }else{
        
        if(fitCov.immR){
          if(rCov.idx){
            for(t in 1:(Tmax+1)){
              log(immR[t]) <- log(Mu.immR) + betaR.immR[RodentIndex2[t]] + epsilon.immR[t]
            }
          }else{
            log(immR[1:(Tmax+1)]) <- log(Mu.immR) + betaR.immR*RodentAbundance2[1:(Tmax+1)] + epsilon.immR[1:(Tmax+1)]
          }
        }else{
          log(immR[1:(Tmax+1)]) <- log(Mu.immR) + epsilon.immR[1:(Tmax+1)]
        }
      }
      
      for(t in 1:(Tmax+1)){ 
        Imm[t] ~ dpois(sum(R[2:Amax, t])*immR[t])
      }
      
      
    }else{
      
      ## Discrete uniform prior for immigrant numbers
      Imm[1] <- 0 # (Immigration in the first year cannot be disentangled from reproduction)
      #ImmT[1] <- 0 
      
      for(t in 2:(Tmax+1)){
        Imm[t] ~ dcat(DU.prior.Imm[1:uLim.Imm]) 
        #Imm[t] ~ dpois(ImmT[t])
        #ImmT[t] ~ T(dnorm(Mu.Imm, sd = sigma.Imm), 0, uLim.Imm)
      }
      
      DU.prior.Imm[1:uLim.Imm] <- 1/uLim.Imm
      
      ## Derivation of immigration rates
      immR[1] <- 0
      for(t in 2:(Tmax+1)){
        immR[t] <- Imm[t] / sum(R[2:Amax, t])
      }

    }
    
    ## Prior for rodent effect
    if(fitCov.immR){
      if(rCov.idx){
        betaR.immR[1] <- 0 # --> Lowest level corresponds to intercept
        for(x in 2:nLevels.rCov){
          betaR.immR[x] ~ dunif(-5, 5)
        }
      }else{
        betaR.immR ~ dunif(-5, 5)
      }
    }
    
    #---------------------------------------------------------------------------------------------
    
    
    ## Initial population size (discrete uniform prior) 

    N[1:Amax, 1] <- initN[1:Amax]
    
    for(a in 1:Amax){
      initN[a] ~ dcat(DU.prior.N[1:uLim.N]) 
    }
    
    DU.prior.N[1:uLim.N] <- 1/uLim.N
    
    #---------------------------------------------------------------------------------------------
    
    
    ## Random year variation
    for(t in 1:Tmax){  
      epsilon.mH[t] ~ dnorm(0, sd = sigma.mH)
      epsilon.mO[t] ~ dnorm(0, sd = sigma.mO)
    }
    
    for(t in 1:(Tmax+1)){
      epsilon.Psi[t] ~ dnorm(0, sd = sigma.Psi)
      epsilon.rho[t] ~ dnorm(0, sd = sigma.rho) 
      # epsilon.m0[t] ~ dnorm(0, sd = sigma.m0)
    }
    
    sigma.mH ~ dunif(0, 5)
    sigma.Psi ~ dunif(0, 5)
    sigma.rho ~ dunif(0, 5)
    #sigma.m0 ~ dunif(0, 5)
    
    if(mO.varT){
      sigma.mO ~ dunif(0, 5)
    }else{
      sigma.mO <- 0
    }
    
    if(imm.asRate & poolYrs.genData){
      for(t in 1:(Tmax+1)){
        epsilon.immR[t] ~ dnorm(0, sd = sigma.immR)
      }
      sigma.immR ~ dunif(0, 10)
    }
    
    #===============================================================================================
    
    
    
    ##############################
    #### COVARIATE IMPUTATION ####
    ##############################
    
    ## Missing covariate value(s) in number of successful hunters
    if(fitCov.mH){
      for(t in 1:Tmax){
        HarvestEffort[t] ~ dnorm(0, sd = 1)
      }
    }
    
    ## Missing covariate value(s) in rodent abundance
    if(rCov.idx){
      
      for(t in 1:Tmax+1){
        RodentIndex[t] ~ dcat(DU.prior.rCov[1:nLevels.rCov]) 
        RodentIndex2[t] ~ dcat(DU.prior.rCov[1:nLevels.rCov]) 
      }
      DU.prior.rCov[1:nLevels.rCov] <- 1/nLevels.rCov
      
    }else{
      
      for(t in 1:Tmax+1){
        RodentAbundance[t] ~ dnorm(0, sd = 1)
        RodentAbundance2[t] ~ dnorm(0, sd = 1)
      }
    }

    ## Missing covariate values in reindeer information
    for(t in 1:Tmax+1){
      Reindeer[t] ~ dnorm(0, sd = 1)
    }
     
  })
  
  return(redfox.code)
}

