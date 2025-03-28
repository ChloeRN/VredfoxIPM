#' Write NIMBLE code for red fox IPM-PVA
#'
#' @param indLikelihood.genData logical. If TRUE, writes a model that includes 
#' an individual-level likelihood for genetic data which treats and individual's
#' true immigration status as the outcome of a Bernouilli trial with probability
#' corresponding to a "probability of belonging" derived from Geneclass 2 output.
#' If FALSE (default), writes a model that includes a group-level likelihood for
#' genetic data which treats the number of immigrants (identified using 
#' Geneclass 2 and a pre-defined p threshold) as the outcome of a Poisson trial
#' with lambda equal to the number of resident individuals times immigration rate.
#' For models that do not make use of genetic data, this choice is irrelevant.
#' 
#' @return an R call object specifying the model structure for the red fox IPM. 
#' @export
#'
#' @examples

writeCode_redfoxIPM_PVA <- function(indLikelihood.genData = FALSE){
  
  ## Check for incompatible toggles
  if(fitCov.mO & rCov.idx){
    stop("Incompatible model settings. Rodent effects on natural mortality (fitCov.mO = TRUE) are only implemented with continuous covariates (rCov.idx = FALSE).")
  }
  
  if(fitCov.mO & !mO.varT){
    warning("Attempting to fit a model containing environmental covariates but no random year variation for natural mortality. This is not recommended as it may result in inflated effect size/precision due to due to pseudo-replication.")
  }
  
  if(rCov.idx & (DD.mO | DD.immR)){
    stop("Incompatible model settings. The model is not set up to combine categorical rodent effects alongside density-dependence.")
  }
  
  if(!poolYrs.genData & imm.asRate & (comp.immR | DD.immR)){
    stop("Incompatible model settings. The year-specific genetic data likelihood (poolYrs.genData = TRUE) is currently not adapted to work with density and/or compensation effects.")
  }
  
  ## Write model code (individual-level likelihood for genetic data)
  if(indLikelihood.genData){
    
    if(fitCov.immR & !imm.asRate){
      stop("Incompatible model settings. 
         Covariate effects on immigration numbers can only be fit (fitCov.immR = TRUE) when the pooled genetic data likelihood is used (indLikelihood.genData = FALSE).")
    }
    
    if(useData.gen & !imm.asRate){
      warning("You have specified a model using genetic data with an individual likelihood but modelling immigration as a number, not a rate (imm.asRate = FALSE).
              This is currently not implemented. Your model will be run without genetic data.")
    }
    
    redfox.code <- nimbleCode({
      
      
      ##########################  
      #### POPULATION MODEL ####
      ##########################
      
      ### Likelihood (age classes: 1, 2, 3+)
      
      ## Survival
      
      for(t in 1:(Tmax+Tmax_sim)){
        # Age class 0 (index = 1): local reproduction & immigration
        N[1, t+1] <- sum(R[2:Amax, t+1]) + Imm[t+1]
        
        # Age classes 1 to 3 (indeces = 2, 3, 4): age classes 0, 1, and 2 survivors    
        for(a in 1:(Amax-2)){
          N[a+1, t+1] ~ dbin(S[a, t], N[a, t])
        }			
        
        # Age class 4+ (index = Amax = 5): age class 4 and 5+ survivors
        N[Amax, t+1] ~ dbin(S[Amax, t], N[Amax-1, t] + N[Amax, t])
      }
      
      
      ## Reproduction
      
      # Age class 0 (young of the year --> do not reproduce in year of birth)
      B[1, 1:(Tmax+Tmax_sim+1)] <- 0
      L[1, 1:(Tmax+Tmax_sim+1)] <- 0
      R[1, 1:(Tmax+Tmax_sim+1)] <- 0
      
      # First year (reproduction not modelled separately)
      B[2:Amax, 1] <- 0
      L[2:Amax, 1] <- 0
      R[2:Amax, 1] <- 0
      
      # Age classes 1 to 3+    	    
      for(t in 2:(Tmax+Tmax_sim+1)){        				
        
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
      
      for(t in 1:(Tmax+Tmax_sim+1)){
        N.tot[t] <- sum(N[1:Amax, t])
        R.tot[t] <- sum(R[1:Amax, t])		
        B.tot[t] <- sum(B[1:Amax, t])
        
        localN.tot[t] <- sum(R[2:Amax, t]) + sum(N[2:Amax, t])
      }
      
      #===============================================================================================
      
      
      
      ######################################
      #### WINTER AGE-AT-HARVEST MODULE ####
      ######################################
      
      ### Parameters:
      # N = number of individuals in a given age class at a given time (start of October)
      # h = age- and time-dependent probability of dying from winter hunting
      
      ### Data:
      # C_w = winter age-at-harvest matrix
      # pData_w = annual proportion of winter harvests with (complete) carcass data
      
      ### Likelihood
      
      for(t in 1:Tmax){
        for(a in 1:Amax){
          C_w[a, t] ~ dbin(h[a, t]*pData_w[t], N[a, t])
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
          if(imm.asRate){
            immR_pre[1:Tmax_Gen_pre] <- calculateImmR(ImmData = ImmData_pre[1:Xgen_pre], 
                                                      yearIdx = pImm_yrs_pre[1:Xgen_pre],
                                                      Tmax = Tmax_Gen_pre, skip_t1 = FALSE)
          }
        }
      }
      
      
      #===============================================================================================
      
      
      ################################
      #### PUP OBSERVATION MODULE ####
      ################################
      
      ## Likelihood (Number of pups)
      if(useData.pup){
        for(x in 1:X3){
          NoPups[x] ~ dpois(meanLS[NoPups_year[x]])
        }
        
        meanLS[1] <- 0
        for(t in 2:Tmax){
          meanLS[t] <- (sum(R[2:Amax, t])*2)/sum(B[2:Amax, t])
        }
      }
      
      #===============================================================================================
      
      
      ################################
      #### PRIORS AND CONSTRAINTS ####
      ################################
      
      ## Survival and mortality
      
      # Winter harvest mortality hazard rate
      for(t in 1:(Tmax+Tmax_sim+1)){
        if(fitCov.mH){
          mH[1:Amax, t] <- exp(log(Mu.mH[1:Amax]) + betaHE.mH*HarvestEffort[t] + epsilon.mH[t])*pertFac.mH[t]*pertFac.mH.flex[t]
        }else{
          mH[1:Amax, t] <- exp(log(Mu.mH[1:Amax]) + epsilon.mH[t])*pertFac.mH[t]*pertFac.mH.flex[t]
        }
        
        # Log mean deviation of harvest mortality hazard rate
        logDev.mH[t] <- log(mH[1, t]) - log(Mu.mH[1])
      }
      
      for(t in 1:(Tmax+Tmax_sim)){ 

        # Other (natural) mortality hazard rate
        if(fitCov.mO){
          
          # First age class
          mO[1, t] <- exp(log(Mu.mO[1]) + 
            betaR.mO*RodentAbundance[t+1] + 
            betaD.mO*(log(localN.tot[t]) - log(normN)) + 
            betaRxD.mO*RodentAbundance[t+1]*(log(localN.tot[t]) - log(normN)) + 
            gamma.mO*logDev.mH[t] + 
            epsilon.mO[t])*pertFac.mO[t]
          
          # Other age classes
          mO[2:Amax, t] <- exp(log(Mu.mO[2:Amax]) + 
            betaR.mO*RodentAbundance[t+1] +
            gamma.mO*logDev.mH[t] +
            epsilon.mO[t])*pertFac.mO[t]
          
          # All age classes
          # mO[1:Amax, t] <- exp(log(Mu.mO[1:Amax]) + 
          #   betaR.mO*RodentAbundance[t+1] + 
          #   betaD.mO*(log(localN.tot[t]) - log(normN)) + 
          #   betaRxD.mO*RodentAbundance[t+1]*(log(localN.tot[t]) - log(normN)) + 
          #   gamma.mO*logDev.mH[t] + 
          #   epsilon.mO[t])*pertFac.mO[t]
          
        }else{
          mO[1:Amax, t] <- exp(log(Mu.mO[1:Amax]) + epsilon.mO[t])*pertFac.mO[t]
        }
        
        # Survival probability
        S[1:Amax, t] <- exp(-(mH[1:Amax, t] + mO[1:Amax,t]))
        
        # Proportion winter harvest mortality
        alpha[1:Amax, t] <- mH[1:Amax, t]/(mH[1:Amax, t] + mO[1:Amax, t])
        
        # Winter harvest rate
        h[1:Amax, t] <- (1-S[1:Amax, t])*alpha[1:Amax, t]
        
      }
      
      # Median harvest mortality hazard rates
      
      # Age-dependent
      for(a in 1:2){
        Mu.mH[a] ~ dunif(0, 5)
      }
      
      Mu.mH[3:Amax] <- Mu.mH[2]
      
      # Age-independent   
      #Mu.mH.all ~ dunif(0, 5) 
      #Mu.mH[1:Amax] <- Mu.mH.all
      
      # Median other (natural) cause mortality hazard rates
      #* INFORMATIVE PRIOR REQUIRED: LITERATURE VALUE / HOENIG MODEL CALCULATION
      
      if(HoenigPrior){
        # Using prior distributions calculated with Hoenig model
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
      
      if(fitCov.mO){
        betaR.mO ~ dunif(-5, 5)
        
        if(DD.mO){
          betaD.mO ~ dunif(-5, 5)
          if(DDxRodent){
            betaRxD.mO ~ dunif(-5, 5)
          }else{
            betaRxD.mO <- 0
          }
        }else{
          betaD.mO <- 0
          betaRxD.mO <- 0
        }
        
        if(comp.mO & !comp.RE){
          gamma.mO ~ dunif(-5, 5)
        }else{
          gamma.mO <- 0
        }
        
      }
      
      
      #---------------------------------------------------------------------------------------------
      
      
      ## Pregnancy rate
      
      for(t in 1:(Tmax+Tmax_sim+1)){
        Psi[1, t] <- 0
        
        if(fitCov.Psi){
          if(rCov.idx){
            logit(Psi[2:Amax,t]) <- logit(Mu.Psi[2:Amax]) + betaR.Psi[RodentIndex[t]] + epsilon.Psi[t] # Reindeer.rodent interaction not (yet) written in
          }else{
            logit(Psi[2:Amax,t]) <- logit(Mu.Psi[2:Amax]) + betaR.Psi*RodentAbundance_pert[t] + epsilon.Psi[t]
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
      
      
      for(t in 1:(Tmax+Tmax_sim+1)){
        rho[1, t] <- 0
        
        if(fitCov.rho){
          if(rCov.idx){
            log(rho[2:Amax, t]) <- log(Mu.rho[2:Amax]) + betaR.rho[RodentIndex[t]] + epsilon.rho[t]
          }else{
            log(rho[2:Amax, t]) <- log(Mu.rho[2:Amax]) + betaR.rho*RodentAbundance_pert[t] + epsilon.rho[t]
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
      
      for(t in 1:(Tmax+Tmax_sim+1)){ 
        S0[t] <- Mu.S0*pertFac.S0

        #S0[t] <- exp(-m0[t])
        #log(m0[t]) <- log(-log(Mu.S0)) + epsilon.m0[t]
      }
      
      if(useInfPrior.S0){
        Mu.S0 ~ T(dnorm(S0.mean, sd = S0.sd), 0, 1)
      }else{
        Mu.S0 ~ dunif(0, 1)
      }
      
      #---------------------------------------------------------------------------------------------
      
      
      ## Immigration
      
      if(imm.asRate){
        
        if(useData.gen & !poolYrs.genData){
          
          ## Extraction of log mean and sd for immigration rate
          log(Mu.immR) <- log(mean(c(immR[1:Tmax_Gen], immR_pre[1:Tmax_Gen_pre]))) 
          sigma.immR <- calculateLogSD(immR = c(immR[1:Tmax_Gen], immR_pre[1:Tmax_Gen_pre]),
                                       replace0 = 0.01)
          
          ## Projection of immigration rates beyond genetic data coverage
          for(t in (Tmax_Gen+1):(Tmax+Tmax_sim)){
            immR[t] <- exp(log(Mu.immR) + epsilon.immR[t])*pertFac.immR[t]
          }
          
        }else{
          
          if(fitCov.immR){
            if(rCov.idx){
              for(t in 1:(Tmax+Tmax_sim+1)){
                immR[t] <- exp(log(Mu.immR) + betaR.immR[RodentIndex2[t]] + epsilon.immR[t])*pertFac.immR[t]
              }
            }else{
              
              for(t in 1:(Tmax+Tmax_sim+1)){
                immR[t] <- exp(log(Mu.immR) + 
                  betaR.immR*RodentAbundance2[t] + 
                  betaD.immR*(log(localN.tot[t]) - log(normN)) + 
                  betaRxD.immR*RodentAbundance2[t]*(log(localN.tot[t]) - log(normN)) + 
                  gamma.immR*logDev.mH[t] +
                  epsilon.immR[t])*pertFac.immR[t]
              }
              
            }
          }else{
            immR[1:(Tmax+Tmax_sim+1)] <- exp(log(Mu.immR) + epsilon.immR[1:(Tmax+Tmax_sim+1)])*pertFac.immR[1:(Tmax+Tmax_sim+1)]
          }
          
          Mu.immR ~ dunif(0, 10)
        }
        
        for(t in 1:(Tmax+Tmax_sim)){ 
          Imm[t] ~ dpois(R.tot[t]*immR[t])
        }
        
        
      }else{
        
        ## Lognormal prior for immigrant numbers
        for(t in 2:(Tmax+Tmax_sim+1)){
          
          Imm[t] <- round(ImmExp2[t]*pertFac.immR[t])
          if(Imm.logNorm){
            ImmExp2[t] <- ImmExp[t]
          }else{
            ImmExp2[t] ~ T(dnorm(mean = ImmExp[t], sd = sigma.Imm), 0, Inf)
          }
          
          if(fitCov.immR){
            log(ImmExp[t]) <- log(Mu.Imm) + 
              betaR.immR*RodentAbundance2[t] + 
              betaD.immR*(log(localN.tot[t]) - log(normN)) + 
              betaRxD.immR*RodentAbundance2[t]*(log(localN.tot[t]) - log(normN)) +
              gamma.immR*logDev.mH[t] +
              epsilon.immR[t] 
          }else{
            Imm[t] <- round(ImmExp[t]*pertFac.immR[t])
            log(ImmExp[t]) <- log(Mu.Imm) + epsilon.immR[t] 
          }

        }
        
        Mu.Imm ~ dunif(1, uLim.Imm)

        ## Derivation of immigration rates
        immR[1] <- 0
        for(t in 2:(Tmax+Tmax_sim+1)){
          immR[t] <- Imm[t] / R.tot[t]
        }
        Mu.immR <- mean(immR[2:(Tmax+1)])
      }
      
      ## Priors for covariate effects
      if(fitCov.immR){
        if(rCov.idx){
          betaR.immR[1] <- 0 # --> Lowest level corresponds to intercept
          for(x in 2:nLevels.rCov){
            betaR.immR[x] ~ dunif(-5, 5)
          }
        }else{
          betaR.immR ~ dunif(-5, 5)
          
          if(DD.immR){
            betaD.immR ~ dunif(-100, 5)
            if(DDxRodent){
              betaRxD.immR ~ dunif(-5, 5)
            }else{
              betaRxD.immR <- 0
            }
          }else{
            betaD.immR <- 0
            betaRxD.immR <- 0
          }
          
          if(comp.immR & !comp.RE){
            gamma.immR ~ dunif(-5, 5)
          }else{
            gamma.immR <- 0
          }
          
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
      
      for(t in 1:(Tmax+Tmax_sim+1)){
        #epsilon.mH[t] ~ dnorm(0, sd = sigma.mH)
        epsilon.mH[t] <- sigma.mH*eta.mH[t]
        eta.mH[t] ~ dnorm(0, sd = 1)
        
        #epsilon.mO[t] ~ dnorm(0, sd = sigma.mO)
        epsilon.mO[t] <- eta.mO[t] + tau.mO*eta.mH[t]
        eta.mO[t] ~ dnorm(0, sd = sigma.mO)
        
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
      
      for(t in 1:(Tmax+Tmax_sim+1)){
        #epsilon.immR[t] ~ dnorm(0, sd = sigma.immR)
        epsilon.immR[t] <- eta.immR[t] + tau.immR*eta.mH[t]
        eta.immR[t] ~ dnorm(0, sd = sigma.immR)
      }
      
      if(!imm.asRate & !Imm.logNorm){
        sigma.immR <- 0
        sigma.Imm ~ dunif(0, uLim.Imm)
      }else{
        sigma.immR ~ dunif(0, 10)
        sigma.Imm <- 0
      }
      
      if(comp.mO & comp.RE){
        tau.mO ~ dnorm(0, sd = 2.25)
      }else{
        tau.mO <- 0
      }
      
      if(comp.immR & comp.RE){
        tau.immR ~ dnorm(0, sd = 2.25)
      }else{
        tau.immR <- 0
      }
      
      # Calculation of correlation coefficients
      if(!imm.asRate & !Imm.logNorm){
        C.immR <- 0
      }else{
        C.immR <- tau.immR / sqrt(pow(sigma.immR, 2) + pow(tau.immR, 2)) 
      }
      C.mO <- tau.mO / sqrt(pow(sigma.mO, 2) + pow(tau.mO, 2)) 
      
      
      #===============================================================================================
      
      
      
      ##############################
      #### COVARIATE IMPUTATION ####
      ##############################
      
      ## Missing covariate value(s) in number of successful hunters
      if(fitCov.mH){
        for(t in 1:(Tmax+Tmax_sim)){
          HarvestEffort[t] ~ dnorm(0, sd = 1)
        }
      }
      
      ## Missing covariate value(s) in rodent abundance
      if(rCov.idx){
        
        for(t in 1:(Tmax+Tmax_sim+1)){
          RodentIndex[t] ~ dcat(DU.prior.rCov[1:nLevels.rCov]) 
          RodentIndex2[t] ~ dcat(DU.prior.rCov[1:nLevels.rCov]) 
        }
        DU.prior.rCov[1:nLevels.rCov] <- 1/nLevels.rCov
        
      }else{
        
        # Naive normal distribution for NAs early in time series
        for(t in 1:2){
          RodentAbundance[t] ~ dnorm(0, sd = 1)
        }
        
        # Varanger: Second-order autoregressive & correlation models for remainder of time series
        for(t in 3:(Tmax+Tmax_sim+1)){

          RodentAbundance_pred[t] <- beta.RodMod[1]*RodentAbundance[t-1] + 
            beta.RodMod[2]*RodentAbundance[t-2] + 
            beta.RodMod[3]*RodentAbundance[t-1]*RodentAbundance[t-2]
          RodentAbundance[t] ~ dnorm(RodentAbundance_pred[t], sd = sigmaT.RodAbun)
        }
          
        # Larger area: Correlative model for rodent dynamics in the larger area
        for(t in 1:(Tmax+Tmax_sim)){
          RodentAbundance2[t] ~ dnorm(beta.RodCorr*RodentAbundance[t+1], sd = sigmaT.RodAbun2)
        }
      }
      
      # Add perturbation 
      for(t in 1:(Tmax+Tmax_sim+1)){
        RodentAbundance_pert[t] <- RodentAbundance[t] + (1-pertFac.rodent[t])
      }
      for(t in 1:(Tmax+Tmax_sim)){
        RodentAbundance2_pert[t] <- RodentAbundance2[t] + (1-pertFac.rodent[t])
      }
      
      # Additional priors for rodent model
      for(i in 1:3){
        beta.RodMod[i] ~ dunif(-5, 5)
      }
      beta.RodCorr ~ dunif(0, 1)
      sigmaT.RodAbun ~ dunif(0, 5)
      sigmaT.RodAbun2 ~ dunif(0, 5)

      
      #########################################
      #### PERTURBATION FACTOR CALCULATION ####
      #########################################
      
      ## Rodent-dependent harvest perturbation factor
      pertFac.mH.flex[1:Tmax] <- 1
      
      if(Tmax_sim > 0){
        for(t in (Tmax+1):(Tmax+Tmax_sim+1)){
          pertFac.mH.flex[t] <- calculate_pertFac(pertFactor = factor.mH.rodent,
                                                  covThreshold = threshold.rodent.mH,
                                                  thresholdAbove = thresholdAbove,
                                                  covValue = RodentAbundance_pert[t])
        }
      }
      
      
    })
  }
  
  #=============================================================================
  #=============================================================================
  #=============================================================================
  
  
  ## Write model code (group-level likelihood for genetic data)
  if(!indLikelihood.genData){
    redfox.code <- nimbleCode({
      
      
      ##########################  
      #### POPULATION MODEL ####
      ##########################
      
      ### Likelihood (age classes: 1, 2, 3+)
      
      ## Survival
      
      #---------------------------#
      # OCT - JUN (AUTUMN-SPRING) #
      #---------------------------#
      
      for(t in 1:(Tmax+Tmax_sim)){
        # Age class 0 (index = 1): local reproduction & immigration
        N[1, t+1] <- sum(R[2:Amax, t+1]) + Imm[t+1]
        
        # Age classes 1 to 3 (indeces = 2, 3, 4): age classes 0, 1, and 2 survivors    
        for(a in 1:(Amax-2)){
          N[a+1, t+1] ~ dbin(S[a, t], N[a, t])
        }			
        
        # Age class 4+ (index = Amax = 5): age class 4 and 5+ survivors
        N[Amax, t+1] ~ dbin(S[Amax, t], N[Amax-1, t] + N[Amax, t])
      }
      
      
      ## Reproduction
      
      # Age class 0 (young of the year --> do not reproduce in year of birth)
      B[1, 1:(Tmax+Tmax_sim+1)] <- 0
      L[1, 1:(Tmax+Tmax_sim+1)] <- 0
      R[1, 1:(Tmax+Tmax_sim+1)] <- 0
      
      # First year (reproduction not modelled separately)
      B[2:Amax, 1] <- 0
      L[2:Amax, 1] <- 0
      R[2:Amax, 1] <- 0
      
      # Age classes 1 to 3+    	    
      for(t in 2:(Tmax+Tmax_sim+1)){        				
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
      
      for(t in 1:(Tmax+Tmax_sim+1)){
        N.tot[t] <- sum(N[1:Amax, t])
        R.tot[t] <- sum(R[1:Amax, t])		
        B.tot[t] <- sum(B[1:Amax, t])
        
        localN.tot[t] <- sum(R[2:Amax, t]) + sum(N[2:Amax, t])
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
          C_w[a, t] ~ dbin(h[a, t]*pData_w[t], N[a, t])
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
      if(useData.gen){
        
        if(poolYrs.genData){
          
          genObs_Imm ~ dpois(genObs_Res*Mu.immR)
          
        }else{
          
          # Within study period 
          for(t in 1:Tmax_Gen){
            genObs_Imm[t] ~ dpois(genObs_Res[t]*immR[t])
          }
          
          # Outside study period
          if(imm.asRate){
            for(t in 1:Tmax_Gen_pre){
              genObs_Imm_pre[t] ~ dpois(genObs_Res_pre[t]*immR_pre[t])
            }
          }
        }
      }
      
      #===============================================================================================
      
      
      ################################
      #### PUP OBSERVATION MODULE ####
      ################################
      
      ## Likelihood (Number of pups)
      if(useData.pup){
        for(x in 1:X3){
          NoPups[x] ~ dpois(meanLS[NoPups_year[x]])
        }
        
        meanLS[1] <- 0
        for(t in 2:Tmax){
          meanLS[t] <- (sum(R[2:Amax, t])*2)/sum(B[2:Amax, t])
        }
      }

      
      #===============================================================================================
      
      
      ################################
      #### PRIORS AND CONSTRAINTS ####
      ################################
      
      ## Survival and mortality
      
      # Winter harvest mortality hazard rate
      for(t in 1:(Tmax+Tmax_sim+1)){
        
        if(fitCov.mH){
          mH[1:Amax, t] <- exp(log(Mu.mH[1:Amax]) + betaHE.mH*HarvestEffort[t] + epsilon.mH[t])*pertFac.mH[t]*pertFac.mH.flex[t]
        }else{
          mH[1:Amax, t] <- exp(log(Mu.mH[1:Amax]) + epsilon.mH[t])*pertFac.mH[t]*pertFac.mH.flex[t]
        }
        
        # Log mean deviation of harvest mortality hazard rate
        logDev.mH[t] <- log(mH[1, t]) - log(Mu.mH[1])
      }
      
      for(t in 1:(Tmax+Tmax_sim)){ 
        
        # Other (natural) mortality hazard rate
        if(fitCov.mO){
          # First age class
          mO[1, t] <- exp(log(Mu.mO[1]) + 
            betaR.mO*RodentAbundance[t+1] + 
            betaD.mO*(log(localN.tot[t]) - log(normN)) + 
            betaRxD.mO*RodentAbundance[t+1]*(log(localN.tot[t]) - log(normN)) + 
            gamma.mO*logDev.mH[t] + 
            epsilon.mO[t])*pertFac.mO[t]
          
          # Other age classes
          mO[2:Amax, t] <- exp(log(Mu.mO[2:Amax]) + 
            betaR.mO*RodentAbundance[t+1] +
            gamma.mO*logDev.mH[t] +
            epsilon.mO[t])*pertFac.mO[t]
          
          # All age classes
          # mO[1:Amax, t] <- exp(log(Mu.mO[1:Amax]) + 
          #   betaR.mO*RodentAbundance[t+1] + 
          #   betaD.mO*(log(localN.tot[t]) - log(normN)) + 
          #   betaRxD.mO*RodentAbundance[t+1]*(log(localN.tot[t]) - log(normN)) + 
          #   gamma.mO*logDev.mH[t] + 
          #   epsilon.mO[t])*pertFac.mO[t]
          
        }else{
          mO[1:Amax, t] <- exp(log(Mu.mO[1:Amax]) + epsilon.mO[t])*pertFac.mO[t]
        }
        
        # Survival probability
        S[1:Amax, t] <- exp(-(mH[1:Amax, t] + mO[1:Amax,t]))
        
        # Proportion winter harvest mortality
        alpha[1:Amax, t] <- mH[1:Amax, t]/(mH[1:Amax, t] + mO[1:Amax, t])
        
        # Winter harvest rate
        h[1:Amax, t] <- (1-S[1:Amax, t])*alpha[1:Amax, t]
        
      }
      
      # Median harvest mortality hazard rates
      
      # Age-dependent
      for(a in 1:2){
        Mu.mH[a] ~ dunif(0, 5)
      }
      
      Mu.mH[3:Amax] <- Mu.mH[2]
      
      # Age-independent   
      #Mu.mH.all ~ dunif(0, 5) 
      #Mu.mH[1:Amax] <- Mu.mH.all
      
      # Median other (natural) cause mortality hazard rates
      #* INFORMATIVE PRIOR REQUIRED: LITERATURE VALUE / HOENIG MODEL CALCULATION
      
      if(HoenigPrior){
        # Using prior distributions calculated with Hoenig model
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
      
      if(fitCov.mO){
        betaR.mO ~ dunif(-5, 5)
        
        if(DD.mO){
          betaD.mO ~ dunif(-5, 5)
          if(DDxRodent){
            betaRxD.mO ~ dunif(-5, 5)
          }else{
            betaRxD.mO <- 0
          }
        }else{
          betaD.mO <- 0
          betaRxD.mO <- 0
        }
        
        if(comp.mO & !comp.RE){
          gamma.mO ~ dunif(-5, 5)
        }else{
          gamma.mO <- 0
        }
        
      }
      
      
      #---------------------------------------------------------------------------------------------
      
      
      ## Pregnancy rate
      
      for(t in 1:(Tmax+Tmax_sim+1)){
        Psi[1, t] <- 0
        
        if(fitCov.Psi){
          if(rCov.idx){
            logit(Psi[2:Amax,t]) <- logit(Mu.Psi[2:Amax]) + betaR.Psi[RodentIndex[t]] + epsilon.Psi[t] # Reindeer.rodent interaction not (yet) written in
          }else{
            logit(Psi[2:Amax,t]) <- logit(Mu.Psi[2:Amax]) + betaR.Psi*RodentAbundance_pert[t] + epsilon.Psi[t]
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
      
      
      for(t in 1:(Tmax+Tmax_sim+1)){
        rho[1, t] <- 0
        
        if(fitCov.rho){
          if(rCov.idx){
            log(rho[2:Amax, t]) <- log(Mu.rho[2:Amax]) + betaR.rho[RodentIndex[t]] + epsilon.rho[t]
          }else{
            log(rho[2:Amax, t]) <- log(Mu.rho[2:Amax]) + betaR.rho*RodentAbundance_pert[t] + epsilon.rho[t]
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
      
      for(t in 1:(Tmax+Tmax_sim+1)){ 
        S0[t] <- Mu.S0*pertFac.S0[t]

        #S0[t] <- exp(-m0[t])
        #log(m0[t]) <- log(-log(Mu.S0)) + epsilon.m0[t]
      }
      
      if(useInfPrior.S0){
        Mu.S0 ~ T(dnorm(S0.mean, sd = S0.sd), 0, 1)
      }else{
        Mu.S0 ~ dunif(0, 1)
      }
      
      #---------------------------------------------------------------------------------------------
      
      
      ## Immigration (within the study period)
      
      if(imm.asRate){
        
        if(fitCov.immR){
          if(rCov.idx){
            for(t in 1:(Tmax+Tmax_sim+1)){
              immR[t] <- exp(log(Mu.immR) + betaR.immR[RodentIndex2[t]] + epsilon.immR[t])*pertFac.immR[t]
            }
          }else{
            
            for(t in 1:(Tmax+Tmax_sim+1)){
              immR[t] <- exp(log(Mu.immR) + 
                               betaR.immR*RodentAbundance2[t] + 
                               betaD.immR*(log(localN.tot[t]) - log(normN)) + 
                               betaRxD.immR*RodentAbundance2[t]*(log(localN.tot[t]) - log(normN)) + 
                               gamma.immR*logDev.mH[t] +
                               epsilon.immR[t])*pertFac.immR[t]
            }
          }
          
        }else{
          immR[1:(Tmax+Tmax_sim+1)] <- exp(log(Mu.immR) + epsilon.immR[1:(Tmax+Tmax_sim+1)])*pertFac.immR[1:(Tmax+Tmax_sim+1)]
        }
        
        Mu.immR ~ dunif(0, 10)  
          
        for(t in 1:(Tmax+Tmax_sim+1)){ 
          Imm[t] ~ dpois(R.tot[t]*immR[t])
        }
        
      }else{
        
        ## Lognormal prior for immigrant numbers
        
        for(t in 2:(Tmax+Tmax_sim+1)){
          
          Imm[t] <- round(ImmExp2[t]*pertFac.immR[t])
          if(Imm.logNorm){
            ImmExp2[t] <- ImmExp[t]
          }else{
            ImmExp2[t] ~ T(dnorm(mean = ImmExp[t], sd = sigma.Imm), 0, Inf)
            #ImmExp2[t] ~ dnorm(mean = ImmExp[t], sd = sigma.Imm)
          }
          
          if(fitCov.immR){
            log(ImmExp[t]) <- log(Mu.Imm) + 
              betaR.immR*RodentAbundance2[t] + 
              betaD.immR*(log(localN.tot[t]) - log(normN)) + 
              betaRxD.immR*RodentAbundance2[t]*(log(localN.tot[t]) - log(normN)) +
              gamma.immR*logDev.mH[t] +
              epsilon.immR[t]
          }else{
            log(ImmExp[t]) <- log(Mu.Imm) + epsilon.immR[t]
          }
          
        }

        Mu.Imm ~ dunif(1, uLim.Imm)
        
        ## Derivation of immigration rates
        immR[1] <- 0
        for(t in 2:(Tmax+Tmax_sim+1)){
          immR[t] <- Imm[t] / R.tot[t]
        }
        Mu.immR <- mean(immR[2:(Tmax+1)])
        
      }
      
      ## Prior for covariate effects
      if(fitCov.immR){
        if(rCov.idx){
          betaR.immR[1] <- 0 # --> Lowest level corresponds to intercept
          for(x in 2:nLevels.rCov){
            betaR.immR[x] ~ dunif(-5, 5)
          }
        }else{
          betaR.immR ~ dunif(-5, 5)
          
          if(DD.immR){
            betaD.immR ~ dunif(-100, 5)
            if(DDxRodent){
              betaRxD.immR ~ dunif(-10, 10)
            }else{
              betaRxD.immR <- 0
            }
          }else{
            betaD.immR <- 0
            betaRxD.immR <- 0
          }
          
          if(comp.immR & !comp.RE){
            gamma.immR ~ dunif(-5, 5)
          }else{
            gamma.immR <- 0
          }
        }
      }
      
      
      ## Immigration outside the study period
      if(useData.gen & !poolYrs.genData){
        
        for(t in 1:Tmax_Gen_pre){
          
          if(fitCov.immR){
            if(rCov.idx){
              immR_pre[t] ~ dlnorm(log(Mu.immR) + betaR.immR[RodentIndex2_pre[t]], sdlog = sigma.immR)
            }else{
              immR_pre[t] ~ dlnorm(log(Mu.immR) + betaR.immR*RodentAbundance2_pre[t], sdlog = sigma.immR)
            }
          }else{
            immR_pre[t] ~ dlnorm(log(Mu.immR), sdlog = sigma.immR)
          }
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

      for(t in 1:(Tmax+Tmax_sim+1)){
        
        #epsilon.mH[t] ~ dnorm(0, sd = sigma.mH)
        epsilon.mH[t] <- sigma.mH*eta.mH[t]
        eta.mH[t] ~ dnorm(0, sd = 1)
        
        #epsilon.mO[t] ~ dnorm(0, sd = sigma.mO)
        epsilon.mO[t] <- eta.mO[t] + tau.mO*eta.mH[t]
        eta.mO[t] ~ dnorm(0, sd = sigma.mO)
        
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
      
      for(t in 1:(Tmax+Tmax_sim+1)){
        
        epsilon.immR[t] <- eta.immR[t] + tau.immR*eta.mH[t]
        
        if(Imm.logNorm){
          #epsilon.immR[t] ~ dnorm(0, sd = sigma.immR)
          eta.immR[t] ~ dnorm(0, sd = sigma.immR)
        }else{
          eta.immR[t] <- 0
        }
        
      }
      
      if(!imm.asRate & !Imm.logNorm){
        sigma.immR <- 0
        sigma.Imm ~ dunif(0, uLim.Imm)
      }else{
        sigma.immR ~ dunif(0, 10)
        sigma.Imm <- 0
      }
      
      if(comp.mO & comp.RE){
        tau.mO ~ dnorm(0, sd = 2.25)
      }else{
        tau.mO <- 0
      }
      
      if(comp.immR & comp.RE){
        tau.immR ~ dnorm(0, sd = 2.25)
      }else{
        tau.immR <- 0
      }
      
      # Calculation of correlation coefficients
      if(!imm.asRate & !Imm.logNorm){
        C.immR <- 0
      }else{
        C.immR <- tau.immR / sqrt(pow(sigma.immR, 2) + pow(tau.immR, 2)) 
      }
      C.mO <- tau.mO / sqrt(pow(sigma.mO, 2) + pow(tau.mO, 2)) 
      
      #===============================================================================================
      
      
      
      ##############################
      #### COVARIATE IMPUTATION ####
      ##############################
      
      ## Missing covariate value(s) in number of successful hunters
      if(fitCov.mH){
        for(t in 1:(Tmax+Tmax_sim)){
          HarvestEffort[t] ~ dnorm(0, sd = 1)
        }
      }
      
      ## Missing covariate value(s) in rodent abundance
      if(rCov.idx){
        
        for(t in 1:(Tmax+Tmax_sim+1)){
          RodentIndex[t] ~ dcat(DU.prior.rCov[1:nLevels.rCov]) 
          RodentIndex2[t] ~ dcat(DU.prior.rCov[1:nLevels.rCov]) 
        }
        DU.prior.rCov[1:nLevels.rCov] <- 1/nLevels.rCov
        
      }else{
        
        # Naive normal distribution for NAs early in time series
        for(t in 1:2){
          RodentAbundance[t] ~ dnorm(0, sd = 1)
        }
        
        # Varanger: Second-order autoregressive & correlation models for remainder of time series
        for(t in 3:(Tmax+Tmax_sim+1)){
          
          RodentAbundance_pred[t] <- beta.RodMod[1]*RodentAbundance[t-1] + 
            beta.RodMod[2]*RodentAbundance[t-2] + 
            beta.RodMod[3]*RodentAbundance[t-1]*RodentAbundance[t-2]
          RodentAbundance[t] ~ dnorm(RodentAbundance_pred[t], sd = sigmaT.RodAbun)
        }
        
        # Larger area: Correlative model for rodent dynamics in the larger area
        for(t in 1:(Tmax+Tmax_sim)){
          RodentAbundance2[t] ~ dnorm(beta.RodCorr*RodentAbundance[t+1], sd = sigmaT.RodAbun2)
        }
      }
      
      # Add perturbation 
      for(t in 1:(Tmax+Tmax_sim+1)){
        RodentAbundance_pert[t] <- RodentAbundance[t] + (1-pertFac.rodent[t])
      }
      for(t in 1:(Tmax+Tmax_sim)){
        RodentAbundance2_pert[t] <- RodentAbundance2[t] + (1-pertFac.rodent[t])
      }
      
      # Additional priors for rodent model
      for(i in 1:3){
        beta.RodMod[i] ~ dunif(-5, 5)
      }
      beta.RodCorr ~ dunif(0, 1)
      sigmaT.RodAbun ~ dunif(0, 5)
      sigmaT.RodAbun2 ~ dunif(0, 5)

      
      #########################################
      #### PERTURBATION FACTOR CALCULATION ####
      #########################################
      
      ## Rodent-dependent harvest perturbation factor
      pertFac.mH.flex[1:Tmax] <- 1
      
      if(Tmax_sim > 0){
        for(t in (Tmax+1):(Tmax+Tmax_sim+1)){
          pertFac.mH.flex[t] <- calculate_pertFac(pertFactor = factor.mH.rodent,
                                                  covThreshold = threshold.rodent.mH,
                                                  thresholdAbove = thresholdAbove,
                                                  covValue = RodentAbundance[t])
        }
      }
      
      
    })
  }
  
  return(redfox.code)
}

