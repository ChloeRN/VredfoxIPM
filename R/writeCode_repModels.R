#' Write NIMBLE code for independent reproduction models
#'
#' @return an R call object specifying the model structure for the red fox IPM. 
#' @export
#'
#' @examples

writeCode_repModels <- function(){
  
  redfox.code <- nimbleCode({
    
    
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
    
    
    
    ################################
    #### PRIORS AND CONSTRAINTS ####
    ################################
    
 
    ## Pregnancy rate
    
    for(t in 1:(Tmax+1)){
      Psi[1, t] <- 0
      
      if(fitCov.Psi){
        if(rCov.idx){
          logit(Psi[2:Amax,t]) <- logit(Mu.Psi[2:Amax]) + betaR.Psi[RodentIndex[t]+1] + epsilon.Psi[t]
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
          log(rho[2:Amax,t]) <- log(Mu.rho[2:Amax]) + betaR.rho[RodentIndex[t]+1] + epsilon.rho[t]
        }else{
          log(rho[2:Amax,t]) <- log(Mu.rho[2:Amax]) + betaR.rho*RodentAbundance[t] + epsilon.rho[t]
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
    
    
    ## Random year variation
    
    for(t in 1:Tmax){  
      epsilon.Psi[t] ~ dnorm(0, sd = sigma.Psi)
      epsilon.rho[t] ~ dnorm(0, sd = sigma.rho) 
    }
    
    sigma.Psi ~ dunif(0, 5)
    sigma.rho ~ dunif(0, 5)

  })
  
  return(redfox.code)
}

