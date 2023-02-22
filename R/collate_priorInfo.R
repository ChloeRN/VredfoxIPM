collate_priorInfo <- function(datafile, mu.t.max, maxAge, nsim, plot.mO = FALSE){
  
  
  ## Simulate / load prior distribution parameters from Hoening model
  if(file.exists("mO_prior_Parameters.rds")){
    mO.prior <- readRDS("mO_prior_Parameters.rds")
  }else{
    mO.prior <- predict_mO_HoeningMod(datafile, mu.t.max, maxAge, nsim = nsim, plot = plot.mO)
  }
  
  
  ## Assemble prior data
  prior.data <- list(
    
    # Early survival (denning period) - from arctic fox IPM
    earlySurv = list(
      S0.mean = 0.74, # From arctic fox IPM: 0.7428115
      S0.sd = 0.06 #  IPM: 0.05983706
    ),

    
    # Annual natural survival
    annSurv = list(
      #* Literature, Bristol (non-hunted)
      Bristol = list(
        Snat.mean = c(0.48, 0.54, 0.53, 0.51, 0.51),
        Snat.sd = c(0.02, 0.03, 0.03, 0.03, 0.03)
      ),
      
      
      #* Literature, North Sweden (lightly hunted)
      NSweden = list(
        Snat.mean = c(0.33, 0.71, 0.50, 0.59, 0.59),
        Snat.sd = c(0.02, 0.04, 0.05, 0.04, 0.04)
      ),
      
      
      #*Literature meta-analysis (all)
      metaAll = list(
        Snat.mean = c(0.38, 0.53, 0.57, 0.50, 0.52),
        Snat.sd = c(0.10, 0.16, 0.18, 0.20, 0.21)
      ),
      
      
      #*Literature meta-analysis (non/lightly hunted)
      metaSub = list(
        Snat.mean <- c(0.38, 0.53, 0.54, 0.46, 0.55),
        Snat.sd <- c(0.09, 0.10, 0.18, 0.25, 0.25)
      )
      # TODO: Check with Matt if we can re-do these meta-analyses properly
    ),
    
    
    # Natural mortality (Base parameters from Hoening model, ratio from arctic fox IPM)
    natMort = list(
      mnat.logmean = mO.prior$mnat.logmean,
      mnat.logsd = mO.prior$mnat.logsd,
      ratioJA.logmean = 0.4807439,
      ratioJA.logsd = 0.355224
    )
  )
  
  return(prior.data)
}