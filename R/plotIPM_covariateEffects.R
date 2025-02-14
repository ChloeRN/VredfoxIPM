#' Plot vital rates as functions of fitted covariates
#'
#' @param MCMC.samples an mcmc.list containing the output of the fitted IPM. 
#' @param rCov.idx logical. If TRUE, stops function execution as this function 
#' has been written for the case of continuous covariates (rCov.idx = FALSE).
#' @param rodentMIN numeric. Minimum value of z-standardized rodent covariate to 
#' use for predictions. 
#' @param rodentMAX numeric. Maximum value of z-standardized rodent covariate to 
#' use for predictions. 
#' @param AgeClass integer. Age class for which to make predictions. Implemented
#' as "years of age" (range 0-4), not as model age index (range 1-5). 
#'
#' @return character vector of plot names. The plots themselves are saved
#' as pdf's in the subfolder "Plots".
#' @export
#'
#' @examples

plotIPM_covariateEffects <- function(MCMC.samples, rCov.idx, rodentMIN, rodentMAX, AgeClass){
  
  ## Check if continuous rodent covariates have been used (categorical not supported yet)
  if(rCov.idx){
    stop("Plotting of covariate effect relationships not (yet) supported for models fit with categorical covariates.")
  }
  
  ## Create plot folder if not present
  if(!dir.exists("Plots")){dir.create("Plots")}
  
  #-----------------------------------#
  # Make rodent covariate predictions #
  #-----------------------------------#
  
  ## Make vector of values for covariates for which to predict for
  rodentCov <- seq(rodentMIN, rodentMAX, length.out = 100)
  
  ## Set up empty dataframes to store results
  covPred.data <- data.frame()
  
  ## Convert MCMC samples to matrix
  out.mat <- as.matrix(MCMC.samples)
  
  ## Make posterior predictions for rodent covariate effects
  for(x in 1:length(rodentCov)){
    
    mO.pred <- Psi.pred <- rho.pred <- imm.pred <- rep(NA, nrow(out.mat))

    for(i in 1:nrow(out.mat)){
      
      if(AgeClass > 0){
        # Pregnancy rate
        Psi.pred[i] <- plogis(qlogis(out.mat[i, paste0("Mu.Psi[", AgeClass + 1, "]")]) + 
                                out.mat[i, "betaR.Psi"]*rodentCov[x])
        
        # Fetus number
        rho.pred[i] <- exp(log(out.mat[i, paste0("Mu.rho[", AgeClass + 1, "]")]) + 
                             out.mat[i, "betaR.rho"]*rodentCov[x])
      }else{
        Psi.pred[i] <- 0
        rho.pred[i] <- 0
      }
      
      # Immigration rate
      imm.pred[i] <- exp(log(out.mat[i, "Mu.immR"]) + 
                           out.mat[i, "betaR.immR"]*rodentCov[x])
      
      # Natural morality 
      mO.pred[i] <- exp(log(out.mat[i, paste0("Mu.mO[", AgeClass + 1, "]")]) + 
                          out.mat[i, "betaR.mO"]*rodentCov[x])
    }
    
    # Assemble in temporary data frame
    data.temp <- rbind(quantile(Psi.pred, probs = c(0.025, 0.5, 0.975)),
                       quantile(rho.pred, probs = c(0.025, 0.5, 0.975)),
                       quantile(imm.pred, probs = c(0.025, 0.5, 0.975)),
                       quantile(mO.pred, probs = c(0.025, 0.5, 0.975))) %>%
      as.data.frame() %>%
      dplyr::rename(lCI = 1, median = 2, uCI = 3) %>%
      dplyr::mutate(VitalRate = c("Pregnancy rate", "Fetus number", "Immigration rate", "Natural mortality"),
                    RodentAbundance = rodentCov[x])
    
    # Combine in main data frame
    covPred.data <- rbind(covPred.data, data.temp)
  }
  
  ## Re-arrange factor levels
  covPred.data$VitalRate <- factor(covPred.data$VitalRate, levels = c("Natural mortality", "Pregnancy rate", "Fetus number", "Immigration rate"))
  
  #----------------------------#
  # Plot covariate predictions #
  #----------------------------#
  
  ## Rodent effects on pregnancy rate, fetus number, and immigration rate
  p.VRs <- covPred.data %>%
    dplyr::mutate(VitalRate = factor(dplyr::case_when(VitalRate == "Natural mortality" ~ paste0("Natural mortality (age ", AgeClass, ")"),
                                                      VitalRate == "Pregnancy rate" ~ paste0("Pregnancy rate (age ", AgeClass, ")"),
                                                      VitalRate == "Fetus number" ~ paste0("Fetus number (age ", AgeClass, ")"),
                                                      TRUE ~ VitalRate), 
                                     levels = c(paste0("Natural mortality (age ", AgeClass, ")"), 
                                                paste0("Pregnancy rate (age ", AgeClass, ")"), 
                                                paste0("Fetus number (age ", AgeClass, ")"), 
                                                "Immigration rate"))) %>%
    ggplot(aes(x = RodentAbundance)) + 
    geom_line(aes(y = median, color = VitalRate)) + 
    geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = VitalRate), alpha = 0.5) + 
    scale_color_manual(values = c("#089392FF", "#52BA88FF", "#BED68AFF", "#E79069FF")) + 
    scale_fill_manual(values = c("#089392FF", "#52BA88FF", "#BED68AFF", "#E79069FF")) + 
    ylab("Predicted value") + 
    facet_wrap(~VitalRate, nrow = 2, scales = "free_y") + 
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
  
  pdf("Plots/RedfoxIPM_rodentEff_VitalRates.pdf", width = 8, height = 6)
  print(p.VRs)
  dev.off()
  
  ## Return list of plots
  plotList <- c("Plots/RedfoxIPM_rodentEff_VitalRates.pdf")
  
  return(plotList)
} 
