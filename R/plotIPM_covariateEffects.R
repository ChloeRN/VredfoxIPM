#' Plot vital rates as functions of fitted covariates
#'
#' @param MCMC.samples an mcmc.list containing the output of the fitted IPM. 
#' @param rCov.idx logical. If TRUE, stops function execution as this function 
#' has been written for the case of continuous covariates (rCov.idx = FALSE).
#' @param rodentMIN numeric. Minimum value of z-standardized rodent covariate to 
#' use for predictions. 
#' @param rodentMAX numeric. Maximum value of z-standardized rodent covariate to 
#' use for predictions. 
#' @param mHdevMIN numeric. Minimum value of log deviation of harvest mortality to 
#' use for predictions of compensation. 
#' @param mHdevMAX numeric. Maximum value of log deviation of harvest mortality to 
#' use for predictions of compensation. 
#' @param densityMIN numeric. Minimum value of centered log local population size to 
#' use for predictions. 
#' @param densityMAX numeric. Maximum value of centered log local population size to 
#' use for predictions. 
#' @param normN numeric. A priori set average of local population size (used for centering).
#' @param AgeClass integer. Age class for which to make predictions. Implemented
#' as "years of age" (range 0-4), not as model age index (range 1-5). 
#'
#' @return character vector of plot names. The plots themselves are saved
#' as pdf's in the subfolder "Plots".
#' @export
#'
#' @examples

plotIPM_covariateEffects <- function(MCMC.samples, rCov.idx, 
                                     rodentMIN, rodentMAX, 
                                     mHdevMIN, mHdevMAX,
                                     densityMIN, densityMAX,
                                     normN, AgeClass){
  
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
  
  
  #---------------------------------------------------#
  # Make density and harvest compensation predictions #
  #---------------------------------------------------#
  
  ## Make vector of values for covariates for which to predict for
  mHdevCov <- seq(mHdevMIN, mHdevMAX, length.out = 100)
  densityCov <- seq(densityMIN, densityMAX, length.out = 100)
  
  ## Set up empty dataframes to store results
  compPred.data <- data.frame()
  densityPred.data <- data.frame()
  
  ## Extract effect sizes if available
  if("gamma.immR" %in% colnames(out.mat)){
    gamma.immR <- out.mat[, "gamma.immR"]
  }else{
    gamma.immR <- rep(0, nrow(out.mat))
  }
  
  if("gamma.mO" %in% colnames(out.mat)){
    gamma.mO <- out.mat[, "gamma.mO"]
  }else{
    gamma.mO <- rep(0, nrow(out.mat))
  }
  
  if("betaD.immR" %in% colnames(out.mat)){
    betaD.immR <- out.mat[, "betaD.immR"]
  }else{
    betaD.immR <- rep(0, nrow(out.mat))
  }
  
  if("betaD.mO" %in% colnames(out.mat)){
    betaD.mO <- out.mat[, "betaD.mO"]
  }else{
    betaD.mO <- rep(0, nrow(out.mat))
  }
  
  ## Make posterior predictions for effects
  for(x in 1:length(mHdevCov)){
    
    mO.pred.comp <- imm.pred.comp <- rep(NA, nrow(out.mat))
    mO.pred.dens <- imm.pred.dens <- rep(NA, nrow(out.mat))
    
    for(i in 1:nrow(out.mat)){
      
      # Immigration rate
      imm.pred.comp[i] <- exp(log(out.mat[i, "Mu.immR"]) + 
                           gamma.immR[i]*mHdevCov[x])
      imm.pred.dens[i] <- exp(log(out.mat[i, "Mu.immR"]) + 
                                betaD.immR[i]*densityCov[x])
      
      # Natural morality 
      mO.pred.comp[i] <- exp(log(out.mat[i, paste0("Mu.mO[", AgeClass + 1, "]")]) + 
                               gamma.mO[i]*mHdevCov[x])
      mO.pred.dens[i] <- exp(log(out.mat[i, paste0("Mu.mO[", AgeClass + 1, "]")]) + 
                               betaD.mO[i]*densityCov[x])
    }
    
    # Assemble in temporary data frame
    data.temp.comp <- rbind(quantile(imm.pred.comp, probs = c(0.025, 0.5, 0.975)),
                       quantile(mO.pred.comp, probs = c(0.025, 0.5, 0.975))) %>%
      as.data.frame() %>%
      dplyr::rename(lCI = 1, median = 2, uCI = 3) %>%
      dplyr::mutate(VitalRate = c("Immigration rate", "Natural mortality"),
                    mHlogDev = mHdevCov[x])
    
    data.temp.dens <- rbind(quantile(imm.pred.dens, probs = c(0.025, 0.5, 0.975)),
                            quantile(mO.pred.dens, probs = c(0.025, 0.5, 0.975))) %>%
      as.data.frame() %>%
      dplyr::rename(lCI = 1, median = 2, uCI = 3) %>%
      dplyr::mutate(VitalRate = c("Immigration rate", "Natural mortality"),
                    PopDensity = densityCov[x])
    
    # Combine in main data frame
    compPred.data <- rbind(compPred.data, data.temp.comp)
    densityPred.data <- rbind(densityPred.data, data.temp.dens)
  }
  
  ## Re-arrange factor levels
  compPred.data$VitalRate <- factor(compPred.data$VitalRate, levels = c("Natural mortality", "Immigration rate"))
  densityPred.data$VitalRate <- factor(densityPred.data$VitalRate, levels = c("Natural mortality", "Immigration rate"))
  
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
  
  ## Density and compensatory effects on natural mortality and immigration
  p.comp <- compPred.data %>%
    dplyr::mutate(VitalRate = factor(dplyr::case_when(VitalRate == "Natural mortality" ~ paste0("Natural mortality (age ", AgeClass, ")"),
                                                      VitalRate == "Pregnancy rate" ~ paste0("Pregnancy rate (age ", AgeClass, ")"),
                                                      TRUE ~ VitalRate), 
                                     levels = c(paste0("Natural mortality (age ", AgeClass, ")"), 
                                                "Immigration rate"))) %>%
    ggplot(aes(x = mHlogDev)) + 
    geom_line(aes(y = median, color = VitalRate)) + 
    geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = VitalRate), alpha = 0.5) + 
    scale_color_manual(values = c("#089392FF", "#E79069FF")) + 
    scale_fill_manual(values = c("#089392FF", "#E79069FF")) + 
    ylab("Predicted value") + 
    xlab("Log deviation of mH[t]") + 
    ggtitle("a) Compensation") + 
    facet_wrap(~VitalRate, nrow = 1, scales = "free_y") + 
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
  
  p.dens <- densityPred.data %>%
    dplyr::mutate(VitalRate = factor(dplyr::case_when(VitalRate == "Natural mortality" ~ paste0("Natural mortality (age ", AgeClass, ")"),
                                                      VitalRate == "Pregnancy rate" ~ paste0("Pregnancy rate (age ", AgeClass, ")"),
                                                      TRUE ~ VitalRate), 
                                     levels = c(paste0("Natural mortality (age ", AgeClass, ")"), 
                                                "Immigration rate")),
                  PopDensity_t = exp(PopDensity + log(normN))) %>%
    ggplot(aes(x = PopDensity_t)) + 
    geom_line(aes(y = median, color = VitalRate)) + 
    geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = VitalRate), alpha = 0.5) + 
    scale_color_manual(values = c("#089392FF", "#E79069FF")) + 
    scale_fill_manual(values = c("#089392FF", "#E79069FF")) + 
    ylab("Predicted value") + 
    xlab("Local population size") + 
    ggtitle("b) Density dependence") + 
    facet_wrap(~VitalRate, nrow = 1, scales = "free_y") + 
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
  
  
  pdf("Plots/RedfoxIPM_DD&Compensation_VitalRates.pdf", width = 8, height = 6)
  print(p.comp / p.dens)
  dev.off()
  
  p.comp2 <- compPred.data %>%
    dplyr::filter(VitalRate == "Natural mortality") %>%
    dplyr::mutate(VitalRate = factor(dplyr::case_when(VitalRate == "Natural mortality" ~ paste0("Natural mortality (age ", AgeClass, ")"),
                                                      VitalRate == "Pregnancy rate" ~ paste0("Pregnancy rate (age ", AgeClass, ")"),
                                                      TRUE ~ VitalRate), 
                                     levels = c(paste0("Natural mortality (age ", AgeClass, ")"), 
                                                "Immigration rate"))) %>%
    ggplot(aes(x = mHlogDev)) + 
    geom_line(aes(y = median, color = VitalRate)) + 
    geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = VitalRate), alpha = 0.5) + 
    scale_color_manual(values = c("#089392FF")) + 
    scale_fill_manual(values = c("#089392FF")) + 
    ylab("Natural mortality") + 
    xlab("Log deviation of mH[t]") + 
    ggtitle("a) Compensatory mortality") + 
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
  
  p.dens2.mO <- densityPred.data %>%
    dplyr::filter(VitalRate == "Natural mortality") %>%
    dplyr::mutate(VitalRate = factor(dplyr::case_when(VitalRate == "Natural mortality" ~ paste0("Natural mortality (age ", AgeClass, ")"),
                                                      VitalRate == "Pregnancy rate" ~ paste0("Pregnancy rate (age ", AgeClass, ")"),
                                                      TRUE ~ VitalRate), 
                                     levels = c(paste0("Natural mortality (age ", AgeClass, ")"), 
                                                "Immigration rate")),
                  PopDensity_t = exp(PopDensity + log(normN)))%>%
    ggplot(aes(x = PopDensity_t)) + 
    geom_line(aes(y = median, color = VitalRate)) + 
    geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = VitalRate), alpha = 0.5) + 
    scale_color_manual(values = c("#089392FF")) + 
    scale_fill_manual(values = c("#089392FF")) + 
    ylab("Natural mortality") + 
    xlab("Local population size") + 
    ggtitle("b) Density-dependent mortality") + 
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
  
  p.dens2.immR <- densityPred.data %>%
    dplyr::filter(VitalRate == "Immigration rate") %>%
    dplyr::mutate(VitalRate = factor(dplyr::case_when(VitalRate == "Natural mortality" ~ paste0("Natural mortality (age ", AgeClass, ")"),
                                                      VitalRate == "Pregnancy rate" ~ paste0("Pregnancy rate (age ", AgeClass, ")"),
                                                      TRUE ~ VitalRate), 
                                     levels = c(paste0("Natural mortality (age ", AgeClass, ")"), 
                                                "Immigration rate")),
                  PopDensity_t = exp(PopDensity + log(normN))) %>%
    ggplot(aes(x = PopDensity_t)) + 
    geom_line(aes(y = median, color = VitalRate)) + 
    geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = VitalRate), alpha = 0.5) + 
    scale_color_manual(values = c("#E79069FF")) + 
    scale_fill_manual(values = c("#E79069FF")) + 
    ylab("Immigration rate") + 
    xlab("Local population size") + 
    ggtitle("c) Density-dependent immigration") + 
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
  
  
  pdf("Plots/RedfoxIPM_DD&Compensation_VitalRates_reduced.pdf", width = 4, height = 9)
  print(p.comp2 / p.dens2.mO / p.dens2.immR)
  dev.off()
  
  ## Return list of plots
  plotList <- c("Plots/RedfoxIPM_rodentEff_VitalRates.pdf",
                "Plots/RedfoxIPM_DD&Compensation_VitalRates.pdf",
                "Plots/RedfoxIPM_DD&Compensation_VitalRates_reduced.pdf")
  
  return(plotList)
} 
