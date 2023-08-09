
compareModels <- function(Amax, Tmax, minYear, post.filepaths, post.list, model.names, plotFolder){
  
  ## Check models are specified correctly
  if((missing(post.filepaths) & missing(post.list)) |
     (missing(post.filepaths) & missing(post.list))){
    stop("Models have to be specified either via file paths (post.filepaths) or 
         using object names (post.objects).")
  }
  
  ## Make plotting directory if it does not exist already
  if(!dir.exists(plotFolder)){
    dir.create(plotFolder)
  }
  
  ## Count number of models
  nModels <- length(model.names)
  
  ## Reformat posterior samples
  post.data <- data.frame()
  for(i in 1:nModels){
    
    # Extract samples for relevant model
    if(!missing(post.list)){
      samples <- post.list[[i]]
    }else{
      samples <- readRDS(post.filepaths[i])
    }
    
    # Change format and add to list
    model.data <- reshape2::melt(as.matrix(samples))
    colnames(model.data) <- c("Sample", "Parameter", "Value")
    model.data$Model <- model.names[i]
    post.data <- rbind(post.data, model.data)
  }
  
  ## Summarize posterior samples into median + 95% CI
  sum.data <- post.data %>%
    
    dplyr::group_by(Parameter, Model) %>%
    
    dplyr::summarise(median = median(Value, na.rm = TRUE),
                     lCI = quantile(Value, probs = 0.025, na.rm = TRUE),
                     uCI = quantile(Value, probs = 0.975, na.rm = TRUE),
                     .groups = "keep") %>%
    
    dplyr::mutate(Parameter = as.character(Parameter)) %>%
    
    dplyr::ungroup()
  
  ## Extract and add age and year information
  idx.data <- data.frame(cbind(unique(sum.data$Parameter), stringr::str_extract_all(unique(sum.data$Parameter), pattern = "\\d+", simplify = TRUE))) %>%
    dplyr::rename("Parameter" = "X1",
                  "Idx1" = "X2", 
                  "Idx2" = "X3") %>%
    
    dplyr::mutate(Idx1 = as.numeric(ifelse(Idx1 %in% c("", 0), NA, Idx1)),
                  Idx2 = as.numeric(ifelse(Idx2 %in% c("", 0), NA, Idx2)),
                  YearIdx = dplyr::case_when(!is.na(Idx2) ~ Idx2,
                                             !is.na(Idx1) & !grepl("Mu", Parameter) ~ Idx1),
                  AgeIdx = dplyr::case_when(!is.na(Idx2) ~ Idx1, 
                                            !is.na(Idx1) & grepl("Mu", Parameter) ~ Idx1),
                  Year = YearIdx + minYear - 1,
                  Age = ifelse(AgeIdx == Amax, paste0(Amax-1, "+"), AgeIdx-1),
                  ParamName = stringr::word(Parameter, 1, sep = "\\["))

    
  sum.data <- sum.data %>%
    dplyr::left_join(idx.data, by = "Parameter")
    
  ## Set parameter groups for plotting posterior density overlaps
  plot.params <- list(
    VRmeans = c(paste0("Mu.mH[", 1:Amax, "]"), 
                paste0("Mu.mO[", 1:Amax, "]"), 
                paste0("Mu.Psi[", 2:Amax, "]"), 
                paste0("Mu.rho[", 2:Amax, "]"),  
                "Mu.S0", "Mu.immR"),
    
    VReffects = c("sigma.mH", "sigma.mO", "sigma.Psi", "sigma.rho", "sigma.immR",
                  "betaHE.mH", 
                  "betaR.Psi", paste0("betaR.Psi[", 2:3, "]"),
                  "betaR.rho", paste0("betaR.rho[", 2:3, "]"),
                  "betaRd.mO", "betaR.mO", "betaRxRd.mO"),
    
    Imm = paste0("Imm[", 1:Tmax, "]"),
    
    Ntot = paste0("N.tot[", 1:Tmax, "]"),
    
    Btot = paste0("B.tot[", 1:Tmax, "]"),
    
    Rtot = paste0("R.tot[", 1:Tmax, "]")
  )
  
  ## Set parameters plotting time series of posterior summaries
  plotTS.params <- list(
    ParamNames = c("N.tot", "B.tot", "R.tot", "Imm",
                   "mO", "S", #"S0",
                   "mH", "Psi", "rho", "immR"),
    ParamLabels = c("Female population size", "# breeding females", "# female recruits", "# female immigrants",
                    "Natural mortality", "Survival", #"Early survival",
                    "Harvest mortality", "Pregnancy rate", "# fetuses/female", "Immigration rate")
  )

  ## Set plotting colors
  plot.cols <- paletteer::paletteer_c("grDevices::Temps", length(model.names))

  ## Plot posterior overlaps
  pdf(paste0(plotFolder, "/PosteriorDensities.pdf"), width = 9, height = 6)
  for(x in 1:length(plot.params)){
    
    print(
      ggplot(subset(post.data, Parameter %in% plot.params[[x]]), aes(x = Value, color = Model, fill = Model)) + 
        geom_density(alpha = 1/nModels) + 
        facet_wrap(~Parameter, scales = "free") + 
        #scale_fill_viridis_d() + scale_color_viridis_d() + 
        scale_fill_manual(values = plot.cols) + scale_color_manual(values = plot.cols) + 
        theme_bw() + theme(panel.grid = element_blank())
    )
  }
  dev.off()
  
  ## Plot posterior summary time series
  pdf(paste0(plotFolder, "/PosteriorSummaries_TimeSeries.pdf"), width = 8, height = 8)
  for(x in 1:length(plotTS.params$ParamNames)){
    
    #sum.data.sub <- subset(sum.data, ParamName == plotTS.params$ParamNames[x])
    
    print(
      ggplot(subset(sum.data, ParamName == plotTS.params$ParamNames[x]), aes(group = Model)) + 
        geom_line(aes(x = Year, y = median, color = Model)) + 
        geom_ribbon(aes(x = Year, ymin = lCI, ymax = uCI, fill = Model), alpha = 1/nModels) + 
        #scale_fill_viridis_d() + scale_color_viridis_d() + 
        scale_fill_manual(values = plot.cols) + scale_color_manual(values = plot.cols) + 
        facet_wrap(~ Age, ncol = 1, scales = "free_y") + 
        ggtitle(plotTS.params$ParamLabels[x]) +  
        theme_bw() + theme(panel.grid = element_blank())
    )
    
    if(plotTS.params$ParamNames[x] %in% c("Imm", "immR")){
      
      subset.years <- (min(sum.data$Year, na.rm = TRUE) + 1):(max(sum.data$Year, na.rm = TRUE) - 1)
      print(
        ggplot(subset(sum.data, ParamName == plotTS.params$ParamNames[x] & Year %in% subset.years), aes(group = Model)) + 
          geom_line(aes(x = Year, y = median, color = Model)) + 
          geom_ribbon(aes(x = Year, ymin = lCI, ymax = uCI, fill = Model), alpha = 1/nModels) + 
          #scale_fill_viridis_d() + scale_color_viridis_d() + 
          scale_fill_manual(values = plot.cols) + scale_color_manual(values = plot.cols) + 
          facet_wrap(~ Age, ncol = 1, scales = "free_y") + 
          ggtitle(paste0(plotTS.params$ParamLabels[x], " (without first/last year)")) +  
          theme_bw() + theme(panel.grid = element_blank())
      )
      
    }
  }
  dev.off()
  
}
