plotSensitivities <- function(sensitivities, Amax){
  
  
  for(i in 1:2){
    
    ## Select relevant part of data
    params <- sensitivities[[i]]$samples
    
    ## Drop pre-fix for generalising
    names(params) <- stringr::str_split_fixed(names(params), pattern = "_", n = 2)[,2]
    
    ## Extract number of samples
    nosamples <- dim(params[[1]])[1]
    
    #--------------------------#
    # Assemble summarised data #
    #--------------------------#
    
    sum.data <- data.frame(
      type = rep(c("Annual survival", "Harvest mortality", "Natural mortality",
                   "Pregnancy rate", "Fetus number",
                   "Denning survival", "Denning mortality",
                   "Immigration rate", "Population structure"), each = nosamples),
      estimate = c(rowSums(params$S), rowSums(params$mH), rowSums(params$mO),
                   rowSums(params$Psi), rowSums(params$rho),
                   params$S0, params$m0,
                   params$immR, rowSums(params$n))
      
    )
    
    ## Re-order factor levels
    sum.data$type <- factor(sum.data$type, levels = c("Annual survival", "Harvest mortality", "Natural mortality",
                                                      "Pregnancy rate", "Fetus number",
                                                      "Denning survival", "Denning mortality",
                                                      "Immigration rate", "Population structure"))
    
    
    #----------------------------#
    # Assemble age-specific data #
    #----------------------------#
    
    ## Bind all data into a data frame
    age.data <- data.frame(rlist::list.cbind(params))
    
    ## Change column names
    colnames(age.data) <- c(
      paste0("S_", 1:Amax), paste0("mH_", 1:Amax), paste0("mO_", 1:Amax),
      paste0("Psi_", 1:Amax), paste0("rho_", 1:Amax),
      "S0", "m0", 
      "immR", paste0("n_", 1:Amax), paste0("N_", 1:Amax)
    )
    
    ## Convert to longitudinal format
    age.data <- reshape2::melt(age.data)
    
    ## Remove unnecessary data
    age.data <- age.data %>%
      dplyr::filter(!(variable %in% c("Psi_1", "rho_1", paste0("N_", 1:Amax)))) %>%
      dplyr::rename(Parameter = variable, 
                    Estimate = value)
    
    
    #---------------------------------#
    # Plot sensitivities/elasticities #
    #---------------------------------#
    
    ## Plot colors
    plot.colors <- paletteer::paletteer_c("grDevices::Temps", length(unique(sum.data$type)) - 3)
    plot.colors2 <- c(rep(plot.colors[1], 3), plot.colors[2:3], rep(plot.colors[4], 2), plot.colors[5:length(plot.colors)])
    
    ## Summed estimates for all parameters
    addline_format <- function(x,...){
      gsub('\\s','\n',x)
    }
    
    p.sum <- ggplot(sum.data, aes(x = type, y = estimate, group = type)) + 
      geom_violin(aes(fill = type, color = type), alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
      geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") + 
      ylab(ifelse(i == 1, "Sensitivity", "Elasticity")) + 
      xlab('') + 
      scale_fill_manual(values = plot.colors2) + 
      scale_color_manual(values = plot.colors2) + 
      scale_x_discrete(labels = addline_format(c("Annual survival", "Harvest mortality", "Natural mortality",
                                                 "Pregnancy rate", "Fetus number",
                                                 "Denning survival", "Denning mortality",
                                                 "Immigration rate", "Population structure"))) + 
      theme_bw() + 
      theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))
    
    p.sum
    
    ## Survival panel
    p.S <- ggplot(subset(age.data, Parameter %in% paste0("S_", 1:Amax)), aes(x = Parameter, y = Estimate, group = Parameter)) + 
      geom_violin(fill = plot.colors[1], color = plot.colors[1], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
      ylab(ifelse(i == 1, "Sensitivity", "Elasticity")) + 
      xlab('') + 
      scale_x_discrete(labels = expression(S[1], S[2], S[3], S[4], S[5])) + 
      theme_bw() + 
      theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))
    p.S
    
    ## Harvest mortality panel
    p.mH <- ggplot(subset(age.data, Parameter %in% paste0("mH_", 1:Amax)), aes(x = Parameter, y = Estimate, group = Parameter)) + 
      geom_violin(fill = plot.colors[1], color = plot.colors[1], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
      ylab(ifelse(i == 1, "Sensitivity", "Elasticity")) + 
      xlab('') + 
      scale_x_discrete(labels = expression(m[1]^H, m[2]^H, m[3]^H, m[4]^H, m[5]^H)) + 
      theme_bw() + 
      theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))
    p.mH
    
    ## Natural mortality panel
    p.mO <- ggplot(subset(age.data, Parameter %in% paste0("mO_", 1:Amax)), aes(x = Parameter, y = Estimate, group = Parameter)) + 
      geom_violin(fill = plot.colors[1], color = plot.colors[1], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
      ylab(ifelse(i == 1, "Sensitivity", "Elasticity")) + 
      xlab('') + 
      scale_x_discrete(labels = expression(m[1]^O, m[2]^O, m[3]^O, m[4]^O, m[5]^O)) + 
      theme_bw() + 
      theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))
    p.mO
    
    ## Pregnancy rate panel
    p.Psi <- ggplot(subset(age.data, Parameter %in% paste0("Psi_", 1:Amax)), aes(x = Parameter, y = Estimate, group = Parameter)) + 
      geom_violin(fill = plot.colors[2], color = plot.colors[2], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
      ylab(ifelse(i == 1, "Sensitivity", "Elasticity")) + 
      xlab('') + 
      scale_x_discrete(labels = expression(Psi[2], Psi[3], Psi[4], Psi[5])) + 
      theme_bw() + 
      theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))
    p.Psi
    
    ## Fetus number panel
    p.rho <- ggplot(subset(age.data, Parameter %in% paste0("rho_", 1:Amax)), aes(x = Parameter, y = Estimate, group = Parameter)) + 
      geom_violin(fill = plot.colors[3], color = plot.colors[3], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
      ylab(ifelse(i == 1, "Sensitivity", "Elasticity")) + 
      xlab('') + 
      scale_x_discrete(labels = expression(rho[2], rho[3], rho[4], rho[5])) + 
      theme_bw() + 
      theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))
    p.rho
    
    ## Population structure panel
    p.n <- ggplot(subset(age.data, Parameter %in% paste0("n_", 1:Amax)), aes(x = Parameter, y = Estimate, group = Parameter)) + 
      geom_violin(fill = plot.colors[6], color = plot.colors[6], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
      ylab(ifelse(i == 1, "Sensitivity", "Elasticity")) + 
      xlab('') + 
      scale_x_discrete(labels = expression(n[1], n[2], n[3], n[4], n[5])) + 
      theme_bw() + 
      theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))
    p.n
    
    ## Combine panels and save to pdf
    pdf(paste0("Plots/RedFoxIPM_", ifelse(i == 1, "Sensitivities", "Elasticities"), "_sum.pdf"), width = 10, height = 6)
    print(
      p.sum
    )
    dev.off()
    
    pdf(paste0("Plots/RedFoxIPM_", ifelse(i == 1, "Sensitivities", "Elasticities"), "_age.pdf"), width = 7, height = 8)
    print(
      (p.mH + labs(tag = 'a)') | p.mO + labs(tag = 'b)')) / (p.S + labs(tag = 'c)')| p.Psi + labs(tag = 'd)')) / (p.rho  + labs(tag = 'e)')| p.n + labs(tag = 'f)'))
    )
    dev.off()
    
  }  
  
  ## Return list of plots
  plotList <- c(paste0("Plots/RedFoxIPM_", c("Sensitivities", "Elasticities"), "_sum.pdf"),
                paste0("Plots/RedFoxIPM_", c("Sensitivities", "Elasticities"), "_age.pdf"))
  
  return(plotList)
  
}