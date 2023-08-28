#' Plots results from a random design transient LTRE
#'
#' @param LTRE_results a list of lists containing results of a random design
#' transient LTRE.
#' @param Amax integer. Number of age classes. 
#' @param HazardRates logical. If TRUE (default), plots results of an LTRE with 
#' mortality hazard rates, else (FALSE) of an LTRE with survival probabilities. 
#' @param PopStructure logical. If TRUE (default), plots results of an LTRE with 
#' population proportions (n), else (FALSE) of an LTRE with age-specific population numbers (N).
#'
#' @return a character vector of plot names. The plots themselves are saved
#' as pdf's in the subfolder "Plots".
#' @export
#'
#' @examples

plotLTRE_randomDesign <- function(LTRE_results, Amax, HazardRates = TRUE, PopStructure = TRUE){
  
  #-------------#
  # Format data #
  #-------------#
  
  ## Select relevant part of data
  contData <- LTRE_results$contData
  
  ## Extract number of samples
  nosamples <- length(LTRE_results$contList$cont[[1]])
  
  ## Split off and format summed data
  contData_sum <- contData %>%
    dplyr::filter(Variable %in% c("S_sum", "mH_sum", "mO_sum",
                                  "Psi_sum", "rho_sum",
                                  "S0", "m0", "immR", 
                                  "n_sum", "N_sum")) %>%
    dplyr::mutate(type = dplyr::case_when(Variable == "S_sum" ~ "Annual survival",
                                          Variable == "mH_sum" ~ "Harvest mortality",
                                          Variable == "mO_sum" ~ "Natural mortality",
                                          Variable == "Psi_sum" ~ "Pregnancy rate",
                                          Variable == "rho_sum" ~ "Fetus number",
                                          Variable == "S0" ~ "Denning survival",
                                          Variable == "m0" ~ "Denning mortality",
                                          Variable == "immR" ~ "Immigration rate",
                                          Variable == "n_sum" ~ "Population structure",
                                          Variable == "N_sum" ~ "Population size/structure"))
  
  ## Make ordered list of present parameter types
  if(HazardRates){
    typeList <- c("Harvest mortality", "Natural mortality",
                  "Pregnancy rate", "Fetus number",
                  "Denning mortality", "Immigration rate")
  }else{
    typeList <- c("Annual survival", 
                  "Pregnancy rate", "Fetus number",
                  "Denning survival", "Immigration rate")
  }
  
  if(PopStructure){
    typeList <- c(typeList, "Population structure")
  }else{
    typeList <- c(typelist, "Population size/structure")
  }
  
  ## Re-order factor levels
  contData_sum$type <- factor(contData_sum$type, levels = typeList)
  
  
  #---------------------------------#
  # Plot sensitivities/elasticities #
  #---------------------------------#
  
  ## Plot colors
  plot.colors <- paletteer::paletteer_c("grDevices::Temps", length(unique(contData_sum$type)))

  ## Summed estimates for all parameters
  addline_format <- function(x,...){
    gsub('\\s','\n',x)
  }
  
  p.sum <- ggplot(contData_sum, aes(x = type, y = Contribution, group = type)) + 
    geom_violin(aes(fill = type, color = type), alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
    geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") + 
    ylab("Contribution") + 
    xlab('') + 
    scale_fill_manual(values = plot.colors) + 
    scale_color_manual(values = plot.colors) + 
    scale_x_discrete(labels = addline_format(typeList)) + 
    theme_bw() + 
    theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))
  
  p.sum
  
  if(!HazardRates){
    ## Survival panel
    p.S <- ggplot(subset(contData, Variable %in% paste0("S_", 1:Amax)), aes(x = Variable, y = Contribution, group = Variable)) + 
      geom_violin(fill = plot.colors[1], color = plot.colors[1], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
      geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") + 
      ylab("Contribution") + 
      xlab('') + 
      scale_x_discrete(labels = expression(S[1], S[2], S[3], S[4], S[5])) + 
      theme_bw() + 
      theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))
    p.S
  }else{
    ## Harvest mortality panel
    p.mH <- ggplot(subset(contData, Variable %in% paste0("mH_", 1:Amax)), aes(x = Variable, y = Contribution, group = Variable)) + 
      geom_violin(fill = plot.colors[1], color = plot.colors[1], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
      geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") + 
      ylab("Contribution") + 
      xlab('') + 
      scale_x_discrete(labels = expression(m[1]^H, m[2]^H, m[3]^H, m[4]^H, m[5]^H)) + 
      theme_bw() + 
      theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))
    p.mH
    
    ## Natural mortality panel
    p.mO <- ggplot(subset(contData, Variable %in% paste0("mO_", 1:Amax)), aes(x = Variable, y = Contribution, group = Variable)) + 
      geom_violin(fill = plot.colors[1], color = plot.colors[1], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
      geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") + 
      ylab("Contribution") + 
      xlab('') + 
      scale_x_discrete(labels = expression(m[1]^O, m[2]^O, m[3]^O, m[4]^O, m[5]^O)) + 
      theme_bw() + 
      theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))
    p.mO
  }
  
  
  ## Pregnancy rate panel
  p.Psi <- ggplot(subset(contData, Variable %in% paste0("Psi_", 2:Amax)), aes(x = Variable, y = Contribution, group = Variable)) + 
    geom_violin(fill = plot.colors[2], color = plot.colors[2], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
    geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") + 
    ylab("Contribution") + 
    xlab('') + 
    scale_x_discrete(labels = expression(Psi[2], Psi[3], Psi[4], Psi[5])) + 
    theme_bw() + 
    theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))
  p.Psi
  
  ## Fetus number panel
  p.rho <- ggplot(subset(contData, Variable %in% paste0("rho_", 2:Amax)), aes(x = Variable, y = Contribution, group = Variable)) + 
    geom_violin(fill = plot.colors[3], color = plot.colors[3], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
    geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") + 
    ylab("Contribution") + 
    xlab('') + 
    scale_x_discrete(labels = expression(rho[2], rho[3], rho[4], rho[5])) + 
    theme_bw() + 
    theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))
  p.rho
  
  ## Population structure panel
  if(PopStructure){
    p.n <- ggplot(subset(contData, Variable %in% paste0("n_", 1:Amax)), aes(x = Variable, y = Contribution, group = Variable)) + 
      geom_violin(fill = plot.colors[6], color = plot.colors[6], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
      geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") + 
      ylab("Contribution") + 
      xlab('') + 
      scale_x_discrete(labels = expression(n[1], n[2], n[3], n[4], n[5])) + 
      theme_bw() + 
      theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))
    p.n
  }else{
    p.n <- ggplot(subset(contData, Variable %in% paste0("N_", 1:Amax)), aes(x = Variable, y = Contribution, group = Variable)) + 
      geom_violin(fill = plot.colors[6], color = plot.colors[6], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
      geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") + 
      ylab("Contribution") + 
      xlab('') + 
      scale_x_discrete(labels = expression(N[1], N[2], N[3], N[4], N[5])) + 
      theme_bw() + 
      theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))
    p.n
  }

  
  ## Combine panels and save to pdf
  pdf(paste0("Plots/RedFoxIPM_", ifelse(HazardRates, "MHR", "SP"), "_sum.pdf"), width = 10, height = 6)
  print(
    p.sum
  )
  dev.off()
  
  if(HazardRates){
    
    pdf(paste0("Plots/RedFoxIPM_randomLTRE_", ifelse(HazardRates, "MHR", "SP"), "_age.pdf"), width = 7, height = 8)
    print(
      (p.mH + labs(tag = 'a)') | p.mO + labs(tag = 'b)')) / (p.Psi + labs(tag = 'c)')| p.rho + labs(tag = 'd)')) / (p.n  + labs(tag = 'e)') | plot_spacer())
    )
    dev.off()
  
  }else{
    pdf(paste0("Plots/RedFoxIPM_randomLTRE_", ifelse(HazardRates, "MHR", "SP"), "_age.pdf"), width = 7, height = 6)
    print(
      (p.S + labs(tag = 'a)')| p.Psi + labs(tag = 'b)')) / (p.rho  + labs(tag = 'c)')| p.n + labs(tag = 'd)'))
    )
    dev.off()
  }
  
  ## Return list of plots
  plotList <- c(paste0("Plots/RedFoxIPM_randomLTRE_", ifelse(HazardRates, "MHR", "SP"), "_sum.pdf"),
                paste0("Plots/RedFoxIPM_randomLTRE_", ifelse(HazardRates, "MHR", "SP"), "_age.pdf"))
  
  return(plotList)

}