
predict_mO_HoeningMod <- function(datafile, mu.t.max, maxAge, nsim, plot){
  
  
  ## Load posterior samples from Hoening model
  postSamples <- read.table(datafile, header = T)
  
  ## Extract posterior samples for relevant parameters
  param.a <- postSamples[,'a_1']
  param.b <- postSamples[,'b_1']
  sigma.z <- postSamples[,'sigma.z_1']

  ## Simulate prior from maximum recorded age
  # Set the number of posterior samples
  nsamples <- length(param.a)

  # Prepare a data frame to store results
  data <- data.frame()
  
  # Loop over posterior samples and simulate mO
  for(i in 1:nsamples){
    
    log_mu.mO <- param.a[i] + param.b[i] * (log(maxAge) - log(mu.t.max))
    
    mO <- rlnorm(nsim, meanlog = log_mu.mO, sdlog = sigma.z[i])
    data <- rbind(data, data.frame(sampleNo = rep(i, nsim), 
                                   simNo = 1:nsim, 
                                   mO_est = mO))
  }
  
  
  ## Extract log-mean and log-sd for prior distribution
  log_mean <- mean(log(data$mO_est)) 
  log_sd <- sd(log(data$mO_est)) 
  
  ## Optional: Visualize simulated distributions
  if(plot){
    mO_sim <- rlnorm(nrow(data), meanlog = log_mean, sdlog = log_sd)
    
    mO_Comp <- data.frame(
      mO_est = c(data$mO_est, mO_sim),
      Source = rep(c('Posterior estimation', 'Simulation'), each = nrow(data))
    )
    
    ggplot(mO_Comp, aes(x = mO_est)) + 
      geom_density(aes(color = Source, fill = Source), alpha = 0.25) + 
      xlab('Natural mortality prior') + 
      xlim(0, 1) + 
      scale_color_manual(values = c("grey40", "#00C7A9")) +
      scale_fill_manual(values = c("grey40", "#00C7A9")) + 
      theme_bw() + theme(panel.grid = element_blank(), legend.position = 'top')
  }

  
  ## Save and return results
  mO.prior <- list(mnat.logmean = log_mean, mnat.logsd = log_sd)
  
  saveRDS(mO.prior, file = "mO_prior_Parameters.rds")
  
  return(mO.prior)
  
}
