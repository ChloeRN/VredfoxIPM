runLTRE_randomDesign <- function(paramSamples, sensitivities, Amax, Tmax, HazardRates = FALSE, PopStructure = TRUE){

  
  #-------#
  # SETUP #
  #-------#
  
  ## Extract number of samples
  nosamples <- length(paramSamples$t_mean$lambda)

  ## List parameters to ignore
  dropParams <- c("N_tot", "B", "B_tot", "L", "L_tot", "R", "R_tot")
  
  if(HazardRates){
    dropParams <- c(dropParams, "S", "S0")
  }else{
    dropParams <- c(dropParams, "mH", "mO", "m0")
  }
  
  if(PopStructure){
    dropParams <- c(dropParams, "N")
  }else{
    dropParams <- c(dropParams, "n")
  }
  
  dropIdx <- which(names(paramSamples$t) %in% dropParams)
  
  paramList <- paramSamples$t[-dropIdx]
  sensList <- sensitivities$sensitivity$samples[-dropIdx]
  
  ## Set up list of arrays for storing calculated LTRE contributions
  contList <- list()
  
  for(x in 1:length(paramList)){
    
    if(names(paramList[x]) == "lambda"){
      next
    }
    
    if(length(dim(paramList[[x]])) == 3){
      
      tempList <- list()
      for(a in 1:Amax){
        tempList <- c(tempList, list(as.numeric(rep(NA, nosamples))))
      }
      names(tempList) <- paste0(names(paramList)[x], "_", 1:Amax)
      
    }else{
      
      tempList <- list(as.numeric(rep(NA, nosamples)))
      names(tempList) <- names(paramList)[x]
    }
    
    contList <- c(contList, tempList)
  }
  
  contCount <- length(contList)
  
  contList$est_var <- matrix(as.numeric(NA), nrow = contCount, ncol = nosamples)
  contList$est_covar <- matrix(as.numeric(NA), nrow = contCount, ncol = nosamples)
  
  
  #---------------------------------------------------------#
  # CALCULATE LTRE CONTRIBUTIONS PER SAMPLE - RANDOM DESIGN #
  #---------------------------------------------------------#
  
  for(i in 1:nosamples){
    
    ## Make lists of vital rates/population structure and of sensitivities
    dp_stoch_list <- list()
    sens_list <- list()
    for(x in 1:length(paramList)){
      
      # Skip lambda
      if(names(paramList)[x] == "lambda"){
        next
      }
      
      # Set time interval based on parameter
      if(names(paramList)[x] %in% c("Psi", "rho", "S0", "m0", "immR")){
        tInt <- 2:Tmax
      }else{
        tInt <- 1:(Tmax-1)
      }
      
      # Expand age class if required and list relevant parameter estimates
      if(length(dim(paramList[[x]])) == 3){
        
        tempList <- list()
        tempListS <- list()
        
        for(a in 1:Amax){
          tempList <- c(tempList, list(paramList[[x]][i, a, tInt]))
          tempListS <- c(tempListS, list(sensList[[x]][i, a]))
        }
        names(tempList) <- paste0(names(paramList)[x], "_", 1:Amax)
        names(tempListS) <- paste0(names(paramList)[x], "_", 1:Amax)
      }else{
        
        tempList <- list(paramList[[x]][i, tInt])
        tempListS <- list(sensList[[x]][i])
        names(tempList) <- names(paramList)[x]
        names(tempListS) <- names(paramList)[x]
      }
      
      dp_stoch_list <- c(dp_stoch_list, tempList)
      sens_list <- c(sens_list, tempListS)
    }
    
    ## Convert parameter list to matrix
    dp_stoch <- as.matrix(dplyr::bind_rows(dp_stoch_list, .id = "column_label"))
    
    ## Derive process variances and covariances
    dp_varcov <- var(dp_stoch)
    
    ## Save total estimated (co)variance per parameter
    contList$est_var[,i] <- diag(dp_varcov)
    contList$est_covar[,i] <- rowSums(dp_varcov, na.rm = T)
    
    ## Convert sensitivity list to vector
    sensvec <- do.call(c, sens_list)
    
    ## Calculate demographic contributions
    # NOTE: Here we multiply sensitivities and (co)variances
    
    cont.mat <- matrix(NA, nrow = length(sensvec), ncol = length(sensvec))
    for(k in 1:length(sensvec)){
      for(m in 1:length(sensvec)){
        cont.mat[k, m] <- dp_varcov[k, m]*sensvec[k]*sensvec[m]
      }
    }
    
    ## Summarise contributions (sum of variances and covariances)
    cont <- rowSums(cont.mat)
    names(cont) <- names(sensvec)
    
    ## Insert contributions into storage list
    for(x in 1:length(cont)){
      contList[[x]][i] <- cont[x]
    }
    
  }
  
  ## Restructure results list
  contList <- list(cont = contList[1:contCount],
                   other = list(est_var = contList$est_var,
                                est_covar = contList$est_covar)
                   )
  
  #-------------------#
  # SUMMARIZE RESULTS #
  #-------------------#
  
  ## Check sum of contributions against variance in lambda
  contList$other$total.contSum <- rowSums(dplyr::bind_rows(contList$cont))
  quantile(contList$other$total.contSum, probs = c(0.05, 0.5, 0.995))
  
  contList$other$tempvar_lambda <- matrixStats::rowVars(paramList$lambda[,1:(Tmax-1)])
  quantile(contList$other$tempvar_lambda, probs = c(0.05, 0.5, 0.995))
  
  
  ## Calculate summed contributions for age-specific parameters
  sumParams <- names(paramList)[which(!(names(paramList) %in% c("S0", "m0", "immR", "lambda")))] 
  
  for(x in 1:length(sumParams)){
    subList <- contList$cont[which(grepl(paste0(sumParams[x], "_"), names(contList$cont)))]
    
    contList$cont$newSum <- rowSums(dplyr::bind_rows(subList))
    names(contList$cont)[which(names(contList$cont) == "newSum")] <- paste0(sumParams[x], "_sum")
  }
  
  
  ## Arrange results as dataframe
  contData <- melt(dplyr::bind_rows(contList$cont, .id = "column_label"))
  contData <- cbind(contData, stringr::str_split_fixed(contData$variable, "_", 2))
  colnames(contData) <- c("Variable", "Contribution", "Parameter", "AgeClass")
  
  contData$AgeClass[which(contData$AgeClass == "")] <- NA
  
  
  ## Make posterior summaries
  contData_summary <- contData %>%
    dplyr::group_by(Variable) %>%
    dplyr::summarise(lCI = quantile(Contribution, 0.025),
                     median = median(Contribution),
                     uCI = quantile(Contribution, 0.975))
  
  ## Collect and return results
  
  results <- list(contList = contList,
                  contData = contData,
                  contData_summary = contData_summary)
  
  saveRDS(results, file = "RedFoxIPM_LTREresults_randomDesign.rds")
  
  return(results)
  
  
  #-----------------------------------#
  # OPTIONAL: PLOT LTRE CONTRIBUTIONS #
  #-----------------------------------#
  
  # TBA
  
  
}
