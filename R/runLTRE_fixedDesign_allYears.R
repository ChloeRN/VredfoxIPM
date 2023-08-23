runLTRE_fixedDesign_allYears <- function(paramSamples, Amax, Tmax, HazardRates = FALSE, PopStructure = TRUE){
  
  ## List all year pairs to run analyses for
  t_pairs <- cbind(1:(Tmax-2), 1:(Tmax-2) + 1)
  
  ## Set up results lists and dataframes
  contList <- list(cont = list(), other = list())
  contData <- data.frame()
  contData_summary <- data.frame()
  
  ## Run fixed design analyses for each pair of years
  for(t in 1:nrow(t_pairs)){
    LTRE_run <- runLTRE_fixedDesign(paramSamples = paramSamples, 
                                    t_pair = t_pairs[t,], 
                                    Amax = Amax,
                                    HazardRates = HazardRates, PopStructure = PopStructure)
    

    contList$cont[[t]] <- LTRE_run$contList$cont
    names(contList$cont)[t] <- paste0("t", t_pairs[t, 1], "t", t_pairs[t, 2])
    
    contList$other[[t]] <- LTRE_run$contList$other
    names(contList$other)[t] <- paste0("t", t_pairs[t, 1], "t", t_pairs[t, 2])
    
    contData_temp <- LTRE_run$contData
    contData_temp$t1 <- t_pairs[t, 1]
    contData_temp$t2 <- t_pairs[t, 2]
    contData <- rbind(contData, contData_temp)
    
    contData_summary_temp <- LTRE_run$contData_summary
    contData_summary_temp$t1 <- t_pairs[t, 1]
    contData_summary_temp$t2 <- t_pairs[t, 2]
    contData_summary <- rbind(contData_summary, contData_summary_temp)
    
  }
  
  ## Arrange all outputs in a new list
  results_allYears <- list(contList = contList, 
                           contData = contData,
                           contData_summary = contData_summary)
  
  
  ## Save and return
  save(results_allYears, file = "RedFoxIPM_LTREresults_fixedDesign.rds")
  return(results_allYears)
}
  
  