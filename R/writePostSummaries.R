#' Summarise posteriors and write to file
#'
#' @param MCMC.samples an mcmc.list containing the output of the fitted IPM. 
#'
#' @return character vector of file names. The files themselves are saved
#' in the working directory. 
#' @export 
#'
#' @examples

writePostSummaries <- function(MCMC.samples){
  
  ## Summarise posterior distributions
  MCMC.sum <- as.matrix(MCMC.samples) %>%
    reshape2::melt() %>%
    dplyr::rename(Parameter = Var2,
                  Value = value) %>%
    dplyr::group_by(Parameter) %>%
    dplyr::summarise(median = median(Value, na.rm = TRUE),
                     lCI = quantile(Value, probs = 0.025, na.rm = TRUE),
                     uCI = quantile(Value, probs = 0.975, na.rm = TRUE))
  
  ## Save to .rds and .csv
  saveRDS(MCMC.sum, file = "RedfoxIPM_PostSummaries.rds")
  write.csv(MCMC.sum, file = "RedfoxIPM_PostSummaries.csv", row.names = FALSE)
  
  ## Return filenames
  fileNames <- c("RedfoxIPM_PostSummaries.rds", "RedfoxIPM_PostSummaries.csv")
  return(fileNames)
  
}