#' Make harvest effort covariate
#' 
#' We have no explicit data on harvest effort, but are here using the number of
#' successful hunters as a proxy.
#' Per now, we are writing the raw data (as provided by Dorothee in an email on
#' 27/10/20) directly in the function. 
#' It remains to clarify where this data is usually stored and to what degree it 
#' is possible to read it in directly from there.
#'
#' @return a dataframe containing original and scaled counts of successful 
#' hunters per year.
#' @export
#'
#' @examples

makeData_hunters <- function(){
  
  ## Write down number of successful hunters
  hunter.data <- data.frame(
    year = c(2004:2019),
    NHunters = c(NA, 41, 39, 44, 38, 38, 42, 51, 72, 35, 44, 47, 33, 45, 43, 44))
  
  # TODO: Double-check with Doro about NA for 2004 / year assignment
  
  ## Standardize and return
  hunter.data$NHunters_std <- (hunter.data$NHunters - mean(hunter.data$NHunters, na.rm = T))/sd(hunter.data$NHunters, na.rm = T)
  return(hunter.data)
  
  # TODO: Double-check where this data is stored and whether it's possible to 
  #       pull it automatically.
}
