#' Assemble demographic and environmental data for running the IPM
#'
#' @param Amax integer. Number of age classes to consider in analyses.
#' @param Tmax integer. The number of years to consider in analyses.
#' @param minYear integer. First year to consider in analyses.
#' @param wAaH.data a list containing an Age-at-Harvest matrix (winterC) and a vector of
#' yearly proportions of individuals aged/included in Age-at-Harvest data (pData).
#' @param rep.data a list containing formatted reproduction data in two data 
#' frames: 'P1' (counts) and 'P2' (presences/absences).
#' @param rodent.data a list containing rodent abundance data as a continuous variable (cont),
#' and categorical variable with two (cat2) and three (cat3) levels.
#' @param hunter.data a dataframe containing original and scaled counts of successful 
#' hunters per year.
#' @param YearInfo a dataframe containing year indices as used in the model with 
#' corresponding reproduction years and winter harvest seasons. 
#'
#' @return a list containing all data necessary for running the IPM. 
#' @export
#'
#' @examples

assemble_inputData <- function(Amax, Tmax, minYear,
                               wAaH.data, rep.data, rodent.data, hunter.data, YearInfo){
  
  ## Select relevant years
  
  # Winter Age-at-Harvest data
  C <- wAaH.data$winterC[,which(colnames(wAaH.data$winterC) == minYear) + 1:Tmax - 1]
  pData <- wAaH.data$pData[which(colnames(wAaH.data$winterC) == minYear) + 1:Tmax - 1]
  
  # Reproduction data
  P1 <- subset(rep.data$P1, repryear %in% c(minYear + 0:Tmax))
  P2 <- subset(rep.data$P2, repryear %in% c(minYear + 0:Tmax))
  
  ## List all relevant data
  inputData <- list(
    Amax = Amax,
    Tmax = Tmax,
    minYear = minYear,
    
    C = C,
    pData = pData,
    
    P1 = P1$P1,
    P1_age = P1$age_adj,
    P1_year = P1$RepYearIndex,
    X1 = nrow(P1),
    
    P2 = P2$P2,
    P2_age = P2$age_adj,
    P2_year = P2$RepYearIndex,
    X2 = nrow(P2),
    
    RodentAbundance = rodent.data$cont,
    RodentIndex2 = rodent.data$cat2,
    RodentIndex3 = rodent.data$cat3,
    
    NHunters = hunter.data$NHunters_std,
    
    YearInfo = YearInfo
  )
  
  return(inputData)
  
}