#' Prepare rodent abundance covariate data
#'
#' @param datafile character string. Path/file name for raw data file.
#' @param minYear integer. First year to consider in analyses.
#' @param adjust logical. Default = FALSE. If TRUE, adjust levels of 3-category
#' covariate in years 2018 and 2019 from 1 to 2 (in accordance with Dorothee's 
#' suggestion). May become redundant once we update data. 
#'
#' @return a list containing rodent abundance data as a continuous variable (cont),
#' and categorical variable with two (cat2) and three (cat3) levels.
#' @export
#'
#' @examples

wrangleData_rodent <- function(datafile, minYear, adjust = FALSE){
  
  ## Load raw data
  RodentData <- read.table(datafile, header = T)
  
  ## Discard earlier years (if present)
  RodentData <- subset(RodentData, year >= minYear)
  
  ## Adjust 3-level indices if necessary
  # TODO: Double-check neccessity of this with Doro
  if(adjust){
    adj_idx <- which(RodentData$year %in% c(2018, 2019))
    RodentData$cat3[adj_idx] <- 2
  }
  
  ## Format rodent abundance data
  # Continuous
  RodentAbundance <- (RodentData$tot-mean(RodentData$tot))/sd(RodentData$tot)
  
  # 2-level categorical
  RodentIndex2 <- RodentData$cat2
  
  # 3-level categorical
  RodentIndex3 <- RodentData$cat3

  ## List and return
  return(list(cont = RodentAbundance,
              cat2 = RodentIndex2,
              cat3 = RodentIndex3))
  
}