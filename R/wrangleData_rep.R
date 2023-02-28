#' Prepare reproduction data
#'
#' This function relies on a helper function adjustAge_repData. 
#' 
#' @param datafiles vector of two character strings. Path/file names for raw 
#' data files containing counts of placental scars/embryos (index = 1) and 
#' presence/absence of placental scars/embryos/pregnancy signs (index = 2).
#' @param Amax integer. Number of age classes to consider in analyses.
#' @param minYear integer. First year to consider in analyses.
#'
#' @return a list containing two formatted dataframes 'P1' (counts) and 'P2'
#' (presences/absences).
#' @export
#'
#' @examples

wrangleData_rep <- function(datafiles, Amax, minYear){
  
  ## Check "datafiles" has two entries
  if(length(datafiles) < 2){
    stop("Not enough file paths/names provided.")
  }
  if(length(datafiles) > 2){
    stop("Too many file paths/names provided.")
  }
  
  ## Load embryo / placental scar count data
  P1.data <- read.table(datafiles[1], header = T)

  ## Load placental scar / pregnancy data (presence/absence)
  P2.data <- read.table(datafiles[2], header = T)

  ## Adjust age indices
  P1.data <- adjustAge_repData(data = P1.data, ageCol = 2, Amax = Amax)
  P2.data <- adjustAge_repData(data = P2.data, ageCol = 2, Amax = Amax)
  
  ## Remove data for age 0 (index 1) individuals
  age0.rep_P1 <- which(P1.data$age_adj == 1)
  if(length(age0.rep_P1) > 0){
    P1.data <- P1.data[-age0.rep_P1, ]
    warning(paste0(length(age0.rep_P1), " instance(s) of age 0 individuals with placental scars removed from data (P1)"))
  }

  age0.rep_P2 <- subset(P2.data, age_adj == 1 & P2 == 1)
  P2.data <- P2.data[-which(P2.data$age_adj == 1), ]
  if(nrow(age0.rep_P2) > 0){
    warning(paste0(nrow(age0.rep_P2), " instance(s) of age 0 individuals with placental scars removed from data (P2)"))
  }

  ## Add year index
  P1.data$RepYearIndex <- P1.data$repryear - minYear + 1
  P2.data$RepYearIndex <- P2.data$repryear - minYear + 1

  ## List and return data
  return(list(P1 = P1.data, P2 = P2.data))
  
}