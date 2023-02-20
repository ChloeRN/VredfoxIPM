#' Prepare winter Age-at-Harvest data
#'
#' @param datafile character string. Path/file name for raw data files.
#' @param Amax integer. Number of age classes to consider in analyses.
#'
#' @return a list containing an Age-at-Harvest matrix (winterC) and a vector of
#' yearly proportions of individuals aged/included in Age-at-Harvest data (pData).
#' @export
#'
#' @examples

wrangleData_winterAaH <- function(datafile, Amax){
  
  ## Load winter Age-at-Harvest data
  winterAaH <- read.table(datafile, header = T)

  ## Extract Age-at-Harvest matrices (C)
  winterC <- t(winterAaH[, paste0("age", (1:Amax)-1)])
  colnames(winterC) <- winterAaH$year
  
  ## Extract proportions aged (pData)
  pData <- as.numeric(winterAaH[, "pData"])

  ## List and return data
  return(list(winterC = winterC, pData = pData))
}