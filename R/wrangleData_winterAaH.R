#' Prepare winter Age-at-Harvest data
#'
#' @param wAaH.datafile Winter age at harvest data file from wAaH.datafile <- carcassData$AaH.matrix
#' @param Amax integer. Number of age classes to consider in analyses.
#'
#' @return a list containing an Age-at-Harvest matrix (winterC) and a vector of
#' yearly proportions of individuals aged/included in Age-at-Harvest data (pData).
#' @export
#'
#' @examples

wrangleData_winterAaH <- function(wAaH.datafile, Amax){

  ## Extract Age-at-Harvest matrices (C)
  winterC <- t(wAaH.datafile [, paste0("age", (1:Amax)-1)])
  colnames(winterC) <- wAaH.datafile$year
  
  ## Extract proportions aged (pData)
  pData <- as.numeric(wAaH.datafile[, "pData"])

  ## List and return data
  return(list(winterC = winterC, pData = pData))
}
