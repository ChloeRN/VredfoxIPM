#' Prepare Age-at-Harvest data
#'
#' @param AaH.datafile Age at harvest data file from carcass.data
#' @param Amax integer. Number of age classes to consider in analyses.
#'
#' @return a list containing an Age-at-Harvest matrix (C) and a vector of
#' yearly proportions of individuals aged/included in Age-at-Harvest data (pData).
#' @export
#'
#' @examples

wrangleData_AaH <- function(AaH.datafile, Amax){

  
  ## Extract Age-at-Harvest matrices (C)
  C <- t(AaH.datafile [, paste0("age", (1:Amax)-1)])
  colnames(C) <- AaH.datafile$year
  
  ## Extract proportions aged (pData)
  pData <- as.numeric(AaH.datafile[, "pAged"]) * as.numeric(AaH.datafile[, "pCarcass"])

  if("pAged_winter" %in% colnames(AaH.datafile)){
    pData_winter <- as.numeric(AaH.datafile[, "pAged_winter"]) * as.numeric(AaH.datafile[, "pCarcass_winter"])
    return(list(C = C, pData = pData, pData_winter = pData_winter))
  }else{
    return(list(C = C, pData = pData))
  }
}
