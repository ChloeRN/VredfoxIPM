#' Prepare harvest count data
#'
#' @param obsH.datafile a vector containing integer counts of individuals harvested
#' each year. The list of years must be consecutive, even if there are years with
#' 0 or NA counts. 
#' @param sAaH.data a list containing an Age-at-Harvest matrix (C) and a vector of
#' yearly proportions of individuals aged/included in gge-at-Harvest data (pData).
#' @param minYear integer. First year to consider in analyses.
#'
#' @return a list containing a vector of annual harvest counts (obsH) and an
#' index at which to start analysis of counts as opposed to AaH data (Tmin_obsH).
#' @export
#'
#' @examples

wrangleData_obsH <- function(obsH.datafile, AaH.data, minYear){
  
  ## Extract years
  obsH.years <- as.numeric(names(obsH.datafile))
  year_range <- minYear:max(obsH.years)
  
  ## Check for gaps
  if(!all(year_range %in% obsH.years)){
    stop("There seem to be years missing in the data series. Double-check input data.")

  }
  
  ## Set value in 2004 to NA (no summer harvest recorded)
  obsH.datafile[1] <- NA
  
  ## Determine starting year index for analysis of count data
  maxYear_AaH <- max(as.numeric(colnames(AaH.data$C)))
  Tmin_obsH <- (maxYear_AaH + 1) - minYear + 1
  
  ## Return data
  return(list(obsH = as.vector(obsH.datafile),
              Tmin_obsH = Tmin_obsH))
}
