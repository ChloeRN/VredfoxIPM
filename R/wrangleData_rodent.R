#' Prepare rodent abundance covariate data
#'
#' @param rodent.datafile Dataframe with rodent abundance data
#' @param minYear integer. First year to consider in analyses.
#'
#' @return a list containing rodent abundance data as continuous variables (cont),
#' and categorical variables with two (cat2). Data are provided for winter 
#' (fall + spring) in Varanger (.wintvar) and for fall for the larger area (.fallstor).
#' Continuous data are provided as total sums of individuals across all species
#' and as sums weighed by species (voles vs. lemmings, .stsp). 
#' Note that the time indices are shifted forward to represent that reproduction
#' is a function of past rodent abundance.
#' @export
#'
#' @examples

wrangleData_rodent <- function(rodent.reform.dat, minYear){
  
  ## Load reformatted data
  RodentData <- rodent.reform.dat
  
  ## Discard earlier years (if present)
  RodentData <- subset(RodentData, year >= minYear)
  
  ## List and return
  return(list(cont.wintvar          =   RodentData$st.tot.wintvar,     #winter varanger continuous, only standardised for seasons
              cont.wintvar.stsp     =   RodentData$st.lemvole.wintvar, #winter varanger continuous, standardised for seasons and species
              cat2.wintvar          =   RodentData$cat2.wintvar,       #winter varanger 2 categories
              
              cont.fallstor         =   RodentData$st.tot.fallstor,    #fall storskala continuous
              cont.fallstor.stsp    =   RodentData$st.lemvole.fallstor,#fall storskala continuous, standardised for species
              cat2.fallstor         =   RodentData$cat2.fallstor))     #fall storskala 2 factors
}
