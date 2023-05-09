#' Prepare rodent abundance covariate data
#'
#' @param rodent.datafile Dataframe with rodent abundance data
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

wrangleData_rodent <- function(rodent.reform.dat, minYear, adjust){
  
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
