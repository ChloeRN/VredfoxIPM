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
  
  ## Format rodent abundance data
  # Continuous winter and fall, standardised for seasons
  RodentAbundance_winter_seas.st <- RodentData$st.tot.winter
  RodentAbundance_fall_seas.st <- RodentData$st.tot.fall
  
  #Continuous winter and fall, standardised for seasons and species
  RodentAbundance_winter_seas.sp.st <- RodentData$st.lemvole.winter
  RodentAbundance_fall_seas.sp.st <- RodentData$st.lemvole.fall
  
  # 2-level categorical, winter and fall
  RodentIndex2_winter <- RodentData$cat2.winter
  RodentIndex2_fall   <- RodentData$cat2.fall

  ## List and return
  return(list(cont.wint          =   RodentAbundance_winter_seas.st,
              cont.fall          =   RodentAbundance_fall_seas.st,
              cont.wint.stsp     =   RodentAbundance_winter_seas.sp.st,
              cont.fall.stsp     =   RodentAbundance_fall_seas.sp.st,
              cat2.wint          =   RodentIndex2_winter,
              cat2.fall          =   RodentIndex2_fall))
  
}
