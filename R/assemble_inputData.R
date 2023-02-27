#' Assemble demographic and environmental data for running the IPM
#'
#' @param Amax integer. Number of age classes to consider in analyses.
#' @param Tmax integer. The number of years to consider in analyses.
#' @param minYear integer. First year to consider in analyses.
#' @param maxPups integer. Upper prior bound for average litter size.
#' @param uLim.N integer. Upper prior bound for initial number of individuals per age class.
#' @param nLevels.rCov integer. Number of levels of categorical rodent abundance to use.
#' @param uLim.Imm integer. Upper prior bound for annual number of immigrants. 
#' @param wAaH.data a list containing an Age-at-Harvest matrix (winterC) and a vector of
#' yearly proportions of individuals aged/included in Age-at-Harvest data (pData).
#' @param rep.data a list containing formatted reproduction data in two data 
#' frames: 'P1' (counts) and 'P2' (presences/absences).
#' @param rodent.data a list containing rodent abundance data as a continuous variable (cont),
#' and categorical variable with two (cat2) and three (cat3) levels.
#' @param hunter.data a dataframe containing original and scaled counts of successful 
#' hunters per year.
#' @param surv.priors a list of lists containing parameters to define informative priors
#' for early survival, age-specific annual survival, and juvenile/adult natural
#' mortality hazard rate.
#' @param survPriorType a list containing information on prior for annual survival.
#' @param YearInfo a dataframe containing year indices as used in the model with 
#' corresponding reproduction years and winter harvest seasons. 
#' @param save logical. If TRUE, saves assembled data as an .rds file in the 
#' working directory. Default = FALSE. 
#'
#' @return a list containing all data necessary for running the IPM. 
#' @export
#'
#' @examples

assemble_inputData <- function(Amax, Tmax, minYear,
                               maxPups, uLim.N, uLim.Imm, nLevels.rCov = NA,
                               wAaH.data, rep.data, rodent.data, hunter.data, 
                               surv.priors, survPriorType, YearInfo,
                               save = FALSE){
  
  ## Select relevant years from observational data
  
  # Winter Age-at-Harvest data
  C <- wAaH.data$winterC[,which(colnames(wAaH.data$winterC) == minYear) + 1:Tmax - 1]
  pData <- wAaH.data$pData[which(colnames(wAaH.data$winterC) == minYear) + 1:Tmax - 1]
  
  # Reproduction data
  P1 <- subset(rep.data$P1, repryear %in% c(minYear + 0:Tmax))
  P2 <- subset(rep.data$P2, repryear %in% c(minYear + 0:Tmax))
  
  
  ## Select relevant categorical rodent covariate
  if(is.na(nLevels.rCov)){
    RodentIndex <- NA
  }else{
    RodentIndex <- ifelse(nLevels.rCov == 2, rodent.data$cat2, rodent.data$cat3)
  }
  
  ## List all relevant data (split into data and constants as used by NIMBLE)
  # Data
  nim.data <- list(
    C = C,
    pData = pData,
    
    P1 = P1$P1,
    
    P2 = P2$P2
  )
  
  # Constants
  nim.constants <- list(
    Amax = Amax,
    Tmax = Tmax,
    minYear = minYear,
    
    maxPups = maxPups,
    uLim.N = uLim.N,
    uLim.Imm = uLim.Imm,
    
    P1_age = P1$age_adj,
    P1_year = P1$RepYearIndex,
    X1 = nrow(P1$P1),
    
    P2_age = P2$age_adj,
    P2_year = P2$RepYearIndex,
    X2 = nrow(P2$P2),
    
    RodentAbundance = rodent.data$cont,
    RodentIndex = RodentIndex,
    nLevels.rCov = nLevels.rCov,
    
    NHunters = hunter.data$NHunters_std,
    
    YearInfo = YearInfo
  )
  
  ## Add relevant prior information
  nim.constants <- c(nim.constants, surv.priors$earlySurv)
  
  if(SurvPriorType$Parameter == "natMort"){
    nim.constants <- c(nim.constants, surv.priors$natMort)
  }else{
    sublistIdx <- which(names(surv.priors$annSurv) == SurvPriorType$Source)
    nim.constants <- c(nim.constants, surv.priors$annSurv[sublistIdx][[1]])
  }
  
  ## Combine data and constants in a list
  inputData <- list(nim.data = nim.data, 
                    nim.constants = nim.constants)
  
  ## Save (optional) and return data
  if(save){
    saveRDS("inputData_formatted.rds")
  }
  
  return(inputData)
  
}