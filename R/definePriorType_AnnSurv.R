#' Defines and stores information on prior for annual survival
#'
#' @param HoeningPrior logical. If TRUE, prior is based on Hoening model. If FALSE, prior is based on literature values for annual survival.
#' @param sPriorSource string. Has to be provided if HoeningPrior = FALSE and specifies which literature source is used to inform parameters for prior on annual survival. 
#' Currently supported options: Bristol (not hunted), NSweden (North Sweden, lightly hunted), metaAll (meta-analysis of all populations), metaSub (meta-analysis of not/lightly hunted populations).
#'
#' @return a list containing information on prior for annual survival.
#' @export
#'
#' @examples

definePriorType_AnnSurv <- function(HoeningPrior, sPriorSource = NULL){
  
  ## Check consistency of information provided
  if(!HoeningPrior & !exists("sPriorSource")){
    stop("When using literature priors for annual survival (HoeningPrior = FALSE), information on prior source (sPriorSource) has to be provided.")
  }
  
  if(!HoeningPrior & !(sPriorSource %in% c("Bristol", "NSweden", "metaAll", "metaSub"))){
    stop("Invalid prior source information provided. The following options are currently supported: Bristol, NSweden, metaAll, metaSub.")
  }
  
  ## List prior type information
  PriorType <- list(Parameter = ifelse(HoeningPrior, "natMort", "annSurv"),
                    Source = ifelse(HoeningPrior, NA, sPriorSource))
  
  return(PriorType)
}