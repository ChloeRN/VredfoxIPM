#' Set up perturbation vectors for Population Viability Analysis (PVA) scenarios
#'
#' @param Tmax integer. The number of years in the analysis.
#' @param Tmax_sim integer. The number of years to consider for simulations 
#' beyond the data collection period. 
#' @param pert.mH logical. Whether to apply a perturbation to harvest mortality
#' hazard rate. Default = FALSE. 
#' @param factor.mH numeric. Relative change to harvest mortality hazard rate to 
#' apply. 1 = no change (default). < 1 = decrease. > 1 = increase. 
#' @param pert.mO logical. Whether to apply a perturbation to natural mortality
#' hazard rate. Default = FALSE. 
#' @param factor.mO numeric. Relative change to natural mortality hazard rate to 
#' apply. 1 = no change (default). < 1 = decrease. > 1 = increase. 
#' @param pert.S0 logical. Whether to apply a perturbation to denning survival. 
#' Default = FALSE. 
#' @param factor.S0 numeric. Relative change to denning survival to apply.
#' 1 = no change (default). < 1 = decrease. > 1 = increase. 
#' @param pert.immR logical. Whether to apply a perturbation to immigration rate. 
#' Default = FALSE. 
#' @param factor.immR numeric. Relative change to immigration rate to apply.
#' 1 = no change (default). < 1 = decrease. > 1 = increase. 
#'
#' @return
#' @export
#'
#' @examples

setupPerturbVecs_PVA <- function(Tmax, Tmax_sim,
                                 pert.mH = FALSE, factor.mH = 1,
                                 pert.mO = FALSE, factor.mO = 1,
                                 pert.S0 = FALSE, factor.S0 = 1,
                                 pert.immR = FALSE, factor.immR = 1){
  
  ## Check that there are no invalid perturbation factors
  if(any(c(factor.mH, factor.mO, factor.S0, factor.immR) < 0)){
    stop("Invalid perturbation factor provided. Perturbation factors have to be
         numerical values >= 0.")
  }
  
  if(any(!is.numeric(c(factor.mH, factor.mO, factor.S0, factor.immR)))){
    stop("Invalid perturbation factor provided. Perturbation factors have to be
         numerical values >= 0.")
  }
  
  ## Set up basics perturbation vectors during study/data period
  pertFac.mH <- pertFac.mO <- rep(1, Tmax)
  pertFac.S0 <- pertFac.immR <- rep(1, Tmax+1)
  
  ## Add factors for perturbation period
  if(Tmax_sim > 0){
    pertFac.mH <- c(pertFac.mH, rep(factor.mH, Tmax_sim))
    pertFac.mO <- c(pertFac.mO, rep(factor.mO, Tmax_sim))
    pertFac.S0 <- c(pertFac.S0, rep(factor.S0, Tmax_sim))
    pertFac.immR <- c(pertFac.immR, rep(factor.immR, Tmax_sim))
  }
  
  ## List and return perturbation vectors
  pertVecs <- list(pertFac.mH = pertFac.mH, pertFac.mO = pertFac.mO, 
                   pertFac.S0 = pertFac.S0, pertFac.immR = pertFac.immR)
  return(pertVecs)
}