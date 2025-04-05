#' Truncate a value on both sides if it exceeds boundaries
#'
#' @param cov_pred numeric. Value to be truncated.
#' @param minValue numeric. Minimum value allowed. If cov_pred < minValue,
#' minValue is returned.
#' @param maxValue numeric. Maximum value allowed. If cov_pred > maxValue,
#' maxValue is returned.
#'
#' @return either cov_pred, minValue, or maxValue depending on whether 
#' minValue < cov_pred < maxValue.
#' @export
#'
#' @examples
#' 
truncateCovPrediction <- nimble::nimbleFunction(
  run = function(cov_pred = double(), 
                 minValue = double(), maxValue = double()){
    
    # Truncate extreme values
    if(cov_pred < minValue){
      cov_pred <- minValue
    }
    
    if(cov_pred > maxValue){
      cov_pred <- maxValue
    }
    
    return(cov_pred)
    
    returnType(double())
  })
