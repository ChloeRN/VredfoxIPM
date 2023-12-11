
reformatData_reindeer <- function(minYear, Tmax){
  
  ## Load data files
  RDcarcass.data <- read.table("Data/Reindeer_carcasses.txt", header = TRUE)

  ## Define year range
  years <- minYear:(minYear + Tmax - 1)
  
  ## Assemble yearly data
  RDcarcass <- RDstomachs <- rep(NA, Tmax+1)
  for(t in 1:Tmax){
    year <- t + minYear - 1
    if(year %in% RDcarcass.data$winter){
      RDcarcass[t] <- RDcarcass.data$Varanger[which(RDcarcass.data$winter == year)]
    }
  }
  
  ## Arrange relevant data in a list
  reindeer.data <- list(RDcarcass = as.vector(scale(RDcarcass)))
  
  ## Return data
  return(reindeer.data)
  
}