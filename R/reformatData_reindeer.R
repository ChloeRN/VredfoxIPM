
reformatData_reindeer <- function(minYear, Tmax){
  
  ## Load data files
  RDcarcass.data <- read.table("Data/Reindeer_carcasses.txt", header = TRUE)
  RDstomachs.data <- read.table("Data/Reindeer_in_fox_stomachs.txt", header = TRUE)
  
  ## Define year range
  years <- minYear:(minYear + Tmax - 1)
  
  ## Assemble yearly data
  RDcarcass <- RDstomachs <- rep(NA, Tmax+1)
  for(t in 1:Tmax){
    year <- t + minYear - 1
    if(year %in% RDcarcass.data$winter){
      RDcarcass[t] <- RDcarcass.data$Varanger[which(RDcarcass.data$winter == year)]
    }
    if(year %in% RDstomachs.data$winter){
      RDstomachs[t] <- RDstomachs.data$prop_rein[which(RDstomachs.data$winter == year)]
    }
  }
  
  ## Arrange relevant data in a list
  reindeer.data <- list(RDcarcass = as.vector(scale(RDcarcass)),
                        RDstomachs = as.vector(scale(RDstomachs)))
  
  ## Return data
  return(reindeer.data)
  
}