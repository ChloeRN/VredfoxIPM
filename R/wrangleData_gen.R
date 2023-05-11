
wrangleData_gen <- function(datapath, minYear, onlyFemales = FALSE, GeneClass.approach, poolYrs.genData){
  
  ## Check that GeneClass.approach is correctly specified
  if(!(GeneClass.approach %in% c(1,2))){
    stop("Invalid GeneClass.approach. Only 1 and 2 are are accepted as inputs.")
  }
  
  ## Read in data (output from GeneClass analyses)
  gen.data <- read.delim(datapath)
  
  ## Optional: keep only females
  if(onlyFemales){
    gen.data <- subset(gen.data, v_sex == "female")
  }
  
  ## Determine variable (= analytical approach) to use
  if(GeneClass.approach == 1){
    gen.data$pvar <- gen.data$pvarV1
  }else{
    gen.data$pvar <- gen.data$pvarV2
  }
  
  ## Remove NAs and transform year variables
  gen.data <- gen.data %>% 
    dplyr::filter(!is.na(pvar)) %>%
    dplyr::mutate(yearIdx_shot = as.numeric(stringr::str_sub(t_hunting_year, 1, 4)) - minYear + 1,
                  yearIdx_birth = v_year_birth - minYear + 1,
                  yearIdx_birth_pre = ifelse(v_year_birth < minYear, v_year_birth - min(v_year_birth) + 1, NA))
  
  ## Extract and list relevant data
  
  if(poolYrs.genData){
    
    # Time-independent analysis (all data pooled)
    gen.data.out <- list(
      pImm = gen.data$pvar,
      Xgen = nrow(gen.data)
    )
    
  }else{
    
    # Time-dependent analysis (data separated into "in study period" vs. "before study period")
    inStudy <- which(gen.data$yearIdx_birth > 0)
    
    gen.data.out <- list(
      
      pImm_in = gen.data$pvar[inStudy],
      pImm_yrsB_in = gen.data$yearIdx_birth[inStudy],
      pImm_yrsH_in = gen.data$yearIdx_shot[inStudy],
      Xgen_in = length(gen.data$pvar[inStudy]),
      
      pImm_pre = gen.data$pvar[-inStudy],
      pImm_yrsB_pre = gen.data$yearIdx_birth_pre[-inStudy],
      pImm_yrsH_pre = gen.data$yearIdx_shot[-inStudy],
      Xgen_pre = length(gen.data$pvar[-inStudy])
    )
  }
  
  ## Return data
  return(gen.data.out)
}
