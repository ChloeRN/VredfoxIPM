
wrangleData_gen <- function(datapath, minYear, onlyFemales = FALSE, poolYears){
  
  ## Read in data (output from GeneClass analyses)
  gen.data <- read.delim(datapath)
  
  ## Optional: keep only females
  if(onlyFemales){
    gen.data <- subset(gen.data, v_sex == "female")
  }
  
  ## Transform year variables
  gen.data <- gen.data %>% 
    dplyr::mutate(yearIdx_shot = as.numeric(stringr::str_sub(t_hunting_year, 1, 4)) - minYear + 1,
                  yearIdx_birth = v_year_birth - minYear + 1,
                  yearIdx_birth_pre = ifelse(v_year_birth < minYear, v_year_birth - min(v_year_birth) + 1, NA))
}
