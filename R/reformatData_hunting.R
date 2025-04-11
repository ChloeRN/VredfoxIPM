#' Reformat the hunting data to make a combined total number foxes hunted for each winter hunting period, summer hunting period, and "yearly" hunting period 
#(from start of summer to end of winter hunt). (all foxes, male and female, in Varanger)
#'
#' @param summer_removal a vector of numerical months to be removed from winter age at harvest data: c(6,7,8,9)
#' @param hunting.dataset list containing the datasets downloaded from the COAT dataportal.
#'
#' @return a dataframe containing the reported numbers of hunted foxes in total and for both seasons
#'
#' @examples


reformatData_hunting <- function (summer_removal, hunting.dataset) {
  
  
  #========= HELPER FUNCTIONS ==============
  '%notin%' <- Negate('%in%')
  
  #========= LOAD DATA ==============
  
  #We have 2 types of data
  #1) foxes hunted by SNO, in 4 columns
  #2) foxes hunted by local hunters, in 16 columns
  
  sno_hunt <- hunting.dataset$dat_4_columns
  sno_hunt <- sno_hunt[sno_hunt$sn_region == "varanger",] #select varanger only
  local_hunt <- hunting.dataset$dat_16_columns
  
  #-----------------Some checks ----------    
  #---check if SNO hunts in summer (allf is from refotmatData_carcass.R)---
  
  #select from all foxes in carcass.examination, those shot in varanger
 # all_varanger <- allf[is.na(allf$sn_region)== F & allf$sn_region == "varanger"  & allf$v_hunter_id != "road_kill" & allf$v_hunter_id != "found_dead", ]
  
 # sno_summer <- all_varanger[all_varanger$v_hunter_id == "sno_varanger" & all_varanger$mnd %in% summer_removal,]
 # sno_winter <- all_varanger[all_varanger$v_hunter_id == "sno_varanger" & all_varanger$mnd %notin% summer_removal,]
  #they do, 12 summer, 718 winter
 
 
 #-----check if all SNO end up in carcass_examination----
 # Step 1: Count observations for each t_hunting_year where v_hunter_id == "sno_varanger"
# all_varanger_counts <- all_varanger %>%
 #  filter(v_hunter_id == "sno_varanger") %>%
 #  group_by(t_hunting_year) %>%
 #  summarise(observation_count = n(), .groups = "drop")
 # Step 2: Merge with sno_hunt dataset to compare with v_abundance
 #comparison <- sno_hunt %>%
#   filter(v_hunter_id == "sno_varanger") %>%
 #  left_join(all_varanger_counts, by = "t_hunting_year") %>%
 #  mutate(difference = v_abundance - observation_count)
# View(comparison)
  

  #-------------------------------------------------------------------------------------------------
  #It doesnt look good but lets continue anyway
 
 #========== PREPARE DATA ===========
  #Oke here we just assume that SNO hunts only in winter
  
  #we want for each t_hunting year, the total, summer, and winter number of foxes hunted
  
  head(local_hunt)
  head(sno_hunt)
  
  # Step 1: Classify observations in local_hunt as summer or winter
  local_hunt$season <- ifelse(local_hunt$t_month %in% summer_removal, "summer", "winter")
  # Step 2: hunted observations in local_hunt
  local_hunted <- aggregate(
    x = list(hunted = local_hunt$v_hunter_id), 
    by = list(t_hunting_year = local_hunt$t_hunting_year, season = local_hunt$season), 
    FUN = length
  )
  # Reshape the data to wide format (separate columns for summer and winter hunted)
  local_hunted_wide <- reshape(
    local_hunted,
    timevar = "season",
    idvar = "t_hunting_year",
    direction = "wide"
  )
  # Replace NA with 0 for missing hunted
  local_hunted_wide$hunted.summer[is.na(local_hunted_wide$hunted.summer)] <- 0
  local_hunted_wide$hunted.winter[is.na(local_hunted_wide$hunted.winter)] <- 0
  # Step 3: Combine with sno_hunt
  # Filter sno_hunt for sno_varanger and rename v_abundance to winter_hunted_sno
  sno_hunt_filtered <- sno_hunt[sno_hunt$v_hunter_id == "sno_varanger", ]
  colnames(sno_hunt_filtered)[colnames(sno_hunt_filtered) == "v_abundance"] <- "winter_hunted_sno"
  # Merge local_hunt hunted with sno_hunt
  result <- merge(sno_hunt_filtered, local_hunted_wide, by = "t_hunting_year", all.x = TRUE)
  # Replace NA with 0 for missing hunted in local_hunt
  result$hunted.summer[is.na(result$hunted.summer)] <- 0
  result$hunted.winter[is.na(result$hunted.winter)] <- 0
  # Step 4: Calculate totals
  result$winter_hunted <- result$winter_hunted_sno + result$hunted.winter
  result$summer_hunted <- result$hunted.summer
  result$total_hunted <- result$winter_hunted + result$summer_hunted
   
  ((result$winter_hunted + result$summer_hunted) - result$total_hunted) #oke so correct
  
  final_result <- result[, c("t_hunting_year", "total_hunted", "winter_hunted", "summer_hunted")]
    #=== return ===
  return(final_result)
  
}
