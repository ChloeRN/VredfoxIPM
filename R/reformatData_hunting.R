#' Reformat the hunting data to make a combined total number foxes hunted for each winter hunting period, summer hunting period, and "yearly" hunting period 
#' (from start of summer to end of winter hunt). Includes all foxes, male and female, in Varanger.
#'
#' @param summer_removal a vector of numerical months to be removed from winter age at harvest data: c(6,7,8,9)
#' @param hunting.dataset list containing the datasets downloaded from the COAT dataportal.
#'
#' @return a dataframe containing the reported numbers of hunted foxes in total and for summer, winter, and unknown season
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
 
  ## Fix error in SNO data (version 3 only)
  if(nrow(sno_hunt) == 20 & all(sno_hunt$t_hunting_year[19:20] == "2022_2023")){
    sno_hunt$t_hunting_year[20] <- "2023_2024"
  }
  
 #========== PREPARE DATA ===========
  #Oke here we just assume that SNO hunts only in winter
  
  #we want for each t_hunting year, the total, summer, and winter number of foxes hunted
  
  # Step 1: Classify observations in local_hunt as summer or winter
  local_hunt$season <- dplyr::case_when(is.na(local_hunt$t_month) ~ "unknown",
                                        local_hunt$t_month %in% summer_removal ~ "summer",
                                        TRUE ~ "winter")

  # Step 2: hunted observations in local_hunt
  local_hunted <- aggregate(
    x = list(hunted = local_hunt$v_hunter_id), 
    by = list(t_hunting_year = local_hunt$t_hunting_year, season = local_hunt$season), 
    FUN = length
  )
  
  # Reshape the data to wide format (separate columns for summer and winter hunted)
  local_hunted_wide <- stats::reshape(
    local_hunted,
    timevar = "season",
    idvar = "t_hunting_year",
    direction = "wide"
  )
  
  # Replace NA with 0 for missing hunted
  local_hunted_wide$hunted.summer[is.na(local_hunted_wide$hunted.summer)] <- 0
  local_hunted_wide$hunted.winter[is.na(local_hunted_wide$hunted.winter)] <- 0
  local_hunted_wide$hunted.unknown[is.na(local_hunted_wide$hunted.unknown)] <- 0
  
  # Step 3: Combine with sno_hunt
  # Filter sno_hunt for sno_varanger and rename v_abundance to winter_hunted_sno
  sno_hunt_filtered <- sno_hunt[sno_hunt$v_hunter_id == "sno_varanger", ]
  colnames(sno_hunt_filtered)[colnames(sno_hunt_filtered) == "v_abundance"] <- "hunted.unknown.sno"
  
  # Merge local_hunt hunted with sno_hunt
  result <- local_hunted_wide %>%
    dplyr::full_join(sno_hunt_filtered[c("t_hunting_year", "hunted.unknown.sno")], by = "t_hunting_year")
  
  # Replace NA with 0 for missing hunted in sno_hunt
  result$hunted.unknown.sno[is.na(result$hunted.unknown.sno)] <- 0

  # Step 4: Calculate totals
  message("Assuming that all SNO kills are from winter hunting season (Oct-May).")
  
  final_result <- result %>%
    dplyr::mutate(winter_hunted = hunted.winter + hunted.unknown.sno,
                  summer_hunted = hunted.summer,
                  unknown_hunted = hunted.unknown, #+ hunted.unknown.sno,
                  total_hunted = hunted.winter + hunted.summer + hunted.unknown + hunted.unknown.sno) %>%
    dplyr::select(t_hunting_year, total_hunted, winter_hunted, summer_hunted, unknown_hunted) %>%
    dplyr::arrange(t_hunting_year)
  
  #=== return ===
  return(final_result)
  
}
