#' Reformat the carcass data to make Age at harvest matrix & P1var (nr of embryos/scars) & P2var (breeding or not breeding)
#'
#' @param Amax  integer. Number of age classes to consider in analyses.
#' @param summer_removal a vector of numerical months to be removed from winter age at harvest data: c(6,7,8,9)
#' @param winter_removal a vector of numerical months to be removed from summer age at harvest data: c(1:5, 10:12)
#' @param area_selection a vector of study-area sub-area names to consider in the analyses: c("Inner", "BB",  "Tana")
#' @param plac_start integer. Julian day (including) from which we use placental scars presence to calculate breeding or not breeding for P2
#' @param plac_end integer. Julian day (not including)until  we use placental scars presence to calculate breeding or not breeding for P2
#' @param embr_start integer. Julian day (including) from which we use embryos presence to calculate breeding or not breeding for P2
#' @param embr_end integer. Julian day (not including) until we use embryos presence to calculate breeding or not breeding for P2
#' @param carcass.dataset dataframe containing the dataset downloaded from the COAT dataportal.
#' @param shapefile.dir string. Directory of shapefiles delineating study areas and sub-areas.
#' @param add.sumr.unaged logical. Add summer harvested individuals as unaged individuals to the total harvested individuals and their proportion aged.
#' @param saAH.years a vector of years for which the summer age at harvest matrix should be constructed
#' @param hunting.data a dataframe containing the reported numbers of hunted foxes in total and for both seasons
#'
#' @return a list containing the age-at-harvest matrix and dataframes with embryo count (P1var) and placental scar presence-absence (P2var) data.
#' @export
#'
#' @examples


reformatData_carcass <- function (Amax, summer_removal, winter_removal, area_selection,
                              plac_start, plac_end , embr_start, embr_end,
                              carcass.dataset, 
                              shapefile.dir,
                              add.sumr.unaged, saAH.years,
                              hunting.data) {

  #========= HELPER FUNCTIONS ==============
  '%notin%' <- Negate('%in%')
  
  #========= LOAD DATA ==============
  allf <- carcass.dataset
  
  shapefile <- sf::st_read(paste(shapefile.dir, sep = "/"))
  shapefile <- sf::st_transform(shapefile, crs = 4326)
  
  #check if loading carcass data worked
  if(!exists("allf") || !length(allf)==42){
    stop("Carcass data not loaded properly or not 42 columns")
  }
  
  #check if loading shapefile data worked
  if(!exists("shapefile")){
    stop("Shapefile not loaded properly")
  }
  
  #========= PREPARE DATA ==============
  #hunting date as date
  allf$t_hunting_date <- as.Date(allf$t_hunting_date)
  
  #we need hunting year as an actual number sometimes so we can add and substract
  allf$start_hunting_year <- substr(allf$t_hunting_year, start = 1, stop = 4)
  allf$start_hunting_year <- as.numeric(allf$start_hunting_year)
  
  #Add julian day
  allf$julian <- as.POSIXlt(allf$t_hunting_date)$yday  
  
  #Numeric months
  allf$mnd <- as.numeric(format(allf$t_hunting_date, "%m"))
  
  #Select only those foxes that are actually shot in Varanger (male and female)
  all_varanger <- allf[is.na(allf$sn_region)== F & allf$sn_region == "varanger"  & allf$v_hunter_id != "road_kill" & allf$v_hunter_id != "found_dead", ]
  
  #select only females
  female_shot_varanger <- all_varanger[is.na(all_varanger$v_sex)== F & all_varanger$v_sex == "female", ]
  
  #Remove females without hunting date
  fvar1 <- female_shot_varanger[is.na(female_shot_varanger$t_hunting_date)==F,]
  
  #Assigning study area - sub area names with shapefile
  fvar1 <- fvar1[is.na(fvar1$e_dd)==FALSE & is.na(fvar1$n_dd)==FALSE,] #only foxes with valid location
  sf_fvar1 <- sf::st_as_sf(fvar1, coords = c("e_dd", "n_dd"), crs="+proj=longlat +datum=WGS84") # change dataset in to sf object with geometry (instead of lat/long)
  sf_fvar1 <- sf::st_join(sf_fvar1, shapefile["name"], join = sf::st_within) #overlap of dataset positions with shapefile of study area - sub areas, taking the name of the subarea
  colnames(sf_fvar1)[colnames(sf_fvar1) == 'name'] <- 'sub_area'
  coords <- data.frame(sf::st_coordinates(sf_fvar1)) #get back coordinates column from geometry
  sf_fvar1$e_dd <- coords$X # getting back lat long columns
  sf_fvar1$n_dd <- coords$Y
  fvar1 <- sf::st_drop_geometry(sf_fvar1) #remove the geometry again and get back to fvar1, #fvar1 now has shape file area names
  
  #check if shapefile subarea assignment worked and subarea names match those given in script
  if(sum(is.na(fvar1$sub_area)) > 0 || !identical(unique(fvar1$sub_area), c("Inner", "Tana",  "BB"))){
    stop("Sub-area assignment with shapefile produced NA's, or sub-area names do not match those in script")
  }
  
  #Selecting for study area - sub areas
  fvar1 <- subset(fvar1,(fvar1$sub_area %in% area_selection))
  
  #Age class selection
  fvar1$alder4 <- fvar1$v_age
  fvar1$alder4[fvar1$alder4 > (Amax-1)] <- (Amax-1)
  
  
  
  #=============== Proportion of foxes that end up in carcass examination from hunting data ========================
  
  ## Count total number of registered carcasses (per hunting year and season)
  
  # Step 1: Classify observations as summer, winter, or unknown
  all_varanger$season <- dplyr::case_when(is.na(all_varanger$mnd) ~ "unknown",
                                          all_varanger$mnd %in% summer_removal ~ "summer",
                                          TRUE ~ "winter")
  
  # Step 2: Count total foxes per year
  total_carcass <- aggregate(v_individual_id ~ t_hunting_year, 
                             data = all_varanger, 
                             FUN = length)
  colnames(total_carcass) <- c("t_hunting_year", "total_carcass")

  # Step 3: Count summer foxes per year
  summer_carcass <- aggregate(v_individual_id ~ t_hunting_year, 
                             data = all_varanger[all_varanger$season == "summer", ], 
                             FUN = length)
  colnames(summer_carcass) <- c("t_hunting_year", "summer_carcass")
  
  # Step 4: Count winter foxes per year
  winter_carcass <- aggregate(v_individual_id ~ t_hunting_year, 
                             data = all_varanger[all_varanger$season == "winter", ], 
                             FUN = length)
  colnames(winter_carcass) <- c("t_hunting_year", "winter_carcass")
  
  # Step 5: Count unknwon season foxes per year
  unknown_carcass <- aggregate(v_individual_id ~ t_hunting_year, 
                              data = all_varanger[all_varanger$season == "unknown", ], 
                              FUN = length)
  colnames(unknown_carcass) <- c("t_hunting_year", "unknown_carcass")
  
  # Step 6: Merge the results
  result <- merge(total_carcass, winter_carcass, by = "t_hunting_year", all.x = TRUE)
  result <- merge(result, summer_carcass, by = "t_hunting_year", all.x = TRUE)
  result <- merge(result, unknown_carcass, by = "t_hunting_year", all.x = TRUE)
  
  # Replace NA with 0 for missing carcass
  result$summer_carcass[is.na(result$summer_carcass)] <- 0
  result$winter_carcass[is.na(result$winter_carcass)] <- 0
  result$unknown_carcass[is.na(result$unknown_carcass)] <- 0
  
  carcass.numbers <- result
  
  #----checking the difference----
  # Ensure the data frames are aligned by `t_hunting_year`
  # Merge the two data frames by `t_hunting_year` to ensure alignment
  merged_data <- merge(carcass.numbers, hunting.data, by = "t_hunting_year", suffixes = c("_carcass", "_hunting"))
  
  # Calculate the differences for each column
  difference_data <- data.frame(
    t_hunting_year = merged_data$t_hunting_year,
    total_count = merged_data$total_hunted - merged_data$total_carcass,
    total_winter = merged_data$winter_hunted - merged_data$winter_carcass,
    total_summer = merged_data$summer_hunted - merged_data$summer_carcass,
    total_unknown = merged_data$unknown_hunted - merged_data$unknown_carcass
  )
  
  #--------Then I want to get to:----------
  #pCarcassData[t] = Number of foxes in carcass data for hunting season t / Number of foxes in hunting data for hunting season t (again for the summer, winter, and combined hunting seasons).
  
  # Calculate the proportions
  proportion_data <- data.frame(
    t_hunting_year = merged_data$t_hunting_year,
    pCarcass_total = merged_data$total_carcass / merged_data$total_hunted,
    pCarcass_winter = merged_data$winter_carcass / merged_data$winter_hunted,
    pCarcass_summer = merged_data$summer_carcass / merged_data$summer_hunted
  )
  
  # Set values > 1 to 1
  proportion_data$pCarcass_total[which(proportion_data$pCarcass_total > 1)] <- 1
  proportion_data$pCarcass_winter[which(proportion_data$pCarcass_winter > 1)] <- 1
  proportion_data$pCarcass_summer[which(proportion_data$pCarcass_summer > 1)] <- 1
  
  # Assign starting year (for matching with age-at-harvest data)
  proportion_data$year <- as.integer(stringr::str_split_fixed(proportion_data$t_hunting_year, pattern = "_", n = 2)[,1])
  
  
  #===============    WINTER AGE AT HARVEST MATRIX BUILDING ==============================
  #here we exclude foxes shot in summer months and foxes with no age info
  
  var.F.C <- as.matrix(table(fvar1$start_hunting_year[!is.na(fvar1$v_age) & fvar1$mnd %notin% summer_removal], 
                             fvar1$alder4[!is.na(fvar1$v_age) & fvar1$mnd %notin% summer_removal]))
  
  varFC2 <- data.frame(cbind(as.numeric(rownames(var.F.C)), unname(var.F.C)))
  colnames(varFC2) <- c("year", paste0("age", (1:Amax)-1))
  
  
  #proportion that was aged
  agedf.ann <- table(fvar1$start_hunting_year[!is.na(fvar1$v_age) & fvar1$mnd %notin% summer_removal]) #where unaged removed, summer also removed
  
  if(add.sumr.unaged){  
    allf.ann <- table(fvar1$start_hunting_year) #where unaged not removed, summer not removed
    allf.win <- table(fvar1$start_hunting_year[fvar1$mnd %notin% summer_removal])
  }else{
    allf.ann <- table(fvar1$start_hunting_year[fvar1$mnd %notin% summer_removal]) #where unaged not removed, summer also removed
  }
  
  prop <- agedf.ann/allf.ann
  
  varFC2$pAged <- prop
  
  if(add.sumr.unaged){
    prop.win <- agedf.ann/allf.win
    varFC2$pAged_winter <- prop.win
  }
  
  # Add proportion of hunted in carcass data
  varFC2 <- varFC2 %>%
    dplyr::left_join(proportion_data[, c("year", "pCarcass_total", "pCarcass_winter")], by = "year")
  
  if(add.sumr.unaged){
    varFC2 <- varFC2 %>%
      dplyr::rename(pCarcass = pCarcass_total)
  }else{
    varFC2 <- varFC2 %>%
      dplyr::rename(pCarcass = pCarcass_winter) %>%
      dplyr::select(-pCarcass_total)
  }
  
  
  #===============    SUMMER AGE AT HARVEST MATRIX BUILDING ==============================
  #here we exclude foxes shot in winter months and foxes with no age info
  
  # NOTE: Years 2004 and 2016 are missing summer harvest. 
  # If data for these years should be included (part of saAH.years) the 0's 
  # need to be coded in manually as the current code drops any years with
  # 0 observations from the summaries (tables, matrices, vectors) 
  
  svar.F.C <- as.matrix(table(fvar1$start_hunting_year[!is.na(fvar1$v_age) & fvar1$mnd %notin% winter_removal & fvar1$start_hunting_year %in% saAH.years], 
                             fvar1$alder4[!is.na(fvar1$v_age) & fvar1$mnd %notin% winter_removal & fvar1$start_hunting_year %in% saAH.years]))
  
  svarFC2 <- data.frame(cbind(as.numeric(rownames(svar.F.C)), unname(svar.F.C)))
  colnames(svarFC2) <- c("year", paste0("age", (1:Amax)-1))
  
  #proportion that was aged
  sagedf.ann <- table(fvar1$start_hunting_year[!is.na(fvar1$v_age) & fvar1$mnd %notin% winter_removal & fvar1$start_hunting_year %in% saAH.years]) #where unaged removed
  sallf.ann  <- table(fvar1$start_hunting_year[fvar1$mnd %notin% winter_removal & fvar1$start_hunting_year %in% saAH.years]) #where unaged not removed
  sprop <- sagedf.ann/sallf.ann
  svarFC2$pAged <- sprop
  
  # Add proportion of hunted in carcass data
  svarFC2 <- svarFC2 %>%
    dplyr::left_join(proportion_data[, c("year", "pCarcass_summer")], by = "year") %>%
    dplyr::rename(pCarcass = pCarcass_summer)
  
  # ========= REPRODUCTION: NR OF EMBRYO'S / PLACENTAL SCARS =========
  #here we exclude foxes with no age info and foxes with no breeding info recorded (no NA, and nr > 0)
  
  P1var.pl <- fvar1[is.na(fvar1$v_nb_placental_scars)==F & fvar1$v_nb_placental_scars > 0 & is.na(fvar1$v_age)==F,
                    c("v_nb_placental_scars", "v_age", "start_hunting_year", "julian")]
  P1var.pl$type <- "pl.scars"
  
  P1var.emb <- fvar1[is.na(fvar1$v_embryos)==F & fvar1$v_embryos > 0 & is.na(fvar1$v_age)==F ,
                     c("v_embryos", "v_age", "start_hunting_year", "julian")]
  P1var.emb$type <- "embryos"
  
  names(P1var.pl)[1:2] <- c("P1", "P1_age")
  names(P1var.emb)[1:2] <- c("P1", "P1_age")
  
  P1var.pl$repryear <- P1var.pl$start_hunting_year
  P1var.emb$repryear <- P1var.emb$start_hunting_year+1
  
  P1var <- rbind(P1var.pl, P1var.emb)
  
  # ========= REPRODUCTION: PROPORTION REPRODUCING  =========
  #here we exclude foxes with no age info & 0 year olds for placental scars & foxes with NA for placental scars & foxes with no breeding info recorded, & select for time periods in year
  
  P2var.pl <- fvar1[is.na(fvar1$v_placental_scars)==F & is.na(fvar1$v_age)==F & fvar1$v_age> 0 & (fvar1$julian < plac_end | fvar1$julian >= plac_start),
                    c("v_placental_scars", "v_age", "start_hunting_year", "julian")]
  
  P2var.pl$type <- "pl.scars"
  P2var.pl$repryear <- P2var.pl$start_hunting_year
  names(P2var.pl)[1:2] <- c("P2", "P2_age")
  
  P2var.preg <- fvar1[is.na(fvar1$v_breeding)==F & is.na(fvar1$v_age)==F & (fvar1$julian >= embr_start & fvar1$julian < embr_end), 
                      c("v_breeding", "v_age", "start_hunting_year", "julian")]
  
  P2var.preg$type <- "pregnant"
  P2var.preg$repryear <- P2var.preg$start_hunting_year + 1
  names(P2var.preg)[1:2] <- c("P2", "P2_age")
  
  P2var <- rbind(P2var.pl, P2var.preg)
  
  
  
  ## Combine age at harvest, nr of embryos, and prop breeding files in a list
  carcassData <- list(WAaH.matrix = varFC2,
                      SAaH.matrix = svarFC2,
                      P1var = P1var,
                      P2var = P2var)
  
  #=== return ===
  return(carcassData)
  
}
