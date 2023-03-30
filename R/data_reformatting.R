library(reshape2)
library(sf)

#' Reformat the carcass data to make Age at harvest matrix & P1var (nr of embryos/scars) & P2var (breeding or not breeding)
#'
#' @param n_ageclasses integer. Number of age classes to consider in analyses.
#' @param summer_removal a vector of numerical months to be removed from age at harvest data: c(6,7,8,9)
#' @param area_selection a vector of study-area sub-area names to consider in the analyses: c("Inner", "BB",  "Tana")
#' @param plac_start integer. Julian day (including) from which we use placental scars presence to calculate breeding or not breeding for P2
#' @param plac_end integer. Julian day (not including)until  we use placental scars presence to calculate breeding or not breeding for P2
#' @param embr_start integer. Julian day (including) from which we use embryos presence to calculate breeding or not breeding for P2
#' @param embr_end integer. Julian day (not including) until we use embryos presence to calculate breeding or not breeding for P2
#' @param carcass.dir directory of carcass data
#' @param shapefile.dir directory of shapefile data with study areas sub-areas
#'
#' @return a list containing the Age at harvest matrix, P1var and P2var
#' @export
#'
#' @examples


data_reformatting <- function (n_ageclasses, summer_removal, area_selection,
                              plac_start, plac_end , embr_start, embr_end,
                              carcass.dir, shapefile.dir) {
  
  # === 2 things that I use later down ====
  '%notin%' <- Negate('%in%')
  max_age <- (n_ageclasses-1)
  
  #========= LOAD DATA ==============

  carcass_filenames<-dir(paste(carcass.dir, sep = "/"))
  mylist<-c()
  for (i in 1:length(carcass_filenames)){
    mylist[[i]]<-read.table(paste(carcass.dir, carcass_filenames[i], sep = "/"), header=T, sep = ";")
  }
  myfile<-do.call(rbind, mylist)  # combine all files
  allf <- myfile
  rm(myfile)
  
  shapefile <- st_read(paste(shapefile.dir, sep = "/"))
  shapefile <- st_transform(shapefile, crs = 4326)
  
  #check if loading carcass data worked
  if(!exists("allf") || !length(allf)==45){
    stop("Carcass data not loaded properly or not 45 columns")
  }
  
  #check if loading shapefile data worked
  if(!exists("shapefile")){
    stop("Shapefile not loaded properly")
  }
  
  #========= PREPARE DATA ==============
  #hunting date as date
  allf$t_hunting_date <- as.Date(allf$t_hunting_date)
  
  #we need hunting year as an actual number sometimes so we can add and substract
  allf$start_hunting_year<-substr(allf$t_hunting_year, start = 1, stop = 4)
  allf$start_hunting_year<- as.numeric(allf$start_hunting_year)
  
  #Add julian day
  allf$julian <- as.POSIXlt(allf$t_hunting_date)$yday  
  
  #Numeric months
  allf$mnd <- as.numeric(format(allf$t_hunting_date, "%m"))
  
  #Remove foxes without hunting date
  fvar1 <- allf[is.na(allf$t_hunting_date)==F,]
  
  #Select Varanger and female only
  fvar1 <- fvar1[is.na(fvar1$sn_region)== F & fvar1$sn_region == "varanger" & is.na(fvar1$v_sex)== F & fvar1$v_sex == "female", ]
  
  #Remove the ones that were not shot
  fvar1 <- fvar1[fvar1$v_hunter_id != "road_kill", ]
  
  #Assigning study area - sub area names with shapefile
  fvar1   <- fvar1[is.na(fvar1$e_dd)==FALSE & is.na(fvar1$n_dd)==FALSE,] #only foxes with valid location
  sf_fvar1<- st_as_sf(fvar1, coords = c("e_dd", "n_dd"), crs="+proj=longlat +datum=WGS84") # change dataset in to sf object with geometry (instead of lat/long)
  sf_fvar1<- st_join(sf_fvar1, shapefile["name"], join = st_within) #overlap of dataset positions with shapefile of study area - sub areas, taking the name of the subarea
  colnames(sf_fvar1)[colnames(sf_fvar1) == 'name'] <- 'sub_area'
  coords  <- data.frame(st_coordinates(sf_fvar1)) #get back coordinates column from geometry
  sf_fvar1$e_dd <- coords$X # getting back lat long columns
  sf_fvar1$n_dd <- coords$Y
  fvar1 <- st_drop_geometry(sf_fvar1) #remove the geometry again and get back to fvar1, #fvar1 now has shape file area names
  
  #check if shapefile subarea assignment worked and subarea names match those given in script
  if(sum(is.na(fvar1$sub_area)) > 0 || !identical(unique(fvar1$sub_area), c("Inner", "Tana",  "BB"))){
    stop("Sub-area assignment with shapefile produced NA's, or sub-area names do not match those in script")
  }
  
  #Selecting for study area - sub areas
  fvar1 <- subset(fvar1,(fvar1$sub_area %in% area_selection))
  
  #Age class selection
  fvar1$alder4 <- fvar1$v_age
  fvar1$alder4[fvar1$alder4 > max_age] <- max_age
  
  #=============== AGE AT HARVEST MATRIX BUILDING ==============================
  #here we exclude foxes shot in summer months and foxes with no age info
  
  var.F.C <- as.matrix(table(fvar1$t_hunting_year[is.na(fvar1$v_age)==F & fvar1$mnd %notin% summer_removal], 
                             fvar1$alder4[is.na(fvar1$v_age)==F & fvar1$mnd %notin% summer_removal]))
  
  varFC <- as.data.frame (var.F.C)
  varFC2 <- dcast(varFC, Var1~Var2)
  
  column_names <- c("year","age0", "age1", "age2", "age3", "age4", "age5", "age6", "age7", "age8", "age9", "age10", "age11", "age12", "age13", "age14", "age15", "age16", "age17", "age18")
  names(varFC2) <- head(column_names, n_ageclasses+1)
  
  #proportion that was aged
  agedf.ann <- table(fvar1$t_hunting_year[is.na(fvar1$v_age)==F & fvar1$mnd %notin% summer_removal]) #where unaged removed
  allf.ann  <- table(fvar1$t_hunting_year[                        fvar1$mnd %notin% summer_removal]) #where unaged not removed
  prop <- agedf.ann/allf.ann
  varFC2$pData <- prop
  
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
  carcassData<- list(AaH.matrix = varFC2, 
                     P1var = P1var,
                     P2var = P2var)
  
  #=== return ===
  return(carcassData)
  
}


carcassData<- data_reformatting (
  n_ageclasses       = 5, 
  summer_removal     = c(6,7,8,9) ,
  area_selection     = c("Inner", "BB",  "Tana"),
  plac_start         = 140,
  plac_end           = 80 ,
  embr_start         = 100 ,
  embr_end           = 140,
  carcass.dir        ="C:\\Users\\sho189\\OneDrive - UiT Office 365\\PhD\\RedfoxIPM\\Data from google disk\\carcass_examination",
  shapefile.dir      ="C:\\Users\\sho189\\OneDrive - UiT Office 365\\PhD\\RedfoxIPM\\Fox areas shapefile\\tana rest"
)
