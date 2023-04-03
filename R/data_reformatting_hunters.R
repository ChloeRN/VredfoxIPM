hunting.dir        <-"C:\\Users\\sho189\\OneDrive - UiT Office 365\\PhD\\RedfoxIPM\\Data from google disk\\hunting"
shapefile.dir      <- "C:\\Users\\sho189\\OneDrive - UiT Office 365\\PhD\\RedfoxIPM\\Fox areas shapefile\\tana rest"


#========= LOAD DATA ==============

hunting_filenames<-dir(paste(hunting.dir, sep = "/"))
mylist<-c()
for (i in 1:length(hunting_filenames)){
  mylist[[i]]<-read.table(paste(hunting.dir, hunting_filenames[i], sep = "/"), header=T, sep = ";")
}
myfile<-do.call(rbind, mylist)  # combine all files
allhunt <- myfile
rm(myfile)

shapefile <- sf::st_read(paste(shapefile.dir, sep = "/"))
shapefile <- sf::st_transform(shapefile, crs = 4326)

#check if loading carcass data worked
if(!exists("allhunt") || !length(allhunt)==16){
  stop("Hunting data not loaded properly or not 16 columns")
}

#check if loading shapefile data worked
if(!exists("shapefile")){
  stop("Shapefile not loaded properly")
}


#========= PREPARE DATA ==============
#we need hunting year as an actual number sometimes so we can add and substract
allhunt$start_hunting_year<-substr(allhunt$t_hunting_year, start = 1, stop = 4)
allhunt$start_hunting_year<- as.numeric(allhunt$start_hunting_year)

#Select Varanger only
hvar1 <- allhunt[is.na(allhunt$sn_region)== F & allhunt$sn_region == "varanger" ,]

#Remove the ones that were not shot
hvar1 <- hvar1[hvar1$v_hunter_id != "road_kill", ]
hvar1 <- hvar1[hvar1$v_hunting_method != "road_kill" & hvar1$v_hunting_method != "found_dead", ]

#hmmm here found_dead and road_kill in v_hunting method instead of v_hunter_id.. check if I didnt miss any in carcass.reformatting script. how many road kill there?
check<- allhunt[is.na(allhunt$v_hunting_method)==F & (allhunt$v_hunting_method == "road_kill" | allhunt$v_hunting_method == "found_dead"), ]
check2<- allf[is.na(allf$v_hunter_id)==F & (allf$v_hunter_id == "road_kill" | allf$v_hunter_id == "found_dead"), ]
#116 -83 non shot foxes "went missing". Possibly 33 not shot foxes registered as hunted? or were all v_hunting method non-shot transferred to v_hunting. ID?


#Assigning study area - sub area names with shapefile
#we have municipality or coordinates here, but not both. So use coordinates to get municipality, then use municipality for sub-areas. Although they do not overlap entirely with shapefile borders


hvar1   <- hvar1[is.na(hvar1$e_dd)==FALSE & is.na(hvar1$n_dd)==FALSE,] #only foxes with valid location
sf_hvar1<- sf::st_as_sf(hvar1, coords = c("e_dd", "n_dd"), crs="+proj=longlat +datum=WGS84") # change dataset in to sf object with geometry (instead of lat/long)
sf_hvar1<- sf::st_join(sf_hvar1, shapefile["name"], join = sf::st_within) #overlap of dataset positions with shapefile of study area - sub areas, taking the name of the subarea
colnames(sf_hvar1)[colnames(sf_hvar1) == 'name'] <- 'sub_area'
coords  <- data.frame(sf::st_coordinates(sf_hvar1)) #get back coordinates column from geometry
sf_hvar1$e_dd <- coords$X # getting back lat long columns
sf_hvar1$n_dd <- coords$Y
hvar1 <- sf::st_drop_geometry(sf_hvar1) #remove the geometry again and get back to hvar1, #hvar1 now has shape file area names

#check if shapefile subarea assignment worked and subarea names match those given in script
if(sum(is.na(hvar1$sub_area)) > 0 || !identical(unique(hvar1$sub_area), c("Inner", "Tana",  "BB"))){
  stop("Sub-area assignment with shapefile produced NA's, or sub-area names do not match those in script")
}

#Selecting for study area - sub areas
hvar1 <- subset(hvar1,(hvar1$sub_area %in% area_selection))

