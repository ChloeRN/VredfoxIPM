#======= OPTION 1 : use hunting data files =========

hunting.dir        <-"C:\\Users\\sho189\\OneDrive - UiT Office 365\\PhD\\RedfoxIPM\\Data from google disk\\hunting"
shapefile.dir      <-"C:\\Users\\sho189\\OneDrive - UiT Office 365\\PhD\\RedfoxIPM\\Fox areas shapefile\\tana rest"


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
hvar1 <- hvar1[hvar1$v_hunting_method != "road_kill" & hvar1$v_hunting_method != "found_dead" | is.na(hvar1$v_hunting_method), ]

#hmmm here found_dead and road_kill in v_hunting method instead of v_hunter_id.. check if I didnt miss any in carcass.reformatting script. how many road kill there?
check<- allhunt[is.na(allhunt$v_hunting_method)==F & (allhunt$v_hunting_method == "road_kill" | allhunt$v_hunting_method == "found_dead"), ]
#check2<- allf[is.na(allf$v_hunter_id)==F & (allf$v_hunter_id == "road_kill" | allf$v_hunter_id == "found_dead"), ]
#116 -83 = 33 non shot foxes "went missing". Possibly 33 not shot foxes registered as hunted? or were all v_hunting method non-shot transferred to v_hunter_ID?


#Assigning study area - sub area names with shapefile
#we have municipality or coordinates here, but not both. So use coordinates to get municipality, then use municipality for sub-areas. Although they do not overlap entirely with shapefile borders

#library(ggplot2)
#ggplot(hvar1, aes(x = v_municipality)) +
#  geom_bar()


#option 1: just use all areas
hunterfile.data<- data.frame()

y <- unique(hvar1$start_hunting_year[!is.na(hvar1$start_hunting_year)])
for(i in y){
  newdata <- data.frame(
    Year =i,
    NHunters = length(unique(hvar1$v_hunter_id[hvar1$start_hunting_year==i ]))
  )
  hunterfile.data<-rbind(hunterfile.data, newdata)
}
#set 2004 to NA?
hunterfile.data$NHunters[hunterfile.data$Year==2004]<-NA

#these are different numbers than the manual ones...?
hunterfile.data$doro.manual<- c(NA, 41, 39, 44, 38, 38, 42, 51, 72, 35, 44, 47, 33, 45, 43, 44,NA)


#select some sub area
#hvar1   <- hvar1[is.na(hvar1$e_dd)==FALSE & is.na(hvar1$n_dd)==FALSE,] #only foxes with valid location
#sf_hvar1<- sf::st_as_sf(hvar1, coords = c("e_dd", "n_dd"), crs="+proj=longlat +datum=WGS84") # change dataset in to sf object with geometry (instead of lat/long)
#sf_hvar1<- sf::st_join(sf_hvar1, shapefile["name"], join = sf::st_within) #overlap of dataset positions with shapefile of study area - sub areas, taking the name of the subarea
#colnames(sf_hvar1)[colnames(sf_hvar1) == 'name'] <- 'sub_area'
#coords  <- data.frame(sf::st_coordinates(sf_hvar1)) #get back coordinates column from geometry
#sf_hvar1$e_dd <- coords$X # getting back lat long columns
#sf_hvar1$n_dd <- coords$Y
#hvar1 <- sf::st_drop_geometry(sf_hvar1) #remove the geometry again and get back to hvar1, #hvar1 now has shape file area names

#check if shapefile subarea assignment worked and subarea names match those given in script
#if(sum(is.na(hvar1$sub_area)) > 0 || !identical(unique(hvar1$sub_area), c("Inner", "Tana",  "BB"))){
#  stop("Sub-area assignment with shapefile produced NA's, or sub-area names do not match those in script")
#}

#Selecting for study area - sub areas
#hvar1 <- subset(hvar1,(hvar1$sub_area %in% area_selection))




#======= OPTION 2: why not use the allf carcass datafile with the hunter ID column, that one includes coordinates ===========0

library(reshape2)
library(sf)
'%notin%' <- Negate('%in%')

# choosing varanger sub area ("Inner" / "BB" / "Tana)
area_selection<- c("Inner")
## set directories
carcass.dir        <-"C:\\Users\\sho189\\OneDrive - UiT Office 365\\PhD\\RedfoxIPM\\Data from google disk\\carcass_examination"
shapefile.dir      <-"C:\\Users\\sho189\\OneDrive - UiT Office 365\\PhD\\RedfoxIPM\\Fox areas shapefile\\tana rest"

#========= LOAD DATA ==============
## load data
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

#========= PREPARE DATA ==============
#we need hunting year as an actual number sometimes so we can add and substract
allf$start_hunting_year<-substr(allf$t_hunting_year, start = 1, stop = 4)
allf$start_hunting_year<- as.numeric(allf$start_hunting_year)

#Select Varanger only
fvar1 <- allf[is.na(allf$sn_region)== F & allf$sn_region == "varanger" , ]

#Remove the ones that were not shot
fvar1 <- fvar1[fvar1$v_hunter_id != "road_kill" & fvar1$v_hunter_id != "found_dead", ]

#Assigning study area - sub area names with shapefile
fvar1   <- fvar1[is.na(fvar1$e_dd)==FALSE & is.na(fvar1$n_dd)==FALSE,] #only foxes with valid location
sf_fvar1<- st_as_sf(fvar1, coords = c("e_dd", "n_dd"), crs="+proj=longlat +datum=WGS84") # change dataset in to sf object with geometry (instead of lat/long)
sf_fvar1<- st_join(sf_fvar1, shapefile["name"], join = st_within) #overlap of dataset positions with shapefile of study area - sub areas, taking the name of the subarea
colnames(sf_fvar1)[colnames(sf_fvar1) == 'name'] <- 'sub_area'
coords  <- data.frame(st_coordinates(sf_fvar1)) #get back coordinates column from geometry
sf_fvar1$e_dd <- coords$X # getting back lat long columns
sf_fvar1$n_dd <- coords$Y
fvar1 <- st_drop_geometry(sf_fvar1) #remove the geometry again and get back to fvar1, #fvar1 now has shape file area names

#Selecting for study area - sub areas
fvar1 <- subset(fvar1,(fvar1$sub_area %in% area_selection))

unique(fvar1$sub_area)
unique(fvar1$v_hunter_id)

length(unique(fvar1$v_hunter_id)) #7 hunters less than in the hunting file
length(fvar1[fvar1$v_hunter_id=="sno_varanger",]) #the carcass file somehow has more foxes than the hunting file


#nr of successfull hunters per year
hunter.data<- data.frame()

y <- unique(fvar1$start_hunting_year[!is.na(fvar1$start_hunting_year)])
for(i in y){
  newdata <- data.frame(
    Year =i,
    NHunters = length(unique(fvar1$v_hunter_id[fvar1$start_hunting_year==i ]))
  )
  hunter.data<-rbind(hunter.data, newdata)
}
#set 2004 to NA?
hunter.data$NHunters[hunter.data$Year==2004]<-NA

#differences other data sources
hunter.data$from.huntingfiles <- c(hunterfile.data$NHunters,NA)
hunter.data$doro.manual<- c(NA, 41, 39, 44, 38, 38, 42, 51, 72, 35, 44, 47, 33, 45, 43, 44,NA, NA)

#make a plot over time of all the different data options
plot(hunter.data$Year, hunter.data$NHunters, type = "l", ylim =c(10, 80))
lines(hunter.data$Year, hunter.data$from.huntingfiles, col = "red")
lines(hunter.data$Year, hunter.data$doro.manual, col = "green")
#rerun script with different area selection
lines(hunter.data$Year, hunter.data$NHunters, col = "blue") #-BB
#rerun script with different area selection
lines(hunter.data$Year, hunter.data$NHunters, col = "cyan") #-Tana
#rerun script with different area selection
lines(hunter.data$Year, hunter.data$NHunters, col = "orange") #-BB - Tana

legend(2004,80,c("doro.manual","hunting.files", "carc.all", "carc-BB", "carc-Tana", "carc-BB+Tana"),  lwd=1, col=c("green","red", "black", "blue", "cyan", "orange"))

#lessons: if hunting pressure in BB has increased over time, this measure does not capture hunting pressure at all because succesfull hunters per year = same over years in BB
#