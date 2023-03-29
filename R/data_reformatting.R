rm(list=ls())
options(max.print = 1000000)

library(reshape2)
library(ggplot2)
library(sf)
library(tmap)
'%notin%' <- Negate('%in%')

#====== OPTIONS ===========
#number of age classes:
n_ageclasses <- 5 
max_age <- (n_ageclasses-1)
#removal of summer months:
summer_removal <- c(6,7,8,9) #numerical months to be removed from age at harvest data
# choosing varanger sub area ("Inner" / "BB" / "Tana)
area_selection<- c("Inner", "BB",  "Tana")
# start and end of placental scars and embryo sample periods (julian day)
plac_start <- 140 #including
plac_end   <- 80  #until, not including
embr_start <- 100 #including
embr_end   <- 140 #until, not including

#========= LOAD DATA ==============
## set directories
in.dir        <-"C:\\Users\\sho189\\OneDrive - UiT Office 365\\PhD\\RedfoxIPM\\Data from google disk"
shapefile.dir <-"C:\\Users\\sho189\\OneDrive - UiT Office 365\\PhD\\RedfoxIPM\\Fox areas shapefile\\tana rest"

## dataset names
dataset<-"V_redfox_carcass_examination"
shapefile_name <-"tana-rest.shp"

## load data
filenames<-dir(paste(in.dir, "carcass_examination", sep = "/"))
mylist<-c()
for (i in 1:length(filenames)){
  mylist[[i]]<-read.table(paste(in.dir,"carcass_examination", filenames[i], sep = "/"), header=T, sep = ";")
}
myfile<-do.call(rbind, mylist)  # combine all files
allf <- myfile
rm(myfile)

shapefile <- st_read(paste(shapefile.dir,shapefile_name, sep = "/"))
shapefile <- st_transform(shapefile, crs = 4326)

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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

