## This is a script I made to assemble the data Chloe needed to start with the integrated population model
library(reshape2)
library(ggplot2)
library(sf)
library(tmap)

rm(list=ls())
options(max.print = 1000000)

#====== dataset guide ==========

# allf  = all foxes
# fvar1 = allf [hunting date registered, female only, varanger only, no road kill, area within varanger selection]
# fvar2 = fvar1 [without months 6,7,8,9]      --> TO CALCULATE % CHECKED FOR AGE (fvar3 / fvar2)
# fvar3 = fvar2 [only foxes with age info]    --> AGE AT HARVEST MATRIX
# fvar4 = fvar1 [only foxes with age info]    --> REPRODUCTIVE INFO

#====== OPTIONS ===========
#number of age classes:
n_ageclasses <- 5
max_age <- (n_ageclasses-1)
#removal of summer months:
summer_removal <- c(6,7,8,9)
# ------choosing varanger sub area------------------
# https://www.norgeskart.no/#!?project=norgeskart&layers=1002&zoom=8&lat=7862448.02&lon=998682.39&markerLat=7847177.458430003&markerLon=1062026.4302111468&p=searchOptionsPanel&drawing=W9XDE4cBWPXm_Pf_0Gxo&sok=Austertanaveien
# hist(fvar1$e_dd, breaks = 80)
# can also use utm33 / longlat in combination to get different angled lines
north_max <- 7894000
east_min <- 28.4


#========= LOAD DATA ==============
#First load that data that can be applied to all info needed : Age at harvest matrix, breeding pop %, nr of fetuses.
## Load the files that are formated for the database

## set directories
in.dir<-"C:\\Users\\sho189\\OneDrive - UiT Office 365\\PhD\\RedfoxIPM\\Data from google disk" 
# path to folder where formated files are saved

## choose dataset
dataset<-"V_redfox_carcass_examination"

## load data
filenames<-dir(paste(in.dir, "carcass_examination", sep = "/"))

mylist<-c()
for (i in 1:length(filenames)){
  mylist[[i]]<-read.table(paste(in.dir,"carcass_examination", filenames[i], sep = "/"), header=T, sep = ";")
}

myfile<-do.call(rbind, mylist)  # combine all files

allf <- myfile
rm(myfile)

#----hunting date as date---
allf$t_hunting_date <- as.Date(allf$t_hunting_date)

#  -------- we need hunting year as an actual number sometimes so we can add and substract --------
allf$num_hunting_year<-substr(allf$t_hunting_year, start = 1, stop = 4)
allf$num_hunting_year<- as.numeric(allf$num_hunting_year)

# -------julian day----------
allf$julian <- as.POSIXlt(allf$t_hunting_date)$yday  

#  ------   remove foxes without hunting date?------------
fvar1 <- allf[is.na(allf$t_hunting_date)==F,]

#--------Select Varanger and female only----------
fvar1 <- fvar1[is.na(fvar1$sn_region)== F & fvar1$sn_region == "varanger" & is.na(fvar1$v_sex)== F & fvar1$v_sex == "female", ]
dim(fvar1) # 1721

#-----Remove the ones that were not shot? For the file with age at harvest-----
fvar1 <- fvar1[fvar1$v_hunter_id != "road_kill", ]

# ------choosing varanger sub area------------------
fvar1 <- fvar1[is.na(fvar1$n_utm33)==F & fvar1$n_utm33 < north_max & is.na(fvar1$e_dd)==F & fvar1$e_dd > east_min, ]

# ------ Alternative with shapefile ------------
# ??
# ??
# ??
# ??

# --------age class selection process------
fvar1$alder4 <- fvar1$v_age
fvar1$alder4[fvar1$alder4 > max_age] <- max_age


#=============== AGE AT HARVEST MATRIX BUILDING ==============================
fvar2<-fvar1
rm(fvar1)

# ----remove hunting in summer months -------
fvar2$mnd <- as.numeric(format(fvar2$t_hunting_date, "%m"))
fvar2 <- subset(fvar2,!(fvar2$mnd %in% summer_removal))

# ----------Remove foxes with lacking age info------------
fvar3 <- fvar2[is.na(fvar2$v_age)==F, ]  #This removes foxes with no age info

#------Building age at harvest matrix---------
table(fvar3$alder4)
var.F.C <- as.matrix(table(fvar3$t_hunting_year, fvar3$alder4))
varFC <- as.data.frame (var.F.C)
varFC2 <- dcast(varFC, Var1~Var2)

column_names <- c("year","age0", "age1", "age2", "age3", "age4", "age5", "age6", "age7", "age8", "age9", "age10", "age11", "age12", "age13", "age14", "age15", "age16", "age17", "age18")
names(varFC2) <- head(column_names, n_ageclasses+1)

#percentage that were aged
fvar3.ann <- table(fvar3$t_hunting_year) #where unaged removed
fvar2.ann <- table(fvar2$t_hunting_year) #where unaged not removed
prop <- fvar3.ann/fvar2.ann 
varFC2$pData <- prop

rm(fvar2)
rm(fvar3)


# ========= REPRODUCTION: NR OF EMBRYO'S / PLACENTAL SCARS =========
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#here restart with a new dataset that includes the summer months, they have age info, and they are checked for repr. info
# start again with fvar 1 (summer months were removed in fvar2)
rm(fvar2)
rm(fvar3)
(fvar1)

#remove unaged foxes
# ----------Remove foxes with lacking age info------------
fvar4 <- fvar1[is.na(fvar1$v_age)==F, ]  #This removes foxes with no age info

fvar<-fvar1

#==== Nr of scars / embryos

# numbers are all the same, no need to selecto for period when estimating nr ob scars or embryos


#doro stuff:
P1var.pl <- fvar[is.na(fvar$nbplacentalscars)==F & fvar$nbplacentalscars > 0 & is.na(fvar$alder3)==F & is.na(fvar$foxyear)==F,
                 c("nbplacentalscars", "alder3", "foxyear", "julian")]

dim(P1var.pl) #186 tot, 168 south
summary(P1var.pl)
P1var.pl$type <- "pl.scars"

P1var.emb <- fvar[is.na(fvar$nbembryos)==F & fvar$nbembryos > 0 & is.na(fvar$alder3)==F & is.na(fvar$foxyear)==F,
                  c("nbembryos", "alder3", "foxyear", "julian")]

dim(P1var.emb) #118 tot, 112 tot
summary(P1var.emb)
P1var.emb$type <- "embryos"

names(P1var.pl)[1:2] <- c("P1", "P1_age")
names(P1var.emb)[1:2] <- c("P1", "P1_age")
P1var.pl$repryear <- P1var.pl$foxyear
P1var.emb$repryear <- P1var.emb$foxyear+1
head(P1var.pl)
head(P1var.emb)

P1var <- rbind(P1var.pl, P1var.emb)
summary(P1var)
dim(P1var) #tot 304, 280 south

write.table(P1var, file="P1var_tot.txt", sep="\t", quote=F, row.names = F)
write.table(P1var, file="P1var_south.txt", sep="\t", quote=F, row.names = F)










# ---choosing thresholds ------
try<-data.frame()
for (i in unique(fvar$num_hunting_year[!is.na(fvar$num_hunting_year)])) {
  
  data.temp <- data.frame(
    Year = i,
    propPC_outnext_plac= (length(fvar$v_placental_scars[is.na(fvar$v_placental_scars)==F & fvar$v_placental_scars==1 & (fvar$julian < 80 | fvar$julian>140) & fvar$v_age> 0 & is.na(fvar$v_age)==F & is.na(fvar$num_hunting_year)==F & fvar$num_hunting_year == i+1] ) / 
                          length(fvar$v_placental_scars[is.na(fvar$v_placental_scars)==F &                             (fvar$julian < 80 | fvar$julian>140) & fvar$v_age> 0 & is.na(fvar$v_age)==F & is.na(fvar$num_hunting_year)==F & fvar$num_hunting_year == i+1] ))
    ,
    
    n_plac =       (length(fvar$v_placental_scars[is.na(fvar$v_placental_scars)==F &                                   (fvar$julian < 80 | fvar$julian>140) & fvar$v_age> 0 & is.na(fvar$v_age)==F & is.na(fvar$num_hunting_year)==F & fvar$num_hunting_year == i+1] ))
    ,
    propPC_in_embryo = (length(fvar$v_breeding[is.na(fvar$v_breeding)==F & fvar$v_breeding==1 & (fvar$julian > 100 & fvar$julian <140) & is.na(fvar$num_hunting_year)==F & fvar$num_hunting_year == i] ) / 
                        length(fvar$v_breeding[is.na(fvar$v_breeding)==F &                      (fvar$julian > 100 & fvar$julian <140) & is.na(fvar$num_hunting_year)==F & fvar$num_hunting_year == i] ))
    
    ,
    n_embr =       (length(fvar$v_breeding[is.na(fvar$v_breeding)==F &                          (fvar$julian > 100 & fvar$julian <140) & is.na(fvar$num_hunting_year)==F & fvar$num_hunting_year == i] ))
  )
  
  data.temp$Diff <- data.temp$propPC_in_embryo -data.temp$propPC_outnext_plac
  data.temp$sample_size <- data.temp$n_plac  + data.temp$n_embr        
  
  try <- rbind(try, data.temp)
}
(try)
mean(try$propPC_outnext_plac, na.rm = TRUE ) #.67 = approximately correct number
mean(try$propPC_in_embryo   , na.rm = TRUE ) 



#so 100-140 best but same size is quite low

#Conclusions: for % placental scars you can use all periods, but % during breeding period is less, likely because you then have a subsample of foxes that are less likely to breed
#             for % with embryo's its best to use a later period, about 100-130


# numbers are all the same, no need to selecto for period when estimating nr ob scars or embryos







#===========doro below =====================
P2var.pl1 <- fvar[is.na(fvar$v_nb_placental_scars)==F & is.na(fvar$v_age)==F & is.na(fvar$t_hunting_year)==F & is.na(fvar$julian)==F & # 1 jan - day 60
                    fvar$julian & fvar$julian < 61, c("v_nb_placental_scars", "v_age", "t_hunting_year", "julian")]
P2var.pl2 <- fvar[is.na(fvar$v_nb_placental_scars)==F & is.na(fvar$v_age)==F & is.na(fvar$t_hunting_year)==F & is.na(fvar$julian)==F & # from sommer to 31.12
                    fvar$julian & fvar$julian > 179, c("v_nb_placental_scars", "v_age", "t_hunting_year", "julian")]
dim(P2var.pl1) # tot 208; south 198
dim(P2var.pl2) # tot 212; south 165
P2var.pl <- rbind(P2var.pl1, P2var.pl2)
head(P2var.pl)

P2var.pl$type <- "pl.scars"
P2var.pl$repryear <- P2var.pl$t_hunting_year
names(P2var.pl)[1:2] <- c("P2", "P2_age")
dim(P2var.pl) # tot 420; south 363

# starting day 90, we assume that they gave birth in the current year
# but they can also have been pregnant 

# between 61 and 90 use pregnant. placentalscars reflect last years reproduction and repr this years reproduction. 
# But from day 61 we cannot detect last years reproduction any more
# The pregnant fox on day 54 will be added to the pregnant ones in the current year

# pregnant foxes
# between day 60 and day 90 only pregnant is used
P2var.preg <- fvar[is.na(fvar$v_age)==F & is.na(fvar$t_hunting_year)==F & is.na(fvar$julian)==F & is.na(fvar$repr)==F &  # between day 60 and day 89 only pregnant is used
            fvar$julian > 60  & fvar$julian < 90, c("repr", "v_age", "t_hunting_year", "julian")]
P2var.preg.add <- fvar[is.na(fvar$v_age)==F & is.na(fvar$t_hunting_year)==F & is.na(fvar$julian)==F & is.na(fvar$repr)==F & fvar$repr==1 &  # pregnant fox on day 54
                         fvar$julian ==54, c("repr", "v_age", "t_hunting_year", "julian")]
P2var.preg <- rbind(P2var.preg, P2var.preg.add)

P2var.preg$type <- "pregnant"
P2var.preg$repryear <- P2var.preg$t_hunting_year + 1
names(P2var.preg)[1:2] <- c("P2", "P2_age")
head(P2var.preg)
dim(P2var.preg) # tot 192, south 181

#From day 90 to 150 
X <- fvar[is.na(fvar$v_age)==F & is.na(fvar$t_hunting_year)==F & is.na(fvar$julian)==F &  
       fvar$julian <179  & fvar$julian > 89 & fvar$julian & fvar$julian < 150, c("placentalscars", "repr", "v_age", "t_hunting_year", "julian", "Date","comment")]
#X[order(X$julian), ]
# Here we need a combined variable between repr and placental scars.
X$reprcomb <- X$repr
X$reprcomb[is.na(X$repr)==F & X$repr == 0] <- X$placentalscars[is.na(X$repr)==F & X$repr == 0]
#X[order(X$julian), c("placentalscars", "repr", "reprcomb", "v_age", "t_hunting_year", "julian", "Date","comment")]
# this is OK

P2var.comb <- X[is.na(X$reprcomb)==F,c("reprcomb", "v_age", "t_hunting_year", "julian")]
head(P2var.comb)

P2var.comb$type <- "combined"
P2var.comb$repryear <- P2var.comb$t_hunting_year + 1
names(P2var.comb)[1:2] <- c("P2", "P2_age")
head(P2var.comb)
dim(P2var.comb) # tot 161; south 156

head(P2var.pl)
head(P2var.preg)
head(P2var.comb)

P2var.tot <- rbind(P2var.pl, P2var.preg, P2var.comb)

write.table(P2var.tot, file="P2var_tot.txt", row.names=F, sep="\t", quote=F)
write.table(P2var.tot, file="P2var_south.txt", row.names=F, sep="\t", quote=F)




