library(reshape2)
library(ggplot2)

## This is a script I made to assemble the data Chloe needed to start with the integrated population model



#=============== AGE AT HARVEST MATRIX BUILDING ==============================

rm(list=ls())
options(max.print = 1000000)

##*****************************************************************************
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
head(allf, 3)

allf$t_hunting_date <- as.Date(allf$t_hunting_date)


#  --------need hunting year as an actual number here --------

allf$t_hunting_year
allf$num_hunting_year<-substr(allf$t_hunting_year, start = 1, stop = 4)
allf$num_hunting_year<- as.numeric(allf$num_hunting_year)


#--------Select Varanger and female only----------
fvar1 <- allf[is.na(allf$sn_region)== F & allf$sn_region == "varanger" & is.na(allf$v_sex)== F & allf$v_sex == "female", ]
dim(fvar1) # 1598

#-----Remove the ones that were not shot? For the file with age at harvest-----
fvar1 <- fvar1[fvar1$v_hunter_id != "road_kill", ]
dim(fvar1) # 1571 

## Subsets
# Make two groups of foxes
# varanger, but not the northern parts of the peninsula. In the absence of ArcGIS, cut at 70.57
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#----------- Add option here, select south or total--------

#(now south)

# ------Remove the ones from B?tsfjord and Berlev?g, ie north of 70.58------
# as one of the variants of the file
fvar1[is.na(fvar1$n_dd)==F & fvar1$n_dd > 70.58, ]
dim(fvar1[is.na(fvar1$n_dd)==F & fvar1$n_dd > 70.58, ]) # 158
fvar1 <- fvar1[is.na(fvar1$n_dd)==F & fvar1$n_dd < 70.58, ]
dim(fvar1) # 1415

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ----remove hunting in summer months -------
fvar1$mnd <- as.numeric(format(fvar1$t_hunting_date, "%m"))

table(fvar1$mnd, useNA="ifany") 
fvar1[is.na(fvar1$mnd), ] # Number of foxes that dont have a month assignmed to when they are shot
fvar1[is.na(fvar1$mnd)==F & fvar1$mnd %in% c(5,6), ] # This checks which foxes where shot in june or may
#fvar1 <- fvar1[fvar1$mnd != 6, ] #This removes foxes shot in june
table(fvar1$t_hunting_year[is.na(fvar1$mnd)==F & fvar1$mnd %in% c(7:9)]) # yearly variation in foxes shot july-sept
fvar1[is.na(fvar1$mnd)==F & fvar1$mnd %in% c(7:9), ]  # This shows foxes shot july-sept
fvar1 <- fvar1[fvar1$mnd != 7, ] #This removes foxes shot in july
fvar1 <- fvar1[fvar1$mnd != 8, ] #This removes foxes shot in august
fvar1 <- fvar1[fvar1$mnd != 9, ] #This removes foxes shot in september 

# So use only from October to May for the age distribution
dim(fvar1) # 1378 oct may - 1263 no north

# ----------Remove foxes with lacking age and reproduction info------------
fvar <- fvar1[is.na(fvar1$v_age)==F, ]  #This removes foxes with no age info

# Why here to fvar from fvar1?

dim(fvar) #xx - 681 - xx no north
table(fvar$v_breeding, fvar$v_placental_scars, useNA="ifany") # Breeding is on y axis, plac scar x axis: Shows a table of foxes that have NA in both placental scars and embryo
fvar <- fvar[is.na(fvar$v_breeding)==F, ] #this removes foxes that were not checked for reproductive info (breeding = NA)
dim(fvar) #xxx - 667 - xxx no north
fvar[is.na(fvar$v_breeding)==F & fvar$v_breeding==0 & is.na(fvar$v_placental_scars), ]
# remove one from oct which lacks info about placental scars. WHY? to make it fair with the other ones that were removed? i.e. this one wasnt checked either, but 0 for reproductive because it was october?
# the other 3 are from march, so OK. 
fvar <- fvar[fvar$v_individual_id != "2013_117", ] # make rule here about which ones should be removed, from which months when not checked for breeding
dim(fvar) # 666 - xxx no north


# --------add age class selection process here------


# use 3 age classes
fvar$alder4 <- fvar$v_age
fvar$alder4[fvar$alder4 > 2] <- 2

# use 5 age classes
fvar$alder4 <- fvar$v_age
fvar$alder4[fvar$alder4 > 4] <- 4

#------Building age at harvest matrix---------
table(fvar$alder4)
var.F.C <- as.matrix(table(fvar$t_hunting_year, fvar$alder4))
varFC <- as.data.frame (var.F.C)
varFC2 <- dcast(varFC, Var1~Var2)

#-------Here insert some other way if different amount of age classes-----

names(varFC2) <- c("year","age0", "age1", "age2")

#with 5 classes
names(varFC2) <- c("year","age0", "age1", "age2", "age3", "age4")

fvar1.ann <- table(fvar1$t_hunting_year)
fvar.ann <- table(fvar$t_hunting_year)
prop <- fvar.ann/fvar1.ann 
varFC2$pData <- prop

#------Add option here again for south or total----------

#write.table(varFC2, file="Cvar.tot_oct-mai_5age.txt", sep="\t", quote = F, row.names = F)
#write.table(varFC2, file="Cvar.south_oct-mai_5age.txt", sep="\t", quote = F, row.names = F)



# ========= REPRODUCTION: NR OF EMBRYO'S / PLACENTAL SCARS =========
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fvar1$julian <- as.POSIXlt(fvar1$t_hunting_date)$yday  

fvar <- fvar1
P1var.pl <- fvar[is.na(fvar$v_nb_placental_scars)==F & fvar$v_nb_placental_scars > 0 & is.na(fvar$v_age)==F & is.na(fvar$t_hunting_year)==F,
              c("v_nb_placental_scars", "v_age", "t_hunting_year", "julian")]

dim(P1var.pl) #186 tot, 168 south
summary(P1var.pl)
P1var.pl$type <- "pl.scars"

P1var.emb <- fvar[is.na(fvar$v_embryos)==F & fvar$v_embryos > 0 & is.na(fvar$v_age)==F & is.na(fvar$t_hunting_year)==F,
                 c("v_embryos", "v_age", "t_hunting_year", "julian")]

dim(P1var.emb) #118 tot, 112 tot
summary(P1var.emb)
P1var.emb$type <- "embryos"

names(P1var.pl)[1:2] <- c("P1", "P1_age")
names(P1var.emb)[1:2] <- c("P1", "P1_age")
P1var.pl$repryear <- P1var.pl$t_hunting_year

P1var.emb$repryear <- P1var.emb$t_hunting_year+1
head(P1var.pl)
head(P1var.emb)


##-------- here use the new numeric hunting year --------------


# ------- Also, doesnt chloe do this +1 -1 for embryo en plac scar in her script later?---------
# ===================== end 02.03.2023 =================


P1var <- rbind(P1var.pl, P1var.emb)
summary(P1var)
dim(P1var) #tot 304, 280 south

#write.table(P1var, file="P1var_tot.txt", sep="\t", quote=F, row.names = F)
#write.table(P1var, file="P1var_south.txt", sep="\t", quote=F, row.names = F)

# ========= REPRODUCTION: PROPORTION REPRODUCING  =========
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fvar <- fvar1

########################################################################################################
#                                                                                                      #
#  from day 1-60 use scars                                                                             #
#  from day 61 to 89 use only pregnant (scars here still reflect the previous year)                    #  
#  from day 90 to 150 use a combined variable, scars here reflect reproduction in the current year     #
#  June is excluded (only two pups in the data)                                                        #
#  from day 180 to end of the year use scars                                                           #
#                                                                                                      #  
########################################################################################################

fvar[is.na(fvar$v_nb_placental_scars)==F & is.na(fvar$v_age)==F & is.na(fvar$t_hunting_year)==F,
     c("v_nb_placental_scars", "v_age", "t_hunting_year", "julian")]


#======= DECIDING REPRODUCTION CUT OFF ================
# earliest reproduction
aggregate(fvar$julian[is.na(fvar$v_breeding)==F & fvar$v_breeding==1], by=list(fvar$t_hunting_year[is.na(fvar$v_breeding)==F & fvar$v_breeding==1]), min)
aggregate(fvar$julian[is.na(fvar$v_breeding)==F & fvar$v_breeding==1], by=list(fvar$t_hunting_year[is.na(fvar$v_breeding)==F & fvar$v_breeding==1]), max)

# the earliest pregant is on day 43
# From day 42, we consider that foxes can be pregnant, and thus absence of placental scars cannot be documented.
# what if we take only the ones that were aged?

# ----- Well here I would maybe take the full sample, and not just the ones that were aged? for deciding cut off ---------

aggregate(fvar$julian[is.na(fvar$v_age) ==F & is.na(fvar$v_breeding)==F & fvar$v_breeding==1], 
          by=list(fvar$t_hunting_year[is.na(fvar$v_age) ==F& is.na(fvar$v_breeding)==F & fvar$v_breeding==1]), min)
aggregate(fvar$julian[is.na(fvar$v_age) ==F & is.na(fvar$v_breeding)==F & fvar$v_breeding==1], 
          by=list(fvar$t_hunting_year[is.na(fvar$v_age) ==F& is.na(fvar$v_breeding)==F & fvar$v_breeding==1]), max)
# Then the earliest pregnant is on day 54 and the next on day 61
# Use 1 March as cutoff as earlier, which is day 61 in leap years. 

table(fvar$julian[is.na(fvar$v_age) ==F & is.na(fvar$v_breeding)==F & fvar$v_breeding==1], fvar$t_hunting_year[is.na(fvar$v_age) ==F & is.na(fvar$v_breeding)==F & fvar$v_breeding==1])

hist(fvar$julian[is.na(fvar$v_age) ==F & is.na(fvar$v_breeding)==F & fvar$v_breeding==1])
#including unaged foxes
hist(fvar$julian[is.na(fvar$v_breeding)==F & fvar$v_breeding==1])
hist(fvar$julian[is.na(fvar$v_breeding)==F & fvar$v_breeding==1], breaks=20)
hist(fvar$julian[is.na(fvar$v_breeding)==F & fvar$v_breeding==1], breaks=40)
boxplot(fvar$julian[is.na(fvar$v_breeding)==F & fvar$v_breeding==1])

quantile((fvar$julian[is.na(fvar$v_age) ==F & is.na(fvar$v_breeding)==F & fvar$v_breeding==1]), 0.95)
quantile((fvar$julian[is.na(fvar$v_age) ==F & is.na(fvar$v_breeding)==F & fvar$v_breeding==1]), 0.05)


#does this histogram reflect when they are pregant or when they are hunted though...
hist(fvar$julian, breaks = 40)
hist(fvar$julian[is.na(fvar$v_breeding)==F & fvar$v_breeding==1], col='green', add=TRUE, breaks=10)

#maybe bit of both?

#placental scar info as proportion of hunted
hist(fvar$julian, breaks = 40)
hist(fvar$julian[is.na(fvar$v_placental_scars)==F & fvar$v_placental_scars==1], col='green', add=TRUE, breaks=40)

#so from day 60 -70 they start to get pregnant
# and from day 60- 70 plac scar infor starts to drop too


#plot of when reproduction starts in different years
par(mfrow=c(4,1))
for (i in unique(fvar$t_hunting_year[!is.na(fvar$t_hunting_year)])) {
  hist(fvar$julian[fvar$t_hunting_year== i & fvar$julian > 40 & fvar$julian < 150], xlim = c(40, 150), breaks = c(40,50,60,70,80,90,100,110,120,130,140,150))
   hist(fvar$julian[is.na(fvar$v_breeding)==F & fvar$v_breeding==1 & fvar$t_hunting_year== i], col='green', xlim = c(40, 150), breaks = c(40,50,60,70,80,90,100,110,120,130,140,150), add =TRUE)
}
par(mfrow=c(1,1))

#oke so 70-120 days = 95% of pregnancy data = +- 50 days of pregancy duration
#70-120 = good cut-off
#Alternatives: 90-110 (must be overlap) , 80 - 110 (conservative) and 60-130 (all)
#

#what are the sample sizes of different cut offs
# pregnant foxes
aggregate(fvar$julian[is.na(fvar$v_breeding)==F & fvar$v_breeding==1], by=list(fvar$t_hunting_year[is.na(fvar$v_breeding)==F & fvar$v_breeding==1]), length)
length(fvar$julian[is.na(fvar$v_breeding)==F & fvar$v_breeding==1])
length(fvar$julian[is.na(fvar$v_breeding)==F & fvar$v_breeding==1 & fvar$julian > 70 & fvar$julian < 120])
length(fvar$julian[is.na(fvar$v_breeding)==F & fvar$v_breeding==1 & fvar$julian > 40 & fvar$julian < 130])
length(fvar$julian[is.na(fvar$v_breeding)==F & fvar$v_breeding==1 & fvar$julian > 40 & fvar$julian < 50])
length(fvar$julian[is.na(fvar$v_breeding)==F & fvar$v_breeding==1 & fvar$julian > 50 & fvar$julian < 70])
length(fvar$julian[is.na(fvar$v_breeding)==F & fvar$v_breeding==1 & fvar$julian > 120 & fvar$julian < 130])

# foxes with placental scar info
aggregate(fvar$julian[is.na(fvar$v_placental_scars)==F ], by=list(fvar$t_hunting_year[is.na(fvar$v_placental_scars)==F ]), length)
length(fvar$julian[is.na(fvar$v_placental_scars)==F ])
length(fvar$julian[is.na(fvar$v_placental_scars)==F  & fvar$julian > 70 & fvar$julian < 120])
length(fvar$julian[is.na(fvar$v_placental_scars)==F  & fvar$julian > 40 & fvar$julian < 130])
length(fvar$julian[is.na(fvar$v_placental_scars)==F  & fvar$julian > 40 & fvar$julian < 50])
length(fvar$julian[is.na(fvar$v_placental_scars)==F  & fvar$julian > 50 & fvar$julian < 70])
length(fvar$julian[is.na(fvar$v_placental_scars)==F  & fvar$julian > 120 & fvar$julian < 130])
length(fvar$julian[is.na(fvar$v_placental_scars)==F  & fvar$julian > 90 & fvar$julian < 130])
length(fvar$julian[is.na(fvar$v_placental_scars)==F  & fvar$julian > 50 & fvar$julian < 90])

#===== Placental scar % in pregnancy period compared to plac scar % in autumn period ========
#total percentage foxes with placental scar outside breeding period (n = 480) 
(length(fvar$v_placental_scars[is.na(fvar$v_placental_scars)==F & fvar$v_placental_scars==1 & (fvar$julian < 40 | fvar$julian>130) ] ) / 
 length(fvar$v_placental_scars[is.na(fvar$v_placental_scars)==F &                             (fvar$julian < 40 | fvar$julian>130) ] ))
# = .26 

#total percentage foxes with placental scar inside breeding period (n = 1055)
(length(fvar$v_placental_scars[is.na(fvar$v_placental_scars)==F & fvar$v_placental_scars==1 & (fvar$julian > 40 & fvar$julian < 130)] ) /
 length(fvar$v_placental_scars[is.na(fvar$v_placental_scars)==F                             & (fvar$julian > 40 & fvar$julian < 130)] ))
# = .31

#So this is quite equal

# try again with only 1 and above year olds -----------
#total percentage foxes with placental scar outside breeding period (n = 480) (and exclude autumn pups and unknown age (72 left))
(length(fvar$v_placental_scars[is.na(fvar$v_placental_scars)==F & fvar$v_placental_scars==1 & (fvar$julian < 40 | fvar$julian>130) & fvar$v_age> 0 & is.na(fvar$v_age)==F] ) / 
   length(fvar$v_placental_scars[is.na(fvar$v_placental_scars)==F &                             (fvar$julian < 40 | fvar$julian>130) & fvar$v_age> 0 & is.na(fvar$v_age)==F] ))
# = .71 

#total percentage foxes with placental scar inside breeding period (n = 1055) (and exclude autumn pups and unknown age (215 left))
(length(fvar$v_placental_scars[is.na(fvar$v_placental_scars)==F & fvar$v_placental_scars==1 & (fvar$julian > 40 & fvar$julian < 130) & fvar$v_age> 0 & is.na(fvar$v_age)==F] ) /
    length(fvar$v_placental_scars[is.na(fvar$v_placental_scars)==F                          & (fvar$julian > 40 & fvar$julian < 130) & fvar$v_age> 0 & is.na(fvar$v_age)==F] ))
# = .50

# for all foxes: % placental scars within and outside breeding period the same
# For adult foxes: the ones that were checked for placental scars in breeding period were not breeding, and these were less likely to have bred in the past


#===== Nr. of plac scar in pregnancy period compared to nr in autumn period ========

# all foxes ---------
par(mfrow=c(1,2))
boxplot(fvar$v_nb_placental_scars[is.na(fvar$v_nb_placental_scars)==F  & (fvar$julian < 40 | fvar$julian>130) ], ylim =c(0,11) )  
mean(fvar$v_nb_placental_scars[is.na(fvar$v_nb_placental_scars)==F  & (fvar$julian < 40 | fvar$julian>130) ] )  

#total percentage foxes with placental scar inside breeding period (n = 74)
boxplot(fvar$v_nb_placental_scars[is.na(fvar$v_nb_placental_scars)==F  & (fvar$julian > 40 & fvar$julian < 130) ] , ylim =c(0,11))  
mean(fvar$v_nb_placental_scars[is.na(fvar$v_nb_placental_scars)==F  & (fvar$julian > 40 & fvar$julian < 130) ] )  
par(mfrow=c(1,1))

# adults only---------
#nr of scars outside breeding period (n = 29) (and exclude autumn pups)
boxplot(fvar$v_nb_placental_scars[is.na(fvar$v_nb_placental_scars)==F  & (fvar$julian < 40 | fvar$julian>130) & fvar$v_age> 0 & is.na(fvar$v_age)==F] )  
mean(fvar$v_nb_placental_scars[is.na(fvar$v_nb_placental_scars)==F  & (fvar$julian < 40 | fvar$julian>130) & fvar$v_age> 0 & is.na(fvar$v_age)==F] )  

#nr of placental scar inside breeding period (n = 74)
boxplot(fvar$v_nb_placental_scars[is.na(fvar$v_nb_placental_scars)==F  & (fvar$julian > 40 & fvar$julian < 130) & fvar$v_age> 0 & is.na(fvar$v_age)==F] )  
mean(fvar$v_nb_placental_scars[is.na(fvar$v_nb_placental_scars)==F  & (fvar$julian > 40 & fvar$julian < 130) & fvar$v_age> 0 & is.na(fvar$v_age)==F] )  

#the number of scars are similar


# what about for each hunting year? (now set threshold lower (maybe 100) and a bit later (140) so you dont get the "new" placental scars) --------------
try<-data.frame()
for (i in unique(fvar$t_hunting_year[!is.na(fvar$t_hunting_year)])) {
  
data.temp <- data.frame(
  Year = i,
  propPC_out= (length(fvar$v_placental_scars[is.na(fvar$v_placental_scars)==F & fvar$v_placental_scars==1 & (fvar$julian < 40 | fvar$julian>140) & is.na(fvar$t_hunting_year)==F & fvar$t_hunting_year == i] ) / 
                          length(fvar$v_placental_scars[is.na(fvar$v_placental_scars)==F &                         (fvar$julian < 40 | fvar$julian>140) & is.na(fvar$t_hunting_year)==F & fvar$t_hunting_year == i] ))
  ,
  propPC_in = (length(fvar$v_placental_scars[is.na(fvar$v_placental_scars)==F & fvar$v_placental_scars==1 & (fvar$julian > 40 & fvar$julian <100) & is.na(fvar$t_hunting_year)==F & fvar$t_hunting_year == i] ) / 
                         length(fvar$v_placental_scars[is.na(fvar$v_placental_scars)==F &                         (fvar$julian > 40 & fvar$julian <100) & is.na(fvar$t_hunting_year)==F & fvar$t_hunting_year == i] ))
  
)

data.temp$Diff <- data.temp$propPC_in -data.temp$propPC_out

try <- rbind(try, data.temp)
}
(try)
#probably not enough data to check per year
mean(try$propPC_in)
mean(try$propPC_out) # = similar proportions as previous all year together calculation


# ======== comparing embryos to placental scar next year ==============
# previous histograms with green: from day 60 -70 they start to get pregnant
# and from day 60- 70 plac scar infor starts to drop too
# day 43 earliest, so day 80-90 is earliest they can have plac scar from this year
# but most after 60-70 so we can safely take 90 or even 100
#BUT: adult %plac scar within breeding < than outside because see previous, so include as little as possible within breeding. Start at 70 or 80

# choosing thresholds 
try<-data.frame()
for (i in unique(fvar$num_hunting_year[!is.na(fvar$num_hunting_year)])) {
  
  data.temp <- data.frame(
    Year = i,
    propPC_outnext_plac= (length(fvar$v_placental_scars[is.na(fvar$v_placental_scars)==F & fvar$v_placental_scars==1 & (fvar$julian < 80 | fvar$julian>140) & fvar$v_age> 0 & is.na(fvar$v_age)==F & is.na(fvar$num_hunting_year)==F & fvar$num_hunting_year == i+1] ) / 
                   length(fvar$v_placental_scars[is.na(fvar$v_placental_scars)==F &                                    (fvar$julian < 80 | fvar$julian>140) & fvar$v_age> 0 & is.na(fvar$v_age)==F & is.na(fvar$num_hunting_year)==F & fvar$num_hunting_year == i+1] ))
    ,
    
    n_plac =       (length(fvar$v_placental_scars[is.na(fvar$v_placental_scars)==F &                                    (fvar$julian < 80 | fvar$julian>140) & fvar$v_age> 0 & is.na(fvar$v_age)==F & is.na(fvar$num_hunting_year)==F & fvar$num_hunting_year == i+1] ))
    ,
    propPC_in_embryo = (length(fvar$v_breeding[is.na(fvar$v_breeding)==F & fvar$v_breeding==1 & (fvar$julian > 100 & fvar$julian <130) & is.na(fvar$num_hunting_year)==F & fvar$num_hunting_year == i] ) / 
                   length(fvar$v_breeding[is.na(fvar$v_breeding)==F &                         (fvar$julian > 100 & fvar$julian <130) & is.na(fvar$num_hunting_year)==F & fvar$num_hunting_year == i] ))
    
    ,
    n_embr =       (length(fvar$v_breeding[is.na(fvar$v_breeding)==F &                         (fvar$julian > 100 & fvar$julian <130) & is.na(fvar$num_hunting_year)==F & fvar$num_hunting_year == i] ))
  )
  
  data.temp$Diff <- data.temp$propPC_in_embryo -data.temp$propPC_outnext_plac
  data.temp$sample_size <- data.temp$n_plac  + data.temp$n_embr        
  
  try <- rbind(try, data.temp)
}
(try)
mean(try$propPC_outnext_plac, na.rm = TRUE ) #.67 = approximately correct number
mean(try$propPC_in_embryo   , na.rm = TRUE ) 

#threshold:     40-130    70-120      80-110        90 - 100      100-130   70-130  90-130
#% pregnant:    24%       35%           40%          35%            64%       36%     49%

#so 100-130 best but same size is quite low

#Conclusions: for % placental scars you can use all periods, but % during breeding period is less, likely because you then have a subsample of foxes that are less likely to breed
#             for % with embryo's its best to use a later period, about 100-130


#show in histogram form
#does this histogram reflect when they are pregant or when they are hunted though...
hist(fvar$julian, breaks = 40)
hist(fvar$julian[is.na(fvar$v_breeding)==F], col = "orange", breaks = 40, add = TRUE)
hist(fvar$julian[is.na(fvar$v_breeding)==F & fvar$v_breeding==1], col='green', add=TRUE, breaks=10)

#placental scar info as proportion of hunted
hist(fvar$julian, breaks = 40)
hist(fvar$julian[is.na(fvar$v_placental_scars)==F], col = "orange", breaks = 40, add=TRUE)
hist(fvar$julian[is.na(fvar$v_placental_scars)==F & fvar$v_placental_scars==1], col='green', add=TRUE, breaks=40)

#adults only
hist(fvar$julian[is.na(fvar$v_age)==F & fvar$v_age> 0] , breaks = 40)
hist(fvar$julian[is.na(fvar$v_placental_scars)==F & fvar$v_age> 0 & is.na(fvar$v_age)==F ], col = "orange", breaks = 40, add=TRUE)
hist(fvar$julian[is.na(fvar$v_placental_scars)==F & fvar$v_placental_scars==1& fvar$v_age> 0 & is.na(fvar$v_age)==F ], col='green', add=TRUE, breaks=40)

#AHHHH MAYBE JUST THE PLACENTAL SCARS START TO DISAPPEAR!
#THATS WHY LESS IN BREEDING, NOT BECAUSE ITS A SUBSAMPLE OF NON BREEDERS


#===========HIER WAS IK =========================

# what about nr of placental scar of nr of embryo's? some periods best?

#nr of scars kinda equal?

#===== Nr. of plac scar in pregnancy period compared to nr in autumn period ========

# adults only---------
#nr of scars outside breeding period (n = 29) (and exclude autumn pups)
boxplot(fvar$v_nb_placental_scars[is.na(fvar$v_nb_placental_scars)==F  & (fvar$julian < 60 | fvar$julian>130) & fvar$v_age> 0 & is.na(fvar$v_age)==F] )  
mean(fvar$v_nb_placental_scars[is.na(fvar$v_nb_placental_scars)==F  & (fvar$julian < 60 | fvar$julian>130) & fvar$v_age> 0 & is.na(fvar$v_age)==F] )  

#nr of placental scar inside breeding period (n = 74)
boxplot(fvar$v_nb_placental_scars[is.na(fvar$v_nb_placental_scars)==F  & (fvar$julian > 60 & fvar$julian < 130) & fvar$v_age> 0 & is.na(fvar$v_age)==F] )  
mean(fvar$v_nb_placental_scars[is.na(fvar$v_nb_placental_scars)==F  & (fvar$julian > 60 & fvar$julian < 130) & fvar$v_age> 0 & is.na(fvar$v_age)==F] )  

#only for those with placental scars
#nr of scars outside breeding period (n = 29) (and exclude autumn pups)
boxplot(fvar$v_nb_placental_scars[fvar$v_placental_scars> 0 & is.na(fvar$v_placental_scars)==F & (fvar$julian < 60 | fvar$julian>130) & fvar$v_age> 0 & is.na(fvar$v_age)==F])  
mean   (fvar$v_nb_placental_scars[fvar$v_placental_scars> 0 & is.na(fvar$v_placental_scars)==F & (fvar$julian < 60 | fvar$julian>130) & fvar$v_age> 0 & is.na(fvar$v_age)==F])   

#nr of placental scar inside breeding period (n = 74)
boxplot(fvar$v_nb_placental_scars[fvar$v_placental_scars> 0 & is.na(fvar$v_placental_scars)==F & (fvar$julian > 60 & fvar$julian<130) & fvar$v_age> 0 & is.na(fvar$v_age)==F])  
mean   (fvar$v_nb_placental_scars[fvar$v_placental_scars> 0 & is.na(fvar$v_placental_scars)==F & (fvar$julian > 60 & fvar$julian<130) & fvar$v_age> 0 & is.na(fvar$v_age)==F])   


# ========= Nr. of embryo =================

boxplot(fvar$v_embryos[fvar$v_embryos> 0 & is.na(fvar$v_embryos)==F & (fvar$julian > 60 | fvar$julian<130) & fvar$v_age> 0 & is.na(fvar$v_age)==F], ylim = c(0,10))  
mean   (fvar$v_embryos[fvar$v_embryos> 0 & is.na(fvar$v_embryos)==F & (fvar$julian > 60 | fvar$julian<130) & fvar$v_age> 0 & is.na(fvar$v_age)==F])   

#adult only

boxplot(fvar$v_embryos[fvar$v_embryos> 0 & is.na(fvar$v_embryos)==F & (fvar$julian > 60 | fvar$julian<130) & is.na(fvar$v_age)==F], ylim = c(0,10))  
mean   (fvar$v_embryos[fvar$v_embryos> 0 & is.na(fvar$v_embryos)==F & (fvar$julian > 60 | fvar$julian<130) & is.na(fvar$v_age)==F])   



boxplot(fvar$v_embryos[is.na(fvar$v_embryos)==F  & (fvar$julian < 60 | fvar$julian>130) & fvar$v_age> 0 & is.na(fvar$v_age)==F], ylim= c(0,10) )  
mean(fvar$v_embryos[is.na(fvar$v_embryos)==F  & (fvar$julian < 60 | fvar$julian>130) & fvar$v_age> 0 & is.na(fvar$v_age)==F] )  






#doro below
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





# ============= WHAT ABOUT IESJAVRI? ==================



# look at distribution by month
fies1 <- allf[is.na(allf$reg)== F & allf$reg == "ies" & is.na(allf$sex)== F & allf$sex == "female", ]
dim(fies1) # 366
table(fies1$mnd) # remove July, these are pups
fies1[is.na(fies1$mnd)==F & fies1$mnd==7, ]
fies1 <- fies1[is.na(fies1$mnd)==F & fies1$mnd < 7, ] 
dim(fies1) # 363

# only the ones with age and full repr data
fies <- fies1[is.na(fies1$v_age)==F, ] 
dim(fies)#306
table(fies$v_age) # max 7 years

table(fies$repr, useNA="ifany")
table(fies$placentalscars, useNA="ifany") 
# 5 with missing information about reproduction, removed them as well, ie keep only females with full data
fies[is.na(fies$repr), ]
fies <- fies[is.na(fies$repr)==F, ]
dim(fies) # 301

# use 3 age classes
fies$alder4 <- fies$v_age
fies$alder4[fies$alder4 > 2] <- 2
table(fies$alder4)
ies.F.C <- as.matrix(table(fies$t_hunting_year, fies$alder4))
library(reshape2)
iesFC <- as.data.frame (ies.F.C)
iesFC2 <- dcast(iesFC, Var1~Var2)
names(iesFC2) <- c("year","age0", "age1", "age2")

# annual proportion of aged females
fies1.ann <- table(fies1$t_hunting_year)
fies.ann <- table(fies$t_hunting_year)
prop <- fies.ann/fies1.ann 
iesFC2$pData <- prop

write.table(iesFC2, file="ies_C.txt", sep="\t", quote = F)

table(is.na(fies$nbplacentalscars))

table(fies$repr, useNA="ifany")
table(fies$placentalscars, useNA="ifany")
fies[fies$repr==0 & is.na(fies$placentalscars), ]
firstrepr <- aggregate(fies$julian[fies$repr == 1], by=list(fies$year[fies$repr == 1]), min)

table(fies$repr, fies$placentalscars, useNA="ifany")
table(fies$repr, fies$nbplacentalscars, useNA="ifany") 
table(fies$placentalscars, fies$nbplacentalscars, useNA="ifany")


# how to best solve the problem of overlapping time
# check for each year when the first pregnant female was observed. 
# All those which were shot before that are considered based on placental scars, 
# and all those which were shot after that are considered based on pregnancy

# observed number of placental scars and embryos
iesP1 <- fies[is.na(fies$nbplacentalscars)==F & fies$nbplacentalscars > 0, 
              c("nbplacentalscars", "v_age", "t_hunting_year", "julian")]

fies[is.na(fies$nbplacentalscars)==F & fies$nbplacentalscars > 0 & fies$v_age==0, 
     c("nr", "nbplacentalscars", "v_age", "t_hunting_year", "julian")]      

