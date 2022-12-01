##########################################
#### RED FOX DATA PREPARATION FOR IPM ####
##########################################

#*********************#
# Age-at-Harvest data #
#*********************#

## Load winter Age-at-Harvest data
wAaH.tot <- read.table('Data_fromDoro_200927/Cvar.tot_oct-mai_5age.txt', header = T)
wAaH.south <- read.table('Data_fromDoro_200927/Cvar.south_oct-mai_5age.txt', header = T)

## Extract Age-at-Harvest matrices (C)
wC.tot <- rbind(
	as.numeric(wAaH.tot[,'age0']), 
	as.numeric(wAaH.tot[,'age1']), 
	as.numeric(wAaH.tot[,'age2']),
	as.numeric(wAaH.tot[,'age3']),
	as.numeric(wAaH.tot[,'age4'])
	)

wC.south <- rbind(
	as.numeric(wAaH.south[,'age0']), 
	as.numeric(wAaH.south[,'age1']), 
	as.numeric(wAaH.south[,'age2']),
	as.numeric(wAaH.south[,'age3']),
	as.numeric(wAaH.south[,'age4'])
	)

## Extract proportions aged (pData)
pData.tot <- as.numeric(wAaH.tot[,'pData'])
pData.south <- as.numeric(wAaH.south[,'pData'])

## Collect information on year-indeces
YearInfo <- data.frame(index = 1:ncol(wC.tot), WinterHarvestSeason = paste0('Oct ', wAaH.tot[,'year'], ' - May ', as.numeric(wAaH.tot[,'year'])+1))


#*******************#
# Reproduction data #
#*******************#

## Load placental scar data (counts)
P1.tot <- read.table('Data_fromDoro_200418/P1var_tot.txt', header = T)
P1.south <- read.table('Data_fromDoro_200418/P1var_south.txt', header = T)

## Load placental scar / pregnancy data (presence/absence)
P2.tot <- read.table('Data_fromDoro_200418/P2var_tot.txt', header = T)
P2.south <- read.table('Data_fromDoro_200418/P2var_south.txt', header = T)

## Adjust age indeces
P_age_adjust <- function(data){
	
	data$age_adj <- NA
	for(x in 1:nrow(data)){
		
		if(data$type[x] == 'pl.scars'){
			data$age_adj[x] <- ifelse(data[x,2] < 4, data[x,2], 4) + 1
		}
		
		if(data$type[x] %in% c('embryos', 'pregnant', 'combined')){
			data$age_adj[x] <- ifelse(data[x,2]+1 < 4, data[x,2]+1, 4) + 1
		}
	}
	
	return(data)
}

P1.tot <- P_age_adjust(P1.tot)
P1.south <- P_age_adjust(P1.south)
P2.tot <- P_age_adjust(P2.tot)
P2.south <- P_age_adjust(P2.south)

## Check age distributions
table(P1.tot$age_adj)
table(P1.south$age_adj)
table(P2.tot$age_adj)
table(P2.south$age_adj)

## Remove age 1 individual with placental scars
#  (This should not be possible)
P1.tot <- P1.tot[-which(P1.tot$age_adj == 1),]
P1.south <- P1.south[-which(P1.south$age_adj == 1),]

## Add year index
P1.tot$RepYearIndex <- P1.tot$repryear - 2003
P1.south$RepYearIndex <- P1.south$repryear - 2003
P2.tot$RepYearIndex <- P2.tot$repryear - 2003
P2.south$RepYearIndex <- P2.south$repryear - 2003

## Add information to year-indeces
YearInfo$ReproductionYear <- YearInfo$index + 2003


#********************#
# Environmental data #
#********************#

## Load rodent abundance data
RodentData <- read.table('Data_fromDoro_200927/stor_intensiv_04_20-year-var.txt', header = T)

## Format rodent abundance data
RodentAbundance <- (RodentData$tot-mean(RodentData$tot))/sd(RodentData$tot)
RodentIndex2 <- RodentData$cat2
RodentIndex3 <- RodentData$cat3
RodentIndex3[15:16] <- c(2,2)
# TODO: Double-check both indeces (+adjustments) with Doro

## Make covariate for number of successful hunters
#  (Data from Doro's email on 27/10/20)
NHunters <- c(NA, 41, 39, 44, 38, 38, 42, 51, 72, 35, 44, 47, 33, 45, 43, 44)
NHunters <- (NHunters - mean(NHunters, na.rm = T))/sd(NHunters, na.rm = T)
# TODO: Double-check with Doro about NA for 2004 / year assignment


#******************#
# Data arrangement #
#******************#

## Arrange data into lists
IPM.data <- list(C = wC.tot,
				 A = dim(wC.tot)[1],
				 Tmax = dim(wC.tot)[2],
				 pData = pData.tot,
				 
				 P1 = P1.tot$P1,
				 P1_age = P1.tot$age_adj,
				 P1_year = P1.tot$RepYearIndex,
				 X1 = nrow(P1.tot),
				 
				 P2 = P2.tot$P2,
				 P2_age = P2.tot$age_adj,
				 P2_year = P2.tot$RepYearIndex,
				 X2 = nrow(P2.tot),
				 
				 RodentAbundance = RodentAbundance,
				 RodentIndex2 = RodentIndex2,
				 RodentIndex3 = RodentIndex3,
				 NHunters = NHunters,
				 
				 YearInfo = YearInfo)

IPM.data.south <- list(C = wC.south,
				 	   A = dim(wC.south)[1],
				 	   Tmax = dim(wC.south)[2],
				 	   pData = pData.south,
				 
				 	   P1 = P1.south$P1,
				 	   P1_age = P1.south$age_adj,
				 	   P1_year = P1.south$RepYearIndex,
				 	   X1 = nrow(P1.south),
				 
				 	   P2 = P2.south$P2,
				 	   P2_age = P2.south$age_adj,
				 	   P2_year = P2.south$RepYearIndex,
				 	   X2 = nrow(P2.south),
				 	   
				 	   RodentAbundance = RodentAbundance,
				 	   RodentIndex2 = RodentIndex2,
				 	   RodentIndex3 = RodentIndex3,
				 	   NHunters = NHunters,
				 	   
				 	   YearInfo = YearInfo)

## Save workspace containing data lists
rm(list=setdiff(ls(), c('IPM.data', 'IPM.data.south')))

save.image('210426_RedFoxData_IPM.RData')

