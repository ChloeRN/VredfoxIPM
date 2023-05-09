#' Makes a yearly rodent abundance variable, mean nr of rodents in a plot, for each year
#'
#' @param rodent.dataset dataframe containing the carcass dataset downloaded from the COAT dataportal
#'
#' @return a dataframe containing the mean nr of rodents (in a plot), for each year
#' @export
#'
#' @examples


reformatData_rodent <- function(rodent.dataset) {

#========= LOAD DATA ==============
allrod <- rodent.dataset
#================ Select rodents from storskala areas ==============================

storskala<- reshape2::dcast(allrod, sn_locality + sn_site + t_year + t_season ~ v_species, value.var = "v_abundance", fun.aggregate = sum) 

#remove row where season is NA
storskala <- storskala[!is.na(storskala$t_season) ,]

# remove aves sorex and ranidae and neo_fod (keep "cricetidae","lem_lem", "mic_oec", "myo_ruf", "myo_rut")
storskala<- subset(storskala, select=-c(aves,neo_fod, ranidae, sor_ara, sor_cae, sor_min, sor_sp))
#rename some columns
names(storskala)[1:4] <- c("region", "plot", "year", "season")
names(storskala)[names(storskala) == 'cricetidae'] <- 'rodsp'
names(storskala)[names(storskala) == 'mic_oec'] <- 'Moec'
names(storskala)[names(storskala) == 'myo_ruf'] <- 'Mruf'
names(storskala)[names(storskala) == 'myo_rut'] <- 'Mrut'
names(storskala)[names(storskala) == 'lem_lem'] <- 'Llem'

#================================================================
storskala$session <- as.factor(paste(storskala$year, storskala$season, sep=""))

# aggregate by region, plot, session, season, year
stor <- aggregate(storskala$Llem, by=list(storskala$region, storskala$plot, storskala$session, storskala$season, storskala$year), sum, na.rm=T) 
sum(storskala$Llem, na.rm=T)
sum(stor$x) # OK except for one with lacking plot number
names(stor) <- c("reg", "plot", "session", "season", "year", "Llem") 
stor$Moec <- aggregate(storskala$Moec, by=list(storskala$region, storskala$plot, storskala$session, storskala$season, storskala$year), sum, na.rm=T)$x 
stor$Mruf <- aggregate(storskala$Mruf, by=list(storskala$region, storskala$plot, storskala$session, storskala$season, storskala$year), sum, na.rm=T)$x 
stor$Mrut <- aggregate(storskala$Mrut, by=list(storskala$region, storskala$plot, storskala$session, storskala$season, storskala$year), sum, na.rm=T)$x 
stor$rodsp <- aggregate(storskala$rodsp, by=list(storskala$region, storskala$plot, storskala$session, storskala$season, storskala$year), sum, na.rm=T)$x
stor$vole <- stor$Mruf + stor$Moec + stor$Mrut + stor$rodsp
stor$tot <- stor$vole + stor$Llem

# Aggregate annual means into three regions; Varanger; Nordkyn; Ifjord-gaissene
stor$foxreg <- NA
stor$foxreg[stor$reg %in% c("komagdalen",  "stjernevann", "vestre_jakobselv")] <- "varanger"
stor$foxreg[stor$reg %in% c("nordkynn",  "bekkarfjord")] <- "nordkynn"
stor$foxreg[stor$reg == "ifjordfjellet"] <- "ifjordfjellet"

#check if 2 seasons for each year
nrow(stor[stor$season =="fall",])
nrow(stor[stor$season =="spring",])

check<-aggregate(. ~ season +year, stor, FUN = length)
if( (length(check$season)) %% 2 !=0 ){    #check if even nr of seasons
  stop("uneven number of seasons in rodent dataset")  #this will return an error when using google disc folders because autumn 2022 not included
}

agdat <-  aggregate(cbind(Llem, Moec, Mruf, Mrut, rodsp, vole, tot) ~ year + season + foxreg, stor, mean)  #the mean nr of rodents per plot, for each year and season
agdat<-agdat[agdat$foxreg=="varanger",] # we only use varanger

#--- continuous rodent variables winter and fall ---
agdat$start_hunting_year <- ifelse(agdat$season == "spring", agdat$year-1, agdat$year) # -1 year from spring, so we can add fall to form a hunting year representative rodent variable instead of adding spring and fall from the same year

agdat_fall <- agdat[agdat$season=="fall",]
agdat_fall$st.tot <- (agdat_fall$tot-mean(agdat_fall$tot))/sd(agdat_fall$tot) #vole and lemming actual numbers together
agdat_fall$st.lem <- (agdat_fall$Llem-mean(agdat_fall$Llem))/sd(agdat_fall$Llem)  #standardise lemming
agdat_fall$st.vole <- (agdat_fall$vole-mean(agdat_fall$vole))/sd(agdat_fall$vole) #standardise vole
agdat_fall$st.lemvole<-(agdat_fall$st.lem+ agdat_fall$st.vole) # vole and lemming together after standardised

agdat_spring <- agdat[agdat$season=="spring",]
agdat_spring$st.tot <- (agdat_spring$tot-mean(agdat_spring$tot))/sd(agdat_spring$tot)#vole and lemming actual numbers together
agdat_spring$st.lem <- (agdat_spring$Llem-mean(agdat_spring$Llem))/sd(agdat_spring$Llem) #standardise lemming
agdat_spring$st.vole <- (agdat_spring$vole-mean(agdat_spring$vole))/sd(agdat_spring$vole)#standardise vole
agdat_spring$st.lemvole<-  (agdat_spring$st.lem+ agdat_spring$st.vole) # vole and lemming together after standardised

agdat_winter <- rbind(agdat_fall, agdat_spring)
agdat_winter <- aggregate(cbind(Llem, Moec, Mruf, Mrut, rodsp, vole, tot, st.tot, st.lem, st.vole, st.lemvole) ~ start_hunting_year , agdat_winter, mean)

#remove start_hunting year 2003 (only spring 2004, no rodents from 2003 fall)
agdat_winter[ agdat_winter$start_hunting_year==2003, -1] <-NA

#if spring and fall have dissimilar "start_hunting_year" in the final row, this "start_hunting_year" should also be NA (because spring rodents not collected yet for that year)
if (agdat_spring[nrow(agdat_spring), 11] != agdat_fall[nrow(agdat_fall), 11]) {
 agdat_winter[nrow(agdat_winter), -1] <-NA    
}
  

#--- categorical rodent variable winter---
#if we split the continuous abundance into 2 -> high and low, at threshold -0.34, 9 years are high and 9 are low, 
# and the categories are the same whether we split based on st.tot or st.lemvole (as long as you choose -0.34 as your threshold)
agdat_winter$cat2 <- ifelse(agdat_winter$st.tot>-0.34, 1,0)

#--- categorical rodent variable fall---
#if we split the continuous abundance into 2 -> high and low, at threshold 0, 8 years are high and 11 are low, 
# and the categories are the same whether we split based on st.tot or st.lemvole (as long as you choose 0 as your threshold)
agdat_fall$cat2 <- ifelse(agdat_fall$st.tot>0, 1,0)

#--- Add all rodent variables together
merged_fallwinter<-merge(agdat_winter, agdat_fall, by="start_hunting_year", all = T, suffixes = c(".winter",".fall"))

return(merged_fallwinter)
}


