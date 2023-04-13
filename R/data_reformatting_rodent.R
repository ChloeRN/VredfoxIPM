#' Makes a yearly rodent abundance variable, mean nr of rodents in a plot, for each year
#'
#' @param rodent.dir directory of rodent data
#'
#' @return a dataset with the mean nr of rodents (in a plot), for each year
#' @export
#'
#' @examples


data_reformatting_rodent <- function(rodent.dir) {

  #========= LOAD DATA ==============

rodent_filenames<-dir(paste(rodent.dir, sep = "/"))
mylist<-c()
for (i in 1:length(rodent_filenames)){
  mylist[[i]]<-read.table(paste(rodent.dir, rodent_filenames[i], sep = "/"), header=T, sep = ";")
}
myfile<-do.call(rbind, mylist)  # combine all files
allrod <- myfile
rm(myfile)

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
  stop("uneven number of seasons in rodent dataset")
}

agdat <-  aggregate(cbind(Llem, Moec, Mruf, Mrut, rodsp, vole, tot) ~ year + foxreg, stor, mean)  #the mean nr of rodents per plot, for each year
agdat<-agdat[agdat$foxreg=="varanger",] # we only use varanger

#make categories of rodent abundance
agdat$cat2 <- ifelse(agdat$tot > 1.5,1,0)
agdat$cat3 <- ifelse(agdat$tot > 3,2,agdat$cat2)

return(agdat)
}


