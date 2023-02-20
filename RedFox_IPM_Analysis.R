

#**********#
# 0) SETUP #
#**********#

## Set general parameters
Amax <- 5 # Number of age classes
Tmax <- 15  # Number of years
minYear <- 2004 # First year to consider

## Load all relevant functions


#*********************#
# 1) DATA PREPARATION #
#*********************#

# 1a) Winter Age-at-Harvest data #
#--------------------------------#

## Set data path/filename
wAaH.datafile <- "Data/Cvar.tot_oct-mai_5age.txt"

## Prepare winter AaH data
wAaH.data <- wrangleData_winterAaH(datafile = wAaH.datafile, 
                                   Amax = Amax)


# 1b) Reproduction data #
#-----------------------#

## Set data paths/filenames
rep.datafiles <- c("Data/P1var_tot.txt", # Placental scar/embryo count
                   "Data/P2var_tot.txt") # Presence of placental scars/embryos/pregnancy signs

## Prepare reproduction data
rep.data <- wrangleData_rep(datafiles = rep.datafiles, 
                            Amax = Amax, 
                            minYear = minYear)


# 1c) Environmental data #
#------------------------#



# 1d) Conceptual year information #
#---------------------------------#



