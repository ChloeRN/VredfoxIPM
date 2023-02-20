

#**********#
# 0) SETUP #
#**********#

## Set general parameters
Amax <- 5 # Number of age classes
Tmax <- 15  # Number of years

## Load all relevant functions


#*********************#
# 1) DATA PREPARATION #
#*********************#

# 1a) Winter Age-at-Harvest data #
#--------------------------------#

## Set data path/filename
datafile <- "Data/Cvar.tot_oct-mai_5age.txt"

## Prepare winter AaH data
wAaH.data <- wrangleData_winterAaH(datafile = datafile, Amax = Amax)


# 1b)