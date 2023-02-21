

#**********#
# 0) SETUP #
#**********#

## Set general parameters
Amax <- 5 # Number of age classes
Tmax <- 15  # Number of years
minYear <- 2004 # First year to consider

## Source all functions in "R" folder
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir('R')



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

## Set data paths/filenames
rodent.datafile <- "Data/stor_intensiv_04_20-year-var.txt"

## Prepare rodent abundance data
rodent.data <- wrangleData_rodent(datafile = rodent.datafile,
                                  minYear = minYear,
                                  adjust = TRUE)

## Prepare harvest effort data
hunter.data <- makeData_hunters()


# 1d) Conceptual year information #
#---------------------------------#

YearInfo <- collate_yearInfo(minYear = minYear,
                             Tmax = Tmax)


# 1e) Assemble IPM input data #
#-----------------------------#

input.data <- assemble_inputData(Amax = Amax, 
                                 Tmax = Tmax, 
                                 minYear = minYear,
                                 wAaH.data = wAaH.data, 
                                 rep.data = rep.data, 
                                 rodent.data = rodent.data, 
                                 hunter.data = hunter.data, 
                                 YearInfo = YearInfo)
  

#**********************#
# 2) PRIOR INFORMATION #
#**********************#

surv.priors <- collate_priorInfo
