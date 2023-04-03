rodent.dir        <-"C:\\Users\\sho189\\OneDrive - UiT Office 365\\PhD\\RedfoxIPM\\Data from google disk\\Plot_based_data-database"


#========= LOAD DATA ==============

rodent_filenames<-dir(paste(rodent.dir, sep = "/"))
mylist<-c()
for (i in 1:length(rodent_filenames)){
  mylist[[i]]<-read.table(paste(rodent.dir, rodent_filenames[i], sep = "/"), header=T, sep = ";")
}
myfile<-do.call(rbind, mylist)  # combine all files
allrod <- myfile
rm(myfile)

# I think I need doro's previous script to see how she prepared the data exactly so I don't have to do things double
# do we need shapefile here too or we take rodents from all of varanger?