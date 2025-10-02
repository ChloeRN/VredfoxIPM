
plotMap_Data <- function(carcass.data, rodent.data,
                         genetics.datapath, pups.datapath,
                         ){
  
  #---------carcass data-------------
  
  # Convert data frame to sf object
  carcass.points <- sf::st_as_sf(x = carcass.data$fvar1, 
                                 coords = c("e_dd", "n_dd"),
                                 crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  
  #--------- genetic samples inside Varanger ----------------
  
  # Retrieve raw genetic data
  datapath <- genetics.datapath
  gen.data.raw <- read.delim(datapath)
  
  # Extract the coordinates from carcass.data.raw, where the individual ID from gen.data.raw matches (and no NA in coordinates)
  gened.carcass <- carcass.data.raw[carcass.data.raw$v_individual_id %in% gen.data.raw$v_individual_id &!is.na(carcass.data.raw$e_dd),]
  gen_varanger.points <- sf::st_as_sf(x = gened.carcass, 
                                      coords = c("e_dd", "n_dd"),
                                      crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  
  # ---- genetic samples outside Varanger ---------
  # (manually give them a location)
  # iesjavri is actually 69.657, 24.219 but we move it a bit so we can have a smaller map
  
  
  gensamples <- data.frame(site = c("iesjavri", "sorvaranger", "nordkinn"),
                           lat  = c(69.73, 69.897, 70.72 ),
                           long = c(24.63, 29.072, 27.7),
                           n    = c(158, 8, 28 ),
                           dist = c("", "", ""))
  
  gen_source.points <- sf::st_as_sf(x = gensamples, 
                                    coords = c("long", "lat"),
                                    crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  
  #----Rodent data---------
  
  # Retrieve full raw data frame
  stor <- rodent.data$stor
  
  # Extract plot location ID's from stor, and get the locations from "coordtrap.txt", where the location ID's match
  all.trapcoord <- read.delim("data/coordtrap.txt")
  rodent.used   <- all.trapcoord[all.trapcoord$Site_ID %in% stor$plot,] 
  rodent.points <- sf::st_as_sf(x = rodent.used, 
                                coords = c("Long", "Lat"),
                                crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  
  
  
  # ----- opportunistic litter size observations. Run IPM_Analysis.R to get pup datapath-----------
  
  # Load data file
  pupData <- read.csv(pups.datapath, sep = ";")
  
  #manually add locations
  obs_loc <- data.frame(site = c("Iesjavri", "Hubehiet", "nfi103", "nfi110", "nfi124"),
                        lat  = c(69.73, 70.3, 70.263 - 0.02 , 70.296, 70.140 - 0.00),
                        long = c(24.45, 29.9, 27.157 , 27.204, 27.021),
                        n    = c(1, 3, 2, 5, 1),
                        dist = c("", "", "", "", ""))
  
  pup.points <- sf::st_as_sf(x = obs_loc, 
                             coords = c("long", "lat"),
                             crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  
  
  # ---- Dummy coordinates for map annotation ----#
  
  names_loc <- data.frame(label = c("Nordkinn", "Ifjordfjellet", "Varanger Penninsula"),
                          lat  = c(71.05, 70.54, 70.8),
                          long = c(27.2, 27.5, 30.4),
                          dist = c("", "", ""))
  names.points <- sf::st_as_sf(x = names_loc, 
                               coords = c("long", "lat"),
                               crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  
  
  
  #----Basic plotting-----
  
  ## Get world map
  worldMap <- rnaturalearth::ne_countries(returnclass = "sf", scale = "large")
  
  ## Draw main map
  mainMap <- ggplot() + theme_void() +
    geom_sf(data = worldMap, fill = "grey70") +
    geom_sf(data = carcass.points, alpha = 0.7, size = 5, pch = 21, fill = "#A71B4B") +
    geom_sf(data = gen_varanger.points, alpha = 0.8, size = 5, pch = 1, color = "#0675BD", stroke = 1.5) +
    geom_sf(data = gen_source.points, alpha = 1, size = 30, pch = 1, color = "#0675BD", stroke = 2) +
    geom_sf(data = pup.points, alpha = 1, size = 10, pch = 24, fill = "#BED68A", stroke = 1.5) +
    geom_sf(data = rodent.points, alpha = 1, size = 6, color = "black", fill = "white", shape = 22, stroke = 0.5) +
    geom_sf_text(data = gen_source.points, aes(label = n), color = "#0675BD", fontface ="bold", size = 6) +
    geom_sf_text(data = pup.points, aes(label = n), color = "black", fontface ="bold", size = 5) +
    geom_sf_text(data = names.points, aes(label = label), color = "black", fontface = "italic", size = 5) +
    coord_sf(xlim = c(24, 31.5), ylim = c(69.62, 71.2), expand = F) +
    annotation_scale(location = "tr", text_cex = 1) +
    theme(panel.background = element_rect(fill = "lightblue"),
          panel.ontop = FALSE) 
  
  ## Draw insert map
  worldMap <- rnaturalearth::ne_countries(returnclass = "sf", scale = "large")
  
  insertMap <- ggplot() + 
    theme_void() +
    geom_sf(data = worldMap, fill = "white") +
    coord_sf(xlim = c(4, 42), ylim = c(55, 72), expand = FALSE) +
    theme(panel.background = element_rect(fill = "lightblue"),
          panel.ontop = FALSE) +
    geom_rect(aes(xmin = 24,
                  xmax = 31.6,
                  ymin = 69.62,
                  ymax = 71.2),
              fill = NA,
              colour = "black",
              size = 0.6) + 
    # Draw a border around insertMap
    geom_rect(aes(xmin = 4, xmax = 42, ymin = 55, ymax = 72), 
              color = "black", 
              fill = NA, 
              size = 1.2)
  
  ## Draw combined map
  combinedMap <- ggdraw(mainMap) +
    draw_plot(insertMap, x = -0.04, y = 0.493, width = 0.4, height = 0.4) # Adjust x, y, width, height as needed
  
  ## Save map in different formats
  ggsave(filename = "Plots/Map_DataAvailability.pdf", plot = combinedMap, width = 10, height = 8)
  ggsave(filename = "Plots/Map_DataAvailability.png", plot = combinedMap, width = 10, height = 8, dpi = 600)
  
  ## Retrun vector of plot names
  return(c("Plots/Map_DataAvailability.pdf", "Plots/Map_DataAvailability.png"))
}

