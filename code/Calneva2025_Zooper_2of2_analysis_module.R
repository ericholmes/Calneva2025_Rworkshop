## Zooper module: mapping and zooplankton analysis tutorial -----------------
## Written for the Cal-Neva AFS workshop 2025-05-30 in Lodi
## Author: Eric Holmes, contact: eric.holmes@water.ca.gov

## This is a data exploration script working with data from the Zooper R package created by folks at the IEP.
## More details on the Zooper package can be found here: https://github.com/InteragencyEcologicalProgram/zooper

## ****Power of intuition****
## Scenario: We would like to understand the spatial and temporal patterns in zooplankton community data
## collected over the past 20 years in the Sacramento-San Joaquin Delta

##  Steps: 
##  1) Find and download data (use Calneva2025_Zooper_1of2_data_grabber.R)
##  2) Prepare the data (column formatting, data aggregation/binning)
##  3) Exploratory analysis
##      (Time series plots, map sites, boxplots, species/taxa group barplots)
##  4) Advanced analysis 
##      (NMDS community analysis, spatial join, spatial interpolation, 
##        add covariates, mixed effects modeling) 
##  5) mechanistic understanding
##      (boosted regression trees, flow thresholds)
##  6) Visualizing spatiotemporal variation (create animated interpolation)

# Load libraries ----------------------------------------------------------

## Data manipulation package
library(tidyverse) # umbrella library with many useful packages for manipulating, summarizing, and plotting data

## mapping packages
library(leaflet) # Create interactive maps
library(sf) # simple features spatial package
library(gstat) # spatial modeling, prediction and simulation package
library(raster) # package for working with spatial raster data
library(sp) # package for working with spatial data

## Ecological community analysis
library(vegan) # package for standardizing and analyzing community data

##Load the imagemagick library. This imagemagick toolset is very useful for image operations
library(magick)

# Load data -----------------------------------------------------------------

load("data/Calneva2025_workshop_data.Rdata")

### Aggregate data ----
##Quick check on counts by Class
table(zoop$Class)

plot(rev(sort(table(zoop$Class)[1:20])), las = 3)

##dplyr: summarize
zooply <- zoop %>% group_by(Source, Station, Latitude, Longitude, 
                            Date, jday, year, month, Year, wy, group, SalSurf, SalBott) %>% 
  summarize(sumcatch = sum(CPUE), n = length(Station))

##Convert date formats to be compatible with flow data
zooply$Date <- as.Date(zooply$Date)
##Join zooply and wytype data
zooply <- merge(zooply, wytype, by = "wy", all.x = T)
##Join delta outflow data
zooply <- merge(zooply, dto[,c("Date", "param_val")], by = "Date", all.x = T)

##Calculate unique sampling events per site
siteN <- zooply %>% filter(is.na(Longitude) == F) %>% group_by(Source, Station, Latitude, Longitude, Date) %>% 
  summarize(sumtot = sum(sumcatch, na.rm = T)) %>% 
  group_by(Source, Station,Latitude, Longitude) %>% summarize(n = length(sumtot), 
                                  logzoop = log(mean(sumtot + 1))) %>% 
  filter(n > 50) %>% 
  mutate(sourcecol = ifelse(Source %in% "20mm", "red", ifelse(Source %in% "EMP", "blue", "black")))

## Visually inspect data ----

##Plot histogram of site sample sizes
ggplot(siteN, aes(x = Station, y = n)) + 
  geom_bar(aes(fill = Source), stat = "identity") + 
  coord_flip() +
  theme_bw() + theme(axis.text.x = element_text(angle = 90))

### Plot sites zoop sampling sites on a map ---------------------------------

##Plot the sampling sites
ggplot(siteN, aes(x = Longitude, y = Latitude, color = Source, shape = Source)) + 
  geom_point() + theme_bw()

##Convert sites dataframe to sf object
siteNsf <- st_as_sf(siteN, coords = c("Longitude", "Latitude"), crs = 4269)

##transform CRS to match deltamapr projection
siteNsf <- st_transform(siteNsf, st_crs(R_EDSM_Regions_1718P1))

##Use geom_sf() function in ggplot to create a map of the SFE with zoop sampling sites
(sfemap <- ggplot() + 
  geom_sf(data = WW_Delta, fill = "blue", color = "blue") + 
  geom_sf(data = R_EDSM_Regions_1718P1$geometry, fill = 1:4, alpha = .3) +
  geom_sf(data = R_EDSM_Regions_1718P1$geometry, fill = NA, linewidth = 1, color = 1) +
  geom_sf(data = siteNsf$geometry, size = 3) + 
  geom_sf(data = siteNsf, aes(color = Source)) +
  theme_bw())

##Plot sampling sites in an interactive map
leaflet(siteN) %>% addCircles(lng = ~Longitude, lat = ~Latitude,
                              label = ~Station, color = ~sourcecol,
                              labelOptions = labelOptions(noHide = F, textOnly = F)) %>% addTiles()

##Spruce up the interactive map with persistent labels, subsetted sites, and a basemap toggle
leaflet(siteN[siteN$n>50,]) %>% addCircles(lng = ~Longitude, lat = ~Latitude,
                                           label = ~Station, color = ~sourcecol,
                                           labelOptions = labelOptions(noHide = T, textOnly = F)) %>%
  addTiles(options = providerTileOptions(noWrap = TRUE), group="Base") %>%
  addProviderTiles("Esri.WorldImagery", group="Imagery") %>%
  addLayersControl(baseGroups = c("Base","Imagery"), options = layersControlOptions(collapsed = FALSE))

##Highlight one site at the confluence
sfemap + geom_sf(data = siteNsf[siteNsf$Station == "NZ048", ], color = "chartreuse", size = 7) +
  geom_sf(data = siteNsf[siteNsf$Station == "NZ048", ], size = 4)

##Plot time series of single site
ggplot(zooply[zooply$Station == "NZ048",], aes(x = Date, y = sumcatch)) + geom_point() + theme_bw()

##Plot log transformed CPUE by julien day, faceted by group
ggplot(zooply[zooply$Station == "NZ048",], aes(x = jday, y = log10(sumcatch + 1))) + 
  geom_point() + stat_smooth(se = F) + theme_bw() +
  facet_grid(group ~ ., scales = "fixed")

##Plot log CPUE distribution by month
ggplot(zooply[zooply$Station == "NZ048",], aes(x = month, y = log10(sumcatch + 1))) + 
  geom_boxplot(aes(group = month), fill = "black") +
  facet_grid(group ~ ., scales = "fixed") + stat_smooth(se = F, color = "brown") + theme_bw()

zooptax <- zoop %>% group_by(Taxname, group) %>% summarize(sumtot = sum(CPUE))
zoopgroup <- zoop %>% group_by(group) %>% summarize(sumtot = sum(CPUE))

ggplot(zooptax, aes(x = reorder(Taxname, sumtot), y = sumtot)) + coord_flip() +
  geom_bar(stat = "identity") + facet_wrap(group ~ ., scales = "free_y") + scale_y_log10()

ggplot(zooptax[zooptax$group == "Copepoda",], aes(x = reorder(Taxname, sumtot), y = sumtot)) + coord_flip() +
  geom_bar(stat = "identity") + facet_wrap(group ~ ., scales = "free_y") #+ scale_y_log10()

ggplot(zoopgroup, aes(x = reorder(group, -sumtot), y = sumtot)) + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

##Look at individual taxa

ggplot(zoop[zoop$Station == "NZ048" & 
              zoop$Taxname == "Pseudodiaptomus forbesi",], 
       aes(x = jday, y = log10(CPUE + 1))) + 
  geom_point() + stat_smooth(se = F) + 
  theme_bw() + labs(x = "julien day", title = "Pseudodiaptomus forbesi") +
  facet_wrap(wy ~ ., scales = "fixed")

ggplot(zoop[zoop$Station == "NZ048" & 
                 zoop$Taxname %in% c("Bosmina longirostris" ),], 
       aes(x = jday, y = log10(CPUE + 1))) + 
  geom_point() + stat_smooth(se = F) + theme_bw() +
  labs(x = "julien day", title = "Daphia sp.") +
  facet_wrap(wy ~ ., scales = "fixed")

ggplot(zoop[zoop$Station == "NZ048" & 
              zoop$Taxname %in% c("Daphnia_UnID", "Daphniidae_all_Meso", "Daphniidae_UnID" ),], 
       aes(x = jday, y = log10(CPUE + 1))) + 
  geom_point() + stat_smooth(se = F) + theme_bw() +
  labs(x = "julien day", title = "Daphia sp.") +
  facet_wrap(wy ~ ., scales = "fixed")

# NDMS ----

# Change data structure from long to wide
zoopcast <- pivot_wider(zooply[zooply$Source == "20mm" & zooply$Year %in% c(2015, 2017),#c(2015, 2017),
                         c("Year", "group", "sumcatch", "Station", "Latitude", "Longitude", "Date")],
                  names_from = group, values_from = "sumcatch") %>% data.frame()

hist(zoopcast$Malacostraca)

# zoopcast[is.na(zoopcast)] <- 0
## Community data transformation: using hellinger (see Legendre & Gallagher 2001)
## other common transformations include: “chi.square”, “log”, “normalize”, “range”
deco <- decostand(zoopcast[, c("Branchiopoda", "Cirripedia", "Copepoda", 
                               "Rotifera", "Ostracoda")], 
                  method = "hellinger") 

##run the nmds analysis with bray distance matrix
vare.mds <- metaMDS(deco[,], distance = "bray", trymax = 10)
vare.mds

## Basic NMDS visualization
plot(vare.mds, type = "t")

##Extract scores from NMDS analysis into a data.frame for a better visualization
data.scores <- as.data.frame(scores(vare.mds, "site"))

## append metadata columns from zoop community data.frame
data.scores <- cbind(data.scores, 
                     zoopcast[, c("Station", "Year", "Latitude", "Longitude", "Date")])

data.scores$Year <- as.factor(data.scores$Year)

#Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores <- as.data.frame(scores(vare.mds, "species")) 
# create a column of species, from the rownames of species.scores 
species.scores$species <- rownames(species.scores)  

##Plot nmds data in ggplot
ggplot() +
  geom_point(data = data.scores, aes(x = NMDS1, y = NMDS2, color = Year), size = 1.5) + theme_bw() +
  geom_segment(data = species.scores, aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               colour="black", alpha = .2, linewidth = 2) +
  stat_ellipse(data = data.scores, aes(x = NMDS1, y = NMDS2, color = Year))+
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), alpha=0.8, size = 2) +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank())

ggplot() +
  geom_point(data = data.scores, aes(x = NMDS1, y = NMDS2, color = Year), size = 1.5) + theme_bw() +
  geom_segment(data = species.scores, aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               colour="black", alpha = .2, linewidth = 2) +
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), alpha=0.8, size = 2) +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank())

### Spatial join delta regions -----
sfemap
##convert nmds data output to sf object
data.scoressf <- st_as_sf(data.scores, 
                          coords = c("Longitude", "Latitude"),
                          crs = 4269)

##transform CRS to match deltamapr spatial data
data.scoressf <- st_transform(data.scoressf, st_crs(R_EDSM_Regions_1718P1))

##Spatial join
dsjoin <- st_join(data.scoressf, R_EDSM_Regions_1718P1, left = T)
dsjoin <- dsjoin[is.na(dsjoin$Region) == F, ]

##add a season info
dsjoin$Month <- as.numeric(format(dsjoin$Date, format = "%m"))
dsjoin$Season <- ifelse(dsjoin$Month %in% c(1:5), "early", "late")

dsjoin$Regionfac <- factor(dsjoin$Region, levels = c("Far West", "North", "West", "South"))


ggplot() +
  geom_point(data = dsjoin[,], aes(x = NMDS1, y = NMDS2, color = Year), size = 1.5) + theme_bw() +
  stat_ellipse(data = dsjoin[,], aes(x = NMDS1, y = NMDS2, color = Year))+
  geom_segment(data = species.scores, aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               colour="black", alpha = .2, linewidth = 2) +
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), alpha=1, size = 3) +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) + 
  facet_wrap(Regionfac ~ .)

##Subset to early season points
ggplot() +
  geom_point(data = dsjoin[dsjoin$Season == "early",], aes(x = NMDS1, y = NMDS2, color = Year), size = 1.5) + theme_bw() +
  stat_ellipse(data = dsjoin[dsjoin$Season == "early",], aes(x = NMDS1, y = NMDS2, color = Year))+
  geom_segment(data = species.scores, aes(x=0, xend=NMDS1, y=0, yend=NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               colour="black", alpha = .2, linewidth = 2) +
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), alpha=1, size = 3) +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank()) + 
  facet_wrap(Regionfac ~ .)

# Plot water year type ------------------------------------------------------

##Plot histogram of Sac River water year runoff index
ggplot(wytype, aes(x = sac_index)) + 
  geom_histogram(aes(fill = sac_yr_typefac)) + 
  scale_fill_brewer(palette = "RdYlBu", direction = -1) +
  theme_bw()

##Plot single station response to water year type
ggplot(zooply[zooply$Station == "NZ048",], aes(x = sac_yr_typefac, y = log10(sumcatch + 1))) + geom_boxplot() +
  facet_wrap(group ~ ., scales = "free") + theme_bw()

##Plot log cladocera CPUE versus Oct-mar flow index for a single station
ggplot(zooply[zooply$Station == "NZ048" &
                     zooply$group %in% "Branchiopoda" & zooply$month %in% c(10:12, 1:3),], 
       aes(x = sac_oct_mar, y = log10(sumcatch + 1))) + geom_point() +
  facet_grid(group ~ ., scales = "free") + stat_smooth(se = F) + labs(x = "Oct-Mar flow index")

##Plot log cladocera CPUE versus Oct-mar flow index faceted by station
ggplot(zooply[zooply$Station %in% unlist(siteN[siteN$n > 50, "Station"]) & 
                grepl(zooply$Source, pattern = "EMP") &
                zooply$group %in% "Branchiopoda" & 
                zooply$month %in% c(10:12, 1:3),], aes(x = sac_oct_mar, y = log10(sumcatch + 1))) + geom_point() +
  facet_wrap(Station ~ ., scales = "free") + stat_smooth(se = F, color = "red") + 
  labs(x = "Oct-Mar flow index", y = "Log cladocera CPUE")

# Plot flow data ----------------------------------------------------------

##Calculate mean flow for the jan-may period for each wy
dtoply <- dto[dto$month %in% 1:3,] %>% group_by(wy) %>% summarize(meanflow = mean(param_val))

##Plot log Cladocera CPUE response to delta outflow
ggplot(zooply[zooply$Station %in% unlist(siteN[siteN$n > 50, "Station"]) & 
                grepl(zooply$Source, pattern = "20mm") &
                zooply$group %in% "Branchiopoda",], aes(x = param_val, y = log10(sumcatch + 1))) + geom_point() +
  facet_wrap(Station ~ .) + stat_smooth(se = F, color = "red") + theme_bw()

ggplot(zooply[zooply$Station %in% "342" & 
                zooply$group %in% "Branchiopoda" &
                zooply$month %in% c(1:3),], aes(x = param_val, y = log10(sumcatch + 1))) + geom_point() +
  facet_wrap(Station ~ .) + stat_smooth(se = F, color = "red") + theme_bw() + 
  labs(x = "Delta outflow (cfs)", y = "log(Branchiopoda CPUE)")

ggplot(zooply[zooply$Station %in% unlist(siteN[siteN$n > 50, "Station"]) & 
                grepl(zooply$Source, pattern = "EMP") &
                zooply$group %in% "Branchiopoda" & 
                zooply$month %in% c(10:12, 1:3),], aes(x = param_val, y = log10(sumcatch + 1))) + geom_point() +
  facet_wrap(Station ~ .) + stat_smooth(se = F, color = "red") + theme_bw() + 
  labs(x = "Delta outflow (cfs)", y = "log(Branchiopoda CPUE)")

ggplot(zooply[zooply$Station %in% "NZ048" & 
                zooply$group %in% "Branchiopoda" &
                zooply$month %in% c(1:3),], aes(x = param_val, y = log10(sumcatch + 1))) + geom_point() +
  facet_wrap(Station ~ .) + stat_smooth(se = F, color = "red") + theme_bw() + 
  labs(x = "Delta outflow (cfs)", y = "log(Branchiopoda CPUE)")

##Look at Branchiopoda response to salinity

ggplot(zooply[zooply$Station %in% unlist(siteN[siteN$n > 50, "Station"]) & 
                grepl(zooply$Source, pattern = "EMP") &
                zooply$group %in% "Branchiopoda",], aes(x = SalSurf, y = log10(sumcatch + 1))) + geom_point() +
  facet_wrap(Station ~ .) + stat_smooth(se = F, color = "red") + theme_bw() +
  labs(x = "Salinity (PSU)", y = "log(Branchiopoda CPUE)")

ggplot(zooply[zooply$Station %in% "NZ048" & 
                zooply$group %in% "Branchiopoda" &
                zooply$month %in% c(1:3),], aes(x = SalSurf, y = log10(sumcatch + 1))) + geom_point() +
  facet_wrap(Station ~ .) + stat_smooth(se = F, color = "red") + theme_bw() + 
  labs(x = "Salinity (PSU)", y = "log(Branchiopoda CPUE)")

# Spatial interpolation ---------------------------------------------------
##adapted from this example: https://mgimond.github.io/Spatial/interpolation-in-r.html

##Transform delta regions to WGS84 geographic coordinate system
deltawgs84 <- st_transform(R_EDSM_Regions_1718P1, CRS("+ellps=WGS84 +proj=longlat +datum=WGS84 +no_defs"))

##Convert deltawgs84 object from sf to sp - spatial points dataframe
deltasp <- as(deltawgs84, "Spatial")

##Convert siteN from data.frame to sf object
zoopsf <- st_as_sf(siteN, coords = c("Longitude", "Latitude"))
##Set coordinate reference system to wgs84
zoopsf <- st_set_crs(zoopsf, CRS("+ellps=WGS84 +proj=longlat +datum=WGS84 +no_defs"))
##Convert object from sf to sp
zoopsp <- as(zoopsf, "Spatial")
##Create an empty grid where n is the total number of cells
grd <- as.data.frame(spsample(zoopsp, "regular", n=100000))
names(grd) <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd) <- TRUE  # Create SpatialPixel object
fullgrid(grd) <- TRUE  # Create SpatialGrid object

##Add pojection information to the empty grid
proj4string(grd) <- proj4string(zoopsp)
zoopsp@bbox <- deltasp@bbox

##Interpolate the grid cells using a power value of 2 (idp=2.0)
zoopidw <- gstat::idw(logzoop ~ 1, zoopsp, newdata = grd, idp=2.0)

##Convert spatial grid dataframe to raster object
zoopras <- raster(zoopidw)
plot(zoopras)
plot(deltasp, col = NULL, add = T)

##clip raster to delta open water geometry polygon
zooprasdelta <- mask(zoopras, deltasp)
plot(zooprasdelta)

##Spruce it up with a better color ramp
plot(zooprasdelta,
     breaks = seq(zoopras@data@min, zoopras@data@max, length.out = 10), 
     col = hcl.colors(10, "Lajolla"))

# Create a loop to visualize annual Branchiopoda spatial distributions --------

##Subset to Branchiopoda class and convert dataframe to simple features spatial object
cladsf <- st_as_sf(zooply[zooply$group %in% "Branchiopoda" & is.na(zooply$Longitude) == F,], 
                   coords = c("Longitude", "Latitude"))

##Set CRS to WGS geographic coordinate system
cladsf <- st_set_crs(cladsf, CRS("+ellps=WGS84 +proj=longlat +datum=WGS84 +no_defs"))

##Calculate log total CPUE
cladsf$logzoop <- log(cladsf$sumcatch + 1)

##Convert sf object to spatial points dataframe
cladsp <- as(cladsf, "Spatial")

##Set standard breaks for symbolizing the interpolated CPUE values
brks <- seq(0, max(cladsf$logzoop), .5)

##Loop
for(i in unique(cladsp$wy)){
  ##Print the water year the console for progress update
  print(i)
  ##Interpolation
  cladidw <- gstat::idw(logzoop ~ 1, cladsp[cladsp$wy == i & cladsp$month %in% c(1:5), ], newdata=grd, idp=2.0)
  ##Mask by delta open water polygon
  cladrasdelta <- mask(raster(cladidw), deltasp)
  
  ##Save plot of interpolated values for each water year
  ##NOTE: will fail if you do not have an output/Branchiopoda_annual folder in your working directory
  png(paste("output/Branchiopoda_annual/Branchiopoda_", i, ".png", sep = ""), 
      height = 7, width = 9, unit = "in", res = 300)
  
  plot(cladrasdelta, breaks = brks, col = hcl.colors(length(brks), "Lajolla"), 
       main = paste(i, ": Jan-May Branchiopoda log CPUE", sep = ""))
  points(cladsp[cladsp$year == i & cladsp$month %in% c(1:5), ])
  
  # Add an inset plot with a timeline of the sac wy index
  # Define inset area
  par(fig = c(0.12, 0.45, 0.20, 0.41), new = TRUE)
  # highlight rule for selected wy
  bar_colors <- ifelse(unique(cladsp$wy) == i, "salmon3", "gray95")
  # create the barplot
  barplot(sac_index~wy, data = wytype[wytype$wy %in% unique(cladsp$wy),], 
          col = bar_colors, cex.axis = .8, cex.names = .7, cex.lab = .8, las = 2, 
          xlab = NA, ylab = "Sac R. Index")

  dev.off()
}

# Create animation --------------------------------------------------------

##Use the magick package to load, scale, and join images in an animated gif
list.files(path='output/Branchiopoda_annual/', pattern = '*.png', full.names = TRUE) %>% 
  image_read() %>% # reads each path file
  image_scale("1000") %>% # resize image
  image_join() %>% # joins image
  image_animate(fps=2) %>% # animates, can opt for number of loops
  image_write("output/deltazoop_animation.gif")

print("End of module, thank you!")
