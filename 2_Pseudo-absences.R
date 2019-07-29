#### Generate pseudoabsences

# Load required packages
library(raster)
library(rgdal)
library(rgeos)
library(plyr)
library(dismo)

##### Within HH elevational range
#
#
# Load points and reproject
hh.pt <- readOGR(dsn = "Shapefiles/Hen Harrier data", layer = "HH_all_data")
proj.crs <- CRS("+proj=merc +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs
                +ellps=WGS84 +towgs84=0,0,0")
hh.pt <- spTransform(hh.pt, proj.crs)

# Subset shapefile by year and set year as factor to remove levels
#
# E.g. sub.fac(species.occurrences, 2015)
sub.fac <- function(sp, yr){
  o <- subset(sp, Year == yr)
  o$Year <- factor(o$Year)
  return(o)
}

hh.pt <- sub.fac(hh.pt, c(2010, 2015))
hh.pt$Year <- factor(hh.pt$Year)

# Remove duplicated rows
hh.pt_t <- remove.duplicates(hh.pt)

# Create buffer around points,dissolve borders and generate random points for pseudo-absence
x <- gBuffer(hh.pt_t, width=20000, byid=TRUE )
x$ID <- "pseudo"
y <- gUnaryUnion(x, id = x@data$ID)

samp1 <- spsample(y, 500, type='random', iter=25)
plot(samp1)

# Create coordinats dataframe and merge to create SpatialPointsDataFrame
df <- data.frame(lat = coordinates(samp1)[,1], lon =  coordinates(samp1)[,2])
pseudo <- SpatialPointsDataFrame(samp1, df)

# Write shapefile
writeOGR(pseudo, "Shapefiles", "pseudoabsences", driver="ESRI Shapefile")
#
#
##############


####### Rest of Ireland
#
#
# Load Ireland shapefile
eire <- readOGR(dsn = "Shapefiles", layer = "Ireland")

# Load points and reproject
hh.pt <- readOGR(dsn = "Shapefiles/Hen Harrier data", layer = "HH_all_data")
proj.crs <- CRS("+proj=merc +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs
                +ellps=WGS84 +towgs84=0,0,0")
hh.pt <- spTransform(hh.pt, proj.crs)

# Subset shapefile by year and set year as factor to remove levels
#
# E.g. sub.fac(species.occurrences, 2015)
sub.fac <- function(sp, yr){
  o <- subset(sp, Year == yr)
  o$Year <- factor(o$Year)
  return(o)
}

hh.pt <- sub.fac(hh.pt, c(2010, 2015))
hh.pt$Year <- factor(hh.pt$Year)

# Remove duplicated rows
hh.pt_t <- remove.duplicates(hh.pt)

# Create buffer around points,dissolve borders and generate random points for pseudo-absence
x <- gBuffer(hh.pt_t, width=20000, byid=TRUE )
x$ID <- "pseudo"
y <- gUnaryUnion(x, id = x@data$ID)

# Erase buffers from Ireland
t <- gDifference(eire, y)

# Generate random samples in Ireland, outside buffers
samp1 <- spsample(t, 500, type='random', iter=30)

# Create coordinats dataframe and merge to create SpatialPointsDataFrame
df <- data.frame(lat = coordinates(samp1)[,1], lon =  coordinates(samp1)[,2])
pseudo2 <- SpatialPointsDataFrame(samp1, df)

# Write shapefile
writeOGR(pseudo2, "Shapefiles", "pseudoabsences_Eire", driver="ESRI Shapefile")
#
#
#######################################
