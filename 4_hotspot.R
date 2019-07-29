###### Moran's I and Hotspot analyses
#
# Load required packages
library(raster)
library(rgdal)
library(rgeos)
library(ggplot2)
library(adehabitatHR)
library(ape)

# Load points and reproject
hh.pt <- readOGR(dsn = "../../shapefiles/Hen harrier data", layer = "HH_all_data")
proj.crs <- CRS("+proj=merc +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs
                +ellps=WGS84 +towgs84=0,0,0")
hh.pt <- spTransform(hh.pt, proj.crs)
hh.pt$Fledged[is.na(hh.pt$Fledged)] <- 0
hh.pt$Fledged <- as.character(hh.pt$Fledged)
hh.pts <- subset(hh.pt, Year %in% c("2010", "2015")) # Subset to focal years
hh.pts$Year <- factor(hh.pts$Year)
hh.pts$Territor_1[is.na(hh.pts$Territor_1)] <- "Confirmed"
hh.pts <- hh.pts[!(hh.pts$Territor_1 == "Possible"),] # Confirmed territories only
hh.pts$Territor_1 <- factor(hh.pts$Territor_1)
hh.pts <- hh.pts[, -c(2:16)]
hh.pts <- hh.pts[, -c(3:5)]
hh.pts$Fledged[hh.pts$Fledged %in% c("No", "Unknown", "no", "unknown")] <- 0
hh.pts$Fledged[hh.pts$Fledged == "yes"] <- 1
hh.pts$Fledged <- as.numeric(hh.pts$Fledged)

ps.pt <- readOGR(dsn = "../../shapefiles", layer = "pseudoabsences")
ps.pt <- spTransform(ps.pt, proj.crs)
ps.pt$Year <- "2017"
ps.pt$Fledged <- 0
ps.pt <- ps.pt[, -c(1:3)]

pts <- rbind(hh.pts, ps.pt)

eire <- readOGR(dsn = "../../shapefiles", layer = "Ireland")
hh.pt <- spTransform(eire, proj.crs)


### Hotspot
#
#
Extent <- extent(eire) #this is the geographic extent of the grid.

#Here we specify the size of each grid cell in metres (since those are the units our data are projected in).
resolution<- 1000

#Create an empty grid for the hotspot
x <- seq(Extent[1],Extent[2],by=resolution)  # where resolution is the pixel size you desire 
y <- seq(Extent[3],Extent[4],by=resolution) 
xy <- expand.grid(x=x,y=y) 
coordinates(xy) <- ~x+y 
gridded(xy) <- TRUE 

#You can see the grid here (this may appear solid black if the cells are small)
plot(xy)
plot(eire, border="red", add=T)

dens.map <- raster(kernelUD(hh.pts, h="LSCV", grid = xy)) #Note we are running two functions here - first KernelUD then converting the result to a raster object.

# Rescale cell values between 0 and 1 using the following function.
rasterRescale<-function(r){
  ((r-cellStats(r,"min"))/(cellStats(r,"max")-cellStats(r,"min")))
}

dens.map <- rasterRescale(dens.map)

## crop and mask
dens.map2 <- crop(dens.map, Extent)
dens.map2 <- mask(dens.map2, eire)

# Save map
png(filename="../../figures/autocorrelation/20181119_hotspot.png", type = "cairo", units="in", width=5.5, height=8, res=300)
par(family = "serif")
plot(dens.map)
plot(eire, lwd = 1, add=T)
dev.off()

### Moran's I
#
#
# To calculate Moran's I, we will need to generate a matrix of inverse distance weights.  
# In the matrix, entries for pairs of points that are close together are higher than for pairs of points that are far apart.  
# We can first generate a distance matrix, then take inverse of the matrix values and replace the diagonal entries with zero:
mat <- as.matrix(dist(cbind(coordinates(pts)[,1], coordinates(pts)[,2])))
mat.inv <- 1/mat
diag(mat.inv) <- 0

mat.inv[1:5, 1:5]

# We can now calculate Moran's I using the command Moran.I.
Moran.I(pts$Fledged, mat.inv)
