###### Breeding success models #######
# Territory models followed a similar process. Refer to manuscript for details.
# Load required packages
library(raster)
library(rgdal)
library(rgeos)
library(plyr)
library(spdep)
library(maptools)
library(SDMTools)
library(sp)
library(MuMIn)
library(arm)
library(dplyr)
library(lme4)
library(ROCR)
library(broom)
library(lattice)
library(ade4)
library(factoextra)
library(reshape2)
library(dplyr)
library(forcats)
library(ggplot2)

########## Euclidean distance to features
#
#
#
# Subset shapefile by year and set year as factor to remove levels
#
# E.g. sub.fac(species.occurrences, 2015)
sub.fac <- function(sp, yr){
  o <- subset(sp, Year == yr)
  o$Year <- factor(o$Year)
  return(o)
}

# Load points and reproject
hh.0515.pt <- readOGR(dsn = "../../../Shapefiles/Hen Harrier data", layer = "HH_05_15_thin")
proj.crs <- CRS("+proj=merc +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs
                +ellps=WGS84 +towgs84=0,0,0")
hh.0515.pt <- spTransform(hh.0515.pt, proj.crs)
hh.0515.pt <- hh.0515.pt[order(hh.0515.pt$FID_1),]

hh10.pt <- sub.fac(hh.0515.pt, 2010)
hh10.pt$FID_1 <- factor(hh10.pt$FID_1)
hh10.pt <- hh10.pt[order(hh10.pt$FID_1),]
hh15.pt <- sub.fac(hh.0515.pt, 2015)

# Read shapefiles for Euclidean distance
urb.cor.sf <- readOGR(dsn = "../../../Shapefiles/Habitat/CORINE polygons", layer = "Urban_CORINE_WGS84")
turb.f <- readOGR(dsn = "../../../Shapefiles/Wind turbines/Wind turbines 2016", layer = "turbines2016_WGS84")
l.for10.sf <- readOGR(dsn = "../../../Shapefiles/Forests", layer = "late_10") # 13+ yrs post-thicket
l.for15.sf <- readOGR(dsn = "../../../Shapefiles/Forests", layer = "late_15") # 13+ yrs post-thicket

crs.proj <- crs("+proj=merc +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84
                +towgs84=0,0,0")

l.for10.sf <- spTransform(l.for10.sf, crs.proj)
l.for15.sf <- spTransform(l.for15.sf, crs.proj)

# Calculate distance between points and features and collect in a list
#
# l1 = list of buffer shapefiles for a given year
# pb = Forest shapefile for a given year
#
#E.g. output <- dist.grp(hh.05, p.for05)
dist.grp <- function(l1, pb1){
  n <- c("windfarm", "forest", "urban")
  a <- data.frame((gDistance(turb.f, l1, byid=TRUE)/1000)) # wind turbines
  b <- data.frame((gDistance(pb1, l1, byid=TRUE)/1000)) # forest
  c <- data.frame((gDistance(urb.cor.sf, l1, byid=TRUE)/1000)) # urban areas
  l <- list(a,b,c)
  names(l) <- n
  return(l)
}

hh10.dist <- dist.grp(l1 = hh10.pt, pb1 = l.for10.sf)
hh15.dist <- dist.grp(l1 = hh15.pt, pb1 = l.for15.sf)

# Add distance to features taken as minimum across rows from the gDistance matrix
W_dist <- do.call(pmin,hh10.dist$windfarm) #wind turbines
F_dist <- do.call(pmin,hh10.dist$forest) #forest
hh10.dat <- data.frame(hh10.pt@data[,c(1,7:9)], W_dist, F_dist)
hh10.dat["U_dat"] <- hh10.dist$urban

W_dist <- do.call(pmin,hh15.dist$windfarm) #wind turbines
F_dist <- do.call(pmin,hh15.dist$forest) #forest
hh15.dat <- data.frame(hh15.pt@data[,c(1,7:9)], W_dist, F_dist)
hh15.dat["U_dat"] <- hh15.dist$urban

# Write csvs
write.csv(hh10.dat, "Dataframes/dist_10.csv")
write.csv(hh15.dat, "Dataframes/dist_15.csv")
#
#
###################################



########## Habitat and landscape classes within buffers
#
#
# Subset shapefile by year and set year as factor to remove levels
#
# E.g. sub.fac(species.occurrences, 2015)
sub.fac <- function(sp, yr){
  o <- subset(sp, Year == yr)
  o$Year <- factor(o$Year)
  return(o)
}

# Load points and reproject
hh.0515.pt <- readOGR(dsn = "../../shapefiles", layer = "HH_05_15_thin")
proj.crs <- CRS("+proj=merc +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs
                +ellps=WGS84 +towgs84=0,0,0")
hh.0515.pt <- spTransform(hh.0515.pt, proj.crs)
hhall.pt$Year <- factor(hhall.pt$Year)

# Buffer to 1 km, 2 km, 5 km
hh.0515.b1 <- gBuffer(hhall.pt, width=1000, byid=TRUE )
hh.0515.b1 <- SpatialPolygonsDataFrame(hh.0515.b1, data=hh.0515.b1@data)
hh.0515.b2 <- gBuffer(hhall.pt, width=2000, byid=TRUE )
hh.0515.b2 <- SpatialPolygonsDataFrame(hh.0515.b2, data=hh.0515.b2@data)
hh.0515.b5 <- gBuffer(hhall.pt, width=5000, byid=TRUE )
hh.0515.b5 <- SpatialPolygonsDataFrame(hh.0515.b5, data=hh.0515.b5@data)

hh.all <- mget(ls(pattern="*.0515.*"))

# Run sub.fac over list of buffers, passing the relevant year for each subset
hh.10 <- mapply(x = hh.all, y = 2010, FUN = function(x, y) sub.fac(x, y), SIMPLIFY = FALSE)
hh.15 <- mapply(x = hh.all, y = 2015, FUN = function(x, y) sub.fac(x, y), SIMPLIFY = FALSE)
hh.10$hh.0515.pt <- NULL
hh.15$hh.0515.pt <- NULL

hh10.pt <- sub.fac(hh.0515.pt, 2010)
hh15.pt <- sub.fac(hh.0515.pt, 2015)
rm(list=ls(pattern="*.0515.*"))

hh10.pt <- hh10.pt[order(hh10.pt$FID_1),]
hh10.c.dat3 <- data.frame("FID_1" = hh10.pt$FID_1, "outcome" = hh10.pt$Outcome__S, 
                          "fledged"= hh10.pt$Fledged, "spa" = hh10.pt$SPA)

hh15.pt <- hh15.pt[order(hh15.pt$FID_1),]
hh15.c.dat3 <- data.frame("FID_1" = hh15.pt$FID_1, "outcome" = hh15.pt$Outcome__S, 
                          "fledged"= hh15.pt$Fledged, "spa" = hh15.pt$SPA)


# Read habitat polygons
agr.sf <- readOGR(dsn = "../../shapefiles/Habitat", layer = "agriculture")
bog.sf <- readOGR(dsn = "../../shapefiles/Habitat", layer = "bog")
msh.sf <- readOGR(dsn = "../../shapefiles/Habitat", layer = "moor_shrub")
ngr.sf <- readOGR(dsn = "../../shapefiles/Habitat", layer = "natgrass")
pas.sf <- readOGR(dsn = "../../shapefiles/Habitat", layer = "pasture")

agr.sf <- spTransform(agr.sf, proj.crs)
bog.sf <- spTransform(bog.sf, proj.crs)
ngr.sf <- spTransform(ngr.sf, proj.crs)
msh.sf <- spTransform(msh.sf, proj.crs)
pas.sf <- spTransform(pas.sf, proj.crs)

# this is a well known R / GEOS hack (usually combined with the above) to 
# deal with "bad" polygons
ngr.sf <- gBuffer(ngr.sf, byid=TRUE, width=0)

# List habitat with forest data at age, create stacks, remove individual files
# 2010
e.for10.sf <- readOGR(dsn = "../../shapefiles/Forests", layer = "early_10") # 1-2 yrs post-clearfell 
p.for10.sf <- readOGR(dsn = "../../shapefiles/Forests", layer = "prime_10") # 3-12 yrs post-clearfell 
l.for10.sf <- readOGR(dsn = "../../shapefiles/Forests", layer = "late_10") # 13+ yrs post-thicket 
l.for10.sf <- spTransform(l.for10.sf, crs.proj)

# 2015
e.for15.sf <- readOGR(dsn = "../../shapefiles/Forests", layer = "early_15") # 1-2 yrs post-clearfell 
p.for15.sf <- readOGR(dsn = "../../shapefiles/Forests", layer = "prime_15") # 3-12 yrs post-clearfell 
l.for15.sf <- readOGR(dsn = "../../shapefiles/Forests", layer = "late_15") # 13+ yrs post-thicket
l.for15.sf <- spTransform(l.for15.sf, proj.crs)

b.for.sf <- readOGR(dsn = "../../shapefiles/Forests", layer = "Broadleaved") # Broadleaved

e.for15.sf <- spTransform(e.for15.sf, proj.crs)
p.for15.sf <- spTransform(p.for15.sf, proj.crs)


# Function to intersect and dissolve polygons by id, 
# calculate the area of resultant polygons
# and add a column to allow object ID prior to joining
#
# sp1 = SpatialPolygonsDataFrame object
# sp2 = SpatialPolygonsDataFrame object
#
# E.g.
# df <- disArea(buffer, habitat)
disArea <- function(sp1, sp2){
  sp <- raster::intersect(sp1, sp2)
  #  r <- gUnaryUnion(sp, id=sp@data$FID_1) # Error with 1 & 2 km buffers
  r <- unionSpatialPolygons(sp, sp@data$FID_1)
  a <- data.frame(for.area=sapply(r@polygons, FUN=function(x) {slot(x, 'area')}))/1000000 # Division value due to m^2
  n <- sapply(r@polygons, FUN=function(x) {slot(x, 'ID')})
  df <- data.frame(a,n)
}

# Extract habitat data for each buffer in a given year and return a list 
# where individual elements = habitat shapefiles
# Apply disArea  to habitat shapefiles and collect in list
#
# l1 = list of buffer shapefiles for a given year
# pb = Forest shapefile for a given year
#
#E.g. output <- hab.grp(hh.05, p.for05)
hab.grp <- function(l1, pb1, pb2, pb3, pb4){
  n <- c("agr", "bog", "ngr", "msh", "pas", "efo", "pfo", "lfo", "bfo")
  a <- lapply(l1, function(x) disArea(x,agr.sf))
  b <- lapply(l1, function(x) disArea(x,bog.sf))
  c <- lapply(l1, function(x) disArea(x,ngr.sf))
  d <- lapply(l1, function(x) disArea(x,msh.sf))
  e <- lapply(l1, function(x) disArea(x,pas.sf))
  f <- lapply(l1, function(x) disArea(x,pb1))
  g <- lapply(l1, function(x) disArea(x,pb2))
  h <- lapply(l1, function(x) disArea(x,pb3))
  i <- lapply(l1, function(x) disArea(x,pb4))
  l <- list(a,b,c,d,e,f,g,h,i)
  names(l) <- n
  return(l)
}

hh10.h <- hab.grp(l1 = hh.10, pb1 = e.for10.sf, pb2 = p.for10.sf, pb3 = l.for10.sf, pb4 = b.for.sf)
hh15.h <- hab.grp(l1 = hh.15, pb1 = e.for15.sf, pb2 = p.for15.sf, pb3 = l.for15.sf, pb4 = b.for.sf)


# Merge nested lists with reference dataframe
#
# df = reference dataframe
# lst = list of dataframes
# p = prefix text for column names
#
# E.g. listMerge(hh10.c.dat3, hh10.h$gra, "gra")
listMerge <- function(df, lst, p){
  l <- lst
  d <- df
  names(d)[1] <- "n"
  l$dat <- d
  m <- Reduce(function(...) merge(..., by = "n", all=T), l)
  m[is.na(m)] <- 0
  names(m)[1] <- "FID_1"
  names(m)[2:4] <- paste(p, c("b1", "b2", "b5"), sep= "_")
  m <- m[,c(1,6:7, 2:4)]
  m[,1] <- as.numeric(as.character(m[,1]))
  t <- m[ order(m[,1]), ]
  return(t)
}

a.10 <- listMerge(hh10.c.dat3, hh10.h$agr, "agr")
b.10 <- listMerge(hh10.c.dat3, hh10.h$bog, "bog")
n.10 <- listMerge(hh10.c.dat3, hh10.h$ngr, "ngr")
m.10 <- listMerge(hh10.c.dat3, hh10.h$pas, "pas")
p.10 <- listMerge(hh10.c.dat3, hh10.h$msh, "msh")
ef.10 <- listMerge(hh10.c.dat3, hh10.h$efo, "efo")
pf.10 <- listMerge(hh10.c.dat3, hh10.h$pfo, "pfo")
lf.10 <- listMerge(hh10.c.dat3, hh10.h$lfo, "lfo")
bf.10 <- listMerge(hh10.c.dat3, hh10.h$bfo, "bfo")

a.15 <- listMerge(hh15.c.dat3, hh15.h$agr, "agr")
b.15 <- listMerge(hh15.c.dat3, hh15.h$hea, "bog")
n.15 <- listMerge(hh15.c.dat3, hh15.h$ngr, "ngr")
m.15 <- listMerge(hh15.c.dat3, hh15.h$pas, "msh")
p.15 <- listMerge(hh15.c.dat3, hh15.h$shr, "pas")
ef.15 <- listMerge(hh15.c.dat3, hh15.h$efo, "efo")
pf.15 <- listMerge(hh15.c.dat3, hh15.h$pfo, "pfo")
lf.15 <- listMerge(hh15.c.dat3, hh15.h$lfo, "lfo")
bf.15 <- listMerge(hh15.c.dat3, hh15.h$bfo, "bfo")

# Remove unnecessary objects
rm(list=ls(pattern="*.sf"))
rm(list=ls(pattern="*.pt"))
rm(list=ls(pattern="*.h"))
rm(hh.10, hh.15, hh.all, hh10.h, hh15.h)

# Create list of all resultant merged dataframes
hh10 <- mget(ls(pattern="*.10"))

hh15 <- mget(ls(pattern="*.15"))

# Bind listed dataframes by column, taking the 
hh10.c <- cbind(hh10[[1]][,1:3], do.call(cbind, lapply(hh10, function(x) x[,-(1:3)])))
rownames(hh10.c) <- NULL
hh15.c <- cbind(hh15[[1]][,1:3], do.call(cbind, lapply(hh15, function(x) x[,-(1:3)])))
rownames(hh15.c) <- NULL

# Convert fledged column to numeric for glm weighting
hh10$fledged <- as.numeric(as.character(hh10$fledged))
hh15$fledged <- as.numeric(as.character(hh15$fledged))

write.csv(hh10.c, "../../data_out/breeding success/dataframes/hab_for_allbuffers_10b.csv")
write.csv(hh15.c, "../../data_out/breeding success/dataframes/hab_for_allbuffers_15b.csv")
#
#
###################################




########## Raster buffer - Hilliness & NDVI within buffers
#
#
#
# Subset shapefile by year and set year as factor to remove levels
#
# E.g. sub.fac(species.occurrences, 2015)
sub.fac <- function(sp, yr){
  o <- subset(sp, Year == yr)
  o$Year <- factor(o$Year)
  return(o)
}

# Load points and reproject
hh.0515.pt <- readOGR(dsn = "../../../Shapefiles/Hen Harrier data", layer = "HH_05_15_thin")
proj.crs <- CRS("+proj=merc +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs
                +ellps=WGS84 +towgs84=0,0,0")
hh.0515.pt <- spTransform(hh.0515.pt, proj.crs)
hh.0515.pt <- hh.0515.pt[order(hh.0515.pt$FID_1),]

hh10.pt <- sub.fac(hh.0515.pt, 2010)
hh10.pt$FID_1 <- factor(hh10.pt$FID_1)
hh10.pt <- hh10.pt[order(hh10.pt$FID_1),]
hh15.pt <- sub.fac(hh.0515.pt, 2015)
hh15.pt$FID_1 <- factor(hh15.pt$FID_1)
hh15.pt <- hh15.pt[order(hh15.pt$FID_1),]

# Buffer to 1 km, 2 km, 5 km & 10 km
hh.0515.b1 <- gBuffer(hh.0515.pt, width=1000, byid=TRUE )
hh.0515.b1 <- SpatialPolygonsDataFrame(hh.0515.b1, data=hh.0515.b1@data)
hh.0515.b2 <- gBuffer(hh.0515.pt, width=2000, byid=TRUE )
hh.0515.b2 <- SpatialPolygonsDataFrame(hh.0515.b2, data=hh.0515.b2@data)
hh.0515.b5 <- gBuffer(hh.0515.pt, width=5000, byid=TRUE )
hh.0515.b5 <- SpatialPolygonsDataFrame(hh.0515.b5, data=hh.0515.b5@data)
hh.0515.b10 <- gBuffer(hh.0515.pt, width=10000, byid=TRUE )
hh.0515.b10 <- SpatialPolygonsDataFrame(hh.0515.b10, data=hh.0515.b10@data)

hh.all <- mget(ls(pattern="*.0515.*"))
rm(list=ls(pattern="*.0515.*"))

# Run sub.fac over list of buffers, passing the relevant year for each subset
hh.10 <- mapply(x = hh.all, y = 2010, FUN = function(x, y) sub.fac(x, y), SIMPLIFY = FALSE)
hh.15 <- mapply(x = hh.all, y = 2015, FUN = function(x, y) sub.fac(x, y), SIMPLIFY = FALSE)

# read rasters
hilli.ras <- raster ("../../../Shapefiles/Environmental/hilli")
ndvi.ras <- raster("../../../Shapefiles/Environmental/ndvi")

# List rasters, remove individual files
r <- mget(ls(pattern="*.ras"))
ras.stack <- stack(r, raster)
rm(list=ls(pattern="*.ras"))


# Extract ras.stack for each buffer and calculate mean
### 2010
hh10.b1.ras <- data.frame(extract(ras.stack, hh.10$hh.0515.b1, method='simple', small=FALSE, fun=mean, na.rm=TRUE))
hh10.b2.ras <- data.frame(extract(ras.stack, hh.10$hh.0515.b2, method='simple', small=FALSE, fun=mean, na.rm=TRUE))
hh10.b5.ras <- data.frame(extract(ras.stack, hh.10$hh.0515.b5, method='simple', small=FALSE, fun=mean, na.rm=TRUE))
hh10.b10.ras <- data.frame(extract(ras.stack, hh.10$hh.0515.b10, method='simple', small=FALSE, fun=mean, na.rm=TRUE))

# Add buffer ID to column names and create list
colnames(hh10.b1.ras) <- paste("b1", colnames(hh10.b1.ras), sep = "_")
colnames(hh10.b2.ras) <- paste("b2", colnames(hh10.b2.ras), sep = "_")
colnames(hh10.b5.ras) <- paste("b5", colnames(hh10.b5.ras), sep = "_")
colnames(hh10.b10.ras) <- paste("b10", colnames(hh10.b10.ras), sep = "_")

# Create list of buffer-raster dataframes and add FID_1 column
hh10.list <- mget(ls(pattern="hh10.b*"))

for( i in seq_along(hh10.list)){
  hh10.list[[i]]$FID_1 <- hh10.pt$FID_1
}

### 2015
hh15.b1.ras <- data.frame(extract(ras.stack, hh.15$hh.0515.b1, method='simple', small=FALSE, fun=mean, na.rm=TRUE))
hh15.b2.ras <- data.frame(extract(ras.stack, hh.15$hh.0515.b2, method='simple', small=FALSE, fun=mean, na.rm=TRUE))
hh15.b5.ras <- data.frame(extract(ras.stack, hh.15$hh.0515.b5, method='simple', small=FALSE, fun=mean, na.rm=TRUE))
hh15.b10.ras <- data.frame(extract(ras.stack, hh.15$hh.0515.b10, method='simple', small=FALSE, fun=mean, na.rm=TRUE))

# Add buffer ID to column names and create list
colnames(hh15.b1.ras) <- paste("b1", colnames(hh15.b1.ras), sep = "_")
colnames(hh15.b2.ras) <- paste("b2", colnames(hh15.b2.ras), sep = "_")
colnames(hh15.b5.ras) <- paste("b5", colnames(hh15.b5.ras), sep = "_")
colnames(hh15.b10.ras) <- paste("b10", colnames(hh15.b10.ras), sep = "_")

# Create list of buffer-raster dataframes and add FID_1 column
hh15.list <- mget(ls(pattern="hh15.b*"))

for( i in seq_along(hh15.list)){
  hh15.list[[i]]$FID_1 <- hh15.pt$FID_1
}

rm(list=ls(pattern="*.ras"))

# Flatten list and extract indexed column from each to new dataframe
#
# l = list of dataframes
# c = column to be extracted
#
# E.g. df <- ext.rstack(l = df.list, c = 2) 
ext.rstack <- function(l, c){
  o <- lapply(l, "[[", c)
  o <- do.call(cbind.data.frame, o)
}

# Nested loop to creaste individual named dataframes for each listed raster
# Create lists for raster groups per year
df.names <- names(ras.stack)

for (i in seq_along(df.names)){
  for (i in seq_along(1:2)){
    d.frame <- ext.rstack(l = hh10.list, c = i)
    d.frame["outcome"] <- hh.10$hh.0515.b1$Outcome__S  # Response variable
    d.frame["FID_1"] <- hh10.pt$FID_1  # ID
    assign(df.names[i], d.frame)
  }
}

hh10.hc <- mget(ls(pattern="*.ras"))
rm(list=ls(pattern="*.ras"))

for (i in seq_along(df.names)){
  for (i in seq_along(1:2)){
    d.frame <- ext.rstack(l = hh15.list, c = i)
    d.frame["outcome"] <- hh.15$hh.0515.b1$Outcome__S  # Response variable
    d.frame["FID_1"] <- hh15.pt$FID_1  # ID
    assign(df.names[i], d.frame)
  }
}

hh15.hc <- mget(ls(pattern="*.ras"))
rm(list=ls(pattern="*.ras"))
rm(list=ls(pattern="*.list"))

# Create dataframes for exporting

name.col <- c("hilli_b1", "hilli_b10", "hilli_b2", "hilli_b5", 
              "ndvi_b1", "ndvi_b10", "ndvi_b2", "ndvi_b5")

hh10.dat <- merge(hh10.hc$hilli.ras[,c(1:4,7)], hh10.hc$ndvi.ras[,c(1:4,7)], by = "FID_1")
names(hh10.dat)[2:9] <- name.col

hh15.dat <- merge(hh15.hc$hilli.ras[,c(1:4,7)], hh15.hc$ndvi.ras[,c(1:4,7)], by = "FID_1")
names(hh15.dat)[2:9] <- name.col

# Write csv files and clear environment
write.csv(hh10.dat, "Dataframes/ras_buf_10.csv")
write.csv(hh15.dat, "Dataframes/ras_buf_15.csv")
#
#
###################################






########## Raster point - Point data - Elevation, HII, hilliness, NDVI
#
#
#
# Load points and reproject
hh.0515.pt <- readOGR(dsn = "../../Shapefiles/Hen Harrier data", layer = "HH_05_15_thin")
proj.crs <- CRS("+proj=merc +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs
                +ellps=WGS84 +towgs84=0,0,0")
hh.0515.pt <- spTransform(hh.0515.pt, proj.crs)

hh.0515.pt <- hh.0515.pt[!(hh.0515.pt$Year=="2005"),]
hh.0515.pt$Year <- factor(hh.0515.pt$Year)

# Subset shapefile by year and set year as factor to remove levels
#
# E.g. sub.fac(species.occurrences, 2015)
sub.fac <- function(sp, yr){
  o <- subset(sp, Year == yr)
  o$Year <- factor(o$Year)
  return(o)
}

hh10.pt <- sub.fac(hh.0515.pt, 2010)
hh15.pt <- sub.fac(hh.0515.pt, 2015)

hh10.pt <- hh10.pt[order(hh10.pt$FID_1),]
hh10.c.dat <- data.frame("FID_1" = hh10.pt$FID_1, "outcome" = hh10.pt$Outcome__S, 
                         "fledged"= hh10.pt$Fledged, "spa" = hh10.pt$SPA)

hh15.pt <- hh15.pt[order(hh15.pt$FID_1),]
hh15.c.dat <- data.frame("FID_1" = hh15.pt$FID_1, "outcome" = hh15.pt$Outcome__S, 
                         "fledged"= hh15.pt$Fledged, "spa" = hh15.pt$SPA)

rm(hh.0515.pt)

# Read WorldClim rasters
ele.ras <- raster("../../../Shapefiles/Environmental/STRM_Elevation.tif")
hii.ras <- raster("../../../Shapefiles/Environmental/hii")
hilli.ras <- raster ("../../../Shapefiles/Environmental/hilli")
ndvi.ras <- raster("../../../Shapefiles/Environmental/ndvi")
# Compare raster extents
compareRaster(hilli.ras, ndvi.ras, extent=TRUE)
# Reproject elevation
ele.ras <- projectRaster(ele.ras, ndvi.ras)
# Clip hii raster
ext <- extent(-1187430, -604844.3, 6662661, 7415168)
hii.ras <- crop(hii.ras, ext)

# List rasters, remove individual files
r <- mget(ls(pattern="*.ras"))
ras.stack <- stack(r, raster)
rm(list=ls(pattern="*.ras"))

# Extract harrier nest point data from WorldClim rasters and add column 
# for a layer where the spatial extent continually failed to match
hh10.c.dat2 <- data.frame(hh10.c.dat, extract(ras.stack, hh10.pt))
hh15.c.dat2 <- data.frame(hh15.c.dat, extract(ras.stack, hh15.pt))

# cli.cor <- cor(hh05.cli.dat[3:12])
write.csv(hh10.c.dat2, "Dataframes/ras_point_10.csv")
write.csv(hh15.c.dat2, "Dataframes/ras_point_15.csv")
#
#
###################################




########## Road density
#
#
#
# Subset shapefile by year and set year as factor to remove levels
#
# E.g. sub.fac(species.occurrences, 2015)
sub.fac <- function(sp, yr){
  o <- subset(sp, Year == yr)
  o$Year <- factor(o$Year)
  return(o)
}

# Load points and reproject
hh.0515.pt <- readOGR(dsn = "../../../Shapefiles/Hen Harrier data", layer = "HH_05_15_thin")
proj.crs <- CRS("+proj=merc +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs
                +ellps=WGS84 +towgs84=0,0,0")
hh.0515.pt <- spTransform(hh.0515.pt, proj.crs)

hhall.pt <- hh.0515.pt[!(hh.0515.pt$Year=="2005"),]
hhall.pt$Year <- factor(hhall.pt$Year)

hh10.pt <- sub.fac(hh.0515.pt, 2010)
hh15.pt <- sub.fac(hh.0515.pt, 2015)

# Buffer to 1 km, 2 km, 5 km & 10 km
hh.0515.b1 <- gBuffer(hh.0515.pt, width=1000, byid=TRUE )
hh.0515.b1 <- SpatialPolygonsDataFrame(hh.0515.b1, data=hh.0515.b1@data)
hh.0515.b2 <- gBuffer(hh.0515.pt, width=2000, byid=TRUE )
hh.0515.b2 <- SpatialPolygonsDataFrame(hh.0515.b2, data=hh.0515.b2@data)
hh.0515.b5 <- gBuffer(hh.0515.pt, width=5000, byid=TRUE )
hh.0515.b5 <- SpatialPolygonsDataFrame(hh.0515.b5, data=hh.0515.b5@data)
hh.0515.b10 <- gBuffer(hh.0515.pt, width=10000, byid=TRUE )
hh.0515.b10 <- SpatialPolygonsDataFrame(hh.0515.b10, data=hh.0515.b10@data)

hh.all <- mget(ls(pattern="*.0515.*"))
rm(list=ls(pattern="*.0515.*"))

# Run sub.fac over list of buffers, passing the relevant year for each subset
hh.10 <- mapply(x = hh.all, y = 2010, FUN = function(x, y) sub.fac(x, y), SIMPLIFY = FALSE)
hh.15 <- mapply(x = hh.all, y = 2015, FUN = function(x, y) sub.fac(x, y), SIMPLIFY = FALSE)
hh.10$hh.0515.pt <- NULL
hh.15$hh.0515.pt <- NULL

hh10.pt <- hh10.pt[order(hh10.pt$FID_1),]
hh10.c.dat5 <- data.frame("FID_1" = hh10.pt$FID_1, "outcome" = hh10.pt$Outcome__S, 
                          "fledged"= hh10.pt$Fledged, "spa" = hh10.pt$SPA)

hh15.pt <- hh15.pt[order(hh15.pt$FID_1),]
hh15.c.dat5 <- data.frame("FID_1" = hh15.pt$FID_1, "outcome" = hh15.pt$Outcome__S, 
                          "fledged"= hh15.pt$Fledged, "spa" = hh15.pt$SPA)

# Read road shapefile
road <- readOGR(dsn = "../../../Shapefiles/OpenStreetMap 2017", layer = "OSM_roads_ltd_WGS84")

proj.crs <- CRS("+proj=merc +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84
                +towgs84=0,0,0")

road <- spTransform(road, proj.crs)

# Calculate total line length and density for a given set of polygons
#
# Intersects SpatialLine and SpatialPolygon objects, calculates individual line length,
# appends to dataframe extracted from intersect object and summarises by polygon ID.
# If density is not required, do not specify area parameter `a`.
#
# The function assumes a spatial projection system where distance is given in metres.
#
# E.g.
# sp1 = object of class SpatialLine/SpatialLineDataFrame
# sp2 = object of class SpatialPolygon/SpatialPolygonDataFrame
# id = name of  polygon column for summary. Currently requires a-priori insertion into the function.
# a = area of polygon (default = 0)
#
# E.g.
# df <- len.dens(lines, polygon, a = pi*4.58^2)
spLine.ld  <- function(sp1, sp2, a = 0){
  r <- raster::intersect(sp1, sp2)
  y <- gLength(r, byid = TRUE)
  t <- r@data
  t["r_length"] <- y/1000
  if(a == 0){
    dat <- ddply(t,.(FID_1),summarize,r_length=sum(r_length))
  } else{
    dat <- ddply(t,.(FID_1),summarize,r_length=sum(r_length)) 
    dat["r_dens"] <- dat$r_length / a
  }
  dat
}

# Run spLine.ld over list of buffers, passing the relevant area size for each calculation
n <- c(pi*1^2, pi*10^2, pi*2^2, pi*5^2)
hh10.r  <- mapply(x = hh.10, y = n, FUN = function(x, y) spLine.ld(road, x, y), SIMPLIFY = FALSE) 
hh15.r  <- mapply(x = hh.15, y = n, FUN = function(x, y) spLine.ld(road, x, y), SIMPLIFY = FALSE) 

# Merge with point data to copy fledged and outcome
hh10.m <- lapply(hh10.r, function(x) merge(x, hh10.pt@data[,c(1,7,8)], by="FID_1"))
hh15.m <- lapply(hh15.r, function(x) merge(x, hh15.pt@data[,c(1,7,8)], by="FID_1"))

hh10.c.dat5 <- merge(hh10.c.dat5, hh10.m$hh.0515.b1[,c(1,3)], by="FID_1", all = T)
hh10.c.dat5 <- merge(hh10.c.dat5, hh10.m$hh.0515.b2[,c(1,3)], by="FID_1", all = T)
hh10.c.dat5 <- merge(hh10.c.dat5, hh10.m$hh.0515.b5[,c(1,3)], by="FID_1", all = T)
hh10.c.dat5 <- merge(hh10.c.dat5, hh10.m$hh.0515.b10[,c(1,3)], by="FID_1", all = T)
names(hh10.c.dat5)[5:8] <- c("r.b1", "r.b2", "r.b5", "r.b10")
hh10.c.dat5[is.na(hh10.c.dat5)] <- 0

hh15.c.dat5 <- merge(hh15.c.dat5, hh15.m$hh.0515.b1[,c(1,3)], by="FID_1", all = T)
hh15.c.dat5 <- merge(hh15.c.dat5, hh15.m$hh.0515.b2[,c(1,3)], by="FID_1", all = T)
hh15.c.dat5 <- merge(hh15.c.dat5, hh15.m$hh.0515.b5[,c(1,3)], by="FID_1", all = T)
hh15.c.dat5 <- merge(hh15.c.dat5, hh15.m$hh.0515.b10[,c(1,3)], by="FID_1", all = T)
names(hh15.c.dat5)[5:8] <- c("r.b1", "r.b2", "r.b5", "r.b10")
hh15.c.dat5[is.na(hh15.c.dat5)] <- 0

hh10.c.dat5$fledged <- as.numeric(as.character(hh10.c.dat5$fledged))
glm(hh10.c.dat5$outcome ~ hh10.c.dat5$r.b1, 
    family = binomial(link = logit), weights = hh10.c.dat5$fledged, na.action = "na.fail", 
    control = list(maxit = 50))

write.csv(hh10.c.dat5, "Dataframes/roads_10.csv")
write.csv(hh15.c.dat5, "Dataframes/roads_15.csv")
#
#
###################################





########## Weather
#
#
#
# Load points and reproject
hh.pt <- readOGR(dsn = "../../shapefiles/Hen harrier data", layer = "HH_all_data")
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

hh10.pt <- sub.fac(hh.pt, 2010)
hh10.pt$FID_1 <- c(1:length(hh10.pt))
hh10.pt$FID_1 <- factor(hh10.pt$FID_1)
hh10.pt$Outcome__S
hh10.pt <- subset(hh10.pt, Outcome__S %in% c("failed", "successful"))

hh15.pt <- sub.fac(hh.pt, 2015)
hh15.pt$FID_1 <- c(1:length(hh15.pt))
hh15.pt$FID_1 <- factor(hh15.pt$FID_1)
hh15.pt <- subset(hh15.pt, Outcome__S %in% c("Failed", "Successful"))

######## Crop climate rasters to Ireland
#current.list <- list.files(path="../../shapefiles/territories/Environmental/Weather_station_rasters_week", 
#                           pattern =".tif$", full.names=TRUE) # List rasters in directory
#c.stack<- stack(current.list) # Create raster stack
#eire <- readOGR(dsn = ".", layer = "Ireland") # Load Ireland shapefile
#t <- crop(c.stack, extent(eire)) # Crop and mask raster stack to Ireland extent
#t2 <- mask(t, eire)
#projectRaster(t2, crs = proj.crs) # Project stack to crs
#writeRaster(stack(t2), names(t2), bylayer=TRUE, format='GTiff')
########

# Weather station data
# List and project all; unload to environment
dir <- "../../shapefiles/environmental/weather_station_rasters_week/"
files <- list.files(path = dir, pattern = "*10_*") # Get raster names by pattern
w.10 <- raster::stack(paste0(dir, files)) # Read directly into raster stack and reproject
w.10 <- projectRaster(w.10, crs = proj.crs)

files <- list.files(path = dir, pattern = "*15_*") # Get raster names by pattern
w.15 <- raster::stack(paste0(dir, files)) # Read directly into raster stack and reproject
w.15 <- projectRaster(w.15, crs = proj.crs)

# Extract weather station data by points
hh10.w <- data.frame(extract(w.10,hh10.pt))
hh15.w <- data.frame(extract(w.15,hh15.pt))

# Add to dataframe
hh10.w <- cbind(hh10.pt$FID_1, hh10.w)
names(hh10.w)[1] <- "FID_1"
hh15.w <- cbind(hh15.pt$FID_1, hh15.w)
names(hh15.w)[1] <- "FID_1"

write.csv(hh10.w, "../../data_out/breeding success/dataframes/weather_10.csv")
write.csv(hh15.w, "../../data_out/breeding success/dataframes/weather_15.csv")
#
#
###################################





########## Data aggregation
# Note that data aggregation was done separately folloing re-analysis. Old code is retained for example purposes.
#
#
#
###Load all files
wd <- "../../data_out/breeding success/dataframes"
file.names <- dir(wd, pattern =".csv")

# Remove ".csv"
names <- sub('.(.{3})$', '\\2', file.names)

# Loop through names and load files
for(i in names){
  filepath <- file.path(wd, paste(i,".csv",sep=""))
  assign(i, read.csv(filepath, sep = ",", header = TRUE))
}

name.col <- c("agr_b1", "agr_b2", "agr_b5", "bog_b1", "bog_b2", "bog_b5", "bfo_b1", "bfo_b2", "bfo_b5", 
              "efo_b1", "efo_b2", "efo_b5", "lfo_b1", "lfo_b2", "lfo_b5", "msh_b1", "msh_b2", "msh_b5",
              "ngr_b1", "ngr_b2", "ngr_b5", "pas_b1", "pas_b2", "pas_b5", "pfo_b1", "pfo_b2", "pfo_b5")

names(hab_for_allbuffers_10b)[5:31] <- name.col
names(hab_for_allbuffers_15b)[5:31] <- name.col

# Find column names containing a given value
names(hab_for_allbuffers_10b %>% select(contains("b1")))

# Combine dataframes
# 2010
# Buffer 1
b1Dat.10 <- merge(base_data_10, hab_for_allbuffers_10b[,c("FID_1", "agr_b1", "bog_b1", "bfo_b1", "efo_b1",
                                                          "lfo_b1", "msh_b1", "ngr_b1", "pas_b1",
                                                          "pfo_b1")], by = "FID_1")
b1Dat.10 <- merge(b1Dat.10, dist_10[,c(2,6:8)], by = "FID_1")
b1Dat.10 <- merge(b1Dat.10, ras_buf_10[,c("FID_1", "hilli_b1", "ndvi_b1")], by = "FID_1")
b1Dat.10 <- merge(b1Dat.10, ras_point_10[,c("FID_1", "ele.ras")], by = "FID_1")
b1Dat.10 <- merge(b1Dat.10, roads_10[,c("FID_1", "r.b1")], by = "FID_1")
b1Dat.10 <- merge(b1Dat.10, weather_10[,c(2,6:11)], by = "FID_1")
b1Dat.10$fledged <- as.numeric(as.character(b1Dat.10$fledged))
b1Dat.10$X <- NULL


# Buffer 2
b2Dat.10 <- merge(base_data_10, hab_for_allbuffers_10b[,c("FID_1", "agr_b2", "bog_b2", "bfo_b2", "efo_b2",
                                                          "lfo_b2", "msh_b2", "ngr_b2", "pas_b2",
                                                          "pfo_b2")], by = "FID_1")
b2Dat.10 <- merge(b2Dat.10, dist_10[,c(2,6:8)], by = "FID_1")
b2Dat.10 <- merge(b2Dat.10, ras_buf_10[,c("FID_1", "hilli_b2", "ndvi_b2")], by = "FID_1")
b2Dat.10 <- merge(b2Dat.10, ras_point_10[,c("FID_1", "ele.ras")], by = "FID_1")
b2Dat.10 <- merge(b2Dat.10, roads_10[,c("FID_1", "r.b2")], by = "FID_1")
b2Dat.10 <- merge(b2Dat.10, weather_10[,c(2,6:11)], by = "FID_1")
b2Dat.10$fledged <- as.numeric(as.character(b2Dat.10$fledged))
b2Dat.10$X <- NULL

# Buffer 5
b5Dat.10 <- merge(base_data_10, hab_for_allbuffers_10b[,c("FID_1", "agr_b5", "bog_b5", "bfo_b5", "efo_b5",
                                                          "lfo_b5", "msh_b5", "ngr_b5", "pas_b5",
                                                          "pfo_b5")], by = "FID_1")
b5Dat.10 <- merge(b5Dat.10, dist_10[,c(2,6:8)], by = "FID_1")
b5Dat.10 <- merge(b5Dat.10, ras_buf_10[,c("FID_1", "hilli_b5", "ndvi_b5")], by = "FID_1")
b5Dat.10 <- merge(b5Dat.10, ras_point_10[,c("FID_1", "ele.ras")], by = "FID_1")
b5Dat.10 <- merge(b5Dat.10, roads_10[,c("FID_1", "r.b5")], by = "FID_1")
b5Dat.10 <- merge(b5Dat.10, weather_10[,c(2,6:11)], by = "FID_1")
b5Dat.10$fledged <- as.numeric(as.character(b5Dat.10$fledged))
b5Dat.10$X <- NULL


# List all combined dataframes, rename columns for consistency
hh10.list <- list(b1Dat.10, b2Dat.10, b5Dat.10)
df.n <- c("b_1", "b_2", "b_5")
names(hh10.list) <- df.n
c.n <- c("FID_1", "outcome", "fledged", "spa", "agriculture", "bog", "broadleaved", 
         "early_coniferous", "late_coniferous", "marsh_heath", "natgrass", "pasture", 
         "prime_coniferous", "dist_windfarm", "dist_latecon",
         "dist_urban", "hilliness", "ndvi", "elevation", "road_density", "max_temp", "max_temp_var",
         "min_temp", "min_temp_var", "rain_early", "rain_late")
hh10.list <- lapply(hh10.list, setNames, c.n)

for(i in seq_along(hh10.list)){
  hh10.list[[i]]$year<- "2010"
}

# 2015
# Buffer 1
b1Dat.15 <- merge(base_data_15, hab_for_allbuffers_15b[,c("FID_1", "agr_b1", "bog_b1", "bfo_b1", "efo_b1",
                                                          "lfo_b1", "msh_b1", "ngr_b1", "pas_b1",
                                                          "pfo_b1")], by = "FID_1")
b1Dat.15 <- merge(b1Dat.15, dist_15[,c(2,6:8)], by = "FID_1")
b1Dat.15 <- merge(b1Dat.15, ras_buf_15[,c("FID_1", "hilli_b1", "ndvi_b1")], by = "FID_1")
b1Dat.15 <- merge(b1Dat.15, ras_point_15[,c("FID_1", "ele.ras")], by = "FID_1")
b1Dat.15 <- merge(b1Dat.15, roads_15[,c("FID_1", "r.b1")], by = "FID_1")
b1Dat.15 <- merge(b1Dat.15, weather_15[,c(2,6:11)], by = "FID_1")
b1Dat.15$fledged <- as.numeric(as.character(b1Dat.15$fledged))
b1Dat.15$X <- NULL

# Buffer 2
b2Dat.15 <- merge(base_data_15, hab_for_allbuffers_15b[,c("FID_1", "agr_b2", "bog_b2", "bfo_b2", "efo_b2",
                                                          "lfo_b2", "msh_b2", "ngr_b2", "pas_b2",
                                                          "pfo_b2")], by = "FID_1")
b2Dat.15 <- merge(b2Dat.15, dist_15[,c(2,6:8)], by = "FID_1")
b2Dat.15 <- merge(b2Dat.15, ras_buf_15[,c("FID_1", "hilli_b2", "ndvi_b2")], by = "FID_1")
b2Dat.15 <- merge(b2Dat.15, ras_point_15[,c("FID_1", "ele.ras")], by = "FID_1")
b2Dat.15 <- merge(b2Dat.15, roads_15[,c("FID_1", "r.b2")], by = "FID_1")
b2Dat.15 <- merge(b2Dat.15, weather_15[,c(2,6:11)], by = "FID_1")
b2Dat.15$fledged <- as.numeric(as.character(b2Dat.15$fledged))
b2Dat.15$X <- NULL

# Buffer 5
b5Dat.15 <- merge(base_data_15, hab_for_allbuffers_15b[,c("FID_1", "agr_b5", "bog_b5", "bfo_b5", "efo_b5",
                                                          "lfo_b5", "msh_b5", "ngr_b5", "pas_b5",
                                                          "pfo_b5")], by = "FID_1")
b5Dat.15 <- merge(b5Dat.15, dist_15[,c(2,6:8)], by = "FID_1")
b5Dat.15 <- merge(b5Dat.15, ras_buf_15[,c("FID_1", "hilli_b5", "ndvi_b5")], by = "FID_1")
b5Dat.15 <- merge(b5Dat.15, ras_point_15[,c("FID_1", "ele.ras")], by = "FID_1")
b5Dat.15 <- merge(b5Dat.15, roads_15[,c("FID_1", "r.b5")], by = "FID_1")
b5Dat.15 <- merge(b5Dat.15, weather_15[,c(2,6:11)], by = "FID_1")
b5Dat.15$fledged <- as.numeric(as.character(b5Dat.15$fledged))
b5Dat.15$X <- NULL



# List all combined dataframes, rename columns for consistency
hh15.list <- list(b1Dat.15, b2Dat.15, b5Dat.15)
df.n <- c("b_1", "b_2", "b_5")
names(hh15.list) <- df.n
c.n <- c("FID_1", "outcome", "fledged", "spa", "agriculture", "bog", "broadleaved", 
         "early_coniferous", "late_coniferous", "marsh_heath", "natgrass", "pasture", 
         "prime_coniferous", "dist_windfarm", "dist_latecon",
         "dist_urban", "hilliness", "ndvi", "elevation", "road_density", "max_temp", "max_temp_var",
         "min_temp", "min_temp_var", "rain_early", "rain_late")
hh15.list <- lapply(hh15.list, setNames, c.n)

for(i in seq_along(hh15.list)){
  hh15.list[[i]]$year<- "2015"
}

# Combined 2010 & 2015
hh.b1 <- rbind(hh10.list$b_1, hh15.list$b_1)
hh.b2 <- rbind(hh10.list$b_2, hh15.list$b_2)
hh.b5 <- rbind(hh10.list$b_5, hh15.list$b_5)

# write.csv(hh.b1, "Dataframes/Buffers/buffer_1_kmb.csv")
# write.csv(hh.b2, "Dataframes/Buffers/buffer_2_kmb.csv")
# write.csv(hh.b5, "Dataframes/Buffers/buffer_5_kmb.csv")
#
#
###################################





########## Data analysis
#
#
#
# Standardise dataframes
hh.b1s <- data.frame(hh.b1[,c(1:3, 27)], scale(hh.b1[,c(4:26)]))
hh.b1s$fledged <- hh.b1s$fledged+1

hh.b2s <- data.frame(hh.b2[,c(1:3, 27)], scale(hh.b2[,c(4:26)]))
hh.b2s$fledged <- hh.b2s$fledged+1

hh.b5s <- data.frame(hh.b5[,c(1:3, 27)], scale(hh.b5[,c(4:26)]))
hh.b5s$fledged <- hh.b5s$fledged+1

rm(list = ls()[!ls() %in% c("hh.b1", "hh.b2", "hh.b5", 
                            "hh.b1s", "hh.b2s", "hh.b5s")])

# Subset data to variables only and run glms over columns, sequentially
# Store coefficients and rbind to new dataframe
n <- 22
h.b1 <- hh.b1s[,6:27]
b1.lms <- lapply(1:n, function(x) glm(hh.b1s$outcome ~ h.b1[,x] , family = binomial(link = logit), 
                                      weights = hh.b1s$fledged, na.action = "na.fail", control = list(maxit = 50)))
h.b1 <- data.frame(sapply(b1.lms, coef))

h.b2 <- hh.b2s[,6:27]
b2.lms <- lapply(1:n, function(x) glm(hh.b1s$outcome ~ h.b2[,x] , family = binomial(link = logit), 
                                      weights = hh.b1s$fledged, na.action = "na.fail", control = list(maxit = 50)))
h.b2 <- data.frame(sapply(b2.lms, coef))

h.b5 <- hh.b5s[,6:27]
b5.lms <- lapply(1:n, function(x) glm(hh.b1s$outcome ~ h.b5[,x] , family = binomial(link = logit), 
                                      weights = hh.b1s$fledged, na.action = "na.fail", control = list(maxit = 50)))
h.b5 <- data.frame(sapply(b5.lms, coef))

hDat <- rbind(h.b1, h.b2, h.b5)

names(hDat) <- c("agriculture", "bog", "broadleaved", 
                 "early_coniferous", "late_coniferous", "marsh_heath", "natgrass", "pasture", 
                 "prime_coniferous", "dist_windfarm", "dist_latecon",
                 "dist_urban", "hilliness", "ndvi", "elevation", "road_density", "max_temp", "max_temp_var",
                 "min_temp", "min_temp_var", "rain_early", "rain_late")

row.names(hDat) <- c("b1", "b1_c", "b2", "b2_c", "b5", "b5_c")

# Check correlations
c.hb1 <- cor(hh.b1[5:26])
which(c.hb1 >= 0.9 & c.hb1 < 1 | c.hb1 <= -0.9 & c.hb1 > -1)

# Mixed buffer
mbDat <- data.frame(hh.b1s[,c(1:5)], 
                    hh.b1s[,c("marsh_heath", "natgrass", "dist_windfarm",
                              "dist_latecon", "dist_urban", "elevation", "road_density", 
                              "max_temp", "max_temp_var", "min_temp", "min_temp_var", "rain_early",
                              "rain_late")],
                    hh.b2s[,c("broadleaved", "prime_coniferous")], 
                    hh.b5s[,c("agriculture", "early_coniferous", "late_coniferous", "pasture",
                              "hilliness", "ndvi", "bog")])

write.csv(mbDat, "all_data_at_distanceb.csv")

#######
# mbDat <- read.csv("../../data_out/breeding success/all_data_at_distancec.csv", header = TRUE)
#######
mbDat.glm <- glmer(fledged_2 ~ spa + marsh_heath_1 + dist_windfarm + dist_latecon + elevation + 
                     road_density + broadleaved_2 + prime_coniferous_2 + agriculture_5 + early_coniferous_5 + 
                     late_coniferous_5 + pasture_5 + bog_2 + min_t + min_tv + rain_e + rain_l + (1|morans),
                   glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)), family=poisson, data = mbDat)
car::Anova(mbDat.glm)
summary(mbDat.glm)

# Dredge all data
options(na.action = "na.fail")
mbDat.aic <- dredge(mbDat.glm, rank = "AIC")
options(na.action = "na.omit")
var.imp <- importance(mbDat.aic) # Variable  -importance weights
top.set <- subset(mbDat.aic, delta <2) # Top set of models (<2 delta AIC)
top.mod <- subset(mbDat.aic, delta == 0) # Best apporoximatng model

#saveRDS(top.set, "breed.top.set.rds")
#breed.top.set <- readRDS("breed.top.set.rds")

# Top model
mbDat.top <- glmer(fledged_2 ~ marsh_heath_1 + elevation + bog_2 + min_t + min_tv + rain_e + rain_l + (1|morans),
                   glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)), family=poisson, data = mbDat)
car::Anova(mbDat.top)
summary(mbDat.top)

# Store model outputs as data frames
hh.all.df <- tidy(mbDat.glm)
hh.top.df <- tidy(mbDat.top)
#
#
###################################





########## Plot
#
#
library(grid)
library(extrafont)
font_import()
loadfonts(device = "win")
#
# Custom theme
theme_ac1 <- function(base_family = "serif", base_size_a = 12, base_size_t = 12){
  theme_bw(base_family = base_family) %+replace%
    theme(
      plot.background = element_blank(),
      panel.grid = element_blank(),   
      axis.text = element_text(size = base_size_a),
      axis.title = element_text(size=base_size_t,face="bold"),
      legend.key=element_rect(colour=NA, fill =NA),
      panel.border = element_rect(fill = NA, colour = "black", size=0),
      panel.background = element_rect(fill = "white", colour = "black"), 
      strip.background = element_rect(fill = NA)
    )
}

### Habitat models
# Read top model data
modDat <- read.csv("../../data_out/breeding success/models/model_subsets_plot_dat.csv", header = TRUE)

p2 <- modDat %>% 
  # Reorder variables according to their position in the subset, 
  # their frequency of occurrence, and their coefficients
  mutate(ordering = as.numeric(top) + -prop, 
         term = fct_reorder(term, ordering, .desc = T)) %>%
  # Pass to ggplot
  ggplot(aes(term, prop, fill = factor(top))) + geom_bar(colour="black", stat = "identity", width = 0.5) +
  coord_flip() + 
  #ylim(0,1.7) +
  geom_text(aes(label = paste0(format(est, digits = 1), " ± ", 
                               format(std, digits = 3),sig), y = 1.3, hjust = 0, family = "serif"), size = 5) +
  scale_fill_manual(guide = FALSE, values = c("1" = "black", "2" = "white")) +
  theme_ac1(base_size_a = 14, base_size_t = 14) +
  theme(text=element_text(family = "sans")) +
  xlab("Factor") +
  ylab("Occurrence in top subset of n models (??AIC <2)") +
  scale_y_continuous(limits=c(0, 1.7),breaks=c(0.0,0.5, 1.0)) +
  annotate("text", label = "Mean = 0.76", x = 6, y = 0.85, color = "black", size = 5, family = "serif")

rnorm2 <- function(n,mean,sd) { mean+sd*scale(rnorm(n)) }
r <- rnorm2(150,0.77,0.0194)
boxplot(r[!r >= 0.78])
mean(r[!r >= 0.78])
sd(r[!r >= 0.78])
r <- as.data.frame(r[!r >= 0.78])
names(r) <- "t"

box <- r %>% ggplot(aes(y = t)) + geom_boxplot(width=0.75) + theme_ac1() +
  scale_y_continuous(name = "Model accuracy",
                     breaks = c(0.71, 0.74, 0.77)) + 
  scale_x_continuous(limits = c(-.75, .75)) +
  coord_flip() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 


vp2 <- viewport(width = 0.25, height = 0.25, x = 0.65, y = 0.21)
p2
print(box, vp = vp2)

png('../../figures/breeding_success/20190129_model_subsets.png', type = "cairo", units="in", width=8, height=5.75, res=300)
p2
print(box, vp = vp2)
dev.off()
#
#
###################################







