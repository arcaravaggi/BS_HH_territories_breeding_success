library(raster)
library(sp)
library(rgdal)
library(rgeos)
library(maptools)
library(plyr)

proj.crs <- CRS("+proj=merc +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs
                +ellps=WGS84 +towgs84=0,0,0")

pfor <- readOGR(dsn = "../../../Shapefiles/Forests/PrivateForests2016", 
                      layer = "PrivateForests2016_WGS84") # PrivateForests2016

colnames(pfor@data)

# Adjust age accounting for 2 year fallow post-felling
head(pfor@data$X2015_yr)
pfor@data$X2015_yr <- pfor@data$X2015_yr - 2
head(pfor@data$X2015_yr)

head(pfor@data$X2010_yr)
pfor@data$X2010_yr <- pfor@data$X2010_yr - 2
head(pfor@data$X2010_yr)

# Levels and labels for stand age grouping
levels <- c(-Inf, 2, 12, Inf)
labels <- c("1", "2", "3")

# Split between broadleaved and coniferous
bl.pfor <- pfor[ which(pfor@data$sp_grp  =='1'), ]
writeOGR(bl.pfor, "C:/Users/Anthony Caravaggi/Dropbox/Academic/SHINE project/Shapefiles/Forests/temp2",
         "PrivateForests2016_broadleaved", driver="ESRI Shapefile")

co.pfor <- pfor[ which(pfor@data$sp_grp  =='2'), ]
writeOGR(co.pfor, "C:/Users/Anthony Caravaggi/Dropbox/Academic/SHINE project/Shapefiles/Forests/temp2",
         "PrivateForests2016_coniferous", driver="ESRI Shapefile")

# Split data by stage (1 = 1-2 yrs, 2 = 3-12 yrs, 3 = 13+ yrs) - conifers only
# 2015
pfor.15 <- co.pfor[,c(35)]

# Group data
pfor.15@data$grp15 <- cut(pfor.15@data$X2015_yr, levels, labels)

# Remove negative forest ages
pfor.15 <- pfor.15[!(pfor.15@data$X2015_yr < 0),]

# Extract according to group
e.pf15.sf <- pfor.15[ which(pfor.15@data$grp15=='1'), ]
p.pf15.sf <- pfor.15[ which(pfor.15@data$grp15=='2'), ]
l.pf15.sf <- pfor.15[ which(pfor.15@data$grp15=='3'), ]

# 2010
pfor.10 <- co.pfor[,c(36)]

# Group data
pfor.10@data$grp10 <- cut(pfor.10@data$X2010_yr, levels, labels)

# Remove negative forest ages
pfor.10 <- pfor.10[!(pfor.10@data$X2010_yr < 0),]

# Extract according to group
e.pf10.sf <- pfor.10[ which(pfor.10@data$grp10=='1'), ]
p.pf10.sf <- pfor.10[ which(pfor.10@data$grp10=='2'), ]
l.pf10.sf <- pfor.10[ which(pfor.10@data$grp10=='3'), ]


# Coillte forest data

cfor <- readOGR(dsn = "../../../Shapefiles/Forests/Coillte_Forest_Data", 
                layer = "Coillte_Sub_Compartments_WGS84") # Coillte Forest Data

colnames(cfor@data)

# Adjust age accounting for 2 year fallow post-felling
head(cfor@data$X2015_yr)
cfor@data$X2015_yr <- cfor@data$X2015_yr - 2
head(cfor@data$X2015_yr)

head(cfor@data$X2010_yr)
cfor@data$X2010_yr <- cfor@data$X2010_yr - 2
head(cfor@data$X2010_yr)

bl.cfor <- cfor[ which(cfor@data$sp_grp  =='1'), ]
writeOGR(bl.cfor, "C:/Users/Anthony Caravaggi/Dropbox/Academic/SHINE project/Shapefiles/Forests/temp2",
         "CoillteForests_broadleaved", driver="ESRI Shapefile")

co.cfor <- cfor[ which(cfor@data$sp_grp  =='2'), ]
writeOGR(co.cfor, "C:/Users/Anthony Caravaggi/Dropbox/Academic/SHINE project/Shapefiles/Forests/temp2",
         "CoillteForests_coniferous", driver="ESRI Shapefile")

# 2015
cfor.15 <- co.cfor[,c(9)]

# Group data
cfor.15@data$grp15 <- cut(cfor.15@data$X2015_yr, levels, labels)

# Remove negative forest ages
cfor.15 <- cfor.15[!(cfor.15@data$X2015_yr < 0),]

# Extract according to group
e.cf15.sf <- cfor.15[ which(cfor.15@data$grp15=='1'), ]
p.cf15.sf <- cfor.15[ which(cfor.15@data$grp15=='2'), ]
l.cf15.sf <- cfor.15[ which(cfor.15@data$grp15=='3'), ]

# 2010
cfor.10 <- co.cfor[,c(10)]

# Group data
cfor.10@data$grp10 <- cut(cfor.10@data$X2010_yr, levels, labels)

# Remove negative forest ages
cfor.10 <- cfor.10[!(cfor.10@data$X2010_yr < 0),]

# Extract according to group
e.cf10.sf <- cfor.10[ which(cfor.10@data$grp10=='1'), ]
p.cf10.sf <- cfor.10[ which(cfor.10@data$grp10=='2'), ]
l.cf10.sf <- cfor.10[ which(cfor.10@data$grp10=='3'), ]



# NI forest data
nfor <- readOGR(dsn = "../../../Shapefiles/Forests/NI Forestry", 
                layer = "forestServiceSubCompartmentBoundaries2017_WGS84") # NI Forest Data
nfor <- spTransform(nfor, proj.crs)

colnames(nfor@data)

# Group by species
sp <- c("Broadleaf", "Mixed Broadleaf", "Mixed Conifer/Broadleaf", "Conifer", "Mixed Conifer")
grp <- c("1", "1", "1", "2", "2")
nfor@data$sp_grp <- mapvalues(nfor@data$PRIMARY, from = sp, to = grp)

# Adjust age accounting for 2 year fallow post-felling
head(nfor@data$X2015_yr)
nfor@data$X2015_yr <- nfor@data$X2015_yr - 2
head(nfor@data$X2015_yr)

head(nfor@data$X2010_yr)
nfor@data$X2010_yr <- nfor@data$X2010_yr - 2
head(nfor@data$X2010_yr)

bl.nfor <- nfor[ which(nfor@data$sp_grp  =='1'), ]
writeOGR(bl.nfor, "C:/Users/Anthony Caravaggi/Dropbox/Academic/SHINE project/Shapefiles/Forests/temp2",
         "NIForests_broadleaved", driver="ESRI Shapefile")

co.nfor <- nfor[ which(nfor@data$sp_grp  =='2'), ]
writeOGR(co.nfor, "C:/Users/Anthony Caravaggi/Dropbox/Academic/SHINE project/Shapefiles/Forests/temp2",
         "NIForests_coniferous", driver="ESRI Shapefile")

# 2015
nfor.15 <- co.nfor[,c(8)]

# Group data
nfor.15@data$grp15 <- cut(nfor.15@data$X2015_yr, levels, labels)

# Remove negative forest ages
nfor.15 <- nfor.15[!(nfor.15@data$X2015_yr < 0),]

# Extract according to group
e.nf15.sf <- nfor.15[ which(nfor.15@data$grp15=='1'), ]
p.nf15.sf <- nfor.15[ which(nfor.15@data$grp15=='2'), ]
l.nf15.sf <- nfor.15[ which(nfor.15@data$grp15=='3'), ]

# 2010
nfor.10 <- co.nfor[,c(9)]

# Group data
nfor.10@data$grp10 <- cut(nfor.10@data$X2010_yr, levels, labels)

# Remove negative forest ages
nfor.10 <- nfor.10[!(nfor.10@data$X2010_yr < 0),]

# Extract according to group
e.nf10.sf <- nfor.10[ which(nfor.10@data$grp10=='1'), ]
p.nf10.sf <- nfor.10[ which(nfor.10@data$grp10=='2'), ]
l.nf10.sf <- nfor.10[ which(nfor.10@data$grp10=='3'), ]

# Remove unwanted objects
rm(nfor, co.nfor, nfor.10, nfor.15, pfor, co.pfor, pfor.10, pfor.15, cfor, co.cfor, cfor.10, cfor.15, labels, levels, grp, sp)

# Merge early, prime, and late forestry for each year, 
# fixing for duplicate polygon IDs
# 
# shp.lst = list of SpatialPolygon objects
multi.merge <- function(shp.lst){
  p <- lapply(shp.lst, function(x) spChFIDs(x, paste(deparse(substitute(x)), 
                                                     row.names(x), sep = ".")))
  t <- lapply(p, function(x) spTransform(x, proj.crs))
  w <- do.call(bind, t)
  return(w)
}

# Build list of early, prime and late coniferous forests, list lists and apply multi.merge across all
#2010
e.10 <- list(e.pf10.sf, e.nf10.sf, e.cf10.sf)
p.10 <- list(p.pf10.sf, p.nf10.sf, p.cf10.sf)
l.10 <- list(l.pf10.sf, l.nf10.sf, l.cf10.sf)
f.10 <- list(e.10, p.10, l.10)
names(f.10) <- c("e.10", "p.10", "l.10")

con.10 <- lapply(f.10, function(x) multi.merge(x))
names(con.10) <- c("early_10", "prime_10", "late_10")

#2015
e.15 <- list(e.pf15.sf, e.nf15.sf, e.cf15.sf)
p.15 <- list(p.pf15.sf, p.nf15.sf, p.cf15.sf)
l.15 <- list(l.pf15.sf, l.nf15.sf, l.cf15.sf)
f.15 <- list(e.15, p.15, l.15)

con.15 <- lapply(f.15, function(x) multi.merge(x))
names(con.15) <- c("early_15", "prime_15", "late_15")


# Remove individual conifer files
rm(list=ls(pattern="*.sf"))

# Write listed shapefiles to files
mapply(x = con.15, y = c("early_15", "prime_15", "late_15")
       , FUN =  function(x, y)  writeOGR(x, "../../../Shapefiles/Forests/temp2", 
                                         y, driver="ESRI Shapefile"))

mapply(x = con.10, y = c("early_10", "prime_10", "late_10")
       , FUN =  function(x, y)  writeOGR(x, "../../../Shapefiles/Forests/temp2", 
                                         y, driver="ESRI Shapefile"))

# Erasing of early and prime shapefiles from CORINE and merging broadleaved and late shapefiles with CORINE
# were carried out in ArcGIS due to recurrent issues with R
