# Proportional land cover and total number of fledged chicks in each Hen Harrier SPA. See Fig. 5 in Caravaggi et al. (in press),  Factors influencing Hen Harrier, Circus cyaneus, territory site selection and breeding success. Bird Study.

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(raster)
library(rgdal)
library(rgeos)
library(plyr) 
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(gtable)

#install.packages(c("ggplot2", "ggthemes", "ggpubr", "gtable")) 

# Read habitat polygons
agr.sf <- readOGR(dsn = "Shapefiles", layer = "agriculture")
agr.sf@data["var"] <- "agr"
bog.sf <- readOGR(dsn = "Shapefiles", layer = "bog")
bog.sf@data["var"] <- "bog"
msh.sf <- readOGR(dsn = "Shapefiles", layer = "moor_shrub")
msh.sf@data["var"] <- "msh"
ngr.sf <- readOGR(dsn = "Shapefiles", layer = "natgrass")
ngr.sf@data["var"] <- "ngr"
pas.sf <- readOGR(dsn = "Shapefiles", layer = "pasture_erased")
pas.sf@data["var"] <- "pas"
efor.sf <- readOGR(dsn = "Shapefiles", layer = "early_15")
efor.sf@data["var"] <- "efor"
pfor.sf <- readOGR(dsn = "Shapefiles", layer = "prime_15")
pfor.sf@data["var"] <- "pfor"
lfor.sf <- readOGR(dsn = "Shapefiles", layer = "late_15")
lfor.sf@data["var"] <- "lfor"
bfor.sf <- readOGR(dsn = "Shapefiles", layer = "broadleaved")
bfor.sf@data["var"] <- "bfor"

hab <- mget(ls(pattern="*.sf"))

dat <- lapply(hab, function(x) as(x[, c("FID_1", "Site_Name", "size", "var")], "data.frame")) # Extract required columns
hDat <- do.call(rbind.data.frame, dat)
names(hDat) <- c("FID", "name", "area", "var") # Rename columns
hDat <- hDat[!(hDat$name == "Lough Derg (Donegal) SPA"),] # Remove Lough Derg SPA
hDat$name <- factor(hDat$name)

# Create SPA reference datafreame for renaming and merge with above
ref <- data.frame(name = levels(hDat$name),
                  loc = c("MMM", "SAM", "SBe", "SBM", "SSM", "SMW"))
hDat2 <- merge(hDat, ref, by = "name")

hDat$var <- as.factor(hDat$var)
ref2 <- data.frame(var = levels(hDat$var),
                   var2 = c("Arable", "Broadleaved", "Bog", "Conifer (0-2 years)", 
                            "Conifer (13+ years)", "Heath/shrub", "Natural grassland", "Pasture", 
                            "Conifer (3-12 years)"))

hDat3 <- merge(hDat2, ref2, by = "var")
hDat3 <- hDat3[!(hDat3$var2 =="Natural grassland"),]

# write.csv(hDat3, "../../data_out/spa/spa_hab.csv")
hDat3 <- read.csv("../../data_out/spa/spa_hab.csv", header = TRUE)

# Sum area according to each SPA and variable, then calculate relative proportions within each SPA
aDat <- aggregate(area ~ var2 + loc, hDat3, sum, na.rm = TRUE, na.action="na.pass")
aDat <- ddply(aDat, .(loc), transform, prop=area/sum(area))
aDat2 <- dcast(aDat, var2 ~ loc, value.var = "prop")
barplot(as.matrix(aDat2))

# write.csv(aDat, "../../data_out/spa/spa_aggregated.csv")

# Extract rows by condition
dplyr::filter(aDat, grepl('Heath', var2))

# Stacked bar plot of data
library(extrafont)
font_import()
loadfonts(device = "win")
# Load theme
theme_ac1 <- function(base_family = "serif", base_size_a = 18, base_size_t = 18){
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

# Create colour-blind palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
p1 <- ggplot() + 
  geom_bar(aes(y = prop, x = loc, fill = var2), data = aDat, stat="identity") +
  theme_ac1() + 
  theme(text=element_text(family = "sans")) +
  xlab("Special Protection Area (SPA)") +
  ylab("Propotional abundance of land class") +
  scale_fill_manual(values=cbbPalette, name = "Habitat") +
  theme(legend.title=element_text(size=18, face = "bold"), 
        legend.text=element_text(size=16)) 

# Save plot
#png('Figures/SPA_habitats.png', type = "cairo", units="in", width=8, height=5.75, res=300)
#p1
#dev.off()

# Fledged chicks per SPA
bDat <- read.csv("../../data_out/spa/breeding.csv", header = TRUE)

p2 <- ggplot(bDat, aes(y = fledged, x = code, group = year)) + 
  geom_bar(aes(fill = factor(year)), width = 0.4, position = position_dodge(width=0.5), stat="identity") +
  theme_ac1() + 
  theme(text=element_text(family = "sans")) +
  scale_fill_manual(values = c("black", "darkgrey"), name = "Year") +
  xlab("Special Protection Area (SPA)") +
  ylab("Number of fledged chicks") +
  theme(legend.title=element_text(size=18, face = "bold"), 
        legend.text=element_text(size=16))


png('../../figures/spa/20190129_spa.png', type = "cairo", units="in", width=8, height=10, res=300)
ggarrange(p1, p2,
          labels = c("a)", "b)"),
          ncol = 1, nrow = 2,
          align = "v",
          heights = c(2, 1.25))
dev.off()
png('../../figures/spa/20190129_breeding_success.png', type = "cairo", units="in", width=8, height=4, res=300)
p2
dev.off()
