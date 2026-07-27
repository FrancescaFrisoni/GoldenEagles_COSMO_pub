# This code is used to build background points to compare used vs available uplift types
# Martina Scacco - May 8, 2025. Konstanz

library(terra)
library(tidyr)
library(data.table)
library(mapview)
library(lubridate)

setwd("~/ownCloud - mscacco@ab.mpg.de@owncloud.gwdg.de/Martina/ProgettiVari/GoldenEagles_WindMapsFromThermals/COSMO/TesiFrancesca_UpliftClassification/")



#_________________
# Individual MCPs
# Here we use a simple approach. We draw a convex hull for each individual around the soaring points of the entire year
# We then sample randomly a number of background points 5 * the number of their soaring locations
# And we associate these points to the days in which these individuals were flying and the median flight hour

soarPoints <- readRDS("Data/gps_allsoaringpoints_areavol_only2023.rds")
all(soarPoints$segmToKeep==T)

# Extract median flight hour across all individuals, days and soaring segments
summary(hour(soarPoints$timestamp))
medHour <- median(hour(soarPoints$timestamp))
(medHeight_asl <- round(median(soarPoints$height_asl, na.rm=T)))
(medHeight_agl <- round(median(soarPoints$height_agl, na.rm=T)))
(medHeight_ell <- round(median(soarPoints$height_above_ellipsoid, na.rm=T)))

# On average how many segments per day
segmPerDay <- summarise(group_by(soarPoints, individual_local_identifier, date(timestamp)), nSegm=length(unique(unique_segmID)))
summary(segmPerDay$nSegm) # median 33, mean 38. 38 * 5 = 190, so we take on average 200 random locations per MCP

# turn into spatial object
class(soarPoints) <- setdiff(class(soarPoints), c("tbl_df","tbl")) #vect and tbl are not compatible
soarPoints <- terra::vect(soarPoints, geom=c("location_long", "location_lat"), crs="EPSG:4326", keepgeom=T)

# split by individual
ind_ls <- split(soarPoints, as.character(soarPoints$individual_local_identifier))
rm(soarPoints);gc() # free up space

# For each individual, calculate and store minimum convex polygons of all their soaring locations across the year
polyList <- list() # prepare to empty lists to fill
polyInfos <- as.data.frame(matrix(nrow=length(ind_ls), ncol=7, 
                                  dimnames = list(NULL, c("individual_local_identifier","Ndays","minDate","maxDate","medHour","MCParea_km2","NsoarSegm"))))
backPoints_ls <- list()

for(i in 1:length(ind_ls)){
  print(paste0("Individual ",i," of ",length(ind_ls)))
  ind <- ind_ls[[i]]
  # calculate a minimum convex polygon including all soaring locations of that day
  ind_mcp <- terra::convHull(ind) 
  # calculate area, better if in latlong geog coords
  mcpArea <- terra::expanse(ind_mcp)
  # sample 200 random points within the polygon
  random_points <- spatSample(ind_mcp, size = 200, method = "random")
  #plot(ind_mcp); plot(random_points, add=T)
  # store polygons into a list
  polyList[[i]] <- ind_mcp
  # and information about the polygons together with info about soaring time and height for that individual
  polyInfos[i,] <- c(names(ind_ls)[i], length(unique(date(ind$timestamp))),
                     as.character(min(date(ind$timestamp))),as.character(max(date(ind$timestamp))),
                     median(hour(ind$timestamp)),
                     round(mcpArea / 1e6, 2),
                     length(unique(ind$unique_segmID)))
  # associate random points with id and area
  random_points <- as.data.frame(geom(random_points)[,c("x","y")])
  random_points$individual_local_identifier <- names(ind_ls)[i]
  random_points$MCParea_km2 <- round(mcpArea / 1e6, 2)
  # associate random points with all possible timestamps and store in list
  timestamps <- paste0(unique(date(ind$timestamp))," ",medHour,":00:00")
  backPoints_ls[[i]] <- tidyr::expand_grid(random_points, timestamp=timestamps)
}

names(polyList) <- names(ind_ls)
allBackPoints <- rbindlist(backPoints_ls)
polyInfos[,c(2,5:7)] <- lapply(polyInfos[,c(2,5:7)], as.numeric)

# Add median flight height for annotation
allBackPoints$medHeight_asl <- medHeight_asl
allBackPoints$medHeight_agl <- medHeight_agl
allBackPoints$medHeight_ell <- medHeight_ell

save(polyList,polyInfos, file = "individualMCPs.rdata")
save(allBackPoints, file = "backgroundPoints_AtmoAvailability_perInd.rdata")
write.csv(as.data.frame(allBackPoints), "backgroundPoints_AtmoAvailability_perInd.csv", row.names = F)

# Visualise in mapview
class(allBackPoints) <- setdiff(class(allBackPoints), c("tbl_df","tbl")) #vect and tbl are not compatible
backPoints_mapview <- terra::vect(unique(allBackPoints[,c("x","y")]), geom=c("x", "y"), crs="EPSG:4326", keepgeom=T)

myCols <- rainbow(length(polyList))
layerNames <- names(polyList)
maps <- mapply(
  function(poly, col, lname) mapview(poly, col.regions = col, layer.name = lname),
  polyList, myCols, layerNames, SIMPLIFY = FALSE
)
(final_map <- Reduce(`+`, maps))
(final_map <- final_map + mapView(backPoints_mapview, layer.name="Background points", cex=0.8))

# Save to html
htmlwidgets::saveWidget(final_map@map, file = "mapview_MCPs_perIndividual.html")

