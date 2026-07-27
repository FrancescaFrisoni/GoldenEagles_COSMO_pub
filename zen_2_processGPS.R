# This code finds the high resolution GPS bursts
# Extracts the uplift segments and calculates GPS-based movement metrics
# Martina Scacco & Francesca Frisoni - May 16, 2024. Konstanz

setwd("/home/francesca/ownCloud/TesiFrancesca_UpliftClassification/Data")

library(move2)
library(dplyr)
library(data.table)
library(sf)
library(terra)
library(units)
library(ggplot2)
library(raster)


#### 1. FIRST GPS PROCESSING AND FILTERING ####

#_______________________________________________
# Calculate height above ground and GPS metrics 

# Import rasters
geoid <- rast("directory_env_data/EGM96_us_nga_egm96_15.tif") # wgs 84 in metres
dem <- rast("directory_env_data/DEM_Europe_eumap_epsg3035_v03/dtm_elev.lowestmode_gedi.eml_mf_30m_0..0cm_2000..2018_eumap_epsg3035_v0.3.tif") # laea in decimetres

# Crop DEM to the alpine region and reproject it to the data
bbox <- st_bbox(c(xmin = 5.0, ymin = 42.0, xmax = 16.0, ymax = 49.0), crs = st_crs(4326))
bbox_laea <- st_transform(st_as_sfc(bbox), crs = st_crs(3035))
geoid_crop <- crop(geoid, bbox)
dem_crop <- crop(dem, bbox_laea)
writeRaster(geoid_crop, "geoid_EGM96_us_nga_egm96_15_cropAlps_wgs84.tif")
writeRaster(dem_crop, "dtm_elev_eumap_epsg3035_v0.3_cropAlps_laea.tif")
# geoid and dem, cropped in Alps, are saved in the folder Data

# List files
fls <- grep("df", list.files("GPS_allIMUeagles_20240516", full.names = T), invert=T, value=T)

# Calculate vars on each individual file
lapply(fls, function(f){
  gps <- readRDS(f)
  # remove empty locations and sort by time
  gps <- filter(gps, !sf::st_is_empty(gps))
  gps <- gps %>% arrange(mt_track_id(), mt_time())
  # there are duplicated timestamps+locations
  # table(duplicated(gps[,c("timestamp","geometry")]))
  # We create a vector of number of NAs per row to filter out the duplicates with less NAs
  numNA <- rowSums(is.na(gps)) 
  gps <- arrange(gps, numNA) # and we arrange by that as the functions duplicated always takes the second repetition as TRUE duplicate
  gps <- gps[!duplicated(round(as.numeric(mt_time(gps)), 0)),] # remove duplicates in rounded time and position
  gps <- gps %>% arrange(mt_track_id(), mt_time()) # re-arrange by individual id and time
  # min(mt_time_lags(gps, "secs"), na.rm=T) # still there are many values with timelag < 1 s
  # length(which(mt_time_lags(gps, "secs") < set_units(1,"sec")))
  
  # Calculate basic metrics
  gps$timelag <- mt_time_lags(gps, units = "secs")
  gps$turnAngle <- mt_turnangle(gps)
  gps$stepLength <- mt_distance(gps)
  gps$azimuth <- mt_azimuth(gps)
  # Check why we still have NAs in azimuth after removing duplicates
  # these NAs derive from having points in subsequent segments at the exact same location. For now we keep them and we will remove them later on.
  # table(is.na(gps$azimuth))
  # NAtimes <- gps$roundTimestamp_s[is.na(gps$azimuth)]
  # rowIdx <- which(gps$roundTimestamp_s == NAtimes[1])
  # gps[(rowIdx-1):(rowIdx+1),] %>% as.data.frame()
  # reproject data and extract geoid and dem info
  gps$geoid_egm96 <- extract(geoid_crop, st_coordinates(gps), method="bilinear", cells=F)[,1]
  # test <- extract(geoid_crop, st_coordinates(segmdf_areavol$geometry), method="bilinear", cells=F)[,1]
  gps_laea <- st_transform(gps, crs = st_crs(dem_crop)$epsg)
  gps$dem_30m <- extract(dem_crop, st_coordinates(gps_laea), method="bilinear", cells=F)[,1]
  gps$dem_30m <- gps$dem_30m/10 # from dm to m
  # calculate height above sea level and above ground
  gps$height_asl <- gps$height_above_ellipsoid - set_units(gps$geoid_egm96, "m")
  gps$height_agl <- gps$height_asl - set_units(gps$dem_30m, "m")
  # check for NAs in dem and height agl
  # table(is.na(gps$dem_30m))
  # table(is.na(gps$height_agl))
  
  # overwrite existing file with added info
  saveRDS(gps, f)
})

# Some quick checks
fls <- grep("df", list.files("GPS_allIMUeagles_20240516", full.names = T), invert=T, value=T)
lapply(fls, function(f){
  gps <- readRDS(f)
  table(is.na(gps$turnAngle))
  min(gps$timelag,na.rm=T)
})
# Now the timelag look fine, but there are many gps positions that are overlapping causing the azimuth and turning angle to be NA

#________________________
# Subset to HR 1s bursts

# Bind all data and calculate burst ID of 1 sec data based on timelag
fls <- grep("df", list.files("GPS_allIMUeagles_20240516", full.names = T), invert=T, value=T)

lapply(fls, function(f){
  gps <- readRDS(f)
  gps$burstID <- c(0, cumsum(gps$timelag[-nrow(gps)] > set_units(2, "sec")))
  # keep bursts of at least 1 min duration (median 5 min, max 15 min)
  burstsToKeep <- names(table(gps$burstID)[table(gps$burstID) > 60])
  HRgps <- filter(gps, burstID %in% burstsToKeep)
  # recalculate some of the variables per burst
  HRgps <- HRgps %>% group_by(burstID) %>% arrange(timestamp, .by_group = T) %>% 
    mutate(timelag = c(NA, diff(timestamp)), #diff(time) and mt_time_lag(time) do the same
           grSpeed = stepLength/timelag,
           vertDist = c(NA, diff(height_above_ellipsoid)),
           vertSpeed = vertDist/timelag) %>% ungroup()
  # Transform move2 into data.frame
  HRgps <- mt_as_event_attribute(HRgps, c("individual_local_identifier","deployment_id","tag_id","individual_id","study_id","taxon_canonical_name","tag_local_identifier"))
  HRgps$location_long <- st_coordinates(HRgps)[,1]
  HRgps$location_lat <- st_coordinates(HRgps)[,2]
  HRgps <- st_drop_geometry(HRgps)
  HRgps <- data.frame(HRgps)
  # save per individual
  saveRDS(HRgps, sub("/GPS_","/df_1sGPS_", f))
})


#__________________________________________________________
# Keep only GPS that were collected simultaneously to IMU 

imu <- readRDS("IMU_allEagles_20240516.rds") # all imu data
fls <- list.files("GPS_allIMUeagles_20240516", pattern="df_1sGPS_", full.names = T) # gps files per ind

# Keep only GPS that were collected simultaneously to IMU
timeTolerance <- 5 # sec

imu_timestamps_ls <- split(imu, imu$individual_local_identifier)

GPSsub_allInds <- rbindlist(lapply(fls, function(f){
  HRgps <- readRDS(f)
  ind <- as.character(unique(HRgps$individual_local_identifier))
  GPS_burstTimes <- HRgps %>% arrange(burstID, timestamp) %>% group_by(burstID) %>% summarise(min(timestamp))
  IMUtime <- imu_timestamps_ls[[ind]]$timestamp
  #returns logical vector of which GPS bursts (burst starting time) have IMU within time tolerance
  has_imu <- sapply(GPS_burstTimes$`min(timestamp)`, function(g){ 
    #any(IMUtime > g-timeTolerance & IMUtime < g+timeTolerance)
    any(between(IMUtime, g-timeTolerance, g+timeTolerance)) # ~2x faster
  })
  # keep bursts whose starting time is within 5 sec from IMU data
  burstsToKeep <- filter(GPS_burstTimes, has_imu==T)$burstID
  HRgps_sub <- filter(HRgps, burstID %in% burstsToKeep)
  # Make burstID unique across individuals before pulling all individuals together
  HRgps_sub$burstID <- paste0(HRgps_sub$deployment_id,"_",HRgps_sub$burstID)
  return(HRgps_sub)
}))

# save all GPS subsets together
saveRDS(GPSsub_allInds, "GPS_1secBursts_matchingIMUtimes_20240521.rds")

# ____________________________________________________
# Make a subset including only data for 2020 and 2023 
# (these are the years for which we have COSMO downloaded)
table(year(GPSsub_allInds$timestamp))
GPS_2020_2023 <- filter(GPSsub_allInds, year(timestamp) %in% c("2020","2023"))
table(year(GPS_2020_2023$timestamp))

# These is the dataset that will be sent to ETH for the uplift labelling
saveRDS(GPS_2020_2023, "GPS_1secBursts_matchingIMUtimes_2020&2023_4eth.rds")


#### 2. DEFINE CLIMBING SEGMENTS ####

#_______________________________________
# Smoothing window for vertical velocity

df <- readRDS("GPS_1secBursts_matchingIMUtimes_2020&2023_4eth.rds")
# 3219207 obs of 24 variables

# applying the code to a subset of individuals for checking: only 1
# df <- df[df$individual_local_identifier=="Zebru20 (eobs 7554)",] # 3281 observations

# transformed timestamp in the correct date and time format
df$date <- as.Date(df$timestamp) 
df$time <- format(as.POSIXct(df$timestamp), format = "%H:%M:%OS", tz="UTC")

# changing smoothing window of vertical velocity to 21 seconds to be consistent with COSMO labelling
# classification of soaring segments in each burst: grouped by burst for simplicity, not useful later

sub_ls <- split(df, df$burstID) # create a list grouped via burst

df_clf <- as.data.frame(rbindlist(
  lapply(sub_ls, function(b){
    #b <- sub_ls[[1]]
    
    swV = 10 # takes 21 points/observations: a central point + 10 before and after
    b$vertSpeed_smooth_new <- NA
    b$height_asl_smooth <- NA # for visualization purposes
    b$location_long_smooth <- NA # for visualization purposes
    b$location_lat_smooth <- NA # for visualization purposes
    
    for(i in (swV+1):(nrow(b)-swV)){
      b$vertSpeed_smooth_new[i] <- mean(b$vertSpeed[(i-swV):(i+swV)], na.rm=T)
      b$height_asl_smooth[i] <- mean(b$height_asl[(i-swV):(i+swV)], na.rm=T)
      b$location_long_smooth[i] <- mean(b$location_long[(i-swV):(i+swV)], na.rm=T)
      b$location_lat_smooth[i] <- mean(b$location_lat[(i-swV):(i+swV)], na.rm=T)
    }
    # remember that first and last 10 rows are NA
    
    b$new_class <- "glide"
    b$new_class[b$vertSpeed_smooth_new>0] <- "soar"
    
    b$segmID <- c(0, cumsum(b$new_class[-1] != b$new_class[-nrow(b)]))
    b$segmToKeep <- F # to exclude all the gliding segments
    b$segmToKeep[b$new_class == "soar"] <- T
    
    # removed first and last soaring segments to be sure of having the complete sequence: 1st and last segm could be incomplete due to the tag fix
    if(b$new_class[1] == "soar"){b$segmToKeep[b$segmID == b$segmID[1]] <- F}
    if(b$new_class[nrow(b)] == "soar"){b$segmToKeep[b$segmID == b$segmID[nrow(b)]] <- F}
    
    # assign a new column to uniquely identify the GPS point
    b$unique_pointID <- paste0(b$individual_local_identifier,"_", b$event_id)
    # assign a new column to uniquely identify the climbing segment
    b$unique_segmID <- paste0(b$individual_local_identifier,"_", b$burstID,"_",b$new_class,"_", b$segmID)
    
    return(b) # return df with burst newly classified
  })))

saveRDS(df_clf, "GPS_1secBursts_allSegmentsSoarGlide_2020&2023.rds")

# subset with only soaring segments
df_clf <- df_clf[df_clf$segmToKeep==TRUE,] 
saveRDS(df_clf, "GPS_1secBursts_onlySoaringObservations_2020&2023_segmID.rds")

# Sent this to the ETH for the labelling
# Important to have an ID per segment or per GPS point (e.g. the event_id) that we can use as merging column later


#### 3. GPS METRICS ASSOCIATED TO SEGMENTS ####
# The output is still a dataset with one observation per gps point (but only soaring)

# if not yet loaded:
# df_clf <- readRDS("GPS_1secBursts_onlySoaringObservations_2020&2023_segmID.rds")

# Split for each climbing segment
segm_ls <-split(df_clf, df_clf$unique_segmID) 

# First lapply to calculate dist, duration of each segment
# area and volume of the segments was in the end not calculated
segmdf_areavol <- as.data.frame(rbindlist(
  lapply(segm_ls, function(b){
    #b <- segm_ls[[2]]
    
    # add geometry column to calculate space covered (sum of steplength is equivalent only on a linear path)
    b <- st_as_sf(b, coords = c("location_long", "location_lat"), crs = 4326)
    
    # "horizontal displacement", kept the max bc st_distance returns a matrix of all the distances between points
    b$dist <- max(st_distance(b$geometry)) 
    
    # to keep long and lat
    coords <- st_coordinates(b$geometry)
    b$location_long <- coords[,1]
    b$location_lat <- coords[,2]
    
    # calculate segm duration (for subsequent filtering)
    b$duration <- sum(b$timelag, na.rm=T)
    
    return(b) 
  }))
)
# remove segments with less than 3 consecutive observations (bc not able to calculate the turn angle in a meaningful way)
table(segmdf_areavol$duration >= 3)
segmdf_areavol <- segmdf_areavol[segmdf_areavol$duration >= 3,]

# # Therefore we can safely remove rows with empty turning angle and empty height_above_ellipsoid
# table(is.na(segmdf_areavol$turnAngle))
# table(is.na(segmdf_areavol$height_above_ellipsoid))
# segmdf_areavol <- segmdf_areavol[!is.na(segmdf_areavol$turnAngle),]
# segmdf_areavol <- segmdf_areavol[!is.na(segmdf_areavol$height_above_ellipsoid),]
# Add variable with change in turning during the circle, by segment
segmdf_areavol <- segmdf_areavol %>% group_by(unique_segmID) %>% 
  mutate(turnChange = c(0,diff(sign(turnAngle)))) %>% ungroup()
table(is.na(segmdf_areavol$turnChange))
# Checking for NAs in height columns
table(is.na(segmdf_areavol$height_above_ellipsoid))
table(is.na(segmdf_areavol$height_agl)) #why are there NAs in the height agl and asl?
table(is.na(segmdf_areavol$height_asl))
table(is.na(segmdf_areavol$dem_30m))
table(is.na(segmdf_areavol$geoid_egm96)) # because of the geoid
segmsToCheck <- segmdf_areavol$unique_segmID[is.na(segmdf_areavol$height_agl)]
test <- segmdf_areavol[segmdf_areavol$unique_segmID%in%segmsToCheck,]
print(test[,c("unique_segmID","dem_30m","height_above_ellipsoid","height_agl","timelag")], n=50)
# We can easily recalculate the missing geoid and height information
geoid_crop <- rast("geoid_EGM96_us_nga_egm96_15_cropAlps_wgs84.tif")
segmdf_areavol$geoid_egm96 <- extract(geoid_crop, st_coordinates(segmdf_areavol$geometry), method="bilinear", cells=F)[,1]
segmdf_areavol$geoid_egm96 <- set_units(segmdf_areavol$geoid_egm96, "m")
segmdf_areavol$dem_30m <- set_units(segmdf_areavol$dem_30m, "m")
segmdf_areavol$height_asl <- segmdf_areavol$height_above_ellipsoid - segmdf_areavol$geoid_egm96
segmdf_areavol$height_agl <- segmdf_areavol$height_asl - segmdf_areavol$dem_30m
# Checking for NAs values in turning angle, when do they occur?
table(is.na(segmdf_areavol$azimuth))
table(is.na(segmdf_areavol$turnAngle))
segmsToCheck <- segmdf_areavol$unique_segmID[is.na(segmdf_areavol$azimuth)]
test <- segmdf_areavol[segmdf_areavol$unique_segmID%in%segmsToCheck,]
print(test[,c("unique_segmID","turnAngle","azimuth","geometry","turnChange")], n=50)
# Azimuth is NA at the end of each track and for consecutive points (1 sec apart) that have exact same geometry
# Turn change is also NA when turning angle is NA, but since the turnAngle=NA is due to the geometry being the same we can assume no change in turning
# Therefore we assign turnChange=NA to 0
segmdf_areavol$turnChange[is.na(segmdf_areavol$turnChange)] <- 0
# Re-save these dataset that contains all observations of soaring segments with filtered NAs, filteres short segments and additional var calculated
# saveRDS(segmdf_areavol, "GPS_1secBursts_onlySoaringObservations_2020&2023_filterNA_segmID&otherVars.rds") # by Martina
saveRDS(segmdf_areavol, "GPS_allsoaringpoints_areavol.rds") # by Francesca on 1st July 2024 with small changes

