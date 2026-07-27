# This code process the data provided by Tom Carrard ETH with manually labelled uplifts with COSMO
# added topography characterization and modelled DFA classification algorithm for complete dataset
# Francesca Frisoni - June 20, 2024. Konstanz

# bind old labelled dataset (case_studies) with the new one (labelled unique_segmID from 2023)
# change of smoothing window: both 21 sec
# bind to each GPS point the topography parameters
# bind the labelled events to the complete trajectories dataset by event.id/event_id
# summarize useful variables with unique_segmID
# run the DFA on the COSMO and topographical variables: validation (10 runs) and on my complete dataset (1 run)
# from predictions on complete dataset, boxplots to compare parameters distribution
# assumption of correct classification, run random forest in next code
# trajectories for visualization

#  !! topographical variables (coords) are overwritten to keep RAM free enough to work without crashing, be careful if need to rerun some parts
# if code ran on a better computer, probably better to make unique variable names

setwd("/home/francesca/ownCloud/TesiFrancesca_UpliftClassification/Data")

directory <- "/home/francesca/ownCloud/TesiFrancesca_UpliftClassification"

library(data.table)
library(dplyr)
library(units)
library(sf)
library(lubridate)
library(MASS)
library(caret)
library(plotly)
library(terra)
library(sp)
library(circular)
library(ggtern)
library(ggplot2)
library(gridExtra)

# function to calculate wind speed, used below
wind.speed <- function(u, v) {
  return(sqrt(u^2 + v^2))
}

##### 1. MANUALLY LABELLED DATASET: COMBINED WEATHER AND TOPOGRAPHICAL VARIABLES #####
# this chunk of code produces the manually labelled dataset of 209 segments that will be used for the DFA models

# ___________________________________________
#### Tom's first uplift classification 
# 26 case studies, with time start and end, needed to be combined to a unique_segmID
# segmentation was different in the first place, now binded each event.id and put a smoothing window of 21 sec

tomcomp <- read.csv("/home/francesca/ownCloud/ETH_sharedFolder/LabelledEvents_upliftSources/trajectories_case_studies.csv")
tomlab <- read.csv("/home/francesca/ownCloud/ETH_sharedFolder/LabelledEvents_upliftSources/uplifts_case_studies_labelled_new_corrected_annotated_W.csv")

tomlab$date <- as.Date(tomlab$timestamp)
# tomlab$Time <- format(as.POSIXct(tomlab$timestamp), format = "%H:%M:%S")
tomlab$windspeed <- wind.speed(tomlab$U, tomlab$V)

# classification of soaring segments in each burst
complete_ls <- split(tomcomp, tomcomp$burstID) # split each burst: 82 bursts

complete_df <- as.data.frame(rbindlist(
  lapply(complete_ls, function(b) {
    # b <- complete_ls[[3]]
    
    # check if the number of rows is less than 21
    if (nrow(b) < 21) {
      return(NULL) # Return NULL to skip this data frame entirely bc can't calculate the smoothing
    }
    
    swV <- 10 # takes 21 points/observations: 10 before and after the observation point
    b$vertSpeed_smooth_new <- NA
    
    for (i in (swV + 1):(nrow(b) - swV)) {
      b$vertSpeed_smooth_new[i] <- mean(b$vert.speed[(i - swV):(i + swV)], na.rm = T)
    }
    # remember that first 10 and last 10 rows are NA
    
    b$new_class <- "glide"
    # con for loop: if(b$vertSpeed_smooth[i]>0) {b$new_class[i] <- "soar"}
    b$new_class[b$vertSpeed_smooth_new > 0] <- "soar"
    
    b$segmID <- c(0, cumsum(b$new_class[-1] != b$new_class[-nrow(b)]))
    b$segmToKeep <- F
    b$segmToKeep[b$new_class == "soar"] <- T
    
    # dont' keep the soaring segments at the beginning and at the end of the burst bc not sure if it's fully complete
    if (b$new_class[1] == "soar") {
      b$segmToKeep[b$segmID == b$segmID[1]] <- F
    }
    if (b$new_class[nrow(b)] == "soar") {
      b$segmToKeep[b$segmID == b$segmID[nrow(b)]] <- F
    }
    
    # assign a new column to uniquely identify the climbing segment
    b$unique_segmID <- paste0(b$individual.local.identifier, "_", b$burstID, "_", b$new_class, "_", b$segmID)
    
    
    return(b) # return df with burst newly classified
  })
))

df_sub <- complete_df[complete_df$segmToKeep == TRUE, ] # 30937 observations: removed the glide segments and the excluded soaring segments
length(unique(df_sub$unique_segmID)) # 690 unique segments
df_sub$Date <- as.Date(df_sub$timestamp)
df_sub$Time <- format(as.POSIXct(df_sub$timestamp), format = "%H:%M:%S")

# _______________________________
#### Topography parameters 

# Load DEM
dem <- terra::rast("./dtm_elev_eumap_epsg3035_v0.3_cropAlps_laea.tif")
# this dem is in decimeters (dm), not yet in meters
dem <- dem / 10 # convert to m (check summary of the values to make sure it worked)
dem[dem < -1000] <- NA # put values below -1000 as NA to avoid problems with roughness<0
# Calculate, from dem in m, slope, roughness and aspect in degrees
slope <- terra::terrain(dem, v = "slope", neighbors = 8, unit = "degrees", filename = "")
aspect <- terra::terrain(dem, v = "aspect", neighbors = 8, unit = "degrees", filename = "")
roughness <- terra::terrain(dem, v = "roughness", neighbors = 8, unit = "degrees", filename = "")

# topography on points before matching labelled segments
coords <- df_sub[, c("location.long", "location.lat")]
coords <- st_as_sf(coords, coords = c("location.long", "location.lat"), crs = "+proj=longlat +ellps=WGS84")
coords_proj <- st_transform(coords, crs(dem))

df_sub$slope <- terra::extract(slope, coords_proj, method = "bilinear", cells = F, ID = F)[, 1]
df_sub$aspect <- terra::extract(aspect, coords_proj, method = "bilinear", cells = F, ID = F)[, 1] # to have aspect between 0 and 360 needed method "simple"
df_sub$roughness <- terra::extract(roughness, coords_proj, method = "bilinear", cells = F, ID = F)[, 1]
terrainElevation <- terra::extract(dem, coords_proj, method = "bilinear", cells = F, ID = F)[, 1]
df_sub$height_above_ground <- df_sub$height.above.sea.level - terrainElevation

# bind the unique_segmID to tomlab: tomlab is simply a subset of tomcomp where labelled uplift

lab_segmid <- merge(df_sub[, c(2, 50, 53:56)], tomlab, by = "event.id") # columns from df_sub: event.id, unique_segmID, topography variables

lab_segmid$aspect_rad <- circular(lab_segmid$aspect * pi / 180, units = "radians") # convert aspect to radians: circular variable!

# Summarize variables into segments
lab_ls <- split(lab_segmid, lab_segmid$unique_segmID) # 202 segments 

lab_df <- as.data.frame(rbindlist(
  lapply(lab_ls, function(b) {
    # b <- lab_ls[[5]]
    
    summary_caselab <- b %>%
      group_by(date, individual.local.identifier, unique_segmID, event_type, case_studyID) %>%
      summarise(
        time_start = min(timestamp, na.rm = TRUE),
        time_end = max(timestamp, na.rm = TRUE),
        duration = sum(timelag.sec, na.rm = TRUE),
        w_oro_mean = mean(w_oro, na.rm = TRUE),
        N2_max_h_mean = mean(N2_max_h, na.rm = TRUE),
        N2_max_value_mean = mean(N2_max_value, na.rm = TRUE),
        ASHFL_S_mean = mean(ASHFL_S, na.rm = TRUE),
        U_mean = mean(U, na.rm = TRUE),
        V_mean = mean(V, na.rm = TRUE),
        W_mean = mean(abs(W), na.rm = TRUE), # taken as abs value bc in gravity waves andament otherwise would not be detected
        windspeed_mean = mean(windspeed, na.rm = TRUE),
        max_height.ab.gr = max(height_above_ground, na.rm = TRUE),
        mean_slope = mean(slope, na.rm = TRUE), # slope is not a circular variable even if degree
        mean_roughness = mean(roughness, na.rm = TRUE),
        mean_aspect = as.numeric(mean.circular(aspect_rad), na.rm = TRUE) * 180 / pi, # circular variable, calculated in rad and back to degree
        
        .groups = "drop"
      )
    
    return(summary_caselab)
  })
))

# here not applied filter on duration yet, done after: kept lab_df with 202 segments
names(lab_df)[names(lab_df) == "event_type"] <- "uplift_type"
names(lab_df)[names(lab_df) == "individual.local.identifier"] <- "individual_local_identifier"

table(lab_df$uplift_type)
# orog orog/therm  orog/wave      therm       wave 
# 57         19         15         58         53

# _______________________________
#### Tom's second labelled dataset 

# Load and bind lab_ascents with the last 23 segments especially selected from winter 2023
lab_ascents <- read.csv("/home/francesca/ownCloud/ETH_sharedFolder/LabelledEvents_upliftSources/labelled_150ascents_cosmo.csv")
last <- read.csv("/home/francesca/ownCloud/ETH_sharedFolder/LabelledEvents_upliftSources/labelled_23ascents_winter2023_cosmo.csv")
tot_lab <- rbind(lab_ascents, last) # 48917 obs in 41 vars

# quick summary of the 150 ascents
ascent_summary <- lab_ascents %>%
  distinct(ascent_ID, uplift_type)

table(ascent_summary$uplift_type)
# dynamic     orog      thermal thermal/orog thermal/wave      unknown         wave 
# 1            7           57           30           15           32            8 

# quick summary of the 23 ascents
ascent_summary_last <- last %>%
  distinct(ascent_ID, uplift_type)

table(ascent_summary_last$uplift_type)
# orog thermal unknown    wave 
# 6       1       9       7

# file with U, V and W (corrected on 15th july 26)
env_cosmo <- read.csv("/home/francesca/ownCloud/ETH_sharedFolder/GPS_allsoaringpoints_areavol_subset_annotated_fullwind.csv") 

# incomplete trajectories with all COSMO variables, U and V need to be removed before merging bc the correct ones are in env_cosmo
env_cosmo2 <- read.csv("/home/francesca/ownCloud/ETH_sharedFolder/GPS_allsoaringpoints_areavol_subset_annotated_corrected.csv") 
env_cosmo2 <- env_cosmo2[, !names(env_cosmo2) %in% c("U", "V")]

# Merge the wind columns to env_cosmo2 by event_id
env_cosmo <- merge(env_cosmo2, env_cosmo[, c("event_id", "V","U","W")], by = "event_id", all.x = TRUE) 

names(env_cosmo)
# [1] "event_id"                    "X"                           "Unnamed..0.1"                "index"                       "Unnamed..0"                 
# [6] "height_above_ellipsoid"      "deployment_id"               "timestamp"                   "timelag"                     "turnAngle"                  
# [11] "stepLength"                  "azimuth"                     "geoid_egm96"                 "dem_30m"                     "height.above.sea.level"     
# [16] "height_agl"                  "burstID"                     "grSpeed"                     "vertDist"                    "vertSpeed"                  
# [21] "individual_local_identifier" "tag_id"                      "individual_id"               "study_id"                    "taxon_canonical_name"       
# [26] "tag_local_identifier"        "location.long"               "location.lat"                "date"                        "time"                       
# [31] "vertSpeed_smooth"            "location_long_smooth"        "location_lat_smooth_new"     "height_asl_smooth"           "new_class"                  
# [36] "segmID"                      "segmToKeep"                  "unique_pointID"              "unique_segmID"               "area"                       
# [41] "dist"                        "h_iniz"                      "h_fin"                       "deltah"                      "duration"                   
# [46] "turnChange"                  "vol"                         "rlon"                        "rlat"                        "height.above.surf"          
# [51] "N2_max_h"                    "N2_max_value"                "w_oro"                       "ASHFL_S"                     "V"                          
# [56] "U"  

# env_cosmo is the complete dataset of all points already assigned to unique_segmID, 578335 obs in 57 vars
length(unique(env_cosmo$unique_segmID)) # 34059 unique segments, contains both the labelled manually and the behav ones that I will classify later

# extract the segments with a manually labelled uplift
labelled <- env_cosmo[env_cosmo$event_id %in% tot_lab$event.id, ] 
length(unique(labelled$unique_pointID)) # 12884: 11564 from lab_ascents + 1320 from last events-points are contained into the complete traj dataset
length(unique(labelled$unique_segmID)) # 173 segments

# Topography on all points: attention on overwritten, same name of variables. 
coords <- labelled[, c("location.long", "location.lat")]
coords <- st_as_sf(coords, coords = c("location.long", "location.lat"), crs = "+proj=longlat +ellps=WGS84")
coords_proj <- st_transform(coords, crs(dem))

labelled$slope <- terra::extract(slope, coords_proj, method = "bilinear", cells = F, ID = F)[, 1]
labelled$aspect <- terra::extract(aspect, coords_proj, method = "bilinear", cells = F, ID = F)[, 1]
labelled$roughness <- terra::extract(roughness, coords_proj, method = "bilinear", cells = F, ID = F)[, 1]
terrainElevation <- terra::extract(dem, coords_proj, method = "bilinear", cells = F, ID = F)[, 1]
labelled$height_above_ground <- labelled$height.above.sea.level - terrainElevation

labelled$aspect_rad <- circular(labelled$aspect * pi / 180, units = "radians")

tot_lab$windspeed <- wind.speed(tot_lab$U, tot_lab$V)
tot_lab$timestamp <- as.POSIXct(tot_lab$timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
tot_lab$date <- as.Date(tot_lab$timestamp)
# tot_lab$time <- format(as.POSIXct(lab_ascents$timestamp), format = "%H:%M:%S", tz="UTC")

names(tot_lab)
# [1] "X"                             "event.id"                      "height.above.ellipsoid"        "deployment.id"                
# [5] "timestamp"                     "timelag"                       "turnAngle"                     "stepLength"                   
# [9] "azimuth"                       "geoid.egm96"                   "dem.30m"                       "height.above.sea.level"       
# [13] "height.agl"                    "burstID"                       "grSpeed"                       "vertDist"                     
# [17] "vertSpeed"                     "individual.local.identifier"   "tag.id"                        "individual.id"                
# [21] "study.id"                      "taxon.canonical.name"          "tag.local.identifier"          "location.long"                
# [25] "location.lat"                  "location.long_smooth"          "location.lat_smooth"           "height.above.sea.level_smooth"
# [29] "vertSpeed_smooth"              "distance"                      "behavior"                      "ascent_ID"                    
# [33] "iscase"                        "uplift_type"                   "height.above.surf"             "N2_max_h"                     
# [37] "N2_max_value"                  "w_oro"                         "ASHFL_S"                       "U"                            
# [41] "V"                             "windspeed"                     "date"   

# bind tot_lab(all gps points) to get unique_segmID from env_cosmo
names(labelled)[names(labelled) == "event_id"] <- "event.id"
lab_new_segmid <- merge(labelled[, c(1, 39, 57:62)], tot_lab, by = "event.id") # columns from labelled: event.id, unique_segmID, W, topography variables

names(lab_new_segmid)
# [1] "event.id"                      "unique_segmID"                 "W"                             "slope"                        
# [5] "aspect"                        "roughness"                     "height_above_ground"           "aspect_rad"                   
# [9] "X"                             "height.above.ellipsoid"        "deployment.id"                 "timestamp"                    
# [13] "timelag"                       "turnAngle"                     "stepLength"                    "azimuth"                      
# [17] "geoid.egm96"                   "dem.30m"                       "height.above.sea.level"        "height.agl"                   
# [21] "burstID"                       "grSpeed"                       "vertDist"                      "vertSpeed"                    
# [25] "individual.local.identifier"   "tag.id"                        "individual.id"                 "study.id"                     
# [29] "taxon.canonical.name"          "tag.local.identifier"          "location.long"                 "location.lat"                 
# [33] "location.long_smooth"          "location.lat_smooth"           "height.above.sea.level_smooth" "vertSpeed_smooth"             
# [37] "distance"                      "behavior"                      "ascent_ID"                     "iscase"                       
# [41] "uplift_type"                   "height.above.surf"             "N2_max_h"                      "N2_max_value"                 
# [45] "w_oro"                         "ASHFL_S"                       "U"                             "V"                            
# [49] "windspeed"                     "date"  

# summary of each segment: 173 (150 from 2023 and 23 from 2020)
tom_segm_ls <- split(lab_new_segmid, lab_new_segmid$unique_segmID)

sub_env_df <- as.data.frame(rbindlist(
  lapply(tom_segm_ls, function(b) {
    # b <- tom_segm_ls[[5]]
    
    summary_caselab <- b %>%
      group_by(
        date, individual.local.identifier, unique_segmID,
        uplift_type, burstID, ascent_ID
      ) %>%
      summarise(
        time_start = min(timestamp, na.rm = TRUE),
        time_end = max(timestamp, na.rm = TRUE),
        duration = sum(timelag, na.rm = TRUE),
        w_oro_mean = mean(w_oro, na.rm = TRUE),
        N2_max_h_mean = mean(N2_max_h, na.rm = TRUE),
        N2_max_value_mean = mean(N2_max_value, na.rm = TRUE),
        ASHFL_S_mean = mean(ASHFL_S, na.rm = TRUE),
        U_mean = mean(U, na.rm = TRUE),
        V_mean = mean(V, na.rm = TRUE),
        W_mean = mean(abs(W), na.rm = TRUE), # taken as abs value bc in graves wavy andament otherwise would not be detected
        windspeed_mean = mean(windspeed, na.rm = TRUE),
        max_height.ab.gr = max(height_above_ground, na.rm = TRUE),
        mean_slope = mean(slope, na.rm = TRUE), # slope is not a circular variable
        mean_roughness = mean(roughness, na.rm = TRUE),
        mean_aspect = as.numeric(mean.circular(aspect_rad), na.rm = TRUE) * 180 / pi, # circular variable, calculated in rad and back to degree
        
        .groups = "drop"
      )
    
    return(summary_caselab)
  })
))

# change column names to be consistent between the 2 datasets
names(sub_env_df)[names(sub_env_df) == "individual.local.identifier"] <- "individual_local_identifier"
# 173 segments with 21 variables

names(sub_env_df)
# [1] "date"                        "individual_local_identifier" "unique_segmID"               "uplift_type"                 "burstID"                    
# [6] "ascent_ID"                   "time_start"                  "time_end"                    "duration"                    "w_oro_mean"                 
# [11] "N2_max_h_mean"               "N2_max_value_mean"           "ASHFL_S_mean"                "U_mean"                      "V_mean"                     
# [16] "W_mean"                      "windspeed_mean"              "max_height.ab.gr"            "mean_slope"                  "mean_roughness"             
# [21] "mean_aspect"  

# we could possibly match the sub_env_df with the combined_df with all the behavioural variables by unique_segmID but not necessaery at this precise point
# lab_merged <- merge(combined_df, sub_env_df, by = "unique_segmID")

# ____________________________________
#### Build complete labelled dataset
## run dfa on variable COSMO summarizing the final dataset by unique_segmID
# keep only the shared columns of the 2 labelled dataset
# env cosmo + topography are the predictors

new_lab <- sub_env_df[, c(1, 3, 4, 9:21)] # 173 segments in 16 variables
# names(new_lab)
# [1] "date"              "unique_segmID"     "uplift_type"       "duration"          "w_oro_mean"
# [6] "N2_max_h_mean"     "N2_max_value_mean" "ASHFL_S_mean"      "U_mean"            "V_mean"
# [11] "W_mean"            "windspeed_mean"    "max_height.ab.gr" "mean_slope"        "mean_roughness"
# [16] "mean_aspect"

old_lab <- lab_df[, c(1, 3, 4, 8:20)] # 202 segments in 16 variables
# names(old_lab)
# [1] "date"              "unique_segmID"     "uplift_type"       "duration"          "w_oro_mean"        "N2_max_h_mean"     "N2_max_value_mean"
# [8] "ASHFL_S_mean"      "U_mean"            "V_mean"            "W_mean"            "windspeed_mean"    "max_height.ab.gr"  "mean_slope"       
# [15] "mean_roughness"    "mean_aspect" 

# in the old labelling, thermal was defined as therm -> needed to change label before binding
old_lab$uplift_type[old_lab$uplift_type == "therm"] <- "thermal"

# build complete dataset
complete_lab <- rbind(new_lab, old_lab) # 375 segments in 16 variables

# subset complete dataset to surely classfied uplifts
sub_lab <- complete_lab[complete_lab$uplift_type %in% c("orog", "wave", "thermal"), ] # 254 segments with certain cases of orog, wave, thermal
# set order of levels
sub_lab$uplift_type <- factor(sub_lab$uplift_type, levels = c("orog", "thermal", "wave"))

# complete cases because lda skips lines if NAs in the predictors columns: empty wind values where mismatch between height and dem
sub_lab2 <- sub_lab[complete.cases(sub_lab), ] # 243 segments with complete predictors values 

# at this point sub_lab2 has segments with all durations: from 1 to 305 seconds
quantile(sub_lab2$duration, seq(0, 1, 0.01))
# 0%        1%        2%        3%        4%        5%        6%        7%        8%        9%       10%       11%       12%       13%       14%       15% 
#   1.00000   1.41958   2.00000   2.00000   3.00000   3.00000   4.00000   4.00000   5.00000   5.00000   5.00000   5.00062   6.00000   7.00000   7.88000   9.00000 
# 16%       17%       18%       19%       20%       21%       22%       23%       24%       25%       26%       27%       28%       29%       30%       31% 
#   10.72000  12.00000  13.56000  14.00000  14.00000  14.82000  15.00000  15.66000  17.00000  18.50000  19.00000  20.34000  22.00000  23.00000  23.00060  24.02000 
# 32%       33%       34%       35%       36%       37%       38%       39%       40%       41%       42%       43%       44%       45%       46%       47% 
#   25.88056  27.00000  28.28000  29.00000  30.35912  33.54000  35.96000  36.38000  38.00000  39.00000  40.64036  41.06000  43.00000  43.90090  46.32000  47.00000 
# 48%       49%       50%       51%       52%       53%       54%       55%       56%       57%       58%       59%       60%       61%       62%       63% 
#   49.00000  50.57958  52.00000  53.42000  54.84000  57.00000  58.68000  59.10090  60.00000  61.00000  62.00000  62.78022  64.00000  65.00000  66.99904  70.46000 
# 64%       65%       66%       67%       68%       69%       70%       71%       72%       73%       74%       75%       76%       77%       78%       79% 
#   71.87924  73.00000  74.00000  75.28000  78.00000  78.00000  79.80060  81.00082  83.24076  85.32000  87.00000  89.00050  90.92000  92.34000  97.04000 102.18000 
# 80%       81%       82%       83%       84%       85%       86%       87%       88%       89%       90%       91%       92%       93%       94%       95% 
#   104.20000 106.02000 108.44000 109.00000 115.56000 118.99970 129.00000 130.00000 132.95996 135.38000 137.59920 144.42000 156.64064 170.24000 181.76000 190.80180 
# 96%       97%       98%       99%      100% 
# 206.64000 217.48000 238.40000 291.70000 305.00000 


# filter on duration:
sub_lab_20s <- sub_lab2[sub_lab2$duration > 20, ] # 177 segments: 26% segments below 20sec
sub_lab_10s <- sub_lab2[sub_lab2$duration > 10, ] # 204 segments: 15% segments below 10sec
sub_lab_7s <- sub_lab2[sub_lab2$duration > 7, ] # 209 segments: 13% segments below 7 sec

saveRDS(sub_lab_7s, "labelled209segmentsabove7s_predictors_150726.rds")
# THIS IS THE COMPLETE DATASET MANUALLY LABELLED

table(sub_lab_7s$uplift_type)
# orog thermal    wave 
# 57      97      55 

##### 2. DFA MODELS VALIDATION #####
# this chunk of code uses the manually labelled dataset of 209 segments for 10 DFA models as validation step

# ___________________
#### DFA MODELS
# if yet not loaded:
# sub_lab_7s <- readRDS("./labelled209segmentsabove7s_predictors_150726.rds")

# predictors are all except unique_segmID, date, uplift_type, duration (which is only a filter for now). removed U and V from predictors bc same info than horizontal windspeed
predictors <- names(sub_lab_7s)[c(5:8, 11:16)]
# [1] "w_oro_mean"        "N2_max_h_mean"     "N2_max_value_mean" "ASHFL_S_mean"      "W_mean"            "windspeed_mean"    "max_height.ab.gr"
# [8] "mean_slope"        "mean_roughness"    "mean_aspect"

# 1st DFA MODEL for validation : applied to the manually labelled dataset to see if dfa on same predictors automatically can do the job
# 10fold cross-validation

set.seed(12)

# sample on 10 fold so that each segment enters the test dataset only once
folds <- createFolds(sub_lab_7s$uplift_type, k = 10, list = TRUE)

prob_ls <- lapply(1:10, function(i) {
  test_idx <- folds[[i]]
  train <- sub_lab_7s[-test_idx, ]
  test  <- sub_lab_7s[test_idx, ]
  
  # with i=1
  # table(train$uplift_type)
  # orog thermal    wave 
  # 52      88      49 
  
  # table(test$uplift_type)
  # orog thermal    wave 
  # 5       9       6
  
  lda_model <- lda(as.formula(paste0("uplift_type ~", paste(predictors, collapse = "+"))), data = train)
  predicted <- predict(lda_model, newdata = test)
  
  prob <- as.data.frame(round(predicted$posterior, 2))
  prob$pred <- predicted$class
  prob$true <- test$uplift_type
  prob$segID <- test$unique_segmID   # optional, but useful to confirm no duplicates
  prob$accuracy <- sum(predicted$class == test$uplift_type, na.rm = TRUE) / nrow(prob)
  
  return(prob)
})

## Extract accuracies 
accuracies <- sapply(prob_ls, function(x) {
  accuracy <- unique(x$accuracy)
  return(accuracy)
})

mean_accuracy <- mean(accuracies) # 0.8812554

# Calculate min and max accuracies across runs
min_accuracies <- min(accuracies) # 0.7619048
max_accuracies <- max(accuracies) #1

# Calculate standard error
n <- length(accuracies)
se_accuracy <- sd(accuracies) / sqrt(n)

print(paste("Mean accuracy:", round(mean_accuracy, 4)))
print(paste("Standard error of accuracy:", round(se_accuracy, 4)))
print(paste(
  "Accuracy range (mean ± SE):",
  round(mean_accuracy - se_accuracy, 4), "to",
  round(mean_accuracy + se_accuracy, 4)
))

# [1] "Mean accuracy: 0.8813"
# [1] "Standard error of accuracy: 0.0232"
# [1] "Accuracy range (mean ± SE): 0.858 to 0.9045"

## Calculate class-specific accuracies

calculate_class_accuracies <- function(prob_df) {
  class_accuracies <- sapply(c("orog", "thermal", "wave"), function(class) {
    class_rows <- prob_df$true == class
    sum(prob_df$pred[class_rows] == class) / sum(class_rows)
  })
  names(class_accuracies) <- c("orog", "thermal", "wave")
  return(class_accuracies)
}

# Calculate class-specific accuracies for each run
class_accuracies_list <- lapply(prob_ls, calculate_class_accuracies)

class_accuracies_matrix <- do.call(rbind, class_accuracies_list)

# Calculate mean accuracies and standard errors for each class
mean_class_accuracies <- colMeans(class_accuracies_matrix)
se_class_accuracies <- apply(class_accuracies_matrix, 2, function(x) sd(x) / sqrt(length(x)))

print("Class-specific results:")
for (class in c("orog", "thermal", "wave")) {
  print(paste(class, "- Mean accuracy:", round(mean_class_accuracies[class], 4)))
  print(paste(class, "- Standard error:", round(se_class_accuracies[class], 4)))
  print(paste(
    class, "- Accuracy range (mean ± SE):",
    round(mean_class_accuracies[class] - se_class_accuracies[class], 4), "to",
    round(mean_class_accuracies[class] + se_class_accuracies[class], 4)
  ))
  print("---")
}

# [1] "orog - Mean accuracy: 0.83"
# [1] "orog - Standard error: 0.0498"
# [1] "orog - Accuracy range (mean ± SE): 0.7802 to 0.8798"
# [1] "---"
# [1] "thermal - Mean accuracy: 0.99"
# [1] "thermal - Standard error: 0.01"
# [1] "thermal - Accuracy range (mean ± SE): 0.98 to 1"
# [1] "---"
# [1] "wave - Mean accuracy: 0.7433"
# [1] "wave - Standard error: 0.058"
# [1] "wave - Accuracy range (mean ± SE): 0.6854 to 0.8013"
# [1] "---"

# Calculate min and max accuracies for each class
min_class_accuracies <- apply(class_accuracies_matrix, 2, min)
max_class_accuracies <- apply(class_accuracies_matrix, 2, max)

# Print results
print("Class-specific results:")
for (class in c("orog", "thermal", "wave")) {
  print(paste(class, "- Mean accuracy:", round(mean_class_accuracies[class], 4)))
  print(paste(
    class, "- Accuracy range:",
    round(min_class_accuracies[class], 4), "to",
    round(max_class_accuracies[class], 4)
  ))
  print("---")
}

# [1] "orog - Mean accuracy: 0.83"
# [1] "orog - Accuracy range: 0.6667 to 1"
# [1] "---"
# [1] "thermal - Mean accuracy: 0.99"
# [1] "thermal - Accuracy range: 0.9 to 1"
# [1] "---"
# [1] "wave - Mean accuracy: 0.7433"
# [1] "wave - Accuracy range: 0.4 to 1"
# [1] "---"


# ____________________________________________
#### Visualization of pooled across 10 model

# Pool all predictions across the 10 runs
all_predictions <- do.call(rbind, prob_ls)

## 1: Ternary plot
color_palette <- c(
  "thermal" = "#D55E00",
  "orog" = "#CC79A7", 
  "wave" = "#0072B2"
)

q <- ggtern(all_predictions, aes(orog, thermal, wave, color = true)) +  #put prob_ls[[5]] instead of all_predictions if you want the predictions and visualization of a single run
  geom_point(size = 15, alpha = 0.7) +
  scale_color_manual(
    values = color_palette,
    labels = c(
      "orog" = "Orographic",
      "thermal" = "Thermal",
      "wave" = "Wave"
    )
  ) +
  labs( # title = "Ternary Plot of Orographic, Thermal, and Wave Uplifts",
    # subtitle = "cosmo and topography predictors",
    color = "Uplift Type"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 40),
    plot.margin = unit(c(0.6, 0.4, 0.3, 0.3), "cm"), # increase all margins, especially the left margin where written orographic
    tern.axis.title.L = element_text(hjust = 0.5, vjust = 1.5, size = 20), # adjust left label position
    tern.axis.title.R = element_text(hjust = 0.5, vjust = 1.5, size = 20), # right label position
    tern.axis.title.T = element_text(hjust = 0.5, vjust = -0.5, size = 20) # center the top label
  ) +
  Llab("Orographic") +
  Tlab("Thermal") +
  Rlab("Wave") #+
# theme_showarrows()
print(q)
# ignore warning, it is just internal to ggtern but nothign is wrong
# Ignoring unknown labels:
#   • R : "Wave"
# • Rarrow : "Wave"
# • T : "Thermal"
# • Tarrow : "Thermal"
# • L : "Orographic"
# • Larrow : "Orographic"

ggsave(file.path(directory, "figures_july26", "tern_dfa.pdf"), 
       plot = q, width = 297, height = 210, units = "mm", device = "pdf")

# save ggplot obj to reload if some fig tuning is necessary
saveRDS(q, file.path(directory, "figures_july26", "tern_dfa.rds"))
# q <- readRDS(file.path(directory, "figures_july26", "tern_dfa.rds"))
# load libraries ggplot2 and ggtern
# print(q)
# q <- q + theme(text = element_text(size = 20)) # here modify as usual ggplot



## 2: Confusion matrix plot

# cm <- confusionMatrix(prob_ls[[5]]$pred, prob_ls[[5]]$true) # if only 1 models output
# One aggregate confusion matrix based on all runs
cm <- confusionMatrix(all_predictions$pred, all_predictions$true)

# Define the color palette with different Labels for the uplift names
color_palette_lab <- c("Orographic" = "#CC79A7", 
                       "Thermal" =  "#D55E00", 
                       "Wave" = "#0072B2")

# Melt the confusion matrix for ggplot
melted_cm <- as.data.frame(as.table(cm))
names(melted_cm) <- c("Prediction", "Reference", "Value")

# Update labels
melted_cm$Prediction <- factor(melted_cm$Prediction,
                               levels = c("orog", "thermal", "wave"),
                               labels = c("Orographic", "Thermal", "Wave")
)
melted_cm$Reference <- factor(melted_cm$Reference,
                              levels = c("orog", "thermal", "wave"),
                              labels = c("Orographic", "Thermal", "Wave")
)

# Calculate total predictions per reference class
total_per_reference <- aggregate(Value ~ Reference, data = melted_cm, sum)

melted_cm <- merge(melted_cm, total_per_reference, by = "Reference", suffixes = c("", "_Total")) # merge total counts back into melted_cm

# Calculate percentage for color scaling based on correct predictions per column
melted_cm$Percentage <- melted_cm$Value / melted_cm$Value_Total * 100

p <- ggplot(melted_cm, aes(x = Reference, y = Prediction)) +
  geom_tile(aes(fill = Reference, alpha = Percentage), color = "black") +
  geom_text(aes(label = Value), size = 20) + # size for text labels into squares
  scale_fill_manual(values = color_palette_lab) +
  scale_alpha_continuous(range = c(0.1, 1)) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold"),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    # plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 30, face = "bold"),
    panel.grid = element_blank() # Remove background grid
  ) +
  labs(
    # title = "Confusion Matrix",
    x = "References",
    y = "Prediction",
    fill = "Class",
    alpha = "Percentage"
  ) +
  coord_fixed() +
  scale_y_discrete(limits = rev(levels(melted_cm$Prediction))) # reverse y-axis to correct diagonal

# Add connecting lines
p <- p +
  # all borders
  geom_segment(aes(x = 0.5, xend = 3.5, y = 3.5, yend = 3.5), color = "black", size = 0.5) + # top border
  geom_segment(aes(x = 0.5, xend = 0.5, y = 0.5, yend = 3.5), color = "black", size = 0.5) + # left border
  geom_segment(aes(x = 3.5, xend = 3.5, y = 0.5, yend = 3.5), color = "black", size = 0.5) + # right border
  geom_segment(aes(x = 0.5, xend = 3.5, y = 0.5, yend = 0.5), color = "black", size = 0.5) + # bottom border
  # internal grid lines
  geom_segment(aes(x = 1.5, xend = 1.5, y = 0.5, yend = 3.5), color = "black", size = 0.5) + # vertical internal lines
  geom_segment(aes(x = 2.5, xend = 2.5, y = 0.5, yend = 3.5), color = "black", size = 0.5) +
  geom_segment(aes(x = 0.5, xend = 3.5, y = 1.5, yend = 1.5), color = "black", size = 0.5) + # horizontal internal lines
  geom_segment(aes(x = 0.5, xend = 3.5, y = 2.5, yend = 2.5), color = "black", size = 0.5)


print(p)

ggsave(file.path(directory, "figures_july26", "cm_dfa.pdf"), 
        plot = p, width = 297, height = 210, units = "mm", device = "pdf")

# save ggplot obj to reload if some fig tuning is necessary
saveRDS(p, file.path(directory, "figures_july26", "cm_dfa.rds"))
p <- readRDS("/home/francesca/ownCloud/TesiFrancesca_UpliftClassification/figures_july26/cm_dfa.rds")

#check which layer index is your geom_text layer
p$layers # layer 2
# Modify the size parameter of that layer directly
p$layers[[2]]$aes_params$size <- 10 # changed to same font size used for rf cm
print(p)

ggsave(file.path(directory, "figures_july26", "cm_dfa_10font.pdf"), 
       plot = p, width = 297, height = 210, units = "mm", device = "pdf")


##### 3. DFA MODEL ON THE COMPLETE BEHAVIOURAL DATASET #####
# this chunk of code uses the manually labelled dataset of 209 segments to run a DFA model to classify the complete behavioural dataset

# ______________________________________________________
#### prepare complete behavioural dataset to run the DFA

## now use the complete sub_lab_7s to predict my complete dataset
# combined_df is the complete dataset with 29775 segments and all behavioral variables, but not used here
# here, taken complete dataset with env data from COSMO and removed the segments that have already the labelling by Tom

# env_cosmo from above

comp_gps <- env_cosmo[!env_cosmo$unique_segmID %in% sub_lab_7s$unique_segmID, ] # 572905
# removed from the complete dataset only the certainly labelled segments, meaning that Tom' dynamic, unknown etc are left to be labelled by the automatized algorithm
# if want to be conservative, g <- env_cosmo[!env_cosmo$unique_segmID %in% new_lab$unique_segmID,]

names(comp_gps)
# [1] "event_id"                    "X"                           "Unnamed..0.1"                "index"                      
# [5] "Unnamed..0"                  "height_above_ellipsoid"      "deployment_id"               "timestamp"                  
# [9] "timelag"                     "turnAngle"                   "stepLength"                  "azimuth"                    
# [13] "geoid_egm96"                 "dem_30m"                     "height.above.sea.level"      "height_agl"                 
# [17] "burstID"                     "grSpeed"                     "vertDist"                    "vertSpeed"                  
# [21] "individual_local_identifier" "tag_id"                      "individual_id"               "study_id"                   
# [25] "taxon_canonical_name"        "tag_local_identifier"        "location.long"               "location.lat"               
# [29] "date"                        "time"                        "vertSpeed_smooth"            "location_long_smooth"       
# [33] "location_lat_smooth_new"     "height_asl_smooth"           "new_class"                   "segmID"                     
# [37] "segmToKeep"                  "unique_pointID"              "unique_segmID"               "area"                       
# [41] "dist"                        "h_iniz"                      "h_fin"                       "deltah"                     
# [45] "duration"                    "turnChange"                  "vol"                         "rlon"                       
# [49] "rlat"                        "height.above.surf"           "N2_max_h"                    "N2_max_value"               
# [53] "w_oro"                       "ASHFL_S"                     "V"                           "U"                          
# [57] "W"   

# prepare the dataset with the variables needed to run the dfa:
# wind speed
comp_gps$windspeed <- wind.speed(comp_gps$U, comp_gps$V)

# calculate topography on the points annotated by Tom, associated to unique_segmID
coords <- comp_gps[, c("location.long", "location.lat")]
coords <- st_as_sf(coords, coords = c("location.long", "location.lat"), crs = "+proj=longlat +ellps=WGS84")
coords_proj <- st_transform(coords, crs(dem))

comp_gps$slope <- terra::extract(slope, coords_proj, method = "bilinear", cells = F, ID = F)[, 1]
comp_gps$aspect <- terra::extract(aspect, coords_proj, method = "bilinear", cells = F, ID = F)[, 1]
comp_gps$roughness <- terra::extract(roughness, coords_proj, method = "bilinear", cells = F, ID = F)[, 1]
terrainElevation <- terra::extract(dem, coords_proj, method = "bilinear", cells = F, ID = F)[, 1]
comp_gps$height_above_ground <- comp_gps$height.above.sea.level - terrainElevation

comp_gps$aspect_rad <- circular(comp_gps$aspect * pi / 180, units = "radians")

comp_gps$timestamp <- as.POSIXct(comp_gps$timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
comp_gps$date <- as.Date(comp_gps$timestamp)

# # filter on duration: trials
# length(unique(comp_gps$unique_segmID)) # 33984
# quantile(comp_gps$duration, seq(0, 1, 0.01)) # from 3 seconds to 776 seconds
# # 10% data points are below 21 seconds, 3% below 8 sec

# duration_df is the dataset I use to check the duration distribution, grouped the points by unique_segmID
duration_df <- comp_gps %>%
  group_by(unique_segmID, duration) %>%
  summarise(time_start = min(timestamp), .groups = "drop")

quantile(duration_df$duration, seq(0, 0.7, 0.01)) # 47% of segments are below 20 sec
# plot(q) if q=quantile()
# put the threshold at 7 sec: 1st quantile
summary(duration_df$duration)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 3.00    7.00   21.00   48.74   61.00  776.00

# filter above 7 sec
comp_gps_7s <- comp_gps[comp_gps$duration > 7, ] # 550045 

# many NAs are in the columns rlong and rlat, even if not used later I can remove them here
comp_gps_7s <- comp_gps_7s[, !names(comp_gps_7s) %in% c("rlon", "rlat")]

# check on one random segm to see if it looks good
check <- comp_gps_7s[comp_gps_7s$unique_segmID=="Adamello20 (eobs 7548)_1187719633_1856_soar_7",]

# there are some segments without U and V because of mismatch with height above surface, so quick check (15 july)
# checked if specific correlation with height above surface
# summary(comp_gps_7s$height.above.surf[is.na(comp_gps_7s$U)]) 
# 
# quantile(
#   comp_gps_7s$height.above.surf[!is.na(comp_gps_7s$U)],
#   probs = seq(0, 1, 0.1),
#   na.rm = TRUE
# )
# 
# quantile(
#   comp_gps_7s$height.above.surf[is.na(comp_gps_7s$U)],
#   probs = seq(0, 1, 0.1),
#   na.rm = TRUE
# )

segm_ls <- split(comp_gps_7s, comp_gps_7s$unique_segmID) # 25173 segments 

comp_env_df <- as.data.frame(rbindlist(
  lapply(segm_ls, function(b) {
    # b <- segm_ls[[5]]
    
    summary_caselab <- b %>%
      group_by(
        date, individual_local_identifier, unique_segmID,
        burstID, duration
      ) %>%
      summarise(
        time_start = min(timestamp, na.rm = TRUE),
        time_end = max(timestamp, na.rm = TRUE),
        w_oro_mean = mean(w_oro, na.rm = TRUE),
        N2_max_h_mean = mean(N2_max_h, na.rm = TRUE),
        N2_max_value_mean = mean(N2_max_value, na.rm = TRUE),
        ASHFL_S_mean = mean(ASHFL_S, na.rm = TRUE),
        U_mean = mean(U, na.rm=T),
        V_mean = mean(V, na.rm = TRUE),
        W_mean = mean(abs(W), na.rm = TRUE),
        windspeed_mean = mean(windspeed, na.rm = TRUE),
        max_height.ab.gr = max(height_above_ground, na.rm = TRUE),
        mean_slope = mean(slope, na.rm = TRUE), # slope is not a circular variable
        mean_roughness = mean(roughness, na.rm = TRUE),
        mean_aspect = as.numeric(mean.circular(aspect_rad), na.rm = TRUE) * 180 / pi, # circular variable, calculated in rad and back to degree
        
        .groups = "drop"
      )
    
    return(summary_caselab)
  })
))

# 5 warnings on height above ground bc height above surface is NA, but only 5 segments
segm_ls_check <- rbindlist(lapply(segm_ls, function(b) {
  b %>%
    group_by(unique_segmID) %>%
    summarise(n_total = n(),
              n_na_height = sum(is.na(height_above_ground)),
              .groups = "drop")
}))

segm_ls_check %>% filter(n_na_height == n_total)
# unique_segmID n_total n_na_height
# <char>   <int>       <int>
# 1:  Ettenberg22 (eobs 10539)_2191886024_14164_soar_1      65          65
# 2: Ettenberg22 (eobs 10539)_2191886024_22178_soar_11      38          38
# 3:  Ettenberg22 (eobs 10539)_2191886024_22178_soar_5      15          15
# 4:  Ettenberg22 (eobs 10539)_2191886024_22178_soar_7      57          57
# 5:  Ettenberg22 (eobs 10539)_2191886024_22178_soar_9       6           6

# for these 5 sgements, also U, V, windspeed cannot be annotated, but we know this issue, already addressed and corrected at the best possible on 15 july 26

# remove obs with NA values
compsegm_df <- comp_env_df[complete.cases(comp_env_df), ] # 25161 segments

# ________________________________________________________________
### run DFA with environmental variables on my complete dataset ###

set.seed(12) # for repeatability

train <- sub_lab_7s
# table(train$uplift_type)
# orog thermal    wave
# 57      97      55
new_data <- compsegm_df

# if not yet selected:
# predictors<- names(sub_lab_7s)[c(5:8,11,12)] # these are only COSMO predictors, see above for correct ones 
# [1] "w_oro_mean"        "N2_max_h_mean"     "N2_max_value_mean" "ASHFL_S_mean"      "W_mean"            "windspeed_mean"   
# [7] "max_height.ab.gr"  "mean_slope"        "mean_roughness"    "mean_aspect"

# Fit a LDA model
lda_model <- lda(as.formula(paste0("uplift_type ~", paste(predictors, collapse = "+"))), data = train)
predicted <- predict(lda_model, newdata = new_data)
# head(predicted)

prob <- as.data.frame(round(predicted$posterior, 2))
prob$pred <- predicted$class
prob$unique_segmID <- new_data$unique_segmID

# table(prob$pred)
# orog thermal    wave 
# 2271   20055    2835 

segm_pred <- left_join(compsegm_df, prob, by = "unique_segmID")

# here all segments associated with dfa prediction
saveRDS(segm_pred, "./uplift_classification/summary_segments_cosmo_topography_behav_allpredictions_25161.rds")

saveRDS(lda_model, "./uplift_classification/lda_classificationmodel.rds")
# lda_model
# Call:
#   lda(as.formula(paste0("uplift_type ~", paste(predictors, collapse = "+"))), 
#       data = train)
# 
# Prior probabilities of groups:
#   orog   thermal      wave 
# 0.2727273 0.4641148 0.2631579 
# 
# Group means:
#   w_oro_mean N2_max_h_mean N2_max_value_mean ASHFL_S_mean    W_mean windspeed_mean max_height.ab.gr mean_slope mean_roughness mean_aspect
# orog     0.7896728      791.1152      0.0004215197     25.67038 0.7952738      14.527913         282.6331   31.72030       49.59892   -44.35351
# thermal  0.1369945     1396.5903      0.0002344330   -169.42927 0.3870771       3.549477         647.9060   27.05461       42.15457    23.65683
# wave     0.1892774      510.6354      0.0004616639    -10.28729 1.1739261      18.116528         608.5976   29.43033       45.54544    13.01089
# 
# Coefficients of linear discriminants:
#   LD1           LD2
# w_oro_mean        -1.245036e+00 -1.700464e+00
# N2_max_h_mean      4.006002e-04 -5.537535e-04
# N2_max_value_mean  3.225509e+02  5.541515e+02
# ASHFL_S_mean      -7.183980e-03 -2.169339e-03
# W_mean             4.536675e-03  4.334402e-01
# windspeed_mean    -1.266609e-01  5.711407e-02
# max_height.ab.gr   2.991572e-04  1.146219e-03
# mean_slope        -8.563111e-02  5.841834e-02
# mean_roughness     2.973083e-02 -3.911124e-02
# mean_aspect        1.663126e-03  1.507245e-03
# 
# Proportion of trace:
#   LD1    LD2 
# 0.8572 0.1428 

# ____________________________________________________________________
# Define surely and mixed-classes predicted uplift classes based on DFA probabilities ####

# if not yet loaded:
# segm_pred <- readRDS("./uplift_classification/summary_segments_cosmo_topography_behav_allpredictions_25161.rds")

# extract month and year
segm_pred$month <- month(segm_pred$date, label = TRUE, abbr = FALSE, locale = "en_US.UTF-8")
segm_pred$year <- year(segm_pred$date)

# order by month as factor
segm_pred$month <- factor(segm_pred$month, levels = month.name, ordered = TRUE)

# Define uplift types with my thresholds
segm_pred$uplift_type <- NA
segm_pred <- segm_pred %>%
  mutate(uplift_type = case_when(
    pred == "wave" & wave >= 0.8 ~ "wave",
    pred == "orog" & orog >= 0.8 ~ "orog",
    pred == "thermal" & thermal >= 0.8 ~ "thermal",
    thermal > 0.4 & orog > 0.4 ~ "thermal/orog",
    thermal > 0.4 & wave > 0.4 ~ "thermal/wave",
    orog > 0.4 & wave > 0.4 ~ "orog/wave",
    TRUE ~ "unknown"
  ))

table(segm_pred$uplift_type)
# orog    orog/wave      thermal thermal/orog thermal/wave      unknown         wave 
# 1019          318        17167          183          322         5303          849

segm_pred_sure <- segm_pred[segm_pred$uplift_type %in% c("orog", "wave", "thermal"), ] # 19035 segments 
table(segm_pred_sure$uplift_type)
# orog thermal    wave 
# 1019   17167     849 

saveRDS(segm_pred_sure, "./uplift_classification/segm_pred_sure_19035.rds")
# segm_pred_sure is my certain classification with which I will run the next analyses

# summary(segm_pred_sure$duration)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 7.001  17.000  36.000  64.862  85.000 776.000 



