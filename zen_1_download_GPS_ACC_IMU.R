# This code downloads all available IMU data from the LifeTrack Golden Eagle Alps study
# And the corresponding ACC, IMU and GPS data for those same individuals
# Martina Scacco & Francesca Frisoni - May 16, 2024. Konstanz

library(move2)
library(dplyr)

#________________________
# Explore dataset

keyring::key_list()
# movebank_remove_credentials(key_name = "")
options("move2_movebank_key_name" = "movebank")

# Study Info:
movebank_download_study_info() %>% 
  select(id, name, number_of_deployed_locations) %>% 
  filter(grepl("LifeTrack Golden Eagle", name)) %>%
  print(n = Inf)
# study name LifeTrack Golden Eagle Alps | study id 282734839

# select individuals that have IMU data
metadata <- movebank_download_deployment(study_id = 282734839)
length(unique(metadata$individual_local_identifier)) # 98 individuals, at the time of our download
length(unique(metadata$individual_local_identifier[metadata$sensor_type_ids == "magnetometer,orientation,acceleration,gps"])) # but only 32 with IMU


#### DATA DOWNLOAD ####

setwd("/home/francesca/ownCloud/TesiFrancesca_UpliftClassification/Data")

#__________________
#### Download IMU 

mag <- movebank_retrieve(study_id = 282734839, sensor_type_id = "magnetometer", 
                         entity_type = "event",  attributes = "all")
length(unique(mag$individual_local_identifier))
mag$individual_local_identifier <- as.character(mag$individual_local_identifier)
mag %>% group_by(individual_local_identifier) %>% reframe(range(timestamp)) %>% print(n = Inf) 

quat <- movebank_retrieve(study_id = 282734839, sensor_type_id = "orientation", 
                          entity_type = "event",  attributes = "all")
length(unique(quat$individual_local_identifier))

# Check that all individuals and all timestamps of mag and quat match
table(unique(mag$individual_local_identifier) %in% unique(quat$individual_local_identifier))
table(mag$timestamp %in% quat$timestamp)

# Put the quat and mag together
imu <- mag[,-grep("sensor",names(mag))] %>% 
  full_join(quat[,-grep("event|sensor", names(quat))], by = c("individual_local_identifier", "tag_local_identifier", "study_id", "tag_id", "deployment_id","individual_taxon_canonical_name", "timestamp",
                                                              "eobs_start_timestamp", "eobs_key_bin_checksum", "data_decoding_software", "individual_id", "import_marked_outlier", "visible")) %>% 
  as.data.frame()
rm(mag, quat)

# Save IMU data
saveRDS(imu, "IMU_allEagles_20240516.rds")

# Extract names of individuals with IMU data
imu_inds <- unique(imu$individual_local_identifier)[!is.na(unique(imu$individual_local_identifier))]
saveRDS(imu_inds, "IMU_indLocIdentifier.rds")

# Check time range per individual (IMU data were collected already from 2020)
imu_range <- imu %>%
  group_by(individual_local_identifier) %>%
  summarise(time_start=min(timestamp), time_end=(max(timestamp)))
write.csv(imu_range, "IMU_timeRanges_perIndividual.csv")

rm(imu); gc() # clean memory

## Download each individual separately if too heavy for computer
# We could download IMU without binding to quat because not necessary for our purposes
imu_inds <- readRDS("IMU_indLocIdentifier.rds")
list_inds <- unlist(strsplit(imu_inds,'"\\s+"'))  

# Create a directory to store the CSV files
dir.create("individual_magn_data", showWarnings = FALSE)

# Use lapply to run movebank_retrieve for each individual
magn_ind <- lapply(list_inds, function(individual) {
  tryCatch({
    # individual <- list_inds[1]
    data <- movebank_retrieve(
      study_id = 282734839,
      sensor_type_id = "magnetometer",
      entity_type = "event",
      attributes = "all",
      individual_local_identifier = individual
    )
    
    # Convert the data to a data frame (if it's not already)
    if (inherits(data, "move")) {
      data_df <- as.data.frame(data)
    } else {
      data_df <- data
    }
    
    # safe filename: unique for every individual and without empty spaces
    safe_filename <- gsub("[^a-zA-Z0-9]", "_", individual)
    
    # Save the data as a CSV
    csv_filename <- file.path("individual_magn_data", paste0(safe_filename, ".csv"))
    write.csv(data_df, file = csv_filename, row.names = FALSE)
    
    message("Data for ", individual, " saved to ", csv_filename)
    
    # Return the data (optional, in case you want to use it further in R)
    # return(data)
  }, error = function(e) {
    message("Error processing individual: ", individual, "\n", e$message)
    return(NULL)
  })
})

#___________________
#### Download ACC 

acc <- movebank_retrieve(study_id = 282734839, 
                         sensor_type_id = "acceleration", 
                         entity_type = "event",  
                         attributes = "all",
                         individual_local_identifier = imu_inds)

# Save ACC data
saveRDS(acc, "ACC_allIMUeagles_20240516.rds")
rm(acc); gc()

#___________________
#### Download GPS 

# Download GPS (only for the individuals having IMU)
mv2 <- movebank_download_study(study_id = 282734839, # Download directly as move2 object
                               sensor_type_id = "gps", 
                               attributes = c("event_id","height_above_ellipsoid"),
                               individual_local_identifier = imu_inds)
mv2 %>% print(width=Inf)

# Save GPS data per individual
dir.create("GPS_allIMUeagles_20240516")
mv2_ls <- split(mv2, mt_track_id(mv2))
lapply(names(mv2_ls), function(ind){
  saveRDS(mv2_ls[[ind]], paste0("GPS_allIMUeagles_20240516/GPS_",ind,".rds"))
})



