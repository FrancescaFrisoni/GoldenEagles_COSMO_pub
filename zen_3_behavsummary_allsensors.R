# This code process the ACC and IMU data combined to GPS points
# and calculates IMU- and ACC-based movement metrics, then summarizes GPS, ACC and IMU per each segment
# Martina Scacco & Francesca Frisoni - June 20, 2024. Konstanz
# IMU functions adapted from Elham Nourani and Kamran Safi

# list of to-dos:
# 1- work with singular individuals 
# 2- list file ACC and IMU for each individual from folder ACC-IMU_allIMUeagles_20240516
# 3- manipulate ACC raw data 
# 4- bind ACC and GPS points (each second)
# 5- manipulate IMU raw data 
# 6- bind IMU and GPS/ACC points (each second)
# 7- summarize segments: GPS, ACC and IMU variables (each unique segment)
# 8- save a final dataset for each individual with GPS-ACC-IMU, and a complete dataset with all individuals together

# Set working directory
setwd("/home/francesca/ownCloud/TesiFrancesca_UpliftClassification/Data")

# Import the customized function libraries 
source("/home/francesca/R_projects/GoldenEagles_COSMO_pub/ACC_functions_EobsUvaBits.R")
source("/home/francesca/R_projects/GoldenEagles_COSMO_pub/IMU_functions_Elham.R")

# Libraries
library(data.table)
library(tidyverse) # actually this packages already contains some of the following ones
library(stringr)
library(purrr)
library(pbapply)
library(dplyr)
library(circular)
library(plyr)
library(doParallel)
library(foreach)
library(units)

# ACC and IMU saved for each individual as RDS in folder ACC-IMU_allIMUeagles_20240516 splitting the files ACC_allIMUeagles_20240516.rds and IMU_allEagles_20240516.rds
list_acc_names <- list.files("ACC-IMU_allIMUeagles_20240516", pattern="ACC")
list_imu_names <- list.files("ACC-IMU_allIMUeagles_20240516", pattern="IMU")
imu_acc_path <- "/home/francesca/ownCloud/TesiFrancesca_UpliftClassification/Data/ACC-IMU_allIMUeagles_20240516"

#### 1. ACC ANALYSIS ####

#____________________
#### ACC raw analysis 
# divide into 3 axis and calculate ODBA, VeDBA, mean and sd

acc_beforesegm <-  lapply(list_acc_names, function(file) {  
  full_path <- file.path(imu_acc_path, file)
  b <- readRDS(full_path)
  acc_vedba_indiv <- extractACCandVeDBA_df(b) # function written in ff_functions
  acc_vedba_indiv$timestamp <- as.POSIXct(as.character(acc_vedba_indiv$timestamp), format="%Y-%m-%d %H:%M:%OS", tz="UTC")
  acc_vedba_indiv$date <- as.Date(acc_vedba_indiv$timestamp) 
  
  # safe filename: unique for every individual and without empty spaces
  individual <- unique(acc_vedba_indiv$individual_local_identifier)
  safe_filename_acc <- gsub("[^a-zA-Z0-9]", "_", individual)
  
  # Save the data as a CSV
  rds_filename <- file.path("ACC-IMU_allIMUeagles_20240516", paste0(safe_filename_acc, "acc_beforesegm", ".rds"))
  saveRDS(acc_vedba_indiv, file = rds_filename)
  
  message("Data for ", individual, " saved to ", rds_filename)
  
})

#___________________________
#### Bind ACC to GPS points 

list_acc_names <- list.files("ACC-IMU_allIMUeagles_20240516", pattern="acc_beforesegm", full.names = TRUE)

gps_complete <- as.data.frame(readRDS("GPS_allsoaringpoints_areavol.rds")) 


lapply(list_acc_names, function(b){
  # b <- list_acc_names[[19]]
  acc <- as.data.frame(readRDS(b))
  gps_b <- gps_complete[gps_complete$individual_local_identifier %in% acc$individual_local_identifier,]
  
  if(nrow(gps_b) > 0){    # there are 8 individuals that have acc and imu but not gps (problem in download?) -> at the end 24 individuals
    # gps_b$timestamp <- as.POSIXct(gps_b$timestamp, format="%Y-%m-%d %H:%M:%S", tz="UTC")
    # acc$timestamp <- as.POSIXct(acc$timestamp, format="%Y-%m-%d %H:%M:%S", tz="UTC")
    # acc$date <- as.Date(acc$timestamp) 
    names(gps_b)[names(gps_b) == "timestamp"] <- "gps_timestamp"
    
    gps_b$utcDate <- NULL
    gps_b$utcTimestamp <- NULL
    # Order both datasets by timestamp
    acc <- acc[order(acc$timestamp),]
    names(acc)[names(acc) == "event_id"] <-"acc_event_id" # need this name to run the ACCtoGPS function
    
    gps_b <- gps_b[order(gps_b$gps_timestamp),]
    
    # Exclude rows with missing acc info
    # acc <- acc[which(complete.cases(acc$VedBA)),]
    
    acc_sub <- acc[, c("individual_local_identifier","timestamp", "acc_event_id","stdvACC_z","VedBA","ODBA")]
    
    acc$burstDuration_sec <- 0.8
    acc_sub$samplFreq_perAxis <- acc$eobs_acceleration_sampling_frequency_per_axis   #20 or 33.3
    acc_sub$acc_burstID <- 1:nrow(acc)
    # If acc is not empty continues with calculating VeDBA and associating it to the GPS dataset
    if(nrow(acc) > 0){  
      # Extract column names to associate to the GPS dataset (you can decide to include only the VeDBA, or all of the following variables)
      accColsToAssociate <- c("VedBA","ODBA","samplFreq_perAxis", "stdvACC_z")#,""burstDuration_sec",pressureInHPA")
      # Define time tolerance (e.g. 5 mins in seconds) to look for the closest ACC information and associate it to the gps data
      timeTolerance <- 30 # chosen 30 seconds
      accGps <- ACCtoGPS(GPSdata = gps_b, ACCdata = acc_sub,
                         timeTolerance = timeTolerance,
                         accEventCol="acc_burstID",
                         ACCtimeCol="timestamp", GPStimeCol="gps_timestamp",
                         ColsToAssociate = accColsToAssociate) #columns of the acc data that you want to associate to the gps data
      
      # safe filename: unique for every individual and without empty spaces
      individual <- unique(acc$individual_local_identifier)
      safe_filename_acc <- gsub("[^a-zA-Z0-9]", "_", individual)
      
      # save the data as a CSV
      rds_filename <- file.path("acc_gps_imu_allindividuals", paste0(safe_filename_acc, "acc_gps", ".rds"))
      saveRDS(accGps, file = rds_filename)
      
      message("Data for ", individual, " saved to ", rds_filename)
      
    }}
})

#### 2. IMU ANALYSIS ####

#_____________________
#### IMU raw analysis 

imu_beforesegm <-  lapply(list_imu_names, function(file) {
  # file <- "IMU_Mals1_20 (eobs 7560)"
  full_path <- file.path(imu_acc_path, file)
  b <- readRDS(full_path)
  
  options(digits.secs=5) # keeps the milliseconds 
  b$imu_timestamp <- as.POSIXct(b$eobs_start_timestamp, format = "%Y-%m-%d %H:%M:%OS")
  
  # Open orientation dataset, subset for data with high sampling frequency 
  or <- b %>%
    # group_by(individual.local.identifier) %>% not needed when run for only 1 individual, and then add .by_group = TRUE into the arrange below
    arrange(imu_timestamp) %>%
    mutate(timelag = c(0, diff(as.numeric(imu_timestamp))),  #calculate time lag
           imu_burst_id = cumsum(timelag > 1)) %>% #assign a unique burst ID every time timelag is larger than 1 second. The IDs will be unique only within one ind's data
    ungroup() %>%
    as.data.frame()
  
  # to get an idea of the quaternion data input
  # table(nchar(or$orientation.quaternions.raw))
  # table(unlist(sapply(strsplit(or$orientation.quaternions.raw," "), length)))
  
  # Subset or by the length of the bursts. Only keep those that are longer than 3 second
  burst_size <- or %>% 
    group_by(imu_burst_id) %>%  # add 'individual.local.identifier,' in the group_by if more than 1 individual
    dplyr::summarize(#burst_length = n(), #This is not very informative bc the length is not in seconds
      burst_duration = sum(timelag[-1]),
      min_timestamp = min(imu_timestamp),
      max_timestamp = max(imu_timestamp))%>% #the first value of timelag in each burst is a large value indicating timediff between this and the previous burst. remove it
    ungroup()
  
  # # to see distribution of burst durations
  # q=quantile(burst_size$burst_duration, seq(0,1,0.01))
  # plot(q)
  # tail(q, 20)
  # table(burst_size$burst_duration > 8)
  
  # SUPERBURST OF IMU CAN LAST UP TO MINUTES, so filter only to remove segments too short
  or_hfreq <- or %>%
    right_join(burst_size %>% filter(burst_duration > 3)) 
  
  # rm(or);gc()
  
  # to see how much time between different imu_bursts: 28 sec at least from eachother (important when matched with gps)
  # burstTL <- difftime(burst_size$max_timestamp[-1],burst_size$min_timestamp[-nrow(burst_size)], units = "secs")
  # summary(as.numeric(burstTL))
  
  #calculate pitch, roll, and yaw for the whole dataset
  # (st_time <- Sys.time()) # to see how long it takes on the computer
  
  deg <- function(x) x * 180 / pi # function to convert angles from radians to degrees
  
  or_angles <- or_hfreq %>%
    dplyr::group_by(imu_timestamp) %>%
    dplyr::summarise(
      across(c(individual_local_identifier, tag_local_identifier, timestamp, tag_id, imu_burst_id, burst_duration), first),  # keep these columns by taking first value per group bc repeated 
      
      orientation_quaternions_raw_combined = str_c(orientation_quaternions_raw, collapse = " "),
      
      pitch_rad = process_quaternions(orientation_quaternions_raw_combined, function(x) get.pitch(x, type = "eobs")),
      yaw_rad   = process_quaternions(orientation_quaternions_raw_combined, function(x) get.yaw(x, type = "eobs")),
      roll_rad  = process_quaternions(orientation_quaternions_raw_combined, function(x) get.roll(x, type = "eobs")),
      
      # Convert radians strings to degree strings *inside summarise()*
      pitch_deg = {
        vals <- as.numeric(str_split(pitch_rad, " ", simplify = TRUE))
        deg_vals <- deg(vals)
        str_c(format(deg_vals, digits = 6), collapse = " ")
      },
      yaw_deg = {
        vals <- as.numeric(str_split(yaw_rad, " ", simplify = TRUE))
        deg_vals <- deg(vals)
        str_c(format(deg_vals, digits = 6), collapse = " ")
      },
      roll_deg = {
        vals <- as.numeric(str_split(roll_rad, " ", simplify = TRUE))
        deg_vals <- deg(vals)
        str_c(format(deg_vals, digits = 6), collapse = " ")
      }
    ) %>%
    ungroup() %>%
    as.data.frame()
  
  # ________________________________
  ## Summarize data over each second
  #modify the data to have one row per second. this step is not absolutely necessary
  #Useful for summarizing the quat angles for each second. This way I can later summarize the values for each burst, depending on the burst length of interest
  
  or_seconds <- or_angles %>% 
    group_by(second= floor_date(imu_timestamp, "second")) %>% # individual_local_identifier in group_by if multiple inds
    #create a long character string containing all values for each second
    dplyr::summarize(across(c(roll_rad, pitch_rad, yaw_rad, roll_deg, pitch_deg, yaw_deg), ~paste(., collapse = " ")), 
                     #retain the first value of other important columns
                     across(c(individual_local_identifier, 
                              tag_local_identifier, timestamp, tag_id, imu_burst_id, burst_duration), ~head(.,1)),
                     .groups = "keep") %>%
    ungroup() %>%
    as.data.frame()
  
  #for each second, calculate the mean, min, max, sd, cumsum for each axis
  
  seconds_summaries <- or_seconds %>%
    dplyr::group_by(second) %>%   # adjust if needed to your 'second' grouping column
    
    dplyr::summarise(
      yaw_summary = list(angle_summaries(strings_to_numeric(yaw_deg))),
      pitch_summary = list(angle_summaries(strings_to_numeric(pitch_deg))),
      roll_summary = list(angle_summaries(strings_to_numeric(roll_deg))),
      
      # If you want, keep one representative row for other columns via first()
      # individual_local_identifier = first(individual_local_identifier),
      # tag_local_identifier = first(tag_local_identifier),
      
      .groups = "drop"
    ) %>%
    unnest_wider(yaw_summary, names_sep = "_yaw") %>%
    unnest_wider(pitch_summary, names_sep = "_pitch") %>%
    unnest_wider(roll_summary, names_sep = "_roll") %>%
    as.data.frame()
  
  
  #______________________________________________
  ## Find nearest GPS fix to each quat/mag burst 
  
  #open acc_gps data: all points but already segmented with segmentID
  ## this snippet guarantees that to the IMU calculation is assigned the right ACC-GPS for binding
  
  individual <- unique(b$individual_local_identifier)
  print(individual) # visual check: when running the code, printed to see that it actually reached this point and it is starting the merge
  
  safe_filename_acc <- gsub("[^a-zA-Z0-9]", "_", individual)
  
  # List all acc_gps files 
  ind_list <- list.files("./acc_gps_imu_allindividuals/", pattern = "acc_gps", full.names = TRUE)
  
  # Filter files that start with safe_filename_acc AND end exactly with "acc_gps.rds"
  candidate_files <- ind_list[
    startsWith(basename(ind_list), safe_filename_acc) & 
      grepl("acc_gps\\.rds$", basename(ind_list), ignore.case = TRUE)
  ]
  
  if (length(candidate_files) == 0) {
    stop(paste("No suitable file ending with 'acc_gps.rds' found for individual:", safe_filename_acc))
  } else if (length(candidate_files) > 1) {
    warning(paste("Multiple matching files found for individual:", safe_filename_acc, 
                  "- using the first one:", candidate_files[1]))
    candidate_files <- candidate_files[1]
  }
  
  file_to_load <- candidate_files
  
  # Load the file accordingly, remember that gps_data is acc and gps already binded together above
  gps_data <- readRDS(file_to_load)
  
  # or_w_6gps <- find_closest_gps(or_seconds, gps_data) # if you want to bind the gps to imu dataset, Elham's function
  
  # this is my new function to connect imu to the acc-gps dataset: gps_with_or is the combination of the 3
  gps_with_or <- find_closest_imu_for_gps(gps_data,or_seconds)
  # for each second, one row of data with gps info, acc values, and imu closest timestamp and angles
  
  #________________________________________
  ### Bind IMU summaries to ACC-GPS dataset
  # bind one-second orientation data with summarized values per second with the matched GPS-Orientation data 
  
  # Merge using gps_data$gps_timestamp and seconds_summaries$second
  joined_data <- merge(
    x = gps_with_or,
    y = seconds_summaries,
    by.x = "gps_timestamp",
    by.y = "second",
    all.x = TRUE   # keep all GPS rows, add matching IMU summaries or NA if no match
  )
  
  joined_data <- joined_data %>%
    dplyr::mutate(time_diff_imu_s = as.numeric(difftime(timestamp_closest_imu, gps_timestamp, units = "secs")))
  
  # safe filename: unique for every individual and without empty spaces
  # individual already loaded before
  safe_filename_imu <- gsub("[^a-zA-Z0-9]", "_", individual)
  
  # Save the data as a RDS
  rds_filename <- file.path("acc_gps_imu_cluster", paste0(safe_filename_imu, "gps_acc_imu_28july", ".rds"))
  saveRDS(joined_data, file = rds_filename)
  
  message("Data for ", individual, " saved to ", rds_filename)
  
})

#### 3. SUMMARIZE SEGMENTS VARIABLES WITH ALL SENSORS ####

#_______________________________________________________________________________
#### Summarize segments variables: at this point all data from GPS, ACC and IMU 

ind_list <- list.files("acc_gps_imu_cluster", pattern="gps_acc_imu_28july", full.names = T)

# read all inds together to see how many GPS locations associated to ACC and IMU in total for the 24 individuals
# Read and combine all tables
# all_inds_df <- do.call(rbind, lapply(ind_list, readRDS))
# Count total number of rows (GPS locations)
# cat("Total number of GPS locations associated with ACC and IMU:", nrow(all_inds_df), "\n")
# Total number of GPS locations associated with ACC and IMU: 1695418 

# flight hours tracked in total:
# Extract one representative duration per unique segment: 34561 segments before filter on duration and segment summaries
complete_dataset <- all_inds_df %>% drop_na()
unique_segment_durations <- complete_dataset %>%
  distinct(unique_segmID, duration)
# Sum the durations from the distinct segments
total_duration <- sum(unique_segment_durations$duration, na.rm = TRUE)/3600 # in seconds: 1694889, in hours: 470.80


gps_acc_imu_segm_complete <-  lapply(ind_list, function(file) { 
  
  # file <- ind_list[[2]]
  b <- readRDS(file)
  
  setDT(b)  # Convert 'b' to a data.table if problem with group_by and summarise
  summary_seg_dt <- b[, {
    startTime <- min(gps_timestamp, na.rm = TRUE)
    endTime <- max(gps_timestamp, na.rm = TRUE)
    .(
      # gps variables
      startTime = startTime,
      endTime = endTime,
      duration = as.numeric(endTime - startTime),
      grSpeed_mean = mean(grSpeed, na.rm = TRUE),
      h_min = min(height_agl, na.rm = TRUE),
      h_max = max(height_agl, na.rm = TRUE), 
      vel_mean = mean(vertSpeed_smooth, na.rm = TRUE),
      vel_min = min(vertSpeed_smooth, na.rm = TRUE),
      vel_max = max(vertSpeed_smooth, na.rm = TRUE),
      turnangle_mean = mean(abs(turnAngle), na.rm = TRUE),
      turnangle_sum = sum(abs(turnAngle), na.rm = TRUE),
      turnangle_var = var(turnAngle, na.rm = TRUE),
      turnangle_sd = sd(turnAngle, na.rm = TRUE),
      n_turnChange = length(which(turnChange != 0)),
      long_centroid = mean(location_long, na.rm = TRUE),
      lat_centroid = mean(location_lat, na.rm = TRUE),
      
      # acc variables
      VedBA_mean = mean(VedBA, na.rm = TRUE),
      VedBA_sd = sd(VedBA, na.rm = TRUE),
      VedBA_max = max(VedBA, na.rm = TRUE),
      VedBA_min = min(VedBA, na.rm = TRUE),
      ODBA_mean = mean(ODBA, na.rm = TRUE),
      sdACCz_mean = mean(stdvACC_z, na.rm = TRUE),
      sdACCz_max = max(stdvACC_z, na.rm = TRUE),
      sdACCz_min = min(stdvACC_z, na.rm = TRUE),
      diff_acc_time_mean = mean(diff_acc_time_s, na.rm = TRUE),
      acc_closest_timestamp_mean = mean(acc_closest_timestamp, na.rm = TRUE),
      
      # imu variables
      diff_imu_time_mean = mean(time_diff_imu_s, na.rm = TRUE),
      imu_closest_timestamp_mean = mean(timestamp_closest_imu, na.rm = TRUE),
      yaw_mean = mean(yaw_summary_yawmean, na.rm = TRUE),
      yaw_mean_abs = mean(yaw_summary_yawmean_abs, na.rm = TRUE),
      yaw_max = max(yaw_summary_yawmax, na.rm = TRUE),
      yaw_max_abs = max(yaw_summary_yawmax_abs, na.rm = TRUE),
      yaw_min = min(yaw_summary_yawmin, na.rm = TRUE),
      yaw_min_abs = min(yaw_summary_yawmin_abs, na.rm = TRUE),
      yaw_sum = mean(yaw_summary_yawsum, na.rm = TRUE),
      yaw_sum_abs = mean(yaw_summary_yawsum_abs, na.rm = TRUE),
      yaw_sd = mean(yaw_summary_yawsd, na.rm = TRUE),
      
      pitch_mean = mean(pitch_summary_pitchmean, na.rm = TRUE),
      pitch_mean_abs = mean(pitch_summary_pitchmean_abs, na.rm = TRUE),
      pitch_max = max(pitch_summary_pitchmax, na.rm = TRUE),
      pitch_max_abs = max(pitch_summary_pitchmax_abs, na.rm = TRUE),
      pitch_min = min(pitch_summary_pitchmin, na.rm = TRUE),
      pitch_min_abs = min(pitch_summary_pitchmin_abs, na.rm = TRUE),
      pitch_sum = mean(pitch_summary_pitchsum, na.rm = TRUE),
      pitch_sum_abs = mean(pitch_summary_pitchsum_abs, na.rm = TRUE),
      pitch_sd = mean(pitch_summary_pitchsd, na.rm = TRUE),
      
      roll_mean = mean(roll_summary_rollmean, na.rm = TRUE),
      roll_mean_abs = mean(roll_summary_rollmean_abs, na.rm = TRUE),
      roll_max = max(roll_summary_rollmax, na.rm = TRUE),
      roll_max_abs = max(roll_summary_rollmax_abs, na.rm = TRUE),
      roll_min = min(roll_summary_rollmin, na.rm = TRUE),
      roll_min_abs = min(roll_summary_rollmin_abs, na.rm = TRUE),
      roll_sum = mean(roll_summary_rollsum, na.rm = TRUE),
      roll_sum_abs = mean(roll_summary_rollsum_abs, na.rm = TRUE),
      roll_sd = mean(roll_summary_rollsd, na.rm = TRUE)
      
    )
  }, by = .(individual_local_identifier, burstID, unique_segmID, samplFreq_perAxis, deltah, h_iniz, h_fin, dist)] # area, vol (if calculated beforehand)
  
  # Convert the result back to a data frame
  summary_seg <- as.data.frame(summary_seg_dt)
  
  # Calculated rates only where it makes sense (in other cases would be speed)
  # summary_seg$area_r <- summary_seg$area/summary_seg$duration
  # summary_seg$vol_r <- summary_seg$vol/summary_seg$duration
  summary_seg$deltah_r <- summary_seg$deltah/summary_seg$duration # calculated on hiniz and hfin over ellipsoid -> another way to calculate vertical speed
  summary_seg$dist_r <- summary_seg$dist/summary_seg$duration # another way to calculate horizontal speed
  summary_seg$turnangle_r <- summary_seg$turnangle_sum/summary_seg$duration
  summary_seg$tilt_rad <- atan(summary_seg$deltah/summary_seg$dist) # in radians, corrected to be the right tilt angle! arctg(deltah/horiz)
  summary_seg$yaw_sum_r <- summary_seg$yaw_sum_abs/summary_seg$duration # sums of angles are returned for the abs sum over time
  summary_seg$pitch_sum_r <- summary_seg$pitch_sum_abs/summary_seg$duration
  summary_seg$roll_sum_r <- summary_seg$roll_sum_abs/summary_seg$duration
  
  # Fill values for some of the columns
  summary_seg$is_circling <- F
  summary_seg$n_circles <- 0
  summary_seg$n_circles_r <- 0
  summary_seg$turnChange_r <- 0
  # remember that change_head could be 0 also during circling, so always consider this variable associated with iscircling==T
  
  # Change the above values only for circling segments (when sum of segment turnAngle > 360 degrees (6.28 rad))
  circlingTF <- split(summary_seg, set_units(summary_seg$turnangle_sum) > set_units(6.28, "rad"))
  # before was circlingTF <- split(summary_seg, summary_seg$turnangle_sum > set_units(6, "rad"))
  names(circlingTF)
  # Check if "FALSE" exists in circlingTF and process accordingly
  if ("FALSE" %in% names(circlingTF)) {
    circlingTF[["FALSE"]]$n_turnChange <- 0
  }
  
  # Check if "TRUE" exists in circlingTF and process accordingly
  if ("TRUE" %in% names(circlingTF)) {
    circlingTF[["TRUE"]]$is_circling <- TRUE
    circlingTF[["TRUE"]]$n_circles <- drop_units(circlingTF[["TRUE"]]$turnangle_sum/set_units(6.28, "rad"))
    circlingTF[["TRUE"]]$n_circles_r <- circlingTF[["TRUE"]]$n_circles/circlingTF[["TRUE"]]$duration
    circlingTF[["TRUE"]]$turnChange_r <- circlingTF[["TRUE"]]$n_turnChange/circlingTF[["TRUE"]]$duration
  }
  
  # put the two back together and save
  summary_seg <- as.data.frame(rbindlist(circlingTF))
  
})

# Convert all 'units' columns to 'numeric' in each data frame in the list
for (i in seq_along(gps_acc_imu_segm_complete)) {
  gps_acc_segm_complete[[i]] <- lapply(gps_acc_imu_segm_complete[[i]], function(x) {
    if (inherits(x, "units")) {
      unclass(x)
    } else {
      x
    }
  })
}


# Combine all data frames in the list into a single data frame
gps_acc_imu_segm_complete <- rbindlist(gps_acc_imu_segm_complete, fill=TRUE)
gps_acc_imu_segm_complete <- gps_acc_imu_segm_complete[gps_acc_imu_segm_complete$duration>3,] # kept segments only when longer than 3 seconds


saveRDS(gps_acc_imu_segm_complete, "GPS_ACC_IMU_summaryvariables_28july.rds") # all summaries of all individuals together, 29775 segments with duration between 3 and 775 sec
# complete dataset with also IMU data on 1 august 2025

# if want to check:
# summaries <- readRDS("GPS_ACC_IMU_summaryvariables_28july.rds") 

# timelag between fix from ACC and IMu compared to GPS
# summary(summaries$diff_acc_time_mean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2065  0.2951  0.2999  0.2998  0.3047  0.3977 
# summary(summaries$diff_imu_time_mean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.0832  0.1832  0.2488  0.2486  0.3143  0.4979 
