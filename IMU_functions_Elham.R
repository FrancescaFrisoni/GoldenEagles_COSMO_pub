
#script for calculating Euler angles from quaternions provided by eobs.
#Apr. 08.2024. Elham Nourani, PhD.


#library(tidyverse)

# STEP 1: write some functions #####

#define a function to process quaternions in eobs format (eobs specific) 
# process_quaternions <- function(quaternion_string, ftn) {
#   #define a function to split up a vector of multiple quaternions into a list with each element as a vector of 4 floats (eobs specific) 
#   vector_to_quat_ls <- function(x) {
#     n <- length(x)
#     split(x, rep(1:(n/4), each = 4))
#   }
#   
#   quaternions <- str_split(quaternion_string, " ")[[1]] %>%
#     as.numeric()
#   
#   result <- vector_to_quat_ls(quaternions) %>%
#     map(ftn) %>%
#     unlist() %>%
#     as.character()
#   
#   return(str_c(result, collapse = " "))
# }


# ff on 28 july 2025 with printing for debug and keep track of the code running
process_quaternions <- function(quaternion_string, ftn) {
  # small inside function
  vector_to_quat_ls <- function(x) {
    n <- length(x)
    split(x, rep(1:(n/4), each = 4))
  }
  
  quaternions <- str_split(quaternion_string, " ")[[1]] %>% as.numeric()
  
  quat_list <- vector_to_quat_ls(quaternions)
  cat("Number of quaternions:", length(quat_list), "\n")
  
  results <- vector("list", length(quat_list))
  for(i in seq_along(quat_list)) {
    quat <- quat_list[[i]]
    cat("Quaternion", i, ":", paste(quat, collapse = " "), "\n")
    angle <- ftn(quat)
    cat("Angle from ftn:", angle, "\n\n")
    results[[i]] <- angle
  }
  
  result_unlisted <- unlist(results)
  return(str_c(result_unlisted, collapse = " "))
}

#function to calculate summary statistics for a numeric vector 
string_to_numeric <- function(x) {str_split(x, " ")[[1]] %>% 
    as.numeric()
}

# #function to convert the character string to numeric vector
strings_to_numeric <- function(angle_strings, ftn) { #input is multiple rows of data
  
  #convert the multiple strings to one (for summarizing)
  numeric_vec <- unlist(map(angle_strings, string_to_numeric))
  
  return(numeric_vec)
  #numeric_vec <- str_split(angle_string, " ")[[1]] %>%
  #  as.numeric()
  
  #result <- numeric_vec %>%
  #  map(ftn) %>%
  #  unlist() %>%
  #  as.character()
  
  #return(str_c(result, collapse = " "))
}

#function to calculate the statistical mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}



#function to calculate summary statistics for a numeric vector 
angle_summaries <- function(x) {
  data.frame(mean = mean(x, na.rm = T),
             max = max(x, na.rm = T),
             min = min(x, na.rm = T),
             sum = sum(x, na.rm = T),
             sd = sd(x, na.rm = T),
             mean_abs = mean(abs(x), na.rm = T),
             max_abs = max(abs(x), na.rm = T),
             min_abs = min(abs(x), na.rm = T),
             sum_abs = sum(abs(x), na.rm = T))
}


#the following functions are from Kami's IMU_conversion.r

#define a function for converting raw quaternion values mathematically from integers to floats (based on eobs manual) 
.convertEobs <- function(x){
  if(x[1]==-32768){
    x[2:4] <- 0 # corresponds to scalar part of Quaternion 1.0 or -1.0
  }
  # STEP 1: Convert raw quaternion values from integers to floats based on eobs manual
  # Calculate r
  r <- sqrt(x[2]^2 + x[3]^2 + x[4]^2)
  # Calculate the scalar value (w)
  qw <- x[1] / 32768  # The denominator 32768 is used for normalization from a 16 bit signed integer memory size
  # Calculate s
  if (r != 0) {
    s <- sqrt(1 - qw^2) / r
  } else {
    s <- 0
  }
  return(as.numeric(c(qw, s * x[2], s * x[3], s * x[4])))
}

#define a function to Calculate pitch angle from quaternion
get.pitch <- function(x, type=c("eobs", "quaternion")) {
  type <- match.arg(type)  # ensures type is one of the choices
  if(length(x) != 4){
    stop("Improper quaternion passed to function")
  }
  if(any(is.na(x))){
    pitchAngle <- NA
  } else {
    if(type == "eobs"){
      quat <- .convertEobs(x)
    } else {
      quat <- x
    }
    pitchAngle <- asin(2 * (quat[1]*quat[2] + quat[3]*quat[4]))
  }
  return(pitchAngle)
}

#define a function to calculate roll angle from quaternion
get.roll <- function(x, type=c("eobs", "quaternion")) {
  type <- match.arg(type)  # ensures type is one of the choices
  if(length(x) != 4){
    stop("Improper quaternion passed to function")
  }
  if(any(is.na(x))){
    rollAngle <- NA
  } else {
    if(type == "eobs"){
      quat <- .convertEobs(x)
    } else {
      quat <- x
    }
    rollAngle <- -atan2(2 * (quat[1] * quat[3] - quat[2] * quat[4]), 1.0 - 2.0 * (quat[2]^2 + quat[3]^2))
  }
  return(rollAngle)
}


#define a function to Calculate yaw angle from quaternion
get.yaw <- function(x, type=c("eobs", "quaternion")) {
  type <- match.arg(type)  # ensures type is one of the choices
  if(length(x) != 4){
    stop("Improper quaternion passed to function")
  }
  if(any(is.na(x))){
    yawAngle <- NA
  } else {
    if(type == "eobs"){
      quat <- .convertEobs(x)
    } else {
      quat <- x
    }
    yawAngle <- -1*atan2(2.0*(quat[2]*quat[3] - quat[1]*quat[4]) , 1-2*(quat[2]^2 + quat[4]^2) )
  }
  return(yawAngle)
}

# # STEP 2: apply to eobs data #####
# #open sample data for one individual
# sample_data <- read.csv("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/matched_GPS_IMU/matched_gps_orientation/D163_696_quat_mag_w_gps.csv")
# 
# #calculate pitch, roll, and yaw
# sample_data_quat <- sample_data %>%
#   mutate(
#     pitch = process_quaternions(orientation_quaternions_raw, ~ get.pitch(.x, type = "eobs")),
#     yaw = process_quaternions(orientation_quaternions_raw, ~ get.yaw(.x, type = "eobs")),
#     roll = process_quaternions(orientation_quaternions_raw, ~ get.roll(.x, type = "eobs"))
#   )

# Link to Elham's github
# https://github.com/mahle68/HB_ontogeny/blob/main/MS1_IMU_classification/01b_imu_processing.r
# Define a function to find the closest GPS information and associate it with orientation data

# adapted to my segmentations 28 july 25
# function to connect gps to imu dataset
find_closest_gps <- function(or_data, gps_data, time_tolerance = 10 * 60) { # 10 times 60 sec, so 10 min
  map_df(1:nrow(or_data), function(h) {
    or_row_time <- or_data[h, "timestamp"]
    gps_sub <- gps_data %>%
      filter(between(timestamp, or_row_time - time_tolerance, or_row_time + time_tolerance))
    
    if (nrow(gps_sub) >= 1) {
      time_diff <- abs(difftime(gps_sub$timestamp, or_row_time, units = "secs"))
      min_diff <- which.min(time_diff)
      or_data[h, c("timestamp_closest_gps", "location_long_closest_gps", "location_lat_closest_gps", "height_above_ellipsoid_closest_gps")] <- 
        gps_sub[min_diff, c("timestamp", "location_long", "location_lat", "height_above_ellipsoid"
                            )]
    } else {
      or_data[h, c("timestamp_closest_gps", "location_long_closest_gps", "location_lat_closest_gps", "height_above_ellipsoid_closest_gps"
                   )] <- NA
    }
    return(or_data[h, ])
  })
}

# ff 28 july 2025, in my case I need to associate the imu to my gps dataset, time tolerance of 30 sec

find_closest_imu_for_gps <- function(gps_data, imu_data, time_tolerance = 30) {
  closest_imu_df <- map_df(1:nrow(gps_data), function(i) {
    gps_row_time <- gps_data$gps_timestamp[i]
    
    imu_sub <- imu_data %>%
      filter(between(timestamp,
                     gps_row_time - seconds(time_tolerance),
                     gps_row_time + seconds(time_tolerance)))
    
    if (nrow(imu_sub) == 0) {
      return(tibble(
        gps_row_index = i,
        timestamp_closest_imu = as.POSIXct(NA),
        roll_closest_imu = NA_real_,
        pitch_closest_imu = NA_real_,
        yaw_closest_imu = NA_real_
      ))
    }
    
    time_diff <- abs(difftime(imu_sub$timestamp, gps_row_time, units = "secs"))
    min_diff_index <- which.min(time_diff)
    
    timestamp_val <- imu_sub$timestamp[min_diff_index]
    if (!inherits(timestamp_val, "POSIXct")) {
      timestamp_val <- as.POSIXct(timestamp_val, origin = "1970-01-01", tz = "UTC")
    }
    # the above chunk is to transform the time difference into a real date and not seconds from the origin in 1970
    
    tibble(
      gps_row_index = i,
      timestamp_closest_imu = timestamp_val,
      roll_closest_imu = imu_sub$roll_deg[min_diff_index],
      pitch_closest_imu = imu_sub$pitch_deg[min_diff_index],
      yaw_closest_imu = imu_sub$yaw_deg[min_diff_index]
    )
  })
  
  # Join IMU data back to original gps_data by row index
  gps_data %>%
    dplyr::mutate(gps_row_index = row_number()) %>%
    left_join(closest_imu_df, by = "gps_row_index") %>%
    select(-gps_row_index)  # Remove helper column if desired
}
