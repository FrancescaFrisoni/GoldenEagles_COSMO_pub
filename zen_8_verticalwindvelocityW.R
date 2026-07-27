# This code estimates the vertical wind speed (W) from eagles' circular soaring
# and compares it with COSMO annotation
# Tom Carrard & Francesca Frisoni - updated 18 July 2026. Konstanz

setwd("/home/francesca/ownCloud/TesiFrancesca_UpliftClassification/Data")
directory <- "/home/francesca/ownCloud/TesiFrancesca_UpliftClassification"

library(data.table)
library(dplyr)
library(units)
library(sf)
library(lubridate)
library(ggplot2)
library(tidyr)

## import Martina Scacco's functions - used later to calculate airspeed from U, V and ground speed 
source("/home/francesca/R_projects/GoldenEagles_COSMO_pub/airspeed_windsupport_crosswind.R")


##### 1. LOAD DATASET AND EXTRACT CIRCULAR SEGMENTS #####

# ___________________________________________________________
#### Segments with associated GPS-ACC-IMU behavioural metrics

# load here summary of segments with gps, acc and imu metrics: all segments, but not all location points
segments <- readRDS("GPS_ACC_IMU_summaryvariables_28july.rds") # 29775 in 78 obs

segments$turning_angles_deg <- segments$turnangle_sum * (180 / pi) # conversion from rad to degrees, important after for wind components

circ <- segments[segments$is_circling==T & (segments$turning_angles_deg/segments$duration)>20,] # 5899 segments
# selection of is_circling=T bc could be that turning_angles_deg<360 but small duration so bias and could be >20 even if not circling

#____________________________
#### All GPS location points

# load all GPS points, not grouped into soaring segments
gps_complete <- as.data.frame(readRDS("./GPS_allsoaringpoints_areavol.rds"))

# retain all the gps points of the circular segments extracted above
gps_segm <- gps_complete[gps_complete$unique_segmID %in% circ$unique_segmID,]



##### 2. ESTIMATE WIND COMPONENTS ON THE CIRCULAR SEGMENTS #####

#___________________________________________________________
#### ERA5 data to estimate wind components on eagles points

# Load data from era5 point by point
era1 <- read.csv("./envdataERA5_2020_2023/annotation_csv/env_new_1_24july.csv-6091120100567852935.csv")
era2 <- read.csv("./envdataERA5_2020_2023/annotation_csv/env_new_2_24july.csv-8816069588445800851.csv")

era <- rbind(era1,era2)
segm_era <- merge(gps_segm, era, by="event_id")  # gps_segm are all the gps points from the GPS dataset that are in circ
# the points with U,V and W are 204057 because result of Tom's interpolation every 2nd or 3rd point

#______________________________________
#### Direction of the soaring segments 

gps_ls <- split(gps_segm, gps_segm$unique_segmID) #5899 circular segms

# calculate direction of the soaring segment from first and last location point of the segment
# direction is important when extracting wind components
direction_segm <- as.data.frame(rbindlist(
  lapply(gps_ls, function(b){
    # b <- gps_ls[[1]]
    sf_points <- st_as_sf(b, coords = c("location_long", "location_lat"), crs = 4326)
    first_last_point <- sf_points[c(1,nrow(b)),] # 1 and last row
    b$eucl_distance[1] <- max(st_distance(first_last_point))
    b$cross_country_speed[1] <- unique(b$eucl_distance)/unique(b$duration)
    # direction of a segment
    b$segm_direction[1] <- lwgeom::st_geod_azimuth(first_last_point)
    
    return(b) 
  }))
)

# direction eagles segment -> on complete GPS dataset bc needed first and last point!
# associate to my interpolated subset
segm_era_dir <- left_join(segm_era, unique(direction_segm[,c("unique_segmID","eucl_distance","cross_country_speed","segm_direction")]), by = "unique_segmID")

# Rotate directions, direction TO to FROM. Input direction is in degrees, either -180-180 or 0-360:
dirFROM <- function(dirTO){
  return((dirTO + 180) %% 360)
}

dirFROM_rad <- function(dirTO_rad) {
  return((dirTO_rad + pi) %% (2 * pi))
}

# segm_era_dir$segm_direction this is direction to the eagle, need to go to direction FROM
segm_era_dir$dirfromrad <- dirFROM_rad(segm_era_dir$segm_direction)


#___________________________________________________
#### Estimate wind support, cross wind and airspeed 

# calculate wind support -> direction of track already in radians
wind.support_rad <- function(u, v, dg_rad){ 
  if(any(dg_rad[!is.na(dg_rad)] < 0 | dg_rad[!is.na(dg_rad)] > 2*pi)){
    stop("The track direction (dg_rad) has invalid values, it must range from 0 to 2*pi radians!")
  }
  wd <- atan2(u, v)                          # wd (wind direction) from -pi to +pi
  wd_2pi <- ifelse(wd < 0, 2*pi + wd, wd)    # wd from 0 to 2*pi
  beta <- wd_2pi - dg_rad                    # beta is the difference between wd and dg (both from 0 to 2*pi)
  return(cos(beta) * sqrt(u * u + v * v))    # ws = cos(beta)*Vw
}

ws <- wind.support_rad(segm_era$ECMWF.ERA5.SL.Wind..10.m.above.Ground.U.Component., segm_era$ECMWF.ERA5.SL.Wind..10.m.above.Ground.V.Component., segm_era_dir$dirfromrad)

#  calculate cross-wind
cross.wind_rad <- function(u, v, dg_rad){
  if(any(dg_rad[!is.na(dg_rad)] < 0 | dg_rad[!is.na(dg_rad)] > 2*pi)){
    stop("The track direction (dg_rad) has invalid values, it must range from 0 to 2*pi radians!")
  }
  wd <- atan2(u, v)                             # wd (wind direction) from -pi to +pi
  wd_2pi <- ifelse(wd < 0, 2*pi + wd, wd)       # wd from 0 to 2*pi
  beta <- wd_2pi - dg_rad                       # beta is the difference between wd and dg (both from 0 to 2*pi)
  return(abs(sin(beta) * sqrt(u * u + v * v)))  # wc = |sin(beta)*Vw|
}

cw <- cross.wind_rad(segm_era$ECMWF.ERA5.SL.Wind..10.m.above.Ground.U.Component., segm_era$ECMWF.ERA5.SL.Wind..10.m.above.Ground.V.Component., segm_era_dir$dirfromrad)

# calculate airspeed - will enter the equation of sink speed
airspeed <- function(Vg, Ws, Cw) {
  return(sqrt((Vg - Ws)^2 + (Cw)^2))
}

segm_era$airspeed <- airspeed(as.numeric(segm_era$grSpeed),ws,cw)


##### 3. COMPUTE VERTICAL SPEED FROM EAGLE ASCENT #####
# Derive vertical wind speed based on eagle vertical speed and aerodynamic equations

#_________________________________________
# Sequences of wing span, area, body masses

# sequences: from literature 'Birds of the world'
# wing span 160-222 cm
# body mass 2.387-6.665 kg
# wing area 5237-6070 cm2

mass <- seq(2.4, 6.7, by = 0.3)
wing_sp <- seq(1.6, 2.22, by = 0.1)
wing_area <- seq(0.52, 0.61, by = 0.025) # converted in m^2

# Define variables and parameters
rho <- 0.957  # density (assume constant at 2500m asl for now)
CDo <- 0.1    # drag coefficient
k <- 1.1      # drag coefficient
phi <- 25     # bank angle (degrees)
g <- 9.81     # gravitational acceleration

# Create plausible combinations of wing size and weight
small_bird_meas <- expand.grid(mass = mass[1:7], span = wing_sp[1:3], area = wing_area[1:2])
large_bird_meas <- expand.grid(mass = mass[8:length(mass)], span = wing_sp[4:length(wing_sp)], area = wing_area[3:length(wing_area)])
bird_meas <- rbind(small_bird_meas, large_bird_meas)


#____________________________
# Define and apply sink_speed 
sink_speed <- function(V, span, area, mass, rho=0.957, CDo=0.1, k=1.1, phi=25, g=9.81) {
  Vz <- 2 * k * mass * g / (pi * rho * span^2 * V) + CDo * rho * area * V^3 / (2 * mass * g)
  return(Vz)
}

# trajectory data
traj_all <- segm_era

# Create date object
traj_all$date <- as.Date(substr(traj_all$timestamp.x, 1, 10))

# Apply sinking_speed calculation to all eagle data points
Vz_arr <- matrix(0, nrow = nrow(traj_all), ncol = nrow(bird_meas))

for (i in 1:nrow(bird_meas)) {
  Vz_arr[, i] <- sink_speed(V = traj_all$airspeed, span = bird_meas$span[i], 
                            area = bird_meas$area[i], mass = bird_meas$mass[i])
}

# Compute circling sink rate
Vzc <- Vz_arr / sqrt(cos(phi * pi / 180)^3)

# Compute density-corrected sink rate
Vs <- Vzc * sqrt(1.225 / rho)

# ___________________________________________________________
# Compute the estimate of the vertical wind at eagle location
W_eagle <- matrix(0, nrow = nrow(Vs), ncol = ncol(Vs))

for (i in 1:ncol(Vs)) {
  W_eagle[, i] <- traj_all$vertSpeed + Vs[, i]
}

# Convert to dataframe
col_names <- paste0("W_eagle_", 1:ncol(W_eagle))

# Create dataframe with wind values
W_df <- as.data.frame(W_eagle)
names(W_df) <- col_names

W_df$mean <- rowMeans(W_eagle, na.rm=T)
# W_df$median <- apply(W_eagle, 1, median, na.rm=T)
# W_df$perc5 <- apply(W_eagle, 1, quantile, probs = 0.05, na.rm=T)
# W_df$perc95 <- apply(W_eagle, 1, quantile, probs = 0.95,na.rm=T)

traj_all$W_eagle <- W_df$mean

print(length(traj_all$W_eagle)) #204057
print(length(unique(traj_all$date))) #334


##### 4. COMPARISON OF W BETWEEN EAGLES AND COSMO #####

# _______________________________________________________
# Check and remove outliers in W calculated from eagles

summary(traj_all$W_eagle)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -15.051   2.601   3.456   3.666   4.466 503.572 

quantile(traj_all$W_eagle, seq(0, 1, 0.005), na.rm = TRUE)
#     0.0%        0.5%        1.0%        1.5%        2.0%        2.5%        3.0%       
# -15.0508523  -0.4259085   0.2070554   0.5142497   0.7167664   0.8813429   1.0088904 
# 98.0%       98.5%       99.0%       99.5%  100.0% 
#  7.8353086   8.3568444   9.1497179  10.8367074 503.5720205 

# the very extreme values sit within 0.05%
quantile(traj_all$W_eagle, seq(0, 1, 0.0005), na.rm = TRUE)
 #     0.00%        0.05%         0.10%         
 #  -15.050852277  -3.614929282  -2.386068165
 #     99.95%   100.00% -- apply tail(quantile) to check these
 # 23.61658 503.57202 

# remove lower and upper 0.05% of data as outliers
lower <- quantile(traj_all$W_eagle, 0.0005, na.rm = TRUE) #-3.61
upper <- quantile(traj_all$W_eagle, 0.9995, na.rm = TRUE) #23.6

traj_sub <- traj_all[traj_all$W_eagle > lower & traj_all$W_eagle < upper, ] # 203851

# _______________________________
# Load and merge ERA5 wind speed
windera <- readRDS("/home/francesca/ownCloud/TesiFrancesca_UpliftClassification/Data/ERA5_wind_per_id.rds") 

traj_sub_era <- left_join(traj_sub, windera, by = "event_id") 
# ERA5_wind_speed -> is the synoptic vertical wind speed W of ERA5, not the horizontal from U and V
# this value is not actively used here, but could be useful for synoptic comparisons at some point

# _______________________________
# Load and compare W from COSMO

w_cosmo <- read.csv("/home/francesca/ownCloud/ETH_sharedFolder/GPS_allsoaringpoints_areavol_subset_annotated_fullwind.csv") # corrected july 26

# comparison with W from COSMO -> merge the W column to the traj dataset by eventid
traj_w <- merge(traj_sub_era, w_cosmo[, c("event_id", "W")], by = "event_id", all.x = TRUE) 
saveRDS(traj_w, "./traj_weagle_cosmo_era5_july26.rds" )

summary(traj_w$W) #from COSMO: careful, the label of the column is exactly 'W'
#          Min. 1st Qu.  Median    Mean       3rd Qu.    Max.    NA's 
#   -6.638535 -0.138051  0.007987 -0.011732  0.165795  4.442057        72 

summary(traj_w$W_eagle) # label of W from eagles is W_eagle
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -3.590   2.602   3.456   3.646   4.465  23.504

# ___________________________________________
# Density plots of W between eagles and COSMO

w <- ggplot(traj_w) +
  geom_density(aes(x = W, fill = "W COSMO"), alpha = 0.5) +
  geom_density(aes(x = W_eagle, fill = "W Eagle"), alpha = 0.5) +
  scale_fill_manual(values = c("W COSMO" = "coral", "W Eagle" = "darkgreen")) +
  labs(
    x = "Vertical Wind Speed [m s-1]",
    y = "Density",
    fill = "Variable"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 30), 
    axis.title = element_text(size = 30),  
    axis.text = element_text(size = 30), 
    legend.title = element_text(size = 30),  
    legend.text = element_text(size = 30),  
    legend.position = c(0.95, 0.95),      # x, y in [0,1], top-right corner. can be changed to legend.position = "bottom",
    legend.justification = c(1, 1)         # anchors the legend box to that corner
    # plot.title = element_text(hjust = 0.5, size = 16)  
  )

print(w)

ggsave(file.path(directory, "figures_july26", "w_cosmo_eagle_solo.pdf"), 
       plot = w, width = 297, height = 210, units = "mm", device = "pdf")

# save ggplot obj to reload if some fig tuning is necessary
saveRDS(w, file.path(directory, "figures_july26", "w_cosmo_eagle_solo.rds"))


##### 5. COMPARISON OF W BETWEEN EAGLES AND COSMO ACROSS UPLIFT CLASSES #####

#________________________
## Load labelled dataset
segm_pred_sure <- readRDS("./uplift_classification/segm_pred_sure_19035.rds")

segm_df <- merge(traj_w, segm_pred_sure[,c("unique_segmID","uplift_type")], by="unique_segmID") # 158579

length(unique(segm_df$unique_segmID)) # 4448 segments with labelled uplift types

# --- correlation plots ---
# Calculate correlations
# cor_orog <- cor(segm_df$W[segm_df$uplift_type == "orog"], 
#                 segm_df$W_eagle[segm_df$uplift_type == "orog"], 
#                 use = "complete.obs")
# cor_thermal <- cor(segm_df$W[segm_df$uplift_type == "thermal"], 
#                    segm_df$W_eagle[segm_df$uplift_type == "thermal"], 
#                    use = "complete.obs")
# cor_wave <- cor(segm_df$W[segm_df$uplift_type == "wave"], 
#                 segm_df$W_eagle[segm_df$uplift_type == "wave"], 
#                 use = "complete.obs")
# cor_coefficients <- c(orog = cor_orog, thermal = cor_thermal, wave = cor_wave)
# 
# color_palette <- c(
#   "thermal" = "#D55E00",
#   "orog" = "#CC79A7", 
#   "wave" = "#0072B2"
# )
# 
# # named vector for facet labels
# facet_labels <- c("orog" = "Orographic", "thermal" = "Thermal", "wave" = "Wave")
# 
# # fixed limits for both axes
# axis_limits <- c(-10, 10)
# ggplot(segm_df, aes(x = W_eagle, y = W, color = uplift_type)) +
#   geom_point(alpha = 0.3) +  
#   geom_smooth(method = "lm", se = FALSE, color = "black", aes(group = 1)) +  # correlation lines
#   facet_wrap(~ uplift_type, ncol = 3, labeller = labeller(uplift_type = facet_labels)) +
#   coord_fixed(ratio = 1, xlim = axis_limits, ylim = axis_limits) +  # limits to -10 and +10
#   scale_color_manual(values = color_palette) +
#   labs(x = "W Eagle",
#        y = "W COSMO") +
#   theme_minimal() +
#   theme(
#     text = element_text(size = 16),  
#     axis.title = element_text(size = 16),  
#     axis.text = element_text(size = 14),  
#     strip.text = element_text(face = "bold", size = 16),  
#     panel.spacing = unit(1, "cm"),
#     legend.position = "none"
#   ) +  
#   # correlation coefficients as text
#   geom_text(data = data.frame(uplift_type = names(facet_labels),
#                               label = paste("r =", round(cor_coefficients, 2)),
#                               x = -Inf, y = Inf),
#             aes(x = x, y = y, label = label),
#             hjust = -0.1, vjust = 1.5, color = "black", size = 7)  

# ______________________________________________________________
# Violin plot for density distribution across the 3 uplift types

traj_long <- tidyr::pivot_longer(segm_df, cols = c(W, W_eagle), names_to = "variable", values_to = "value")

e <- ggplot(traj_long, aes(x = variable, y = value, fill = variable)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  facet_wrap(~ uplift_type, ncol = 3, labeller = labeller(uplift_type = facet_labels)) +
  scale_fill_manual(values = c("W" = "coral", "W_eagle" = "darkgreen"),
                    labels = c("W" = "W COSMO", "W_eagle" = "W Eagle")) +
  scale_x_discrete(labels = c("W" = "W COSMO", "W_eagle" = "W Eagle")) +
  labs(x = NULL,
       y = "Vertical Wind Speed [m s-1]",
       fill = "Variable") +
  theme_minimal() +
  theme(
    text = element_text(size = 30),  
    axis.title = element_text(size = 30),  
    axis.text = element_text(size = 25),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25),  
    strip.text = element_text(face = "bold", size = 30), 
    # legend.title = element_text(size = 30),  
    # legend.text = element_text(size = 30),  
    panel.spacing = unit(1, "cm"),
    legend.position = "none",
    # legend.position = "bottom",
    # plot.title = element_text(hjust = 0.5, size = 18)  
  )
print(e)

ggsave(file.path(directory, "figures_july26", "w_cosmo_eagle_3up.pdf"), 
       plot = e, width = 297, height = 210, units = "mm", device = "pdf")

# save ggplot obj to reload if some fig tuning is necessary
saveRDS(e, file.path(directory, "figures_july26", "w_cosmo_eagle_3up.rds"))

