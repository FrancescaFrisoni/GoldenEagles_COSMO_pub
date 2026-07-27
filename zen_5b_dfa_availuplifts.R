# Golden eagles use of uplift sources versus available
# this script classifies the available uplifts using the DFA algorithm defined in script zen_4
# and compares it with the uplifts used in our tracking dataset

# Francesca Frisoni - updated 21st July 2026, Konstanz

setwd("/home/francesca/ownCloud/TesiFrancesca_UpliftClassification")
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
library(ggplot2)
library(gridExtra)
library(tidyr)
library(scales)


# function to calculate wind speed, used below
wind.speed <- function(u, v) {
  return(sqrt(u^2 + v^2))
}

# load dataset used to run the DFA (lda) model to exactrct predictors
sub_lab_7s <- readRDS("./Data/labelled209segmentsabove7s_predictors_150726.rds")
# from sub_lab_7s: predictors are all except unique_segmID, date, uplift_type, duration (which is only a filter for now). removed U and V from predictors bc same info than horizontal windspeed
# names(sub_lab_7s)
predictors <- names(sub_lab_7s)[c(5:8, 11:16)]
# [1] "w_oro_mean"        "N2_max_h_mean"     "N2_max_value_mean" "ASHFL_S_mean"      "W_mean"           
# [6] "windspeed_mean"    "max_height.ab.gr"  "mean_slope"        "mean_roughness"    "mean_aspect" 

# Load DEM
dem <- terra::rast("./Data/dtm_elev_eumap_epsg3035_v0.3_cropAlps_laea.tif")
# this dem is in decimeters (dm), not yet in meters
dem <- dem / 10 # convert to m (check summary of the values to make sure it worked)
dem[dem < -1000] <- NA # put values below -1000 as NA to avoid problems with roughness<0
# Calculate, from dem in m, slope, roughness and aspect in degrees
slope <- terra::terrain(dem, v = "slope", neighbors = 8, unit = "degrees", filename = "")
aspect <- terra::terrain(dem, v = "aspect", neighbors = 8, unit = "degrees", filename = "")
roughness <- terra::terrain(dem, v = "roughness", neighbors = 8, unit = "degrees", filename = "")

# _________________________________________________________________________
# Annotate Individual MCPs with COSMO and topography variables (predictors) 
# for ind MCP: I took data for each individual and took the convex hull of all soaring locations of each individual. I then checked on average how many soaring segments each individual had per day, between 33 and 38. So I took 5 times the amount of soaring segments (200 points) and spread them randomly within each individual’s polygon. I then associated to each of these random locations the days in which the individual flew and for the time I took the median hour of the day across the entire dataset, that was 12:00

# load annotated dataset: all background points annotated with COSMO variables
up_ind <- read.csv("/home/francesca/ownCloud/ETH_sharedFolder/UseVsAvailable/ForFrancesca_backgroundPoints_AtmoAvailability_perInd_annotated_corrected_interp.csv")

nrow(up_ind) #144800
sum(is.na(up_ind$V)) #0, no NAs present
head(up_ind)

names(up_ind)
# [1] "X"                           "index"                       "x"                          
# [4] "y"                           "individual_local_identifier" "MCParea_km2"                
# [7] "timestamp"                   "height.above.sea.level"      "medHeight_agl"              
# [10] "medHeight_ell"               "height.above.surf"           "N2_max_h"                   
# [13] "N2_max_value"                "w_oro"                       "ASHFL_S"                    
# [16] "U"                           "V"                           "W"                          
# [19] "interp" 
# The new column 'interp' tells whether the values of U, V and W are interpolated, in which case the value is True, or whether it is filled with the lowest model level, in which case the value is False. 

# table(up_ind$interp)
# True 
# 144800 

# Wind speed
up_ind$windspeed <- wind.speed(up_ind$U, up_ind$V)

# Topography on points (dem and layers loaded and manipulated above)
coords <- up_ind[, c("x", "y")] # x is long, y is lat
coords <- st_as_sf(coords, coords = c("x", "y"), crs = "+proj=longlat +ellps=WGS84")
coords_proj <- st_transform(coords, crs(dem))

up_ind$slope <- terra::extract(slope, coords_proj, method = "bilinear", cells = F, ID = F)[, 1]
up_ind$aspect <- terra::extract(aspect, coords_proj, method = "bilinear", cells = F, ID = F)[, 1]
up_ind$roughness <- terra::extract(roughness, coords_proj, method = "bilinear", cells = F, ID = F)[, 1]
terrainElevation <- terra::extract(dem, coords_proj, method = "bilinear", cells = F, ID = F)[, 1]
up_ind$height_above_ground <- up_ind$height.above.sea.level - terrainElevation

up_ind$aspect_rad <- circular(up_ind$aspect * pi / 180, units = "radians")

up_ind$timestamp <- as.POSIXct(up_ind$timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
up_ind$date <- as.Date(up_ind$timestamp)

# could possibly summarize variables into polygons: here defined one polygon for each individual each date
up_ind <- up_ind %>% mutate(polyID = paste(individual_local_identifier, date, sep = "_")) 

# keep all points and run DFA singularly on them, no group by
# each pixel of 1 km represents the availability of uplift
# I just change col names to be consistent for the DFA model, it is not a real summary

up_ind_summary_df <- up_ind %>%
  mutate(
    w_oro_mean           = w_oro,
    N2_max_h_mean        = N2_max_h,
    N2_max_value_mean    = N2_max_value,
    ASHFL_S_mean         = ASHFL_S,
    U_mean               = U,
    V_mean               = V,
    W_mean               = abs(W),                 # abs as for other models
    windspeed_mean       = windspeed,
    max_height.ab.gr     = height_above_ground,    # or pmax(...) if you ever summarise
    mean_slope           = slope,
    mean_roughness       = roughness,
    mean_aspect          = as.numeric(mean.circular(aspect_rad)) * 180 / pi
  )

#___________________________________
# RUN DFA MODEL ON BACKGROUND POINTS 

set.seed(12)

# Fit the LDA model
lda_model <- readRDS ("./Data/uplift_classification/lda_classificationmodel.rds") ## SAVED IN SCRIPT ZEN_4

predicted <- predict(lda_model, newdata = up_ind_summary_df)
# head(predicted)

prob <- as.data.frame(round(predicted$posterior, 2))
prob$pred <- predicted$class

# table(prob$pred) 
# orog thermal    wave 
# 30215   93107   21478 

prob$index <- up_ind_summary_df$index

up_ind_pred <- left_join(up_ind_summary_df, prob, by = "index")
# here up_ind_pred can be saved as RDS file

# _______________________________________________
# Define uplift classes based on DFA predictions 

# extract month and year
up_ind_pred$month <- month(up_ind_pred$date, label = TRUE, abbr = FALSE, locale = "en_US.UTF-8")
up_ind_pred$year <- year(up_ind_pred$date)

# order by month as factor
up_ind_pred$month <- factor(up_ind_pred$month, levels = month.name, ordered = TRUE)

# Define uplift types with my thresholds
up_ind_pred$uplift_type <- NA
up_ind_pred <- up_ind_pred %>%
  mutate(uplift_type = case_when(
    pred == "wave" & wave >= 0.8 ~ "wave",
    pred == "orog" & orog >= 0.8 ~ "orog",
    pred == "thermal" & thermal >= 0.8 ~ "thermal",
    thermal > 0.4 & orog > 0.4 ~ "thermal/orog",
    thermal > 0.4 & wave > 0.4 ~ "thermal/wave",
    orog > 0.4 & wave > 0.4 ~ "orog/wave",
    TRUE ~ "unknown"
  ))

table(up_ind_pred$uplift_type)
# orog    orog/wave      thermal thermal/orog thermal/wave      unknown         wave 
# 13783         4767        78238         1528         1983        37561         6940 

up_ind_pred_sure <- up_ind_pred[up_ind_pred$uplift_type %in% c("orog", "wave", "thermal"), ] # 98961 
table(up_ind_pred_sure$uplift_type)

# orog thermal    wave 
# 13783   78238    6940 

saveRDS(up_ind_pred_sure, "./Data/uplift_classification/up_ind_pred_sure_polyind_98961.rds")


# _______________________________________________
## Plot availability of uplifts over months 2023

color_palette <- c(
  "thermal" = "#D55E00",
  "orog" = "#CC79A7", #"#6A0DAD"
  "wave" = "#0072B2"
)

# vector of all months as ordered factor
all_months <- factor(month.name, levels = month.name, ordered = TRUE)


## 2 categories (thermal vs dynamic)

prepare_data_dynamic <- function(df, dataset_name) {
  df$month <- factor(
    month(df$date, label = TRUE, abbr = FALSE, locale = "en_US.UTF-8"),
    levels = month.name, 
    ordered = TRUE
  )
  df$year <- year(df$date)
  
  df_summary <- df %>%
    filter(year == 2023) %>%
    group_by(month) %>%
    summarize(
      orog_wave = sum(uplift_type %in% c("orog", "wave"), na.rm = TRUE),
      thermal   = sum(uplift_type == "thermal",           na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    # Add all missing months with zeros
    complete(
      month = all_months,
      fill = list(orog_wave = 0L, thermal = 0L)
    ) %>%
    mutate(total = thermal + orog_wave) %>%
    pivot_longer(
      cols = c(thermal, orog_wave),
      names_to = "type",
      values_to = "count"
    ) %>%
    mutate(proportion = ifelse(total > 0, count / total, 0),
           dataset = dataset_name)
  
  return(df_summary)
}

up_ind_pred_data_dyn <- prepare_data_dynamic(up_ind_pred_sure, "Uplift prediction")

# add column with labels (only for dynamic to avoid duplicates)
up_ind_pred_data_dyn <- up_ind_pred_data_dyn %>%
  group_by(dataset, month) %>%
  mutate(label = if_else(type == "orog_wave", as.character(first(total)), NA_character_)) %>%
  ungroup()

# re_order levels to have thermals on top
up_ind_pred_data_dyn$type <- factor(up_ind_pred_data_dyn$type, levels = c("thermal", "orog_wave"))

up2 <- ggplot(up_ind_pred_data_dyn, aes(x = month, y = proportion, fill = type)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = label, y = 1.02), size = 4.5, na.rm = TRUE) + 
  scale_fill_manual(values = c("thermal" = "#F4A261", "orog_wave" ="#6495ED"),
                    labels = c("thermal" = "Thermal", "orog_wave" = "Dynamic")) +
  labs(
    x = "Month",
    y = "Proportion of Available Uplifts",
    fill = "Uplift Type") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    legend.text = element_text(size = 30),    
    legend.title = element_text(size = 30, face = "bold")  
  ) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1.1), breaks = seq(0, 1, by = 0.2)) +
  # Ensure all months displayed even if missing
  scale_x_discrete(drop = FALSE)

up2

ggsave(file.path(directory, "figures_july26", "dyn_th_upliftavail2023.pdf"), 
       plot = up2, width = 297, height = 210, units = "mm", device = "pdf")

# save ggplot obj to reload if some fig tuning is necessary
saveRDS(up2, file.path(directory, "figures_july26", "dyn_th_upliftavail2023.rds"))

# to report the max and min proportion of dynamic uplift sources
dyn_proportions <- up_ind_pred_data_dyn$proportion[up_ind_pred_data_dyn$type=="orog_wave"]
range(dyn_proportions) # 0.07172405 0.80144928



# if you are interested in availability of each type of uplift for 3 classes:
prepare_data <- function(df, dataset_name) {
  df$month <- factor(lubridate::month(df$date), levels = 1:12, labels = month.name)
  df$year <- lubridate::year(df$date)
  
  df_summary <- df %>%
    filter(year == 2023) %>%
    group_by(month) %>%
    summarize(
      thermal = sum(uplift_type == "thermal", na.rm = TRUE),
      orog = sum(uplift_type == "orog", na.rm = TRUE),
      wave = sum(uplift_type == "wave", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(total = thermal + orog + wave) %>%
    pivot_longer(cols = c(thermal, orog, wave),
                 names_to = "type",
                 values_to = "count") %>%
    mutate(proportion = count / total,
           dataset = dataset_name)
  
  return(df_summary)
}

up_ind_pred_sure_data <- prepare_data(up_ind_pred_sure, "up_ind_pred_sure")

# Plot of 2023 in 3 categories
ggplot(up_ind_pred_sure_data, aes(x = month, y = proportion, fill = type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_palette) +
  labs(title = "Proportion of Uplift Types in 2023",
       x = "Month",
       y = "Proportion",
       fill = "Uplift Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::percent_format()) +
  facet_wrap(~dataset, ncol = 2) +
  coord_cartesian(ylim = c(0, 1)) 

wave_proportions <- up_ind_pred_sure_data$proportion[up_ind_pred_sure_data$type=="wave"]
range(wave_proportions) #  0.006981982 0.421197144
therm_proportions <- up_ind_pred_sure_data$proportion[up_ind_pred_sure_data$type=="thermal"]
range(therm_proportions) # 0.1985507 0.9282760
orog_proportions <- up_ind_pred_sure_data$proportion[up_ind_pred_sure_data$type=="orog"]
range(orog_proportions) # 0.06051508 0.45362319

