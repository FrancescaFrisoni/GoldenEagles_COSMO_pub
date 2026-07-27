# This code plot the 19035 labelled uplift events:
# - table with flight parameters across uplift types
# - monthly proportion of uplift types
# - boxplots of COSMO and topographic parameters
# Francesca Frisoni - updated 18 july 2026. Konstanz

library(lubridate)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(scales)
library(plotly)
library(cowplot)
library(webshot2)

setwd("/home/francesca/ownCloud/TesiFrancesca_UpliftClassification/Data")
directory <- "/home/francesca/ownCloud/TesiFrancesca_UpliftClassification"

segm_pred_sure <- readRDS("./uplift_classification/segm_pred_sure_19035.rds")

color_palette <- c(
  "thermal" = "#D55E00",
  "orog" = "#CC79A7", #"#6A0DAD"
  "wave" = "#0072B2"
)


#____________________________________________________
#### Table of flight parameters across uplift types

summary_pred_sure <- readRDS("./uplift_classification/summarypredsure_15july26.rds")

pred_summary <- summary_pred_sure %>%
  group_by(pred) %>%
  summarise(across(where(is.numeric),     # Only numeric columns
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd = ~sd(.x, na.rm = TRUE),
                        min = ~min(.x, na.rm = TRUE),
                        max = ~max(.x, na.rm = TRUE)),
                   .names = "{col}_{fn}"),
            .groups = "drop")


plot(summary_pred_sure$deltah, summary_pred_sure$vel_mean)

saveRDS(pred_summary, file = "./uplift_classification/pred_summary_tablemetricsbehav.rds")

#____________________________________________________________
# Trajectories visualization for orog, thermal and wave ####

gps_df <- readRDS("GPS_allsoaringpoints_areavol.rds")
segm_df <- merge(gps_df, segm_pred_sure[, c("unique_segmID", "uplift_type")], by = "unique_segmID") # 1234636 obs

# summary(segm_df$duration)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 7.001  59.000 121.000 149.655 209.000 776.000

subset_duration <- segm_df[segm_df$duration > 250, ]

case_ls <- split(subset_duration, subset_duration$unique_segmID) # 643 segments above 250 seconds, 100 above 400 secs

lapply(case_ls, function(b) {
  # b <- case_ls[[3]]
  animalID <- unique(b$individual_local_identifier)
  b <- b[order(b$timestamp), ]
  # calculate aspect ratios for plot along 3 axes
  rangeLong <- max(b$location_long) - min(b$location_long)
  rangeLat <- max(b$location_lat) - min(b$location_lat)
  rangeLat_m <- (rangeLat * 111.139) * 1000 # transform range in metres (approximated) to compare with elevation
  rangeElev <- max(b$height_above_ellipsoid) - min(b$height_above_ellipsoid)
  ratioLat <- 1
  ratioLong <- rangeLong / rangeLat
  ratioElev <- rangeElev / rangeLat_m
  aspects <- c(x = ratioLong, y = ratioLat, z = ratioElev)
  aspects <- aspects / aspects[which.max(aspects)]
  
  # Ensure uplift_type is a factor with the correct levels
  b$uplift_type <- factor(b$uplift_type, levels = names(color_palette))
  
  # Create a color vector matching the data points
  point_colors <- color_palette[b$uplift_type] # color_palette with 3 categories
  
  p <- plot_ly(b, x = ~location_long, y = ~location_lat, z = ~height_above_ellipsoid) %>%
    add_trace(
      type = "scatter3d", mode = "lines",
      line = list(color = "darkgrey", width = 1),
      name = ~unique_segmID,
      showlegend = FALSE
    ) %>%
    add_trace(
      type = "scatter3d", mode = "markers",
      marker = list(
        color = point_colors,
        size = 5,
        line = list(width = 0),
        opacity = 0.8
      ),
      text = ~ paste0(
        "vert.speed: ", round(vertSpeed_smooth, 2), "\n",
        "ground.speed: ", round(grSpeed, 2), "\n",
        "height(aboveground): ", round(height_agl, 2), "\n",
        "uplift_type: ", uplift_type
      ),
      hoverinfo = "text",
      showlegend = FALSE
    ) %>%
    layout(
      title = paste0(
        "Animal ", unique(b$individual_local_identifier),
        " - Date: ", unique(b$date),
        " - event: ", unique(b$uplift_type)
      ),
      scene = list(
        xaxis = list(title = "Longitude"),
        yaxis = list(title = "Latitude"),
        zaxis = list(title = "Height above ellipsoid (m)"),
        aspectmode = "manual",
        aspectratio = list(x = aspects["x"], y = aspects["y"], z = aspects["z"])
      )
    )
  
  # Add a custom legend
  for (type in levels(b$uplift_type)) {
    p <- p %>% add_trace(
      x = c(), y = c(), z = c(),
      type = "scatter3d",
      mode = "markers",
      marker = list(color = color_palette[type], size = 10),
      name = type,
      showlegend = TRUE
    )
  }
  
  print(p)
  readline(prompt = "Press [enter] to continue")
})

# save high resolution pics from html
html_folder <- file.path(directory, "figures_july26", "html_traj")
all_files <- list.files(html_folder, pattern = "\\.html$")

for (f in all_files) {
  webshot2::webshot(
    file.path(html_folder, f),
    file = file.path(directory, "figures_july26", gsub("\\.html$", "_final.png", f)),
    vwidth = 1600, vheight = 1200, zoom = 3
  )
}



# ____________________________________
#### Proportion of uplift used in 2023

## proportions with written sample size on top of bars
# 2023 for both mine and Tom's labelled datasets: 3 categories

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

segm_pred_data <- prepare_data(segm_pred_sure, "segm_pred_sure")

# Plot of 2023 in 3 categories
ggplot(segm_pred_data, aes(x = month, y = proportion, fill = type)) +
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

wave_proportions <- segm_pred_data$proportion[segm_pred_data$type=="wave"]
range(wave_proportions) #  0.002715915 0.705128205
therm_proportions <- segm_pred_data$proportion[segm_pred_data$type=="thermal"]
range(therm_proportions) # 0.1474359 0.9679522
orog_proportions <- segm_pred_data$proportion[segm_pred_data$type=="orog"]
range(orog_proportions) # 0.0000000 0.1944954

# Plot 2023: 2 categories (thermal vs dynamic)

prepare_data_dynamic <- function(df, dataset_name) {
  df$month <- factor(lubridate::month(df$date), levels = 1:12, labels = month.name)
  df$year <- year(df$date)
  
  df_summary <- df %>%
    filter(year == 2023) %>%
    group_by(month) %>%
    summarize(
      orog_wave = sum(uplift_type %in% c("orog", "wave"), na.rm = TRUE),
      thermal = sum(uplift_type == "thermal", na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(total = thermal + orog_wave) %>%
    pivot_longer(cols = c(thermal, orog_wave), names_to = "type", values_to = "count") %>%
    mutate(proportion = count / total,
           dataset = dataset_name)
  
  return(df_summary)
}

segm_pred_data_dyn <- prepare_data_dynamic(segm_pred_sure, "Uplift prediction")

# add column with labels (only for thermals to avoid duplicates)
segm_pred_data_dyn <- segm_pred_data_dyn %>%
  group_by(dataset, month) %>%
  mutate(label = if_else(type == "orog_wave", as.character(first(total)), NA_character_)) %>%
  ungroup()

# re_order levels to have thermals on top
segm_pred_data_dyn$type <- factor(segm_pred_data_dyn$type, levels = c("thermal", "orog_wave"))

up <- ggplot(segm_pred_data_dyn, aes(x = month, y = proportion, fill = type)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = label, y = 1.02), size = 4.5, na.rm = TRUE) + 
  scale_fill_manual(values = c("thermal" = "#F4A261", "orog_wave" ="#6495ED"),
                    labels = c("thermal" = "Thermal", "orog_wave" = "Dynamic")) +
  labs(
    x = "Month",
    y = "Proportion of Uplift Types",
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
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1.1), breaks = seq(0, 1, by = 0.2))

print(up)

ggsave(file.path(directory, "figures_july26", "dyn_th_upliftuse2023.pdf"), 
       plot = up, width = 297, height = 210, units = "mm", device = "pdf")

# save ggplot obj to reload if some fig tuning is necessary
saveRDS(up, file.path(directory, "figures_july26", "dyn_th_upliftuse2023.rds"))

# to report the max and min proportion of dynamic uplift sources
dyn_proportions <- segm_pred_data_dyn$proportion[segm_pred_data_dyn$type=="orog_wave"]
range(dyn_proportions) # 0.0320478 0.8525641

#___________________________________________________
### Visual validation of classification: boxplots 

# Ensure correct factor levels and labels
segm_pred_sure$uplift_type <- factor(segm_pred_sure$pred, levels = c("orog", "thermal", "wave"), labels = c("Orographic", "Thermal", "Wave"))

sub_lab_7s <- readRDS("./labelled209segmentsabove7s_predictors_150726.rds")
sub_lab_7s$uplift_type <- factor(sub_lab_7s$uplift_type, levels = c("orog", "thermal", "wave"), labels = c("Orographic", "Thermal", "Wave"))

color_palette_lab <- c(
  "Thermal" = "#D55E00",
  "Orographic" = "#CC79A7", #"#6A0DAD"
  "Wave" = "#0072B2"
)

# Function to create a customized boxplot
create_boxplot <- function(data, x, y, title, y_label, y_limits) {
  # Ensure uplift_type is a factor
  data[[x]] <- factor(data[[x]])
  
  ggplot(data, aes_string(x = x, y = y, fill = x)) +
    geom_boxplot() +
    scale_fill_manual(values = color_palette_lab, drop = FALSE) +
    labs(title = title, x = "", y = y_label) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.title.x = element_text(size = 20, face = "bold"),  
          axis.title.y = element_text(size = 23), #, face = "bold"), removed bold bc issues with Brunt-Vaisala rendering due to Umlaut
          axis.text.x = element_text(size = 20, angle = 45, hjust = 1),  
          axis.text.y = element_text(size = 20), 
          legend.position = "none") +
    coord_cartesian(ylim = y_limits)
}


# Create individual plots - cosmo variables
p1 <- create_boxplot(segm_pred_sure, "uplift_type", "w_oro_mean", 
                     "", 
                     "Orographic Uplift Proxy \n[m s-1]",
                     range(c(segm_pred_sure$w_oro_mean, sub_lab_7s$w_oro_mean), na.rm = TRUE))

p2 <- create_boxplot(sub_lab_7s, "uplift_type", "w_oro_mean", 
                     "", 
                     "Orographic Uplift Proxy \n[m s-1]",
                     range(c(segm_pred_sure$w_oro_mean, sub_lab_7s$w_oro_mean), na.rm = TRUE))

p3 <- create_boxplot(segm_pred_sure, "uplift_type", "ASHFL_S_mean", 
                     "", 
                     "Surface Sensible Heat Flux \n[W m-2]",
                     range(c(segm_pred_sure$ASHFL_S_mean, sub_lab_7s$ASHFL_S_mean), na.rm = TRUE))

p4 <- create_boxplot(sub_lab_7s, "uplift_type", "ASHFL_S_mean", 
                     "", 
                     "Surface Sensible Heat Flux \n[W m-2]",
                     range(c(segm_pred_sure$ASHFL_S_mean, sub_lab_7s$ASHFL_S_mean), na.rm = TRUE))

p5 <- create_boxplot(segm_pred_sure, "uplift_type", "N2_max_h_mean", 
                     "", 
                     "Brunt–Väisälä Frequency \n(height of max stability) [m]",
                     range(c(segm_pred_sure$N2_max_h_mean, sub_lab_7s$N2_max_h_mean), na.rm = TRUE))

p6 <- create_boxplot(sub_lab_7s, "uplift_type", "N2_max_h_mean", 
                     "", 
                     "Brunt–Väisälä Frequency \n(height of max stability) [m]",
                     range(c(segm_pred_sure$N2_max_h_mean, sub_lab_7s$N2_max_h_mean), na.rm = TRUE))

p7 <- create_boxplot(segm_pred_sure, "uplift_type", "N2_max_value_mean", 
                     "", 
                     "Brunt–Väisälä Frequency \n(value of max stability) [K m-1]",
                     range(c(0, sub_lab_7s$N2_max_value_mean), na.rm = TRUE))

p8 <- create_boxplot(sub_lab_7s, "uplift_type", "N2_max_value_mean", 
                     "", 
                     "Brunt–Väisälä Frequency \n(value of max stability) [K m-1]",
                     range(c(0, sub_lab_7s$N2_max_value_mean), na.rm = TRUE))


# left column my  analyses, right 204 labelled by Tom
left_column <- plot_grid(p1, p3, p5, p7, ncol = 1, align = 'v')
right_column <- plot_grid(p2, p4, p6, p8, ncol = 1, align = 'v')

# subtitles
left_subtitle <- ggdraw() + 
  draw_label("DFA classification of complete dataset (this study)", 
             fontface = 'bold', size = 20, hjust = 0.5)

right_subtitle <- ggdraw() + 
  draw_label("Manual classification of 209 segments (Carrard et al., 2025)", 
             fontface = 'bold', size = 20, hjust = 0.5)

# add subtitles to columns
left_column_with_subtitle <- plot_grid(left_subtitle, left_column, 
                                       ncol = 1, rel_heights = c(0.07, 1))
right_column_with_subtitle <- plot_grid(right_subtitle, right_column, 
                                        ncol = 1, rel_heights = c(0.07, 1))

# Combine columns with some space between them
combined_plot <- plot_grid(left_column_with_subtitle, NULL, right_column_with_subtitle, 
                           ncol = 3, rel_widths = c(1, 0.1, 1))

# main title
title <- ggdraw() + 
  draw_label("Distribution of COSMO parameters", 
             fontface = 'bold', size = 40, hjust = 0.5)

# Combine title and plot
final_plot <- plot_grid(title, combined_plot, ncol = 1, rel_heights = c(0.07, 1))

ggsave(file.path(directory, "figures_july26", "boxplot_cosmo.pdf"), 
       plot = final_plot, width = 18, height = 20, units = "in", device = cairo_pdf) # portrait sizes 210 and 297 mm not big enough for this plot

# save ggplot obj to reload if some fig tuning is necessary
saveRDS(final_plot, file.path(directory, "figures_july26", "boxplot_cosmo.rds"))


#  Same with wind parameters: some variables overwritten bc not needed anymore
p9 <- create_boxplot(segm_pred_sure, "uplift_type", "W_mean", 
                     "", 
                     "Abs. Vertical Wind Speed \n[m s-1]",
                     range(c(segm_pred_sure$W_mean, sub_lab_7s$W_mean), na.rm = TRUE))

p10 <- create_boxplot(sub_lab_7s, "uplift_type", "W_mean", 
                      "", 
                      "Abs. Vertical Wind Speed \n[m s-1]",
                      range(c(segm_pred_sure$W_mean, sub_lab_7s$W_mean), na.rm = TRUE))

p11 <- create_boxplot(segm_pred_sure, "uplift_type", "windspeed_mean", 
                      "", 
                      "Horizontal Wind Speed \n[m s-1]",
                      range(c(segm_pred_sure$windspeed_mean, sub_lab_7s$windspeed_mean), na.rm = TRUE))

p12 <- create_boxplot(sub_lab_7s, "uplift_type", "windspeed_mean", 
                      "", 
                      "Horizontal Wind Speed \n[m s-1]",
                      range(c(segm_pred_sure$windspeed_mean, sub_lab_7s$windspeed_mean), na.rm = TRUE))

left_column <- plot_grid(p9,p11 , ncol = 1, align = 'v')
right_column <- plot_grid(p10,p12, ncol = 1, align = 'v')

left_subtitle <- ggdraw() + 
  draw_label("DFA classification of complete dataset (this study)", 
             fontface = 'bold', size = 20, hjust = 0.5)

right_subtitle <- ggdraw() + 
  draw_label("Manual classification of 209 segments (Carrard et al., 2025)", 
             fontface = 'bold', size = 20, hjust = 0.5)

left_column_with_subtitle <- plot_grid(left_subtitle, left_column, 
                                       ncol = 1, rel_heights = c(0.07, 1))
right_column_with_subtitle <- plot_grid(right_subtitle, right_column, 
                                        ncol = 1, rel_heights = c(0.07, 1))

combined_plot <- plot_grid(left_column_with_subtitle, NULL, right_column_with_subtitle, 
                           ncol = 3, rel_widths = c(1, 0.1, 1))

title <- ggdraw() + 
  draw_label("Distribution of Wind parameters", 
             fontface = 'bold', size = 40, hjust = 0.5)

final_plot <- plot_grid(title, combined_plot, ncol = 1, rel_heights = c(0.07, 1))

ggsave(file.path(directory, "figures_july26", "boxplot_wind.pdf"), 
       plot = final_plot, width = 18, height = 10, units = "in", device = cairo_pdf) # height 10 bc only 4 plots

# save ggplot obj to reload if some fig tuning is necessary
saveRDS(final_plot, file.path(directory, "figures_july26", "boxplot_wind.rds"))


# same with topography: rewritten plots bc not needed anymore above + save space in my RAM
p1 <- create_boxplot(segm_pred_sure, "uplift_type", "max_height.ab.gr", 
                     "", 
                     "Max Height Above Ground \n[m]",
                     range(c(segm_pred_sure$max_height.ab.gr, sub_lab_7s$max_height.ab.gr), na.rm = TRUE))

p2 <- create_boxplot(sub_lab_7s, "uplift_type", "max_height.ab.gr", 
                     "", 
                     "Max Height Above Ground \n[m]",
                     range(c(segm_pred_sure$max_height.ab.gr, sub_lab_7s$max_height.ab.gr), na.rm = TRUE))

p3 <- create_boxplot(segm_pred_sure, "uplift_type", "mean_slope", 
                     "", 
                     "Slope \n[degrees]",
                     range(c(segm_pred_sure$mean_slope, sub_lab_7s$mean_slope), na.rm = TRUE))

p4 <- create_boxplot(sub_lab_7s, "uplift_type", "mean_slope", 
                     "", 
                     "Slope \n[degrees]",
                     range(c(segm_pred_sure$mean_slope, sub_lab_7s$mean_slope), na.rm = TRUE))

p5 <- create_boxplot(segm_pred_sure, "uplift_type", "mean_roughness", 
                     "", 
                     "Roughness \n[m]",
                     range(c(segm_pred_sure$mean_roughness, sub_lab_7s$mean_roughness), na.rm = TRUE))

p6 <- create_boxplot(sub_lab_7s, "uplift_type", "mean_roughness", 
                     "", 
                     "Roughness \n[m]",
                     range(c(segm_pred_sure$mean_roughness, sub_lab_7s$mean_roughness), na.rm = TRUE))

p7 <- create_boxplot(segm_pred_sure, "uplift_type", "mean_aspect", 
                     "", 
                     "Aspect \n[degrees]",
                     range(c(segm_pred_sure$mean_aspect, sub_lab_7s$mean_aspect), na.rm = TRUE))

p8 <- create_boxplot(sub_lab_7s, "uplift_type", "mean_aspect", 
                     "", 
                     "Aspect \n[degrees]",
                     range(c(segm_pred_sure$mean_aspect, sub_lab_7s$mean_aspect), na.rm = TRUE))


left_column <- plot_grid(p1, p3, p5, p7, ncol = 1, align = 'v')
right_column <- plot_grid(p2, p4, p6, p8, ncol = 1, align = 'v')

left_subtitle <- ggdraw() + 
  draw_label("DFA classification of complete dataset (this study)", 
             fontface = 'bold', size = 20, hjust = 0.5)

right_subtitle <- ggdraw() + 
  draw_label("Manual classification of 209 segments (Carrard et al., 2025)", 
             fontface = 'bold', size = 20, hjust = 0.5)

left_column_with_subtitle <- plot_grid(left_subtitle, left_column, 
                                       ncol = 1, rel_heights = c(0.07, 1))
right_column_with_subtitle <- plot_grid(right_subtitle, right_column, 
                                        ncol = 1, rel_heights = c(0.07, 1))

combined_plot <- plot_grid(left_column_with_subtitle, NULL, right_column_with_subtitle, 
                           ncol = 3, rel_widths = c(1, 0.1, 1))

title <- ggdraw() + 
  draw_label("Distribution of topographic parameters", 
             fontface = 'bold', size = 40, hjust = 0.5)

# Combine title and plot
final_plot <- plot_grid(title, combined_plot, ncol = 1, rel_heights = c(0.07, 1))

ggsave(file.path(directory, "figures_july26", "boxplot_topo.pdf"), 
       plot = final_plot, width = 18, height = 20, units = "in", device = cairo_pdf) 

# save ggplot obj to reload if some fig tuning is necessary
saveRDS(final_plot, file.path(directory, "figures_july26", "boxplot_topo.rds"))

