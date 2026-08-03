# This code run a random forest model with 19035 labelled uplift events
# only behavioural predictors
# Francesca Frisoni - August 13th, 2024. Konstanz

setwd("/home/francesca/ownCloud/TesiFrancesca_UpliftClassification/Data")
directory <- "/home/francesca/ownCloud/TesiFrancesca_UpliftClassification"

library(randomForest)
library(caret)
library(lubridate)
library(dplyr)
library(ggtern)
library(pROC)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(data.table)
library(FactoMineR)
library(factoextra)
library(scales)

##### 1. BEHAV DATASET: MERGE PREDICTED UPLIFTS TYPES WITH BEHAVIOURAL METRICS #####
# this chunk of code combines uplift predictions with GPS, ACC and IMU variables

#______________________________________________________________________________________________________
# Complete behav dataset, not only segm pred sure: Combine predictions with GPS, ACC and IMU variables 

segm_pred <- readRDS("./uplift_classification/summary_segments_cosmo_topography_behav_allpredictions_25161.rds")
# extract month and year
segm_pred$month <- lubridate :: month(segm_pred$date, label = TRUE, abbr = FALSE, locale = "en_US.UTF-8")
segm_pred$year <- year(segm_pred$date)

# order by month as factor
segm_pred$month <- factor(segm_pred$month, levels = month.name, ordered = TRUE)

# define uplift types with my thresholds
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

# table(segm_pred$uplift_type)
# orog    orog/wave      thermal thermal/orog thermal/wave      unknown         wave 
# 1019          318        17167          183          322         5303          849

# load here summary of segments with gps, acc and imu metrics
combined_df <- readRDS("GPS_ACC_IMU_summaryvariables_28july.rds") # 29775 in 78 obs

names(segm_pred)
# [1] "date"                        "individual_local_identifier" "unique_segmID"               "burstID"                    
# [5] "duration"                    "time_start"                  "time_end"                    "w_oro_mean"                 
# [9] "N2_max_h_mean"               "N2_max_value_mean"           "ASHFL_S_mean"                "U_mean"                     
# [13] "V_mean"                      "W_mean"                      "windspeed_mean"              "max_height.ab.gr"           
# [17] "mean_slope"                  "mean_roughness"              "mean_aspect"                 "orog"                       
# [21] "thermal"                     "wave"                        "pred"                        "month"                      
# [25] "year"                        "uplift_type"   

# merge the predictions by dfa together with eagles behaviour through GPS, ACC and IMU
# merged columns are unique_segmID, and from w_oro_mean to uplift_type
summary_pred <- merge(segm_pred[,c(3,8:26)], combined_df, by = "unique_segmID")

# I didn't load directly segm_pred_sure because I will use also the other classes later on to identify differences with mixed classes: summary_pred has all the 6 classes (I excluded unknown)
summary_pred_sure <- summary_pred[summary_pred$uplift_type %in% c("orog","wave","thermal"),]
# set order of levels
summary_pred_sure$uplift_type <- factor(summary_pred_sure$uplift_type, levels=c("orog","thermal","wave"))

# quick check
# table(summary_pred_sure$uplift_type)
#orog thermal    wave 
# 1019   17167     849

saveRDS(summary_pred_sure, "./uplift_classification/summarypredsure_15july26.rds") 
# summary pred sure is the dataset I will run my rf with only behavioural variables



##### 2. PCA algorirthm to select predictors for Random Forest #####

#___________________________________________
# PCA for 35 selected behavioural variables

# if not yet loaded:
# summary_pred_sure <- readRDS("./uplift_classification/summarypredsure_15july26.rds")

# predictor names explicitly
predictors <- c(
  "deltah", "dist", "duration", "grSpeed_mean", "h_min", "h_max",
  "vel_mean", "vel_min", "vel_max", "turnangle_mean", "turnangle_sum",
  "turnangle_var", "turnangle_sd", "n_turnChange", "VedBA_mean", "VedBA_sd",
  "VedBA_max", "VedBA_min", "ODBA_mean", "sdACCz_mean", "sdACCz_max",
  "sdACCz_min", "yaw_mean", "yaw_mean_abs", "yaw_max", "yaw_max_abs",
  "yaw_min", "yaw_min_abs", "yaw_sum", "yaw_sum_abs", "yaw_sd",
  "pitch_mean", "pitch_mean_abs", "pitch_max", "pitch_max_abs",
  "pitch_min", "pitch_min_abs", "pitch_sum", "pitch_sum_abs", "pitch_sd",
  "roll_mean", "roll_mean_abs", "roll_max", "roll_max_abs", "roll_min",
  "roll_min_abs", "roll_sum", "roll_sum_abs", "roll_sd", "deltah_r",
  "dist_r", "turnangle_r", "tilt_rad", "yaw_sum_r", "pitch_sum_r",
  "roll_sum_r", "n_circles", "n_circles_r", "turnChange_r"
)

# dataset of predictors from summary_pred_sure to make them PC components
predictors_dt <- summary_pred_sure[,c(predictors)] 

pca_result <- PCA(predictors_dt, scale.unit = TRUE, ncp = 10, graph = FALSE) 
# visual check
summary(pca_result) # first 10pc with 59 vars explain 79.007 of variance

# Call:
#   PCA(X = predictors_dt, scale.unit = TRUE, ncp = 10, graph = FALSE) 
# 
# 
# Eigenvalues
# Dim.1  Dim.2  Dim.3  Dim.4  Dim.5  Dim.6  Dim.7  Dim.8  Dim.9 Dim.10
# Variance             14.101 10.097  5.879  3.543  2.823  2.486  2.341  1.966  1.711  1.666
# % of var.            23.901 17.113  9.964  6.006  4.785  4.214  3.968  3.332  2.900  2.824
# Cumulative % of var. 23.901 41.014 50.978 56.984 61.768 65.983 69.951 73.284 76.184 79.007
# 
# Individuals (the 10 first)
# Dist    Dim.1    ctr   cos2    Dim.2    ctr   cos2    Dim.3    ctr   cos2  
# 1              |  4.399 | -0.288  0.000  0.004 | -0.492  0.000  0.013 | -0.916  0.001  0.043 |
#   2              |  6.276 |  1.076  0.000  0.029 |  1.896  0.002  0.091 | -1.978  0.003  0.099 |
#   3              |  7.063 |  1.056  0.000  0.022 |  1.405  0.001  0.040 | -3.708  0.012  0.276 |
#   13             |  7.178 | -0.612  0.000  0.007 | -0.683  0.000  0.009 | -3.765  0.013  0.275 |
#   14             |  7.753 |  4.891  0.009  0.398 |  3.120  0.005  0.162 | -1.120  0.001  0.021 |
#   15             |  7.755 | -5.294  0.010  0.466 |  1.680  0.001  0.047 |  1.314  0.002  0.029 |
#   16             |  4.528 |  2.277  0.002  0.253 | -1.396  0.001  0.095 |  0.954  0.001  0.044 |
#   17             | 13.176 |  5.683  0.012  0.186 |  9.894  0.051  0.564 |  2.760  0.007  0.044 |
#   18             |  6.045 |  0.141  0.000  0.001 | -0.674  0.000  0.012 | -0.596  0.000  0.010 |
#   40             |  4.479 |  0.872  0.000  0.038 |  1.924  0.002  0.184 |  1.354  0.002  0.091 |
#   
#   Variables (the 10 first)
# Dim.1    ctr   cos2    Dim.2    ctr   cos2    Dim.3    ctr   cos2  
# deltah         |  0.587  2.442  0.344 | -0.407  1.637  0.165 |  0.428  3.110  0.183 |
#   dist           |  0.290  0.597  0.084 | -0.350  1.216  0.123 |  0.590  5.922  0.348 |
#   duration       |  0.645  2.953  0.416 | -0.429  1.821  0.184 |  0.339  1.958  0.115 |
#   grSpeed_mean   | -0.365  0.945  0.133 |  0.071  0.050  0.005 |  0.367  2.294  0.135 |
#   h_min          | -0.220  0.344  0.049 | -0.287  0.817  0.082 |  0.184  0.575  0.034 |
#   h_max          | -0.004  0.000  0.000 | -0.408  1.646  0.166 |  0.331  1.860  0.109 |
#   vel_mean       |  0.532  2.007  0.283 | -0.352  1.225  0.124 |  0.464  3.669  0.216 |
#   vel_min        |  0.137  0.132  0.019 |  0.044  0.019  0.002 |  0.129  0.283  0.017 |
#   vel_max        |  0.541  2.073  0.292 | -0.333  1.097  0.111 |  0.506  4.354  0.256 |
#   turnangle_mean |  0.799  4.523  0.638 |  0.038  0.015  0.001 | -0.481  3.930  0.231 |


plot(pca_result)
saveRDS(pca_result, "./uplift_classification/pca_result_15july26.rds") 

# taken the first 10 components
pca_scores <- pca_result$ind$coord

# visualization for supplementary material

p1 <- fviz_eig(pca_result, addlabels = TRUE, 
               barfill = "steelblue", barcolor = "steelblue",
               linecolor = "black") +
  labs(title = NULL, x = "Principal Components", y = "% Variance Explained") +
  theme(text = element_text(size = 20))

# Find which layer is the text/label layer: enlarge the % on top of bars
p1$layers
# it's the last layer: enlarge size
p1$layers[[4]]$aes_params$size <- 7
print(p1)

p2 <- fviz_cos2(pca_result, choice = "var", axes = 1:10) +
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 24, face = "bold")
  )

print(p2)

p3 <- fviz_pca_var(pca_result, col.var = "cos2",
                   gradient.cols = c("black", "orange", "green"),
                   repel = TRUE) +
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 24, face = "bold")
  )
print(p3)

# save pdf figures and rds of ggplot obj
ggsave(file.path(directory, "figures_july26", "pca_scree.pdf"), 
       plot = p1, width = 297, height = 210, units = "mm", device = "pdf")

ggsave(file.path(directory, "figures_july26", "pca_cos2.pdf"), 
       plot = p2, width = 297, height = 210, units = "mm", device = "pdf")

ggsave(file.path(directory, "figures_july26", "pca_var.pdf"), 
       plot = p3, width = 297, height = 210, units = "mm", device = "pdf")

saveRDS(p1, file.path(directory, "figures_july26", "pca_scree.rds"))
saveRDS(p2, file.path(directory, "figures_july26", "pca_cos2.rds"))
saveRDS(p3, file.path(directory, "figures_july26", "pca_var.rds"))


# merge the 10 PC with summary_pred_sure to run the random forest: rf_data is the dataset to run the rf!
rf_data <- cbind(summary_pred_sure, pca_scores)
predictors_pc <- names(rf_data[(ncol(rf_data) - 9):ncol(rf_data)]) # from Dim 1 to Dim 10 are last 10columns

saveRDS(rf_data, "./uplift_classification/rf_data_15july26.rds") 
saveRDS(predictors_pc, "./uplift_classification/predictors_pc_15july26.rds") 

##### 2. RUN RF MODEL ON COMPLETE DATASE, USING BEHAVIOURAL METRICS AS PREDICTORS #####

#________________________________
# RANDOM FOREST, CLASS-BALANCED 

set.seed(123)

folds_rf <- createFolds(rf_data$uplift_type, k = 10, list = TRUE)

prob_ls_balanced_balanced <- lapply(1:10, function(i){
  
  test_idx <- folds_rf[[i]]
  train <- rf_data[-test_idx, ]
  test  <- rf_data[test_idx, ]
  
  # with i=1
  # table(train$uplift_type)
  # orog thermal    wave 
  # 917   15450     764 
  # table(test$uplift_type)
  # orog thermal    wave 
  # 102    1717      85 
  
  # balanced sample size: draw the same n (= size of smallest class in this fold)
  # from EACH class, for EVERY tree's bootstrap sample
  min_n <- min(table(train$uplift_type)) # 764
  samp_size <- rep(min_n, length(levels(train$uplift_type))) # 764 764 764
  
  rf <- randomForest(
    as.formula(paste0("uplift_type ~", paste(predictors_pc, collapse = "+"))),
    data = train,
    strata = train$uplift_type,
    sampsize = samp_size,
    proximity = TRUE,
    importance = TRUE
  )
  print(rf)
  
  test_pred <- predict(rf, newdata = test)
  cm <- confusionMatrix(test_pred, test$uplift_type)
  cm
  
  test_prob <- predict(rf, newdata = test, type = "prob")
  prob <- as.data.frame(round(test_prob, 2))
  prob$true <- test$uplift_type
  prob$predict <- test_pred
  prob$unique_segmID <- test$unique_segmID
  prob$accuracy <- sum(prob$predict == prob$true) / nrow(test)
  
  return(list(rf, prob))
})

# typical print-out of one run
# Call:
#   randomForest(formula = as.formula(paste0("uplift_type ~", paste(predictors_pc,      collapse = "+"))), data = train, strata = train$uplift_type,      sampsize = samp_size, proximity = TRUE, importance = TRUE) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 3
# 
# OOB estimate of  error rate: 29.74%
# Confusion matrix:
#   orog thermal wave class.error
# orog     322     427  168   0.6488550
# thermal 2243   11377 1830   0.2636246
# wave     168     258  338   0.5575916

# save prob_ls bc very heavy (several hours to run) and easier to just re-load
saveRDS(prob_ls_balanced, file = "./uplift_classification/rfbehav_classbalance_15july26.rds")

#______________________________________________________
# Variable importance and accuracies of Random Forest models

# if not yet loaded:
prob_ls_balanced <- readRDS("./uplift_classification/rfbehav_classbalance_15july26.rds")
pca_result <- readRDS("./uplift_classification/pca_result_15july26.rds") 

### Extract variables contribution for each rf model 
var_contribution_list <- lapply(prob_ls_balanced, function(model) {
  rf <- model[[1]]  # Extract the rf model
  
  # Variable importance
  var_imp <- importance(rf)
  gini_importance <- var_imp[, "MeanDecreaseGini"]
  gini_order <- sort(gini_importance, decreasing = TRUE)
  
  # Keep top 10 Gini variables and their values
  top_gini <- head(gini_order, 10)
  
  # PCA contributions
  var_info <- get_pca_var(pca_result)
  contributions <- var_info$contrib
  
  # For each PC, get top 10 contributions
  pcs_contrib_list <- lapply(1:10, function(pc) {
    sorted_contrib <- sort(contributions[, pc], decreasing = TRUE)
    head_contrib <- head(sorted_contrib, 10)
    # Return as a data.frame
    data.frame(PC = paste0("PC", pc),
               Variable = names(head_contrib),
               Contribution = as.numeric(head_contrib))
  })
  
  pcs_contrib_df <- do.call(rbind, pcs_contrib_list)
  
  # Return a list of two data frames
  list(TopGini = data.frame(Variable = names(top_gini), GiniImportance = as.numeric(top_gini)),
       TopContribs = pcs_contrib_df)
})

# returns a list, each element containing two data.frames (top gini, top contribs)
# inspect them separately:

# # var_contribution_list[[1]]$TopGini
# Variable GiniImportance
# 1     Dim.5       216.4114
# 2     Dim.7       165.4294
# 3     Dim.3       157.1738
# 4     Dim.8       154.1057
# 5     Dim.2       148.4453
# 6     Dim.9       148.2698
# 7    Dim.10       141.9940
# 8     Dim.4       133.8468
# 9     Dim.6       133.2194
# 10    Dim.1       129.1044

combined_gini <- do.call(rbind, lapply(seq_along(var_contribution_list), function(i) {
  df <- var_contribution_list[[i]]$TopGini
  df$Model <- paste0("Model_", i)
  df
}))
# dim 5, 7, 3 almost always at the top across models


# you can check the raw variables across each model, for example for the first:
# var_contribution_list[[1]]$TopContribs
# PC       Variable Contribution
# 1    PC1    n_circles_r   4.72524909
# 2    PC1 turnangle_mean   4.52301250
# 3    PC1  pitch_sum_abs   4.29019198
# 4    PC1 pitch_mean_abs   4.28355501
# 5    PC1    turnangle_r   4.22647426
# 6    PC1         dist_r   4.06303496
# 7    PC1      pitch_sum   4.02914457
# 8    PC1     pitch_mean   4.02269913
# 9    PC1   turnangle_sd   3.86817843
# 10   PC1      n_circles   3.61316181
# 11   PC2       pitch_sd   7.49897735
# 12   PC2      ODBA_mean   7.29734640
# 13   PC2     VedBA_mean   7.17847636
# 14   PC2    sdACCz_mean   6.40611763
# 15   PC2       VedBA_sd   6.03492693
# 16   PC2        roll_sd   5.80211913
# 17   PC2      VedBA_min   4.60876411
# 18   PC2     roll_sum_r   4.37241358
# 19   PC2      VedBA_max   4.21452101
# 20   PC2     sdACCz_max   4.00464390
# 21   PC3  pitch_min_abs   7.73608528
# 22   PC3           dist   5.92231322
# 23   PC3      pitch_min   4.96017313
# 24   PC3    pitch_sum_r   4.86550480
# 25   PC3    turnangle_r   4.41214567
# 26   PC3        vel_max   4.35392943
# 27   PC3      pitch_sum   4.02334510
# 28   PC3     pitch_mean   3.98471930
# 29   PC3 turnangle_mean   3.93006762
# 30   PC3  pitch_sum_abs   3.73430844
# 31   PC4   yaw_mean_abs  25.91837070
# 32   PC4    yaw_sum_abs  25.89536689
# 33   PC4    yaw_min_abs  12.79741419
# 34   PC4    yaw_max_abs  11.89218328
# 35   PC4      yaw_sum_r   9.03653718
# 36   PC4         yaw_sd   5.51371707
# 37   PC4        yaw_min   2.52975441
# 38   PC4        yaw_max   1.86993530
# 39   PC4   grSpeed_mean   0.60003667
# 40   PC4  turnangle_var   0.43647663
# 41   PC5      roll_mean  10.09027129
# 42   PC5       roll_sum  10.07309556
# 43   PC5          h_max   8.04363036
# 44   PC5          h_min   6.35283372
# 45   PC5       vel_mean   4.86850459
# 46   PC5       deltah_r   4.70695591
# 47   PC5   grSpeed_mean   4.26921839
# 48   PC5        vel_max   3.77918002
# 49   PC5  pitch_min_abs   3.47081509
# 50   PC5       roll_min   2.88394443
# 51   PC6       yaw_mean  38.84401420
# 52   PC6        yaw_sum  38.83915142
# 53   PC6        yaw_max  11.27727785
# 54   PC6        yaw_min   9.34371883
# 55   PC6     sdACCz_min   0.11367034
# 56   PC6          h_max   0.09248650
# 57   PC6      VedBA_min   0.07920250
# 58   PC6         deltah   0.07266053
# 59   PC6    sdACCz_mean   0.07101901
# 60   PC6   turnChange_r   0.06600700
# 61   PC7      roll_mean  11.73434790
# 62   PC7       roll_sum  11.73148988
# 63   PC7     sdACCz_min   6.41269278
# 64   PC7    sdACCz_mean   6.33927282
# 65   PC7       roll_max   6.26974344
# 66   PC7     VedBA_mean   4.98965923
# 67   PC7      ODBA_mean   4.68287602
# 68   PC7  roll_mean_abs   4.06244989
# 69   PC7   roll_sum_abs   3.98931957
# 70   PC7      VedBA_min   3.94431154
# 71   PC8       roll_sum  14.57892462
# 72   PC8      roll_mean  14.56580846
# 73   PC8          h_min   5.97185754
# 74   PC8          h_max   5.85833289
# 75   PC8     roll_sum_r   4.86143154
# 76   PC8      pitch_max   4.57947137
# 77   PC8  pitch_max_abs   3.43164071
# 78   PC8   roll_sum_abs   3.28028604
# 79   PC8  roll_mean_abs   3.25283667
# 80   PC8       roll_min   2.87887388
# 81   PC9  turnangle_var   8.62232862
# 82   PC9   n_turnChange   6.86560360
# 83   PC9   turnChange_r   6.40902244
# 84   PC9   turnangle_sd   5.12616028
# 85   PC9          h_min   4.91264640
# 86   PC9     sdACCz_min   4.68272256
# 87   PC9  pitch_max_abs   4.39925233
# 88   PC9       roll_min   4.30132857
# 89   PC9      VedBA_min   4.17622596
# 90   PC9   roll_min_abs   3.83931508
# 91  PC10  turnangle_sum  11.41563436
# 92  PC10      n_circles  10.69040468
# 93  PC10       duration   9.18181874
# 94  PC10       deltah_r   8.51804989
# 95  PC10   turnChange_r   8.36790173
# 96  PC10        vel_min   7.37404274
# 97  PC10       vel_mean   5.37742814
# 98  PC10        vel_max   4.49020163
# 99  PC10         deltah   3.13667236
# 100 PC10     roll_sum_r   2.56072196


### Extract accuracies 
# Calculate accuracies using the ratio from the probability dataframe
accuracies <- sapply(prob_ls_balanced, function(x) {
  # Extract the probability dataframe
  prob <- x[[2]]
  
  # Calculate accuracy as the ratio of correctly classified instances
  true_positives <- sum(prob$true == prob$pred & prob$true == "orog") + 
    sum(prob$true == prob$pred & prob$true == "thermal") + 
    sum(prob$true == prob$pred & prob$true == "wave")
  
  total_instances <- nrow(prob)
  
  accuracy <- true_positives / total_instances
  
  return(accuracy)
})

# Calculate mean accuracy
mean_accuracy <- mean(accuracies)

# Calculate standard error
n <- length(accuracies)
se_accuracy <- sd(accuracies) / sqrt(n)

# Calculate min and max accuracies
min_accuracy <- min(accuracies)
max_accuracy <- max(accuracies)

# Print results
print("Overall accuracy results:")
print(paste("Mean accuracy:", round(mean_accuracy, 4)))
print(paste("Standard error of accuracy:", round(se_accuracy, 4)))
print(paste("Accuracy range (mean ± SE):", 
            round(mean_accuracy - se_accuracy, 4), "to", 
            round(mean_accuracy + se_accuracy, 4)))
print(paste("Minimum accuracy:", round(min_accuracy, 4)))
print(paste("Maximum accuracy:", round(max_accuracy, 4)))
print(paste("Accuracy range (min to max):", 
            round(min_accuracy, 4), "to", 
            round(max_accuracy, 4)))

# "Mean accuracy: 0.7018"
# "Standard error of accuracy:  0.0042"
# "Accuracy range (mean ± SE): 0.6975 to 0.706"
# "Minimum accuracy: 0.6758"
# "Maximum accuracy: 0.725"
# "Accuracy range (min to max): 0.6758 to 0.725"

### Calculate class-specific accuracies 
calculate_class_accuracies <- function(prob_df) {
  class_accuracies <- sapply(c("orog", "thermal", "wave"), function(class) {
    class_rows <- prob_df$true == class
    sum(prob_df$pred[class_rows] == class) / sum(class_rows)
  })
  names(class_accuracies) <- c("orog", "thermal", "wave")
  return(class_accuracies)
}

# Extract probabilities from each model in the list and calculate accuracies
class_accuracies_list <- lapply(prob_ls_balanced, function(model) {
  prob <- model[[2]]  # Extract the probability dataframe
  calculate_class_accuracies(prob)
})

class_accuracies_matrix <- do.call(rbind, class_accuracies_list)

# Calculate mean, standard error, min, and max accuracies for each class
mean_class_accuracies <- colMeans(class_accuracies_matrix)
se_class_accuracies <- apply(class_accuracies_matrix, 2, function(x) sd(x) / sqrt(length(x)))
min_class_accuracies <- apply(class_accuracies_matrix, 2, min)
max_class_accuracies <- apply(class_accuracies_matrix, 2, max)

print("Class-specific results:")
for (class in c("orog", "thermal", "wave")) {
  print(paste(class, "- Mean accuracy:", round(mean_class_accuracies[class], 4)))
  print(paste(class, "- Standard error:", round(se_class_accuracies[class], 4)))
  print(paste(class, "- Accuracy range (mean ± SE):", 
              round(mean_class_accuracies[class] - se_class_accuracies[class], 4), "to", 
              round(mean_class_accuracies[class] + se_class_accuracies[class], 4)))
  print(paste(class, "- Accuracy range (min to max):", 
              round(min_class_accuracies[class], 4), "to", 
              round(max_class_accuracies[class], 4)))
  print("---")
}

# [1] "orog - Mean accuracy: 0.3729"
# [1] "orog - Standard error: 0.0187"
# [1] "orog - Accuracy range (mean ± SE): 0.3542 to 0.3916"
# [1] "orog - Accuracy range (min to max): 0.2843 to 0.4902"
# [1] "---"
# [1] "thermal - Mean accuracy: 0.7343"
# [1] "thermal - Standard error: 0.0047"
# [1] "thermal - Accuracy range (mean ± SE): 0.7296 to 0.7389"
# [1] "thermal - Accuracy range (min to max): 0.7047 to 0.7611"
# [1] "---"
# [1] "wave - Mean accuracy: 0.4394"
# [1] "wave - Standard error: 0.0201"
# [1] "wave - Accuracy range (mean ± SE): 0.4193 to 0.4595"
# [1] "wave - Accuracy range (min to max): 0.3412 to 0.5176"
# [1] "---"



# ____________________________________________
#### Visualization of pooled across 10 models

# Pool all predictions across the 10 runs
all_predictions <- do.call(rbind, lapply(prob_ls_balanced, function(x) x[[2]]))

## 1: Ternary plot
color_palette <- c(
  "thermal" = "#D55E00",
  "orog" = "#CC79A7", 
  "wave" = "#0072B2"
)

t <- ggtern(all_predictions, aes(orog, thermal, wave, color = true)) +  #put prob_ls_balanced[[1]][[2]] instead of all_predictions if you want the predictions and visualization of a single run, here example model 1
  geom_point(size = 5, alpha = 0.5) +
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
print(t)
# ignore warning, it is just internal to ggtern but nothign is wrong
# Ignoring unknown labels:
#   • R : "Wave"
# • Rarrow : "Wave"
# • T : "Thermal"
# • Tarrow : "Thermal"
# • L : "Orographic"
# • Larrow : "Orographic"

ggsave(file.path(directory, "figures_july26", "tern_rf_balanced.pdf"), 
       plot = t, width = 297, height = 210, units = "mm", device = "pdf")

# save ggplot obj to reload if some fig tuning is necessary
saveRDS(t, file.path(directory, "figures_july26", "tern_rf_balanced.rds"))
# t <- readRDS(file.path(directory, "figures_july26", "tern_dfa.rds"))
# load libraries ggplot2 and ggtern
# print(t)
# t <- t + theme(text = element_text(size = 20)) # here modify as usual ggplot


## 2: Confusion matrix plot

# cm <- confusionMatrix(prob_ls_balanced[[1]][[2]]$pred, prob_ls_balanced[[1]][[2]]$true) # if only 1 models output
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

z <- ggplot(melted_cm, aes(x = Reference, y = Prediction)) +
  geom_tile(aes(fill = Reference, alpha = Percentage), color = "black") +
  geom_text(aes(label = Value), size = 10) + # size for text labels into squares
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
z <- z +
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


print(z)

ggsave(file.path(directory, "figures_july26", "cm_rf_balanced.pdf"), 
       plot = z, width = 297, height = 210, units = "mm", device = "pdf")

# save ggplot obj to reload if some fig tuning is necessary
saveRDS(z, file.path(directory, "figures_july26", "cm_rf_balanced.rds"))


#__________________________
### ROC curves and AUC 

# Convert true labels and predicted labels to factors with the same levels
all_predictions$true <- factor(all_predictions$true, levels = c("thermal", "orog", "wave"))
all_predictions$predict <- factor(all_predictions$predict, levels = c("thermal", "orog", "wave"))

# Calculate multiclass ROC
roc_multiclass <- multiclass.roc(as.numeric(all_predictions$true), as.numeric(all_predictions$predict))

# Get the AUC for multiclass
auc_multiclass <- auc(roc_multiclass)
print(auc_multiclass) #0.6526

# empty dataframe to store FPR, TPR, and Class information
roc_data <- data.frame(FPR = numeric(), TPR = numeric(), Class = character())

# Vector to store individual AUC values
auc_values <- numeric()

# Loop through each class to extract individual ROC curves
for (class in levels(all_predictions$true)) {
  # One-vs-all ROC for the current class
  roc_class <- roc(ifelse(all_predictions$true == class, 1, 0), ifelse(all_predictions$predict == class, 1, 0))
  # Add the FPR, TPR, and class name to the dataframe
  roc_data <- rbind(roc_data, data.frame(
    FPR = 1 - roc_class$specificities,
    TPR = roc_class$sensitivities,
    Class = class
  ))
  
  # Store the AUC value for this class
  auc_values <- c(auc_values, auc(roc_class))
}

# Create a string with all AUC values
auc_string <- paste(
  "Mean AUC:", round(auc_multiclass, 3), "\n",
  "Thermal:", round(auc_values[1], 3), "\n",
  "Orographic:", round(auc_values[2], 3), "\n",
  "Wave:", round(auc_values[3], 3)
)

# "Mean AUC: 0.653 \n Thermal: 0.667 \n Orographic: 0.611 \n Wave: 0.659"

# ROC plot
roc <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Class)) +
  geom_line(size = 1.5, alpha = 1) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  labs(
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Uplift Type") +
  scale_color_manual(values = color_palette, labels = c("orog" = "Orographic", "thermal" = "Thermal", "wave" = "Wave")) + 
  theme_minimal() +
  theme(
    text = element_text(size = 24),                   
    axis.title = element_text(size = 24),             
    axis.text = element_text(size = 24),              
    legend.title = element_text(size = 24),           
    legend.text = element_text(size = 24),           
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  annotate("text", x = 0.7, y = 0.2,
           label = auc_string,
           color = "black", hjust = 0,
           fontface = "bold", size = 7) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

print(roc)

ggsave(file.path(directory, "figures_july26", "roc_rf_balanced.pdf"), 
       plot = roc, width = 297, height = 210, units = "mm", device = "pdf")

# save ggplot obj to reload if some fig tuning is necessary
saveRDS(roc, file.path(directory, "figures_july26", "roc_rf_balanced.rds"))

roc <- readRDS(file.path(directory, "figures_july26", "roc_rf_balanced.rds"))
roc$layers

# confirm layer order first 
sapply(roc$layers, function(l) class(l$geom)[1])
# "GeomLine" "GeomAbline" "GeomText"  <- annotate() is layer 3

# edit the existing theme elements' sizes
roc$theme$text$size         <- 30
roc$theme$axis.title$size   <- 30
roc$theme$axis.text$size    <- 30
roc$theme$legend.title$size <- 30
roc$theme$legend.text$size  <- 30
roc$layers[[3]]$aes_params$size <- 7

roc
ggsave(file.path(directory, "figures_july26", "roc_rf_balanced.pdf"), 
       plot = roc, width = 297, height = 210, units = "mm", device = "pdf")




