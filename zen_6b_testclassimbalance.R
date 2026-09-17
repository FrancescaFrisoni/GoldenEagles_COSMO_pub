# This code run a random forest model with 19035 labelled uplift events
# only behavioural predictors
# class-inbalanced version for Supplementary Materials
# Francesca Frisoni - July 2026. Konstanz

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


# load rf_data and predictors_pc as defined in zen_6_rf_behavclass
rf_data <- readRDS("./uplift_classification/rf_data_15july26.rds") 
predictors_pc <- readRDS("./uplift_classification/predictors_pc_15july26.rds") 


#______________________________________________________
# RANDOM FOREST, only with BEHAVIOURAL variables (10 PC)

set.seed(123)

# 10fold cross-validation, prob_ls is not only the list of probabilities but also of the 10 rf models

folds_rf <- createFolds(rf_data$uplift_type, k = 10, list = TRUE)
# with 10 folds, each test dataset contains 10% of the data, 90% used for training, and each segments is not repeated across models

# run rf models
prob_ls <- lapply(1:10, function(i){
  
  test_idx <- folds_rf[[i]]
  train <- rf_data[-test_idx, ]
  test  <- rf_data[test_idx, ]
  
  # with i=1
  # table(train$uplift_type)
  # orog thermal    wave 
  #  917   15450     764 
  
  # table(test$uplift_type)
  # orog thermal    wave 
  # 102    1717      85 
  
  # --
  # if random sampling 80/20%, but now used folds bc more accurate and each segment goes into testing only once
  # trainIndex <- createDataPartition(rf_data$uplift_type,
  #                                   p = 0.80,
  #                                   list = FALSE)
  # train <- rf_data[trainIndex,]
  # # table(train$uplift_type)
  # # orog thermal    wave
  # # 816   13734     680 
  # test <- rf_data[-trainIndex,]
  # # table(test$uplift_type)
  # # orog thermal    wave
  # #  203    3433     169 
  # --
  
  rf <- randomForest(as.formula(paste0("uplift_type ~", paste(predictors_pc, collapse="+"))), data=train, proximity=TRUE, importance=TRUE)
  print(rf)  #to check while running
  
  # var_imp <- importance(rf)
  # gini_importance <- var_imp[, "MeanDecreaseGini"]
  # gini_order <- sort(gini_importance, decreasing = TRUE)
  # top_gini <- head(gini_order, 20)
  # print(top_gini)
  
  test_pred <- predict(rf, newdata = test)
  # create the confusion matrix
  cm <- confusionMatrix(test_pred, test$uplift_type)
  cm # to check while running
  
  test_prob <- predict(rf, newdata = test, type = "prob")
  prob <- as.data.frame(round(test_prob, 2))
  prob$true <- test$uplift_type
  prob$predict <- test_pred
  prob$unique_segmID <- test$unique_segmID
  prob$accuracy <- sum(prob$predict==prob$true) / nrow(test)
  
  
  return(list(rf,prob)) # returns both the rf models and the probability matrices
})

# a typical print-out of one model (the 10th in my prob_ls run)
# Call:
#   randomForest(formula = as.formula(paste0("uplift_type ~", paste(predictors_pc,      collapse = "+"))), data = train, proximity = TRUE, importance = TRUE) 
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 3
# 
# OOB estimate of  error rate: 9.79%
# Confusion matrix:
#   orog thermal wave class.error
# orog      11     903    3 0.988004362
# thermal    6   15436    9 0.000970811
# wave       4     753    8 0.989542484


# save prob_ls bc very heavy (several hours to run) and easier to just re-load
saveRDS(prob_ls, file = "./uplift_classification/rfbehav_10fold_10pc_15july26.rds")

#______________________________________________________
# Variable importance and accuracies of Random Forest models

### Extract variables contribution for each rf model 
var_contribution_list <- lapply(prob_ls, function(model) {
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
# 1     Dim.5       366.3522
# 2     Dim.7       324.3937
# 3    Dim.10       322.0500
# 4     Dim.3       319.2028
# 5     Dim.2       309.1176
# 6     Dim.4       304.2803
# 7     Dim.9       299.9695
# 8     Dim.6       293.8353
# 9     Dim.8       290.5965
# 10    Dim.1       285.0371

# combined_gini <- do.call(rbind, lapply(seq_along(var_contribution_list), function(i) {
#   df <- var_contribution_list[[i]]$TopGini
#   df$Model <- paste0("Model_", i)
#   df
# }))
# dim 10, 5, 7 almost always at the top across models


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
# 71   PC8       roll_sum  14.57892465
# 72   PC8      roll_mean  14.56580843
# 73   PC8          h_min   5.97185755
# 74   PC8          h_max   5.85833289
# 75   PC8     roll_sum_r   4.86143154
# 76   PC8      pitch_max   4.57947139
# 77   PC8  pitch_max_abs   3.43164074
# 78   PC8   roll_sum_abs   3.28028603
# 79   PC8  roll_mean_abs   3.25283665
# 80   PC8       roll_min   2.87887389
# 81   PC9  turnangle_var   8.62236410
# 82   PC9   n_turnChange   6.86556651
# 83   PC9   turnChange_r   6.40897449
# 84   PC9   turnangle_sd   5.12614268
# 85   PC9          h_min   4.91264220
# 86   PC9     sdACCz_min   4.68275689
# 87   PC9  pitch_max_abs   4.39909097
# 88   PC9       roll_min   4.30137430
# 89   PC9      VedBA_min   4.17625725
# 90   PC9   roll_min_abs   3.83931101
# 91  PC10  turnangle_sum  11.41563803
# 92  PC10      n_circles  10.69040653
# 93  PC10       duration   9.18180349
# 94  PC10       deltah_r   8.51805104
# 95  PC10   turnChange_r   8.36792767
# 96  PC10        vel_min   7.37398180
# 97  PC10       vel_mean   5.37742048
# 98  PC10        vel_max   4.49020167
# 99  PC10         deltah   3.13667986
# 100 PC10     roll_sum_r   2.56071334


### Extract accuracies 
# Calculate accuracies using the ratio from the probability dataframe
accuracies <- sapply(prob_ls, function(x) {
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

# "Mean accuracy: 0.902"
# "Standard error of accuracy: 4e-04"
# "Accuracy range (mean ± SE): 0.9016 to 0.9024"
# "Minimum accuracy: 0.8997"
# "Maximum accuracy: 0.9043"
# "Accuracy range (min to max): 0.8997 to 0.9043"

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
class_accuracies_list <- lapply(prob_ls, function(model) {
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

# [1] "orog - Mean accuracy: 0.0088"
# [1] "orog - Standard error: 0.0027"
# [1] "orog - Accuracy range (mean ± SE): 0.0061 to 0.0116"
# [1] "orog - Accuracy range (min to max): 0 to 0.0198"
# [1] "---"
# [1] "thermal - Mean accuracy: 0.999"
# [1] "thermal - Standard error: 2e-04"
# [1] "thermal - Accuracy range (mean ± SE): 0.9988 to 0.9993"
# [1] "thermal - Accuracy range (min to max): 0.9977 to 1"
# [1] "---"
# [1] "wave - Mean accuracy: 0.0118"
# [1] "wave - Standard error: 0.0047"
# [1] "wave - Accuracy range (mean ± SE): 0.0071 to 0.0165"
# [1] "wave - Accuracy range (min to max): 0 to 0.0476"
# [1] "---"



# ____________________________________________
#### Visualization of pooled across 10 models

# Pool all predictions across the 10 runs
all_predictions <- do.call(rbind, lapply(prob_ls, function(x) x[[2]]))

## 1: Ternary plot
color_palette <- c(
  "thermal" = "#D55E00",
  "orog" = "#CC79A7", 
  "wave" = "#0072B2"
)

t <- ggtern(all_predictions, aes(orog, thermal, wave, color = true)) +  #put prob_ls[[1]][[2]] instead of all_predictions if you want the predictions and visualization of a single run, here example model 1
  geom_point(size = 8, alpha = 0.7) +
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

ggsave(file.path(directory, "figures_july26", "tern_rf.pdf"), 
       plot = t, width = 297, height = 210, units = "mm", device = "pdf")

# save ggplot obj to reload if some fig tuning is necessary
saveRDS(t, file.path(directory, "figures_july26", "tern_rf.rds"))
# t <- readRDS(file.path(directory, "figures_july26", "tern_dfa.rds"))
# load libraries ggplot2 and ggtern
# print(t)
# t <- t + theme(text = element_text(size = 20)) # here modify as usual ggplot


## 2: Confusion matrix plot

# cm <- confusionMatrix(prob_ls[[1]][[2]]$pred, prob_ls[[1]][[2]]$true) # if only 1 models output
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

ggsave(file.path(directory, "figures_july26", "cm_rf.pdf"), 
       plot = z, width = 297, height = 210, units = "mm", device = "pdf")

# save ggplot obj to reload if some fig tuning is necessary
saveRDS(z, file.path(directory, "figures_july26", "cm_rf.rds"))


#__________________________
### ROC curves and AUC 

# Convert true labels and predicted labels to factors with the same levels
all_predictions$true <- factor(all_predictions$true, levels = c("thermal", "orog", "wave"))
all_predictions$predict <- factor(all_predictions$predict, levels = c("thermal", "orog", "wave"))

# Calculate multiclass ROC
roc_multiclass <- multiclass.roc(as.numeric(all_predictions$true), as.numeric(all_predictions$predict))

# Get the AUC for multiclass
auc_multiclass <- auc(roc_multiclass)
print(auc_multiclass) #0.5067

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

# "Mean AUC: 0.507 \n Thermal: 0.508 \n Orographic: 0.504 \n Wave: 0.505"

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

ggsave(file.path(directory, "figures_july26", "roc_rf.pdf"), 
       plot = roc, width = 297, height = 210, units = "mm", device = "pdf")

# save ggplot obj to reload if some fig tuning is necessary
saveRDS(roc, file.path(directory, "figures_july26", "roc_rf.rds"))

# roc <-loadRDS()