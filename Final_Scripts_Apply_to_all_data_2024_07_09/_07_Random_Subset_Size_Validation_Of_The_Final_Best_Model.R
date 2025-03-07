









# Author: Antoine KD ABALO
# Date of creation: 08 April 2024
# Last modification date: 


# Objective: Run a Bagging algorithm with a recursive feature elimination (aka backward selection) analysis



rm(list=ls())
# Import my functions
source("./Scripts/_00_List_of_functions.R")

# Install the packages
install(c("magrittr", "dplyr", "glmnet", "vip", "caret", "tidymodels", 
          "ROCR", "gridExtra", "randomForest", "MLmetrics"))



# Create a directory to save the output -------
save_dir <- paste0("./Results/Filtered_NF_", gsub(" ", "_", format(Sys.Date(), "%Y %m %d")))
# save_dir <- "./Results/Filtered_NF_2024_04_06_KI"

# Ensure the directory exists, create if not
if (!file.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

# Open the sink to redirect output to the log file
# sink(file = file.path(save_dir, "Simple_treeBagRFE.Rout"), append = FALSE, type = "output")


# Load the dataset
load("./Output files/trainlist.Rdata")

load("./Results/Filtered_NF_2024_04_09_KI/treeBagRFE_Filtered_NF_Size5_station.Rdata")
# Get the best model
BestModel <- treeBagRFE_results_list$All_filtered


# Define the training control functions -----------------

# RFE using ROC Curve Ranking 
many_stats <-
  function(data, lev = levels(data$obs), model = NULL) {
    c(
      twoClassSummary(data = data, lev = levels(data$obs), model),
      prSummary(data = data, lev = levels(data$obs), model),
      mnLogLoss(data = data, lev = levels(data$obs), model),
      defaultSummary(data = data, lev = levels(data$obs), model)
    )
  }

rfe_funcs <- caret::treebagFuncs
rfe_funcs$summary <- many_stats

# Use the ROC AUC values for feature ranking
rfe_funcs$rank <- function(object, x, y) {
  roc_vals <- filterVarImp(x, y)
  roc_vals$var <- rownames(roc_vals)
  names(roc_vals)[1] <- "Overall"
  rownames(roc_vals) <- NULL
  roc_vals$control <- NULL
  roc_vals[order(-roc_vals$Overall),, drop = FALSE]
}


# Outer control for RFE; same resamples as before
rfe_ctrl <- rfeControl(
  method = "repeatedcv",
  number = 10, # Change from 5 to 10
  repeats = 10,
  functions = rfe_funcs,
  returnResamp = "all",
  verbose = TRUE
)

# RFE with Random Forest Ranking 
imp_ctrl <- rfe_ctrl
imp_ctrl$functions$rank <- caret::rfFuncs$rank



## define the recipe object

treebag_filtered_recipe <-
  recipe(POD24 ~ ., data = trainlist$All) %>%
  step_corr(all_predictors(), threshold = 0.75) %>%
  step_center(all_predictors()) %>%
  step_scale(all_predictors()) %>%
  step_zv(all_predictors())

# Get the number of variables in the best model

treebag_basic_pred <-
  prep(treebag_filtered_recipe, 
       training = trainlist$All, 
       retain = TRUE) %>%
  juice(all_predictors(), composition = "data.frame")

ncol(treebag_basic_pred)
# [1] 134
ncol(trainlist$All)
# [1] 339
# Variables retained in the final best model
length(BestModel$optVariables)
# [1] 58




# Random Subset Size Validation ------------------------------------------------





BestModel
print(BestModel)


nsize <- 5 # The best nsize obtained 

(sizes  <-  floor(2:(ncol(trainlist$All)/nsize)))

# Random subset size validation ------------------------------------------------

(subset_size <- length(BestModel$optVariables))
iter <- 100
# vars <- names(treebag_basic_pred)
vars <- names(trainlist$All[, !(names(trainlist$All) %in% c("POD24"))]) # Just for a try

subset_ctrl <- trainControl(
  method = "repeatedcv",
  number = 10, # Change from 5 to 10
  repeats = 10,
  index = imp_ctrl$index,
  indexOut = imp_ctrl$indexOut,
  classProbs = TRUE,
  summaryFunction = many_stats
)





for (i in 1:iter) {
  set.seed(32024 + i)
  rand_subset <- sample(vars, subset_size)
  
  ## Create the recipe object for the index subset data
  subset_rec <-
    recipe(POD24 ~ ., data = trainlist$All[, c("POD24", rand_subset)]) %>%
    step_corr(all_predictors(), threshold = 0.75) %>%
    step_center(all_predictors()) %>%
    step_scale(all_predictors()) %>%
    step_zv(all_predictors())
  
  ## Train the index subset data
  set.seed(32024 + i)
  subset_model <- train(
    subset_rec,
    data = trainlist$All[, c("POD24", rand_subset)],
    sizes = sizes,
    rfeControl = subset_ctrl,
    method = "treebag",
    # metric = "ROC",
    ntree = 5000)
  
  # Get the performances of the trained index subset data
  
  subset_perf <- getTrainPerf(subset_model)
  # subset_perf$ncomp <- subset_model$bestTune$ncomp
  if (i == 1) {
    size_check <- subset_perf
  } else {
    size_check <- bind_rows(size_check, subset_perf)
  }
  rm(rand_subset, subset_rec, subset_model, subset_perf)
}

save.image(file = file.path(save_dir, paste0("Random_Subset_Size_Validation_Station_Total_Predictors_Set.Rdata")))

subset_size
# mean(size_check$TrainROC <= BestModel$results$ROC)
table(size_check$TrainAccuracy*100 <= 78.2)
mean(size_check$TrainAccuracy*100 <= 78.2)

head(size_check)
size_check
max(size_check$TrainAccuracy)
