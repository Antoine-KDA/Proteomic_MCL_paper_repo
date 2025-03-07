


# Author: Antoine KD ABALO
# Date of creation: 09 August 2024
# Last modification date: 


# Objective: Run a random forest algorithm with a recursive feature elimination (aka backward selection) analysis

#' Recursive feature elimination with random forest using LOOCV
#' Using the algorithm 2 of *20 Recursive Feature Elimination*
#' https://topepo.github.io/caret/recursive-feature-elimination.html#search
#' 




rm(list=ls())
# Import my functions
source("./Scripts/_00_List_of_functions.R")

# Install the packages
install(c("magrittr", "dplyr", "glmnet", "vip", "caret", "tidymodels", 
          "ROCR", "gridExtra", "randomForest", "MLmetrics", "foreach", "doParallel"))



# Create a directory to save the output -------
save_dir <- paste0("./Results/Filtered_NF_", gsub(" ", "_", format(Sys.Date(), "%Y %m %d")))

# Ensure the directory exists, create if not
if (!file.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}


# Load the dataset
load("./Output files/trainlist.Rdata")






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

rfe_funcs <- caret::rfFuncs
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


# Outer control for RFE - same resamples as before
rfe_ctrl <- rfeControl(
  method = "LOOCV",
  functions = rfe_funcs,
  returnResamp = "all",
  verbose = TRUE
)





# Define a function that parallelize the computations

custom_rfe_loocv_parallel <- function(data, recipe_obj, sizes, rfe_ctrl, num_cores = detectCores() - 1) {
  # Setup parallel backend
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  results <- foreach(i = 1:nrow(data), .packages = c("caret", "recipes", "dplyr")) %dopar% {
    # Define training and testing sets for LOOCV
    training_data <- data[-i, ]
    testing_data <- data[i, , drop = FALSE]  # Single observation for testing
    
    # Prep the recipe on the training data
    prepped_recipe <- recipe_obj %>%
      prep(training = training_data)
    
    # Bake the training and testing data using the prepped recipe
    processed_training <- bake(prepped_recipe, new_data = training_data)
    processed_testing <- bake(prepped_recipe, new_data = testing_data)
    
    # Extract predictors and outcome
    x_train <- processed_training[, -which(names(processed_training) == "POD24")]
    y_train <- processed_training$POD24
    
    x_test <- processed_testing[, -which(names(processed_testing) == "POD24")]
    y_test <- processed_testing$POD24
    
       
    # Apply RFE on the processed training data
    set.seed(32024)
    rfe_fit <- rfe(
      x = x_train,
      y = y_train,
      sizes = sizes,
      rfeControl = rfe_ctrl
    )
    
    # Store the results for this fold
    list(
      rfe_fit = rfe_fit,
      x_test = x_test,
      y_test = y_test
    )
  }
  
  # Stop parallel backend
  stopCluster(cl)
  
  return(results)
}



# Compute the LOOCV - RF using the own-defined function ----------------------------

# Define the results container
All_results_filtered_LOOCV <- list()

# Take the important datasets
trainlist <- trainlist[names(trainlist) %in% c("Onco Infl", "All",                                 
                                               "All: Significant p-adjusted proteins",
                                               "All: Significant p-value proteins",   
                                               "All: Onco prot for dup",
                                               "All: Median for dup" )]

# Create a loop to compute over all the datasets in the list

for(t in names(trainlist)){
  # Select the correct dataset to use
  train_data <- trainlist[[t]]
  
  # Customize the name to appear in the results
  (export_name <- names(trainlist[t]))
  (export_name <- gsub(":", "", export_name))
  
  # Define the preprocessing objects using recipe 
  
  # Define the recipe for preprocessing
  recipe_obj <- recipe(POD24 ~ ., data = train_data) %>%
    step_corr(all_predictors(), threshold = 0.75) %>%
    step_center(all_predictors()) %>%
    step_scale(all_predictors()) %>%
    step_zv(all_predictors())
  
  
  # Define the sizes of the predictors to subset
  max_sizes <- recipe_obj %>% 
    prep() %>% 
    juice(all_predictors()) %>% 
    ncol()
  
  sizes <- seq(2, max_sizes, by = 10); rm(max_sizes)
  
  # Record the starting time
  start_time <- Sys.time()
  
  # Execute your function
  results_parallel <- custom_rfe_loocv_parallel(data = train_data, 
                                                recipe_obj = recipe_obj, 
                                                sizes = sizes, 
                                                rfe_ctrl = rfe_ctrl)
  
  # Record the ending time
  end_time <- Sys.time()
  
  # Calculate and print the runtime
  runtime <- end_time - start_time
  print(paste("Start time:", start_time))
  print(paste("End time:", end_time))
  print(paste("Total runtime:", runtime))
  
  # Export time
  results_parallel[["start_time"]] <- start_time 
  results_parallel[["end_time"]] <- end_time 
  results_parallel[["runtime"]] <- runtime 
  
  # Create a new object with the desired name
  assign(paste0(export_name, "_results_parallel_LOOCV"), results_parallel)
  
  # Save the model
  to_export <- paste0(export_name, "_results_parallel_LOOCV")
  save(list = to_export, file = file.path(save_dir, paste0("RFrfe_Filtered_", 
                                                           paste0(export_name, "_results_parallel_LOOCV"), 
                                                           "_station.Rdata")))
  
  All_results_filtered_LOOCV[[paste0(export_name, "_results_parallel_LOOCV")]] <- results_parallel
  
  # Remove the old object if no longer needed
  rm(recipe_obj, train_data, end_time, export_name, runtime, sizes, 
     start_time, to_export)
  
  
}