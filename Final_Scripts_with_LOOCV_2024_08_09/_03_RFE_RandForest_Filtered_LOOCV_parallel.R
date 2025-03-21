


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
install(c("magrittr", "dplyr", "glmnet", "vip", "caret", "tidymodels", "gbm",
          "ROCR", "gridExtra", "randomForest", "MLmetrics", "doParallel", "foreach"))


# Create a directory to save the output -------
save_dir <- paste0("./Results/Filtered_Only_", gsub(" ", "_", format(Sys.Date(), "%Y %m %d")))

# Ensure the directory exists, create if not
if (!file.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}


# Load the dataset
load("./Output files/trainlist_12July2024_Ant.Rdata")

assign("trainlist", trainlist_12July2024_Ant)
rm(trainlist_12July2024_Ant)



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
  
  max_sizes  <- ifelse(max_sizes > 150, 150, max_sizes)
  sizes <- seq(2, max_sizes, by = 3); rm(max_sizes)
  
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



# Save the results
save(All_results_filtered_LOOCV, file = file.path(save_dir, paste0("RFrfe_All_results_filtered_LOOCV_list.Rdata")))









# Identify and the best model with the best predictors ----------------------------------------------------




# Load the test data -------------------------
load("./Output files/testlist_12July2024_Ant.Rdata")
assign("testlist", testlist_12July2024_Ant)
rm(testlist_12July2024_Ant)


# Identify the model to take and the corresponding training and test data

## Identify the panel
names(All_results_filtered_LOOCV)

names(trainlist)

# Change the names of the trainlist elements
names(trainlist) <- gsub(":", "", names(trainlist))
# Change the names of the testlist elements
names(testlist) <- gsub(":", "", names(testlist))
# Change the names of the model list elements 
names(All_results_filtered_LOOCV) <- gsub("_results_parallel_LOOCV", "", names(All_results_filtered_LOOCV))


# Setup parallel backend
# num_cores <- detectCores() - 1
# cl <- makeCluster(num_cores)
# registerDoParallel(cl)



# Record the starting time
start_time <- Sys.time()

# Define the list of the final models
Final_Model_list_Filtered <- list()

# Run over the loop
for(running_name in names(All_results_filtered_LOOCV)){
  
  
  
  training_data <- trainlist[[running_name]]
  testing_data <- testlist[[running_name]]
  model_to_take <- All_results_filtered_LOOCV[[running_name]]
  
  message("\n\nThe current running model is: ", names(All_results_filtered_LOOCV[running_name]), "\n",
          "The current training data used is: ", names(trainlist[running_name]), "\n",
          "The current testing data used is: ", names(testlist[running_name]), "\n",
          "The running name is: ", running_name)
  
  model_to_take <- model_to_take[!names(model_to_take) %in% c("start_time","end_time","runtime")]
  
  model_to_take[[1]]
  
  # Choose a final model -------------------------------------
  
  # Initialize a list to store best subset sizes and corresponding metrics
  best_subsets <- list()
  best_subsets_df <- data.frame()
  
  # Loop through each LOOCV iteration result
  for (i in seq_along(model_to_take)) {
    rfe_fit <- model_to_take[[i]]$rfe_fit
    
    # Extract the subset size and accuracy for the selected subset
    best_index <- which(rfe_fit$results$Variables == rfe_fit$optsize)
    
    # Store in the list
    best_subsets[[i]] <- list(
      cbind(size = rfe_fit$optsize,
            rfe_fit$results[best_index , ] %>% select(-Variables), 
            Variables = rfe_fit$optVariables)
    )
    
    # Store in the list
    temp_df <- data.frame(
      cbind(Seq_Num = i, 
            size = rfe_fit$optsize,
            rfe_fit$results[best_index , ] %>% select(-Variables), 
            Variables = rfe_fit$optVariables)
    )
    
    best_subsets_df <- rbind(best_subsets_df, temp_df)
    
    rm(i, temp_df, best_index, rfe_fit)
    
  }
  
  
  # Find the best overall subset size
  best_overall_size <- best_subsets_df %>% 
    select(-Variables) %>% 
    distinct(.keep_all = T) %>% 
    dplyr::rename(Variables = size)
  
  (n_size <- pickSizeBest(best_overall_size, metric = "AUC", maximize = T))
  
  # Find the Seq_Num of the maximum AUC with the n_size
  best_overall_size %>% 
    filter(Variables == n_size) -> temp
  
  # In case there are same number of features with different AUC, select the highest AUC
  temp <- temp[which.max(temp$AUC), "Seq_Num"]
  
  # Get the predictors corresponding to the best overall size
  best_variables <- best_subsets_df %>% 
    filter(Seq_Num == temp)
  rm(temp)
  
  # Print the final best subset of variables
  print(best_variables$Variables)
  
  
  # Compute how many times the best variables were selected during resampling
  temp <- best_subsets_df %>% 
    dplyr::group_by(Variables) %>% 
    dplyr::summarise(Resampled = n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::left_join(best_subsets_df) %>% 
    dplyr::mutate(In_Final_Model = ifelse(Variables %in% best_variables$Variables, "Yes", "No"), 
                  Selected_Final = ifelse(Seq_Num %in% best_variables$Seq_Num, "Yes", "No"))
  
  best_subsets_df <- temp; rm(temp)
  
  toplot <- best_subsets_df %>% 
    filter(Selected_Final == "Yes") %>% 
    mutate(Percent_Resampled = (Resampled/max(best_subsets_df$Seq_Num))*100)
  
  
  # Plot the distribution of the best predictors inside the resampling process
  toplot <- toplot[order(toplot$Percent_Resampled, decreasing = F),]
  barplot(height = toplot$Percent_Resampled, 
          names.arg = toplot$Variables, 
          col = viridis::viridis(length(toplot$Percent_Resampled), option = "mako"), 
          space = 1, 
          horiz = T, 
          border = T,
          las = 1, 
          cex.names = .57, 
          xlim = c(0, 100),
          xlab = "Number of resampling (per 100)")
  
  abline(v = 50, col = "red", lwd = 2)
  
  
  # Save the plot
  dev.copy(png, 
           file = file.path(save_dir, paste0(running_name, "_", "Filtered", "_Predictors_in_the_final_model.png")), 
           units = "cm", 
           res = 300, 
           height = 25, 
           width = 12)
  
  dev.off()
  
  # Export the results
  rio::export(best_subsets_df, file.path(save_dir, paste0(running_name, "_", "Filtered", "_Resamples_Predictors.xlsx")))
  rio::export(best_variables, file.path(save_dir, paste0(running_name, "_", "Filtered", "_Best_Predictors.xlsx")))
  
  
  
  # Create a final model with the best variables  ------------------------
  
  # Select the best predictors only, from the test data
  training_data <- training_data %>% 
    select("POD24", best_variables$Variables)
  
  # Define the recipe for preprocessing
  recipe_obj_best <- recipe(POD24 ~ ., data = training_data) %>%
    step_center(all_predictors()) %>%
    step_scale(all_predictors()) %>%
    step_zv(all_predictors()) %>% 
    step_pca(threshold = .95) 
  
  
  # Define training control
  train_control <- trainControl(method = "LOOCV", 
                                allowParallel = FALSE)  
  
  # Train the Random Forest model
  set.seed(32024)
  BestPredModel <- train(recipe_obj_best,
                         data = training_data,
                         method = "gbm",
                         trControl = train_control)
  
  # Print the model results
  print(BestPredModel)
  
  # Save the final model
  Final_Model_list_Filtered[[paste0(running_name, "_Final_model_for_prediction")]] <- BestPredModel
  
  
  # Get the variable importance in the final model
  varimportant <- varImp(BestPredModel)
  varimportant <- varimportant$importance
  varimportant$Variables <- row.names(varimportant)
  
  # Plot the importnace of th 20th most important variables
  varimportant <- varimportant[order(varimportant$Overall, decreasing = F),]
  varimportant <- varimportant[1:20 , ]
  
  barplot(height = varimportant$Overall, 
          names.arg = varimportant$Variables, 
          col = viridis::viridis(length(varimportant$Overall)), 
          space = 1, 
          horiz = T, 
          border = T,
          las = 1, 
          cex.names = .57, 
          #xlim = c(0, 100),
          xlab = "Variable importance")
  
  
  # Save the plot
  dev.copy(png, 
           file = file.path(save_dir, paste0(running_name, "_", "Filtered", "_Variable_Importance_in_the_final_model.png")), 
           units = "cm", 
           res = 300, 
           height = 25, 
           width = 12)
  
  dev.off()
  
  
  
  # Get performance on the test data using the final model ------------------------------
  
  # Select the best predictors only, from the test data
  testing_data <- testing_data %>% 
    select("POD24", best_variables$Variables)
  
  
  # Predict on the test data
  predictions <- predict(BestPredModel, newdata = testing_data)
  
  
  
  # Import my functions
  source("./Scripts/_00_List_of_functions.R")
  
  
  
  # Evaluate the accuracy
  
  (confmatTest <- confusionMatrix(data = relevel(predictions, ref = 'Yes'),
                                  reference = relevel(testing_data$POD24, ref = 'Yes')))
  
  (Accuracy_test <- paste0(sprintf("%.1f", confmatTest$overall[["Accuracy"]]*100), " (", 
                           sprintf("%.1f", confmatTest$overall[["AccuracyLower"]]*100), "-", 
                           sprintf("%.1f", confmatTest$overall[["AccuracyUpper"]]*100), ")"))
  
  # Get the AUC
  # Use my own defined function
  modAuc_test <- GetAuc(Actual_class = testing_data$POD24, Pred_class = predictions)
  modAuc_test
  
  ## Export model performances 
  
  perfexport_test <- ModPerfExport(ConfMatrix = confmatTest,
                                   Accuracy = Accuracy_test, 
                                   Model_AUC = modAuc_test)
  perfexport_test
  
  
  # Compute AUC metrics for the model
  
  
  # Create a ROC curve
  pred_numeric <- as.numeric(predictions)
  
  roc_curve <- roc(testing_data$POD24, pred_numeric)
  
  # Calculate AUC
  auc_value <- auc(roc_curve)
  cat("AUC:", auc_value*100, "\n")
  modAuc_test
  
  # Plot ROC curves for cv_model1 and cv_model3
  
  plot(x = (1-roc_curve$specificities), 
       y = roc_curve$sensitivities,
       col = colors()[32], 
       type = "l",
       lwd = 2,
       # main = "Random forest - RFE", 
       xlab = "1-Specificity", 
       ylab = "Sensitivity")
  abline(a = 0, b = 1, col = "gray", lwd = 2)
  text(0.7, 0.35, paste0("AUC \n", modAuc_test))
  
  # Save the plot
  dev.copy(png, 
           file = file.path(save_dir, paste0(running_name, "_", "Filtered", "_ROC_Curve_on_test.png")), 
           units = "cm", 
           res = 300, 
           height = 25, 
           width = 12)
  
  dev.off()
  
  # - Export all results --------------
  
  rio::export(perfexport_test, file = file.path(save_dir, paste0(running_name, "_", "Filtered", 
                                                                 "_perfomance_on_test.xlsx")))
  

  rm(list = (setdiff(ls(), c("All_results_filtered_LOOCV", "trainlist", "testlist", "save_dir", "Final_Model_list_Filtered", "start_time"))))
  
  
}

# Save the final models for predictions
save(Final_Model_list_Filtered, file = file.path(save_dir, "Final_Model_list_Filtered.Rdata"))


# Stop parallel backend
# stopCluster(cl)


# Record the ending time
end_time <- Sys.time()

# Calculate and print the runtime
runtime <- end_time - start_time
print(paste("Start time:", start_time))
print(paste("End time:", end_time))
print(paste("Total runtime:", runtime))




