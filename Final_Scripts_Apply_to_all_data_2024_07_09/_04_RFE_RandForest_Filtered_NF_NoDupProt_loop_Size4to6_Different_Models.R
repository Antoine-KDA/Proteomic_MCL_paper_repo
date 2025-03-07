





# Author: Antoine KD ABALO
# Date of creation: 08 April 2024
# Last modification date: 


# Objective: Run a random forest algorithm with a recursive feature elimination (aka backward selection) analysis



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
# sink(file = file.path(save_dir, "Simple_RFrfe.Rout"), append = FALSE, type = "output")


# Load the dataset
load("./Output files/Prot_train_nodup.Rdata")

Prot_train_nodup <- list(All=Prot_train_nodup)

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


# Outer control for RFE; same resamples as before
rfe_ctrl <- rfeControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  functions = rfe_funcs,
  returnResamp = "all",
  verbose = TRUE
)

# RFE with Random Forest Ranking 
imp_ctrl <- rfe_ctrl
imp_ctrl$functions$rank <- caret::rfFuncs$rank








# Create a loop over different sizes ------------------
# This loop runs the same model
# but only changing the maximum number of variables to keep into the model (nsize). 



for(s in 4:6){
  
  nsize = s
  
  
  
  # Run the models on all the data using a for loop ------------------------------------------
  
  
  results_list <- list()
  
  for(i in names(Prot_train_nodup)){
    traindata = Prot_train_nodup[[i]]
    sizes = floor(2:(ncol(traindata)/nsize))
    # Filtered models
    
    ## define the recipe object
    
    i_filtered_recipe <-
      recipe(POD24 ~ ., data = traindata) %>%
      step_corr(all_predictors(), threshold = 0.75) %>%
      step_center(all_predictors()) %>%
      step_scale(all_predictors()) %>%
      step_zv(all_predictors())
    
    
    ## Train the index data
    set.seed(32024)
    i_filtered <- rfe(
      i_filtered_recipe,
      data = traindata,
      sizes = sizes,
      rfeControl = imp_ctrl,
      # metric = "ROC",
      ntree = 5000
    )
    
    
    ## export the results
    
    results_list[[paste0(i, "_filtered_recipe")]] <- i_filtered_recipe
    results_list[[paste0(i, "_filtered")]] <- i_filtered
    rm(i_filtered, i_filtered_recipe)
    
    
    
    # Non-Filtered models
    
    ## define the recipe object
    
    i_recipe <-
      recipe(POD24 ~ ., data = traindata) %>%
      step_center(all_predictors()) %>%
      step_scale(all_predictors()) %>%
      step_zv(all_predictors())
    
    
    ## Train the index data
    set.seed(32024)
    i_nf <- rfe(
      i_recipe,
      data = traindata,
      sizes = sizes,
      rfeControl = imp_ctrl,
      # metric = "ROC",
      ntree = 5000
    )
    
    ## export the results
    
    results_list[[paste0(i, "_recipe_NF")]] <- i_recipe
    results_list[[paste0(i, "_NF")]] <- i_nf
    rm(i_nf, i_recipe, traindata, i)
    
    
    
  }
  
  # Save the model
  RFrfe_results_list <- results_list
  save(RFrfe_results_list, rfe_ctrl, Prot_train_nodup, file = file.path(save_dir, paste0("RFrfe_NoDupProt_Filtered_NF_Size", nsize, "_station.Rdata")))
  
  # Transform - and Get the best tuning parameters in all models ----------------
  
  # load(file = file.path(save_dir, "RFrfe_NoDupProt_Filtered_NF_station.Rdata"))
  # results_list <- RFrfe_results_list; rm(RFrfe_results_list, rfe_ctrl, Prot_train_nodup)
  
  
  
  # Plot the performances over variable subsets in all models ---------------
  
  
  # Plot the performances over variable subsets in all models ---------------
  
  
  # Combine all data and create one data frame
  
  namtoplot <- grep("recipe", names(results_list), value = T, ignore.case = T, invert = T)
  namtoplot
  
  dattoplot <- NULL
  for(n in namtoplot){
    # ntplot = grep("recipe", names(n), value = T, ignore.case = T, invert = T)
    dta = results_list[[n]]
    dta <- 
      dta %>% 
      pluck("results") %>% 
      mutate(Predictors = gsub("_NF|_filtered", "", n)) %>% 
      dplyr::select(ROC, Variables, Predictors, Num_Resamples) %>% 
      mutate(Model = ifelse(length(grep("_filtered", n)),
                            "Correlation filter", "Unfiltered"))
    dattoplot = dattoplot %>% 
      bind_rows(dta)
    
    
  }
  
  dattoplot %<>% 
    mutate(Model = as.factor(Model), 
           Predictors = as.factor(Predictors))
  
  # Plot
  head(dattoplot)
  names(dattoplot)
  
  ggplot(dattoplot, aes(x = Variables, y = ROC, col = Model)) +
    geom_point(size = 0.75) + 
    geom_line() + 
    geom_hline(yintercept = 0.8, lty = "dashed", col = "red")+
    facet_grid("Predictors")  + 
    scale_color_manual(values = c("#6A3D9A", "#CAB2D6"))+
    theme(legend.position = "bottom", 
          text = element_text(size = 14), 
          title = element_text(size = 16, face = "bold")) -> VarPlot
  
  VarPlot
  
  ggsave(VarPlot, 
         file = file.path(save_dir, paste0("RFrfe_NoDupProt_Filtered_NF_Size", nsize, "_Variables.png")), 
         units = "cm", height = 20, width = 15)
  rm(VarPlot)
  
  
  
  
  # Subtract the best predictors from the best model --------
  dattoplot[which.max(dattoplot$ROC), ]
  
  # Max Roc of each model
  tapply(dattoplot[, "ROC"], dattoplot[, c("Predictors", "Model")], max)
  MaxRocTrain <- data.frame(tapply(dattoplot[, "ROC"], dattoplot[, c("Predictors", "Model")], max))
  
  table_grob <- gridExtra::tableGrob(tapply(dattoplot[, "ROC"], dattoplot[, c("Predictors", "Model")], FUN = function(x)round(max(x), 2)))
  
  png(file = file.path(save_dir, paste0("RFrfe_NoDupProt_Filtered_NF_Size", nsize, "_MaxRocTrain.png")), units = "cm", 
      height = 5, width = 10, res = 300)
  
  plot(table_grob)
  
  dev.off()
  
  
  # take the best model automatically 
  dattoplot %>% 
    filter(Predictors == "All") -> dattoplotbest
  
  (varpred <- dattoplotbest[which.max(dattoplotbest$ROC), ]$Predictors)
  (varmod <- dattoplotbest[which.max(dattoplotbest$ROC), ]$Model)
  
  (varbest <- ifelse(grep("Correlation", varmod, value = T) == "Correlation filter", 
                     paste0(varpred, "_filtered"), 
                     paste0(varpred, "_NF")))
  
  bestmodel <- results_list[[varbest]]; rm(varpred, varmod, varbest, dattoplotbest)
  
  bestpredictors <-
    bestmodel %>% 
    pluck("variables") %>% 
    filter(Variables == bestmodel$optsize) %>% 
    group_by(var) %>% 
    count() %>% 
    arrange(desc(n)) %>% 
    mutate(final = ifelse(var %in% bestmodel$optVariables, "Yes", "No")) %>% 
    ungroup()
  
  # Get the accuracy performance of the best model
  bestmodel$fit$confusion
  confusionMatrix.train(bestmodel)
  
  cfmatBestMod <- confusionMatrix.train(bestmodel)
  names(cfmatBestMod)
  cfmatBestMod$table
  cfmatBestMod$norm
  cfmatBestMod$B  
  cfmatBestMod$text  
  cfmatBestMod
  
  
  
  PerfBestModel_on_train <- Get_Perf_RFE_object(RFE_Model = bestmodel)
  PerfBestModel_on_train <- data.frame(PerfBestModel_on_train)
  
  rio::export(PerfBestModel_on_train, file = file.path(save_dir, paste0("RFrfe_NoDupProt_Filtered_NF_Size", nsize, "_PerfBestModel_on_train.xlsx")))
  rm(PerfBestModel_on_train)
  
  
  
  # Get the performances on the test data ---------------------
  # Load the dataset
  load("./Output files/Prot_test_nodup.Rdata")
  Prot_test_nodup <- list(All = Prot_test_nodup)
  
  ## Predict class on the test data using the final model 
  pred_class_test <- predict(bestmodel, Prot_test_nodup$All)
  
  ## Create confusion matrix 
  
  (confmatTest <- confusionMatrix(data = relevel(pred_class_test$pred, ref = 'Yes'),
                                  reference = relevel(Prot_test_nodup$All$POD24, ref = 'Yes')))
  
  
  (Accuracy_test <- paste0(sprintf("%.1f", confmatTest$overall[["Accuracy"]]*100), " (", 
                           sprintf("%.1f", confmatTest$overall[["AccuracyLower"]]*100), "-", 
                           sprintf("%.1f", confmatTest$overall[["AccuracyUpper"]]*100), ")"))
  
  # Get the AUC
  # Use my own defined function
  modAuc_test <- GetAuc(Actual_class = Prot_test_nodup$All$POD24, Pred_class = pred_class_test$pred)
  modAuc_test
  
  ## Export model performances 
  
  perfexport_test <- ModPerfExport(ConfMatrix = confmatTest,
                                   Accuracy = Accuracy_test, 
                                   Model_AUC = modAuc_test)
  perfexport_test
  
  
  
  
  # Compute AUC metrics for the model
  
  
  # Create a ROC curve
  pred_numeric <- as.numeric(pred_class_test$pred)
  
  roc_curve <- roc(Prot_test_nodup$All$POD24, pred_numeric)
  
  # Plot the ROC curve
  plot(roc_curve, main = "ROC Curve", col = "blue", lwd = 2)
  
  
  # Calculate AUC
  auc_value <- auc(roc_curve)
  cat("AUC:", auc_value*100, "\n")
  modAuc_test
  
  # Plot ROC curves for cv_model1 and cv_model3
  png(file = file.path(save_dir, paste0("RFrfe_NoDupProt_Filtered_NF_Size", nsize, "_ROC_Curve_on_test.png")), units = "cm", 
      height = 10, width = 10, res = 300)
  plot(x = (1-roc_curve$specificities), 
       y = roc_curve$sensitivities,
       col = colors()[32], 
       type = "l",
       lwd = 2,
       main = "Random forest - RFE (No dup prot)", 
       xlab = "1-Specificity", 
       ylab = "Sensitivity")
  abline(a = 0, b = 1, col = "gray", lwd = 2)
  text(0.7, 0.35, paste0("AUC \n", modAuc_test))
  
  
  dev.off()
  
  
  # - Export all results --------------
  
  rio::export(perfexport_test, file = file.path(save_dir, paste0("RFrfe_NoDupProt_Filtered_NF_Size", nsize, "_perfomance_on_test.xlsx")))
  rio::export(bestpredictors, file = file.path(save_dir, paste0("RFrfe_NoDupProt_Filtered_NF_Size", nsize, "_bestpredictors.xlsx")))
  
  # Save the model
  RFrfe_results_list <- results_list
  save(RFrfe_results_list, rfe_ctrl, Prot_train_nodup, file = file.path(save_dir, paste0("RFrfe_NoDupProt_Filtered_NF_Size", nsize, "_station.Rdata")))
  
  
  
  
  # Remove all internal objects
  rm(Accuracy_test, auc_value, modAuc_test, n, namtoplot, pred_numeric, 
     table_grob, roc_curve, RFrfe_results_list, results_list, 
     pred_class_test, MaxRocTrain, dta, dattoplot, confmatTest, cfmatBestMod,
     bestpredictors, bestmodel, perfexport_test, Prot_test_nodup)
  
  
  
}