
#--- Function to describe all categorical variables in a dataset --------------

desc.factor <- function(data, # the dataframe to describe
                        var_to_desc, digits=2, collapse=FALSE){ # a vector of variable to describe

  descript.factor <- NULL
  
  for(i in var_to_desc){
    new = as.factor(as.character(data[, i]))
    t = table(data[, i], exclude = NULL)
    new.desc = cbind(Variable = c(i, rep("", length(names(t))-1)),
                     Category = names(t),
                     Frequency = t,
                     Proportion = sprintf(paste0("%.", digits, "f"), prop.table(t)*100))
    descript.factor = data.frame(rbind(descript.factor, "", new.desc))
    descript.factor
    
  }
  
  if(collapse==TRUE) {
    descript.factor[ , "N (Col %)"] = paste0(descript.factor$Frequency, " (", 
                                             descript.factor$Proportion, "%)")
    descript.factor = descript.factor[, c("Variable", "Category", "N (Col %)")]
    descript.factor
  }

  descript.factor
}





#--- Function to describe all date variables in a dataset ---------------------


desc.date <- function(data){ # the dataframe to describe

  descript.date <- NULL
  for (i in names(data)) {
    if(class(data[, i])=="Date"){
      min = paste(min(data[, i], na.rm = T))
      max = paste(max(data[, i], na.rm = T))
      median = paste(median(data[, i], na.rm = T))
      na = length(which(is.na(data[, i])))
      available = length(which(!is.na(data[, i])))
      new.date = cbind(Variable = i,
                       Minimum = min,
                       Maximum = max,
                       Median = median,
                       Not_Available = na,
                       Available = available)
      descript.date = data.frame(rbind(descript.date, new.date))
      descript.date
    }
    
  }
  descript.date
}




#create a function to check for installed packages and install them if they are not installed
install <- function(packages){
  new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new.packages)) 
    install.packages(new.packages, dependencies = TRUE)
  sapply(packages, require, character.only = TRUE)
}



# Fonction to wrap up the table() and prop.table() functions

# Frequency table
freqtabl <- function(x, dig=2) cbind(N=table(x, exclude = NULL), 
                                     Prop=sprintf(paste0("%.", dig, "f"), prop.table(table(x, exclude = NULL))*100))

freqtabl_na <- function(x, dig=2) cbind(N=table(x), 
                                     Prop=sprintf(paste0("%.", dig, "f"), prop.table(table(x))*100))





# Theme for ggplot objects 
# Customize the plot
custom_theme2 <- function() {
  require(survminer)
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5, face = "bold"),
      axis.line = element_line(# size = .75 deprecated 
        linewidth = .75
      ), 
      panel.grid.major.x = element_line(linewidth =.1),
      panel.grid.minor.x = element_line(linewidth =.05) #,
      # panel.grid.major.y = element_line(linewidth =.5), 
      # panel.grid.minor.y = element_line(linewidth =.1)
    )
  
}




# Create a function to compute the AUC from the models
# this function was created and tested in the 'Simple_PLSDA.R' script

# auc_ci <- function(Roc_Curv){
#   auc_value <- auc(roc_curve)
#   ci <- ci.auc(roc_curve)
#   auc_ci_res <- paste0(sprintf("%.1f", auc_value*100), " (",
#                        sprintf("%.1f", ci[1]*100), "-",
#                        sprintf("%.1f", ci[3]*100), ")")
#   auc_ci_res
# }


GetAuc <- function(Actual_class, Pred_class) {
  if (!is.factor(Pred_class) || !is.factor(Actual_class)) {
    stop("Both Actual_class and Pred_class should be factors.")
  }
  
  if (!all(levels(Pred_class) %in% c("Yes", "No"))) {
    stop("Levels of Pred_class should be 'Yes' and 'No'.")
  }
  
  if (!all(levels(Actual_class) %in% c("Yes", "No"))) {
    stop("Levels of Actual_class should be 'Yes' and 'No'.")
  }
  
  require(pROC)
  
  # Compute the ROC curve
  roc_curve <- roc(as.numeric(Actual_class == 'Yes'), as.numeric(Pred_class == 'Yes'))
  
  # Get the AUC
  (auc_value <- auc(roc_curve))
  (ci <- ci.auc(roc_curve))
  
  auc_ci_res <- paste0(sprintf("%.1f", auc_value*100), " (",
                       sprintf("%.1f", ci[1]*100), "-",
                       sprintf("%.1f", ci[3]*100), ")")
  auc_ci_res
  
}

# # Example usage
# Actual_class <- Prot_train$POD24
# Pred_class <- pred_class
# 
# result <- GetAuc(Actual_class, Pred_class)
# cat("AUC:", result, "\n")



# Function to export model performances from a confusion matrix object generated by 
# confusionMatrix function from caret package


ModPerfExport <- function(ConfMatrix, Accuracy, Model_AUC){
  modprfm <- as.data.frame(ConfMatrix$byClass) 
  modprfm <- data.frame(cbind(Performances = rownames(modprfm), 
                              Values_in_Percent = sprintf("%.1f", modprfm[[1]]*100)))
  modprfm
  
  modprfm[nrow(modprfm)+1, ] <- c("Accuracy (95% CI)", Accuracy)
  modprfm[nrow(modprfm)+1, ] <- c("AUC (95% CI)", Model_AUC)
  
  modprfm
  
}









modperformancies <- function(Model, Prediction_data){ 
  
  
  ## Predict class on the training data using the final model 
  pred_class_train <- predict(Model, Prediction_data)
  
  ## Create confusion matrix 
  (confusionMatrix(data = relevel(pred_class_train, ref = 'Yes'),
                   reference = relevel(Prediction_data$POD24, ref = 'Yes')) -> confmatTrain)
  
  
  (Accuracy_train <- paste0(sprintf("%.1f", confmatTrain$overall[["Accuracy"]]*100), " (", 
                            sprintf("%.1f", confmatTrain$overall[["AccuracyLower"]]*100), "-", 
                            sprintf("%.1f", confmatTrain$overall[["AccuracyUpper"]]*100), ")"))
  
  # Get the AUC
  # Use my own defined function
  modAuc_train <- GetAuc(Actual_class = Prediction_data$POD24, Pred_class = pred_class_train)
  modAuc_train
  
  ## Export model performances 
  
  perfexport_train <- ModPerfExport(ConfMatrix = confmatTrain,
                                    Accuracy = Accuracy_train, 
                                    Model_AUC = modAuc_train)
  perfexport_train
  
}





# Function to get performances of an RFE model -------------------


Get_Perf_RFE_object <- function(RFE_Model, digt = 1){
  
  # Extract the confusion matrix
  conf_matrix <- confusionMatrix.train(RFE_Model)
  
  # Extract the confusion matrix data
  conf_matrix_data <- conf_matrix$table
  
  # Calculate metrics
  accuracy <- sum(diag(conf_matrix_data)) / sum(conf_matrix_data)
  # Calculate the standard error of the accuracy
  se_accuracy <- sqrt(accuracy * (1 - accuracy) / sum(conf_matrix_data))
  
  # Calculate the margin of error for the 95% CI
  z <- qnorm(0.975)  # 95% confidence level
  margin_error <- z * se_accuracy
  
  # Calculate the lower and upper bounds of the confidence interval
  ci_lower <- accuracy - margin_error
  ci_upper <- accuracy + margin_error
  
  
  sensitivity <- conf_matrix_data[2, 2] / sum(conf_matrix_data[2, ])
  specificity <- conf_matrix_data[1, 1] / sum(conf_matrix_data[1, ])
  precision <- conf_matrix_data[2, 2] / sum(conf_matrix_data[, 2])
  recall <- sensitivity
  f1_score <- 2 * (precision * recall) / (precision + recall)
  
  
  Perfresult <- data.frame(cbind(Performances = c("Sensitivity", 
                                                  "Specificity", 
                                                  "Precision", 
                                                  "Recall", 
                                                  "F1_score", 
                                                  "Accuracy", 
                                                  "Accuracy (95% CI)"), 
                                 Values_in_Percent = c(sprintf(paste0("%.", digt, "f"), sensitivity*100),
                                                       sprintf(paste0("%.", digt, "f"), specificity*100),
                                                       sprintf(paste0("%.", digt, "f"), precision*100),
                                                       sprintf(paste0("%.", digt, "f"), recall*100),
                                                       sprintf(paste0("%.", digt, "f"), f1_score*100),
                                                       sprintf(paste0("%.", digt, "f"), accuracy*100),
                                                       paste0(sprintf(paste0("%.", digt, "f"), accuracy*100), " (",
                                                              sprintf(paste0("%.", digt, "f"), ci_lower*100), "-", 
                                                              sprintf(paste0("%.", digt, "f"), ci_upper*100), ")"
                                                       ))))
  Perfresult
  
  
}