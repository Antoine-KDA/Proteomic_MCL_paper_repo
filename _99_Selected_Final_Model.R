





# Author: Antoine KD ABALO
# Date of creation: 10 April 2024
# Last modification date: 


# Objective: Run a description of the final selected best model



rm(list=ls())
# Import my functions
source("./Scripts/_00_List_of_functions.R")

# Install the packages
install(c("magrittr", "dplyr", "vip", "caret", "MLmetrics", "recipes", "ggplot2"))



# Create a directory to save the output -------
save_dir <- paste0("./Results/Description_Selected_Final_Model_", gsub(" ", "_", format(Sys.Date(), "%Y %m %d")))

# Ensure the directory exists, create if not
if (!file.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

# Open the sink to redirect output to the log file
# sink(file = file.path(save_dir, "Simple_treeBagRFE.Rout"), append = FALSE, type = "output")


# Load the finally selected best model which is treeBagRFE_NoDupProt_Filtered_NF_Size5
load("./Results/Results_Final_Scripts_2024_04_10/Filtered_NF_2024_04_10/treeBagRFE_NoDupProt_Filtered_NF_Size5_station.Rdata")
# The best model is the filtered one
bestfinalmodel <- treeBagRFE_NoDupProt_results_list$All_filtered



# Check the performances of this best model
Get_Perf_RFE_object(bestfinalmodel)

# Barplot of the best predictors -----------------------

# Have a look at the best predictors of the best model
# Import the predictors
bestpredictors <- readxl::read_excel("./Results/Results_Final_Scripts_2024_04_10/Filtered_NF_2024_04_10/treeBagRFE_NoDupProt_Filtered_NF_Size5_bestpredictors.xlsx")

barplot(height = bestpredictors$n, names.arg = bestpredictors$var)


mycol <- colours()[3:length(bestpredictors$var)+3]

barplot(height = bestpredictors$n, 
        names.arg = bestpredictors$var, 
        col = mycol, 
        space = 1, 
        horiz = T, 
        border = TRUE,
        las = 1, cex.names = .5)

# Plot only predictors in the final model

bartoplot <- bestpredictors[bestpredictors$final == "Yes", ]

nrow(bestfinalmodel$recipe$last_term_info)
# [1] 330
nrow(bestpredictors)
# [1] 180
nrow(bartoplot)
# [1] 56
bestfinalmodel$optsize
# [1] 56

oldpar <- par()
par(mar=c(5, 5, .5, 1))

bartoplot <- bartoplot[order(bartoplot$n, decreasing = F),]
barplot(height = bartoplot$n, 
        names.arg = bartoplot$var, 
        col = viridis::viridis(length(bartoplot$n), option = "mako"), 
        space = 1, 
        horiz = T, 
        border = T,
        las = 1, 
        cex.names = .57, 
        xlab = "Number of resampling (per 100)")


# Save the plot
dev.copy(png, 
         file = file.path(save_dir, "Predictors_in_the_final_model.png"), 
         units = "cm", 
         res = 300, 
         height = 25, 
         width = 12)

dev.off()



# Variable importance plot --------------------

importance <-  data.frame(varImp(bestfinalmodel))
importance <- data.frame(cbind(Predictors = row.names(importance), 
                               importance))
# Transform the importance into percentage
importance$Importance_percent <- (importance$Overall/sum(importance$Overall))*100
sum(importance$Importance_percent)

# Limit to predictors in the final model only
importance <- importance[importance$Predictors %in% bartoplot$var, ]
importance <- importance[order(importance$Importance_percent, decreasing = F), ]


# Plot the barplot
barplot(height = importance$Importance_percent, 
        names.arg = importance$Predictors, 
        col = viridis::viridis(length(importance$Predictors)), 
        space = 1, 
        horiz = T, 
        las = 1, 
        cex.names = .57, 
        xlab = "Importance", 
        xpd = TRUE, 
        xlim = c(0, max(importance$Importance_percent)+0.5))

# Save the plot
dev.copy(png, 
         file = file.path(save_dir, "Importance_Predictors_in_the_final_model.png"), 
         units = "cm", 
         res = 300, 
         height = 25, 
         width = 12)

dev.off()





barplot(height = importance$Overall, 
        names.arg = importance$Predictors, 
        col = viridis::viridis(length(importance$Predictors)), 
        space = 1, 
        horiz = T, 
        las = 1, 
        cex.names = .57, 
        xlab = "Importance", 
        xpd = TRUE, 
        xlim = c(0, max(importance$Overall)+0.5))


par(oldpar)




# Calibration plot -----------------------


# Load the dataset
load("./Output files/Prot_test_nodup.Rdata")
Prot_test_nodup <- list(All = Prot_test_nodup)



## Predict class on the test data using the final model 
pred_class_test <- predict(bestfinalmodel, Prot_test_nodup$All)

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



# The calibration

testProbs <- data.frame(obs = bestfinalmodel$fit$y,
                        rfe = predict(bestfinalmodel, bestfinalmodel$fit$x)[,3])

calibration(obs ~ rfe, data = testProbs, class = "Yes")

calPlotData <- calibration(obs ~ rfe, data = testProbs, class = "Yes")
calPlotData

xyplot(calPlotData, auto.key = list(columns = 2))

# I am not sure these codes on calibration are correct, need to double check


# Create a frequency plot ---------------
# Load the models:
# rm(list=ls())

library(readxl)

# The best model
read_excel("./Results/Results_Final_Scripts_2024_04_10/Filtered_NF_2024_04_10/treeBagRFE_NoDupProt_Filtered_NF_Size5_bestpredictors.xlsx") %>% 
  filter(final == "Yes") %>% 
  mutate(Model = "treeBagRFE_NoDupProt_Filtered_NF_Size5") -> treeBagRFE_NoDupProt_Filtered_NF_Size5

# The other models

read_excel("./Results/Results_Final_Scripts_2024_04_10/Filtered_NF_2024_04_09/RFrfe_Filtered_NF_Size4_bestpredictors.xlsx") %>% 
  filter(final == "Yes") %>% 
  mutate(Model = "RFrfe_Filtered_NF_Size4") %>% 
  filter(var %in% treeBagRFE_NoDupProt_Filtered_NF_Size5$var) -> RFrfe_Filtered_NF_Size4

read_excel("./Results/Results_Final_Scripts_2024_04_10/Filtered_NF_2024_04_09/RFrfe_Filtered_NF_Size5_bestpredictors.xlsx") %>% 
  filter(final == "Yes")%>% 
  mutate(Model = "RFrfe_Filtered_NF_Size5") %>% 
  filter(var %in% treeBagRFE_NoDupProt_Filtered_NF_Size5$var) -> RFrfe_Filtered_NF_Size5

read_excel("./Results/Results_Final_Scripts_2024_04_10/Filtered_NF_2024_04_09/RFrfe_Filtered_NF_Size6_bestpredictors.xlsx") %>% 
  filter(final == "Yes")%>% 
  mutate(Model = "RFrfe_Filtered_NF_Size6") %>% 
  filter(var %in% treeBagRFE_NoDupProt_Filtered_NF_Size5$var) -> RFrfe_Filtered_NF_Size6



read_excel("./Results/Results_Final_Scripts_2024_04_10/Filtered_NF_2024_04_10/RFrfe_NoDupProt_Filtered_NF_Size4_bestpredictors.xlsx") %>% 
  filter(final == "Yes") %>% 
  mutate(Model = "RFrfe_NoDupProt_Filtered_NF_Size4") %>% 
  filter(var %in% treeBagRFE_NoDupProt_Filtered_NF_Size5$var) -> RFrfe_NoDupProt_Filtered_NF_Size4

read_excel("./Results/Results_Final_Scripts_2024_04_10/Filtered_NF_2024_04_10/RFrfe_NoDupProt_Filtered_NF_Size5_bestpredictors.xlsx") %>% 
  filter(final == "Yes") %>% 
  mutate(Model = "RFrfe_NoDupProt_Filtered_NF_Size5") %>% 
  filter(var %in% treeBagRFE_NoDupProt_Filtered_NF_Size5$var) -> RFrfe_NoDupProt_Filtered_NF_Size5

read_excel("./Results/Results_Final_Scripts_2024_04_10/Filtered_NF_2024_04_10/RFrfe_NoDupProt_Filtered_NF_Size6_bestpredictors.xlsx") %>% 
  filter(final == "Yes") %>% 
  mutate(Model = "RFrfe_NoDupProt_Filtered_NF_Size6") %>% 
  filter(var %in% treeBagRFE_NoDupProt_Filtered_NF_Size5$var) -> RFrfe_NoDupProt_Filtered_NF_Size6



read_excel("./Results/Results_Final_Scripts_2024_04_10/Filtered_NF_2024_04_10/treeBagRFE_Filtered_NF_Size4_bestpredictors.xlsx") %>% 
  filter(final == "Yes") %>% 
  mutate(Model = "treeBagRFE_Filtered_NF_Size4") %>% 
  filter(var %in% treeBagRFE_NoDupProt_Filtered_NF_Size5$var) -> treeBagRFE_Filtered_NF_Size4

read_excel("./Results/Results_Final_Scripts_2024_04_10/Filtered_NF_2024_04_10/treeBagRFE_Filtered_NF_Size5_bestpredictors.xlsx") %>% 
  filter(final == "Yes") %>% 
  mutate(Model = "treeBagRFE_Filtered_NF_Size5") %>% 
  filter(var %in% treeBagRFE_NoDupProt_Filtered_NF_Size5$var) -> treeBagRFE_Filtered_NF_Size5

read_excel("./Results/Results_Final_Scripts_2024_04_10/Filtered_NF_2024_04_10/treeBagRFE_Filtered_NF_Size6_bestpredictors.xlsx") %>% 
  filter(final == "Yes") %>% 
  mutate(Model = "treeBagRFE_Filtered_NF_Size6") %>% 
  filter(var %in% treeBagRFE_NoDupProt_Filtered_NF_Size5$var) -> treeBagRFE_Filtered_NF_Size6



read_excel("./Results/Results_Final_Scripts_2024_04_10/Filtered_NF_2024_04_10/treeBagRFE_NoDupProt_Filtered_NF_Size4_bestpredictors.xlsx") %>% 
  filter(final == "Yes") %>% 
  mutate(Model = "treeBagRFE_NoDupProt_Filtered_NF_Size4") %>% 
  filter(var %in% treeBagRFE_NoDupProt_Filtered_NF_Size5$var) -> treeBagRFE_NoDupProt_Filtered_NF_Size4

read_excel("./Results/Results_Final_Scripts_2024_04_10/Filtered_NF_2024_04_10/treeBagRFE_NoDupProt_Filtered_NF_Size6_bestpredictors.xlsx") %>% 
  filter(final == "Yes") %>% 
  mutate(Model = "treeBagRFE_NoDupProt_Filtered_NF_Size6") %>% 
  filter(var %in% treeBagRFE_NoDupProt_Filtered_NF_Size5$var)  -> treeBagRFE_NoDupProt_Filtered_NF_Size6


Varinothermodels <- rbind(treeBagRFE_NoDupProt_Filtered_NF_Size5,
                          treeBagRFE_NoDupProt_Filtered_NF_Size4,
                          treeBagRFE_NoDupProt_Filtered_NF_Size6,
                          treeBagRFE_Filtered_NF_Size4,
                          treeBagRFE_Filtered_NF_Size5,
                          treeBagRFE_Filtered_NF_Size6,
                          RFrfe_NoDupProt_Filtered_NF_Size4,
                          RFrfe_NoDupProt_Filtered_NF_Size5,
                          RFrfe_NoDupProt_Filtered_NF_Size6,
                          RFrfe_Filtered_NF_Size4,
                          RFrfe_Filtered_NF_Size5,
                          RFrfe_Filtered_NF_Size6)

length(unique(Varinothermodels$Model))

data.frame(tapply(Varinothermodels$n, Varinothermodels[, "var"], sum)) -> Freqvar 

Freqvar <- as.data.frame(cbind(Predictors = row.names(Freqvar), 
                 N = Freqvar$tapply.Varinothermodels.n..Varinothermodels....var....sum.))

Freqvar$N <- as.numeric(as.character(Freqvar$N))
Freqvar$Percent_N = (Freqvar$N/(length(unique(Varinothermodels$Model))*100))*100


Freqvar <- Freqvar[order(Freqvar$N, decreasing = F), ]



# Plot the barplot
oldpar <- par()
par(mar=c(5, 5, .5, 1))


barplot(height = Freqvar$Percent_N, 
        names.arg = Freqvar$Predictors, 
        col = viridis::viridis(length(Freqvar$Predictors), option = "rocket"), 
        space = 1, 
        horiz = T, 
        las = 1, 
        cex.names = .57, 
        xlab = "% of resampling", 
        xpd = TRUE, 
        xlim = c(0, max(Freqvar$Percent_N)+0.5))

abline(v=50, lty = "dashed", lwd = 2)


# Save the plot
dev.copy(png, 
         file = file.path(save_dir, "Predictors_Resampled_in_other_models_in_the_final_model.png"), 
         units = "cm", 
         res = 300, 
         height = 25, 
         width = 12)

dev.off()

par(oldpar)

rm(treeBagRFE_NoDupProt_Filtered_NF_Size5,
    treeBagRFE_NoDupProt_Filtered_NF_Size4,
    treeBagRFE_NoDupProt_Filtered_NF_Size6,
    treeBagRFE_Filtered_NF_Size4,
    treeBagRFE_Filtered_NF_Size5,
    treeBagRFE_Filtered_NF_Size6,
    RFrfe_NoDupProt_Filtered_NF_Size4,
    RFrfe_NoDupProt_Filtered_NF_Size5,
    RFrfe_NoDupProt_Filtered_NF_Size6,
    RFrfe_Filtered_NF_Size4,
    RFrfe_Filtered_NF_Size5,
    RFrfe_Filtered_NF_Size6, Freqvar, importance, Varinothermodels)



#  Try a tree bagging on the processed data  --------- 

# Load the finally selected best model which is treeBagRFE_NoDupProt_Filtered_NF_Size5
load("./Results/Results_Final_Scripts_2024_04_10/Filtered_NF_2024_04_10/treeBagRFE_NoDupProt_Filtered_NF_Size5_station.Rdata")
# The best model is the filtered one
bestfinalmodel <- treeBagRFE_NoDupProt_results_list$All_filtered
# Get the variables selected from the final model
bestfinalmodel$bestSubset
bestfinalmodel$optVariables

# Load the training data
load("./Output files/Prot_train_nodup.Rdata")
# Limit the data to the final selected proteins only
Prot_train_nodup_tree <- Prot_train_nodup[, c("POD24", bestfinalmodel$optVariables)]

## Define the control functions : ---------------
# The control function

tr_ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  classProbs = TRUE, 
  summaryFunction = twoClassSummary
)

# The recipe object
treebag_recipe <-
  recipe(POD24 ~ ., data = Prot_train_nodup_tree) %>%
  # step_corr(all_predictors(), threshold = 0.75) %>% # No more needed as the features were processed and selected
  step_center(all_predictors()) %>%
  step_scale(all_predictors()) %>%
  step_zv(all_predictors())

# Fit the model 
set.seed(32024)
treebagFit <- train(POD24 ~ .,
                data = Prot_train_nodup_tree,  
                method = "treebag", 
                trControl = tr_ctrl,
                metric = "ROC")
treebagFit
treebagFit$finalModel

confusionMatrix.train(treebagFit)

rm(Prot_train_nodup_tree, rfe_ctrl, treeBagRFE_NoDupProt_results_list)

# look at the performance on the test data set
# Load the training data
load("./Output files/Prot_test_nodup.Rdata")
# Limit the data to the final selected proteins only
Prot_test_nodup_tree <- Prot_test_nodup[, c("POD24", bestfinalmodel$optVariables)]
rm(Prot_test_nodup)


pred_treebag <- bind_cols(
  predict(treebagFit, newdata = Prot_test_nodup_tree, type = "prob"),
  Predicted = predict(treebagFit, newdata = Prot_test_nodup_tree, type = "raw"),
  Actual = Prot_test_nodup_tree$POD24
)

# The confusion matrix

ConfMat_treebag <- confusionMatrix(data = relevel(pred_treebag$Predicted, ref= "Yes"), 
                                   reference = relevel(pred_treebag$Actual, ref ="Yes"))
ConfMat_treebag


(Accuracy_train <- paste0(sprintf("%.1f", ConfMat_treebag$overall[["Accuracy"]]*100), " (", 
                          sprintf("%.1f", ConfMat_treebag$overall[["AccuracyLower"]]*100), "-", 
                          sprintf("%.1f", ConfMat_treebag$overall[["AccuracyUpper"]]*100), ")"))


# Get the AUC
# Use my own defined function
modAuc_test <- GetAuc(Actual_class = Prot_test_nodup_tree$POD24, Pred_class = pred_treebag$Predicted)
modAuc_test

## Export model performances -------------

perfexport_train <- ModPerfExport(ConfMatrix = ConfMat_treebag,
                                  Accuracy = Accuracy_train, 
                                  Model_AUC = modAuc_test)
perfexport_train


# Similar to what I had previously


mdl_auc <- Metrics::auc(actual = pred_treebag$Actual == "Yes", pred_treebag$Yes)
yardstick::roc_curve(pred_treebag, Actual, No) %>%
  autoplot() +
  labs(
    title = "OJ Bagging ROC Curve",
    subtitle = paste0("AUC = ", round(mdl_auc, 4))
  )

yardstick::gain_curve(pred_treebag, Actual, Yes) %>%
  autoplot() +
  labs(title = "OJ Bagging Gain Curve")

plot(varImp(treebagFit), main="Variable Importance with Bagging")



rm(testProbs, treebag_recipe, tr_ctrl, pred_treebag, pred_class_test, 
   perfexport_test, ConfMat_treebag, confmatTest, calPlotData, bartoplot)


# How





# Random subset size validation ------------------------------------------------


subset_size <- length(bestfinalmodel$optVariables)
iter <- 100
vars <- names(Prot_train_nodup)
vars <- vars[vars != "POD24"]

# The control function
subset_ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  classProbs = TRUE, 
  summaryFunction = twoClassSummary )



for (i in 1:iter) {
  set.seed(24 + 2*i)
  rand_subset <- sample(vars, subset_size)
  
  # Create a subdata
  SubData <- Prot_train_nodup[, c("POD24", rand_subset)]
  
  ## Create the recipe object for the index subset data
  subset_rec <-
    recipe(POD24 ~ ., data = SubData) %>%
    # step_corr(all_predictors(), threshold = 0.75) %>% # No more needed as the features were processed and selected
    step_center(all_predictors()) %>%
    step_scale(all_predictors()) %>%
    step_zv(all_predictors())
  
  ## Train the index subset data
  set.seed(32024 + i)
  subset_model <- train(
    subset_rec,
    data = SubData,  
    method = "treebag", 
    trControl = subset_ctrl,
    metric = "ROC")
  
  
  # Get the performances of the trained index subset data
  
  subset_perf <- getTrainPerf(subset_model)
  # subset_perf$ncomp <- subset_model$bestTune$ncomp
  if (i == 1) {
    size_check <- subset_perf
  } else {
    size_check <- bind_rows(size_check, subset_perf)
  }
  rm(rand_subset, subset_rec, subset_model, subset_perf, SubData)
}

# Save the results first
save(size_check, file = file.path(save_dir, "Simulated_Random_Subset_size_check.xlsx"))



subset_size
size_check
getTrainPerf(treebagFit)
treebagFit$results$ROC
# [1] 0.801875

mean(size_check$TrainROC <= treebagFit$results$ROC)
# [1] 0.94
table(size_check$TrainROC <= treebagFit$results$ROC)
# FALSE  TRUE 
# 6    94 
mean(size_check$TrainROC > treebagFit$results$ROC)
# [1] 0.06
mean(size_check$TrainROC > treebagFit$results$ROC)*100


# Plot
barplot(size_check$TrainROC, ylim = c(0, 1), 
        xlab = paste0("Performance of 56 random predictors from a total of 330 predictors \n
                Simulation over 100 iterations"),
        ylab = "ROC")
abline(h=treebagFit$results$ROC)


# Create the barplot
(ggplot(size_check, aes(x = 1:100, y = TrainROC)) +
    geom_point(fill = "lightblue") +
    ylim(0, 1) +
    labs(x = "Performance of 56 random predictors\nFrom a total of 330 predictors\nSimulation over 100 iterations",
         y = "ROC") +
    geom_hline(yintercept = treebagFit$results$ROC, linetype = "dashed", color = "red", lwd = 1) +
    theme_minimal() -> Simulated_Random_Subset)




custom_theme2 <- function() {
  require(survminer)
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5, face = "bold"),
      axis.line = element_line(# size = .75 deprecated 
        linewidth = .75
      ), 
      panel.grid.major.x = element_line(linewidth =.5),
      panel.grid.minor.x = element_line(linewidth =.2),
      panel.grid.major.y = element_line(linewidth =.5),
      panel.grid.minor.y = element_line(linewidth =.2)
    )
  
}

(ggplot(size_check, aes(x = 1:100, y = TrainROC*100)) +
    geom_point(fill = "lightblue") +
    labs(x = "Performance of 56 random predictors\nFrom a total of 330 predictors\nSimulation over 100 iterations",
         y = "ROC") +
    geom_hline(yintercept = treebagFit$results$ROC*100, linetype = "dashed", color = "red", lwd = 1) +
    scale_y_continuous(limits = c(50, 100))+
    scale_x_continuous(n.breaks = 6)+
    custom_theme2() -> Simulated_Random_Subset)

Simulated_Random_Subset +
  annotate("text", x=50, y=95, 
           label = paste("The selected model, ROC =", round(treebagFit$results$ROC*100, 1), "is >= than 94% of simulated models"), 
           hjust = 0.5, vjust = -0.5, color = "darkblue") -> Simulated_Random_Subset_annotate
  


# Save the plot
ggsave(Simulated_Random_Subset, 
       file = file.path(save_dir, "Simulated_Random_Subset.png"), 
       units = "cm", 
       height = 15, 
       width = 20)

ggsave(Simulated_Random_Subset_annotate, 
       file = file.path(save_dir, "Simulated_Random_Subset_annotate.png"), 
       units = "cm", 
       height = 15, 
       width = 20)

dev.off()

