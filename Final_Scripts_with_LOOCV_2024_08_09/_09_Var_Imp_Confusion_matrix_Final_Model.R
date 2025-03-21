



# Author: Antoine KD ABALO
# Date of creation: 20 August 2024
# Last modification date: 


# Objective: Investigating the distribution of the best feature set in the final model





# Import my functions
source("./Scripts/_00_List_of_functions.R")
rm(list=setdiff(ls(), "install"))


# Install the packages
install(c("magrittr", "dplyr", "ggplot2", "gbm", "glmnet", "caret"))


# Load the final models' list
load("./Results/Results_LOOCV_Parallelized_KI_Computer_2024_08_11/Non_Filtered_Only_2024_08_11/Final_Model_list_NonFiltered.Rdata")

# Take the corresponding best model
#' Note: form the summary excel sheet, this is *All Med For Dup NOn Filtered*

selected_mod <- Final_Model_list_NonFiltered$`All Median for dup_Final_model_for_prediction`

selected_mod

# Get the variable importance in the final model
varimportant <- varImp(selected_mod)
varimportant <- varimportant$importance
varimportant$Variables <- gsub(pattern = "__Neur|__Infl|__Onco|__Card", "", row.names(varimportant))

# Plot the importnace of th 20th most important variables
varimportant <- varimportant %>% arrange(desc(Overall))
#varimportant <- varimportant[1:20 , ]

# using ggplot2 

names(varimportant)

p <- ggplot(varimportant, aes(x = Overall, y = reorder(Variables, Overall))) +
  geom_bar(stat = "identity", fill = "darkred") +
  scale_x_continuous(limits = c(0, max(varimportant$Overall) + 0.5), name = "Variable importance (%)") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_blank(),  # No y-axis title
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(size = 10, colour = "black"),
    panel.grid.major.y = element_blank()  # Remove grid lines for y-axis
  )

p

# Save the plot
ggsave(plot = p,  
       file = file.path("./Results/Results_LOOCV_Parallelized_KI_Computer_2024_08_11/Variable_importance_final_model.png"), 
       units = "cm", 
       bg = "white",
       height = 25, 
       width = 15)

dev.off()


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

rm(p, varimportant)



# Plot the confusion matrix
## load the test data
load("./Output files/testlist_12July2024_Ant.Rdata")

# Get the final set of predictors
names(testlist_12July2024_Ant)

testing_data <- testlist_12July2024_Ant[["All: Median for dup"]]; rm(testlist_12July2024_Ant)

testing_data <- testing_data %>% 
  select("POD24", rownames(varImp(selected_mod)$importance))


# Predict on the test data
predictions <- predict(selected_mod, newdata = testing_data)


# Plot the confusion matrice

(confmatTest <- confusionMatrix(data = relevel(predictions, ref = 'Yes'),
                                reference = relevel(testing_data$POD24, ref = 'Yes')))


confmatTest$table

cfm <- data.frame(confmatTest$table)
names(cfm)
names(cfm) <- c("Prediction","Target","N")
names(cfm)
cfm

install("cvms")

plot_confusion_matrix(cfm)



dev.copy(png, 
         file = file.path("./Results/Results_LOOCV_Parallelized_KI_Computer_2024_08_11/Confusion_Matrix_final_model.png"), 
         units = "cm", 
         res = 300,
         bg = "white",
         height = 12, 
         width = 12)

dev.off()


# Plot the roc curve 

install("pROC")

roc_curve <- roc(testing_data$POD24, as.numeric(predictions))

# Extract the AUC (Area Under the Curve)
auc_value <- auc(roc_curve)



# Create a data frame from the ROC curve
roc_data <- data.frame(
  Specificity = rev(roc_curve$specificities),  # Reverse for correct order
  Sensitivity = rev(roc_curve$sensitivities)
)

# Plot using ggplot2
(ggplot(roc_data, aes(x = 1 - Specificity, y = Sensitivity)) +
  geom_line(color = "red4", size = 1) +
  geom_abline(linetype = "dashed", color = "gray") +  # Diagonal line for random guessing
  labs(
    #title = paste("ROC Curve (AUC =", round(auc_value, 2), ")"),
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  ) +
  theme_light() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 13, colour = "black"),
    axis.title.y = element_text(size = 13, colour = "black"),
    axis.text = element_text(size = 12, colour = "black")
  ) +
  geom_text(
    aes(x = 0.6, y = 0.5, label = paste("AUC =", sprintf("%.2f (%.2f - %.2f)", auc_value, 0.583, 1))), 
    size = 5, 
    color = "black"
  ) -> Roc_Curve_plot)




ggsave(plot = Roc_Curve_plot, 
         file = file.path("./Results/Results_LOOCV_Parallelized_KI_Computer_2024_08_11/Roc_Curve_final_model.png"), 
         units = "cm", 
         bg = "white",
         height = 12, 
         width = 12)

dev.off()
