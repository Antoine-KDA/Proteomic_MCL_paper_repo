

# Author: Antoine KD ABALO
# Date of creation: 20 August 2024
# Last modification date: 


# Objective: Investigating the distribution of the best feature set in the final model





# Import my functions
source("./Scripts/_00_List_of_functions.R")
rm(list=setdiff(ls(), "install"))


# Install the packages
install(c("magrittr", "dplyr", "ggplot2"))

# Load the resampling models
load("./Results/Results_LOOCV_Parallelized_KI_Computer_2024_08_11/Non_Filtered_Only_2024_08_11/RFrfe_ll_results_NonFiltered_LOOCV_list.Rdata")

# Check the names in the list
ls()
names(All_results_NonFiltered_LOOCV)


# Take the corresponding best model
#' Note: form the summary excel sheet, this is *All Med For Dup NOn Filtered*

bestmodel <- All_results_NonFiltered_LOOCV$`All Median for dup_results_parallel_LOOCV`
names(bestmodel)

# Take the final best model
bestfinalmodel <- bestmodel[[45]]
bestfinalmodel

rm(All_results_NonFiltered_LOOCV, bestmodel)


# plot the distribution of variable selection

plot(bestfinalmodel$rfe_fit$results$Variables, bestfinalmodel$rfe_fit$results$ROC,
     xlim = c(0, 150), ylim = c(0.5, 1), xlab = "Number of features", ylab = "ROC")
lines(bestfinalmodel$rfe_fit$results$Variables, bestfinalmodel$rfe_fit$results$ROC, col = "red", lwd = 2)

x = bestfinalmodel$rfe_fit$optsize
y = bestfinalmodel$rfe_fit$results$ROC[bestfinalmodel$rfe_fit$results$Variables == bestfinalmodel$rfe_fit$optsize]

points(x=x, y=y, col = "blue", lwd = 2)
points(x=x, y=y, col = "blue", cex = 3, lwd = 2)


# Plot using ggplot 
# Assuming bestfinalmodel$rfe_fit is available
rfe_results <- bestfinalmodel$rfe_fit$results
opt_size <- bestfinalmodel$rfe_fit$optsize
opt_roc <- rfe_results$ROC[rfe_results$Variables == opt_size]

dtseg <- data.frame(cbind(opt_roc, opt_size))


# Adjust the data to create a simulated break
rfe_results$adjusted_variables <- ifelse(rfe_results$Variables > 150, 
                                         rfe_results$Variables - 1246, 
                                         rfe_results$Variables)

# Create the plot
p <- ggplot(rfe_results, aes(x = adjusted_variables, y = ROC)) +
  geom_line(color = "red", size = 1.5) +
  geom_point(color = "black", size = 2) +
  geom_point(data = dtseg, aes(x = ifelse(opt_size > 150, opt_size - 1246, opt_size), y = opt_roc), 
             color = "blue", size = 5) +
  geom_segment(data = dtseg, aes(x = ifelse(opt_size > 150, opt_size - 1246, opt_size), y = 0.6, 
                                 xend = ifelse(opt_size > 150, opt_size - 1246, opt_size), yend = opt_roc), 
               linetype = "dashed", color = "blue") +
  geom_segment(data = dtseg, aes(x = 0, y = opt_roc, 
                                 xend = ifelse(opt_size > 150, opt_size - 1246, opt_size), yend = opt_roc), 
               linetype = "dashed", color = "blue") +
  scale_x_continuous(name = "Number of Features", 
                     breaks = c(seq(0, 150, by = 50), 200, 32),
                     labels = c(seq(0, 150, by = 50), "1446", 32)) +
  scale_y_continuous(limits = c(0.6, 1), name = "ROC (%)", 
                     breaks = c(seq(0.6, 1, by = 0.1), round(opt_roc,2)), 
                     labels = c(seq(0.6, 1, by = 0.1), round(opt_roc,2))*100) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) 


p

ggsave(plot = p, filename = "./Results/Results_LOOCV_Parallelized_KI_Computer_2024_08_11/Best_variables_from_the_best_model_AllmednonFilt.png", 
       units = "cm", height = 10, width = 15)

