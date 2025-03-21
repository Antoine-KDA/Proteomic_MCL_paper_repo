


# Author: Antoine KD ABALO
# Date of creation: 20 August 2024
# Last modification date: 


# Objective: Investigating the distribution of the best models from the resampling process





# Import my functions
source("./Scripts/_00_List_of_functions.R")
rm(list=setdiff(ls(), "install"))


# Install the packages
install(c("magrittr", "dplyr"))

# Load the resampling models
load("./Results/Results_LOOCV_Parallelized_KI_Computer_2024_08_11/Non_Filtered_Only_2024_08_11/RFrfe_ll_results_NonFiltered_LOOCV_list.Rdata")

# Check the names in the list
ls()
names(All_results_NonFiltered_LOOCV)


# Take the corresponding best model
#' Note: form the summary excel sheet, this is *All Med For Dup NOn Filtered*

bestmodel <- All_results_NonFiltered_LOOCV$`All Median for dup_results_parallel_LOOCV`
names(bestmodel)
bestmodel <- bestmodel[!names(bestmodel) %in% c("start_time","end_time","runtime")]

# Plot the rsampling of one model

bestmodel[[1]]$rfe_fit$results%>% View

plot(bestmodel[[1]]$rfe_fit$results$Variables, bestmodel[[1]]$rfe_fit$results$ROC,
     xlim = c(0, 150), ylim = c(0.5, 1), type = "n")
lines(bestmodel[[1]]$rfe_fit$results$Variables, bestmodel[[1]]$rfe_fit$results$ROC, col = "red", lwd = 1)

x = bestmodel[[1]]$rfe_fit$optsize
y = bestmodel[[1]]$rfe_fit$results$ROC[bestmodel[[1]]$rfe_fit$results$Variables == bestmodel[[1]]$rfe_fit$optsize]

points(x=x, y=y, col = "blue", lwd = 2)
#points(x=x, y=y, col = "blue", cex = 3, lwd = 2)

# add a second plot

lines(bestmodel[[2]]$rfe_fit$results$Variables, bestmodel[[2]]$rfe_fit$results$ROC, col = "red", lwd = 1)

x = bestmodel[[2]]$rfe_fit$optsize
y = bestmodel[[2]]$rfe_fit$results$ROC[bestmodel[[2]]$rfe_fit$results$Variables == bestmodel[[2]]$rfe_fit$optsize]

points(x=x, y=y, col = "blue", lwd = 2)
#points(x=x, y=y, col = "blue", cex = 3, lwd = 2)

length(bestmodel)

# Plot all together

Data_output <- data.frame()

plot(bestmodel[[1]]$rfe_fit$results$Variables, bestmodel[[1]]$rfe_fit$results$ROC,
     xlim = c(0, 150), ylim = c(0.5, 1), type = "n")

for(n in 1:length(bestmodel)){
  lines(bestmodel[[n]]$rfe_fit$results$Variables, bestmodel[[n]]$rfe_fit$results$ROC, col = "red", lwd = 1)
  
  Size = bestmodel[[n]]$rfe_fit$optsize
  ROC = bestmodel[[n]]$rfe_fit$results$ROC[bestmodel[[n]]$rfe_fit$results$Variables == bestmodel[[n]]$rfe_fit$optsize]

  points(x=Size, y=ROC, col = "blue", lwd = 2)
  
  # Export 
  data.frame(
    cbind(
      Size, ROC, Resample = n
    )) %>% 
    rbind(Data_output) -> Data_output
  
  
  rm(Size, ROC)
}


## Add a point to the final model picked, Nm_Seq = 45 
x = bestmodel[[45]]$rfe_fit$optsize
y = bestmodel[[45]]$rfe_fit$results$ROC[bestmodel[[45]]$rfe_fit$results$Variables == bestmodel[[45]]$rfe_fit$optsize]

points(x=x, y=y, col = "darkgreen", cex = 3, lwd = 2)

plot(Data_output$Resample, Data_output$Size)
hist(Data_output$Size, xlim = c(0, 100))


