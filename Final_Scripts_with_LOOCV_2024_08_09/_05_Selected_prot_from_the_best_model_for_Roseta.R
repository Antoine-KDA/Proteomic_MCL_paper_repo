


# Author: Antoine KD ABALO
# Date of creation: 12 August 2024
# Last modification date: 


# Objective: Select the proteins from the best models for R.Rosetta





# Import my functions
source("./Scripts/_00_List_of_functions.R")
rm(list=setdiff(ls(), "install"))


# Install the packages
install(c("magrittr", "dplyr"))


# Load the dataset
load("./Output files/trainlist_12July2024_Ant.Rdata")
load("./Output files/testlist_12July2024_Ant.Rdata")

# Take the data corresponding to the final best model
#' Note: form the summary excel sheet, this is *All Med For Dup NOn Filtered*

# Import the best predictors from the final best model
best_predictors <- rio::import("./Results/Results_LOOCV_Parallelized_KI_Computer_2024_08_11/Non_Filtered_Only_2024_08_11/All Median for dup_NonFiltered_Best_Predictors.xlsx")

best_predictors <- best_predictors$Variables


# Take the data corresponding to the final best model

trainlist_12July2024_Ant[["All: Median for dup"]] %>% 
  bind_rows(testlist_12July2024_Ant[["All: Median for dup"]]) %>% 
  select("POD24", all_of(best_predictors)) -> Data_for_Rosetta



rm(trainlist_12July2024_Ant, testlist_12July2024_Ant, best_predictors)


# Export the results
rio::export(Data_for_Rosetta, file = "./Results/Results_LOOCV_Parallelized_KI_Computer_2024_08_11/Data_for_Rosetta.xlsx")




