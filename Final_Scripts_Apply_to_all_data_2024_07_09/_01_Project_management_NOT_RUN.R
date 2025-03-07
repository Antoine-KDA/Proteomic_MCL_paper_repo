
# Author: Antoine KD ABALO
# Date of creation: 20 February 2024
# Last modification date: 




# This project assesses the risk of infections after Mantle Cell Lymphoma diagnose
# We used a population-based registers data to identify MCL patients and 1:10 matched 
# comparators on sex and age


## 1- Project management

rm(list=ls())

# Run only one time

dir.create("Datasets")
dir.create("Results")
dir.create("Output files")
dir.create("Scripts")


## 2- Data import 


# Save the data for further analyses

renv::snapshot()
