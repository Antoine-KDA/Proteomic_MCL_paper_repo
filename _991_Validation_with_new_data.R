



# Validation protein data 

rm(list=ls())
# Import my functions
source("./Scripts/_00_List_of_functions.R")
rm(list=setdiff(ls(), "install"))

# Install the packages
install(c("magrittr", "dplyr", "vip", "caret", "MLmetrics", "recipes", "ggplot2"))



## Load the proteomic data
proteom <- rio::import("./Datasets/Proteomic_Patrick/protein expression data_2024_02_07.xlsx")

## Load the new data 

newprot <- rio::import("./Datasets/New_data/New_POD24_2024_04_10.xlsx")

names(newprot)
table(newprot$`sample-id` %in% proteom$SampleID)

proteom %>% 
  filter(SampleID %in% newprot$`sample-id`) %>% 
  filter(!duplicated(.data$SampleID)) %>% 
  nrow()


# Import the original prot data, take patients with Missing POD24 previously 
absentid <- rio::import("./Datasets/MCL Olink/mcl_metadata.xlsx") %>% 
  filter(POD24 %in% grep("na|miss", .data$POD24, ignore.case = T, value = T))


# Manual copy from excel
# absentid <- c(4070401964,
#                 4069735534,
#                 4069733036,
#                 4059675189,
#                 4059671148,
#                 4050412472,
#                 4033753770,
#                 4031498433,
#                 4030924493,
#                 4017484688,
#                 4016555864,
#                 4014406432,
#                 4014400111,
#                 4014399745,
#                 4000612244,
#                 4000612118,
#                 4000612065,
#                 4000612062,
#                 4000611153,
#                 4000611052,
#                 4000610837,
#                 3500609823)


table(newprot$`sample-id` %in% absentid$SampleID)  

newprot %>% 
  filter(`sample-id` %in% absentid$SampleID) -> newprotid
table(newprotid$POD24, exclude = NULL)

newprotid %<>% 
  filter(POD24 != "NA") %>% 
  rename(SampleID = `sample-id`)

nrow(newprotid)
# [1] 17
rm(newprot, proteom)

# Import the dataset that identify specifically each protein ----------------


olink_mcl <- readxl::read_excel("Datasets/olink_mcl_final_EXTRACTION.xlsx")

names(olink_mcl)

# NOW: select variables that identify uniquely the proteins

table(olink_mcl$SampleID %in% newprotid$SampleID, exclude = NULL)
table(newprotid$SampleID %in% olink_mcl$SampleID, exclude = NULL)

olink_mcl %>% 
  dplyr::select("SampleID", "Assay", "Panel", "NPX") %>% 
  filter(SampleID %in% newprotid$SampleID) %>% 
  mutate(Assay_panel = paste0(Assay,"__",  substr(Panel, 1, 4))) %>%
  right_join(newprotid) %>% 
  dplyr::select(!c("Assay", "Panel")) -> newprotid2

table(is.na(newprotid2$NPX))
# # 
# FALSE  TRUE 
# 24735   289 

newprotid2 %>% 
  filter(!is.na(NPX)) -> newprotid; rm(newprotid2)


# and re-Pivot all wider to a rectangular dataset format:
newprotid %>% 
  filter(Assay_panel %in% c("CXCL8__Card", "CXCL8__Infl", "CXCL8__Neur", "CXCL8__Onco")) %>% 
  dplyr::summarise(n = dplyr::n(), .by = c(SampleID, Assay_panel)) |>
  dplyr::filter(n > 1L) |> View()

newprotid %>% 
  tidyr::pivot_wider(names_from=Assay_panel, values_from=NPX) -> dataprotValidation


# Check if these sampleID were already used in the previously data
load("./Output files/clinicProtData.Rdata")

table(dataprotValidation$SampleID %in% clinicProtData$SampleID, exclude = NULL)
# FALSE 
# 17 
# None of these samples were previously used. 
rm(clinicProtData)

# Keep only proteins dat were present in the training data
load("./Output files/Prot_train_nodup.Rdata")
table(names(Prot_train_nodup) %in% names(dataprotValidation), exclude = NULL)

dataprotValidation2 <- dataprotValidation[, names(dataprotValidation) %in% names(Prot_train_nodup)]

table(names(dataprotValidation2) %in% names(Prot_train_nodup), exclude = NULL)
table(names(Prot_train_nodup) %in% names(dataprotValidation2), exclude = NULL)
cbind(table(Prot_train_nodup$POD24, exclude = NULL), prop.table(table(Prot_train_nodup$POD24, exclude = NULL))*100)
cbind(table(dataprotValidation$POD24, exclude = NULL), prop.table(table(dataprotValidation$POD24, exclude = NULL))*100)

dataprotValidation <- dataprotValidation2; rm(dataprotValidation2, Prot_train_nodup)

# reformat the data 

dataprotValidation %<>% 
  # dplyr::select(!SampleID) %>% 
  mutate(across(contains(c("_Card", "_Onco", "_Infl", "_Neur")), as.numeric)) %>%
  mutate(across(contains(c("female", "POD24")), as.factor))

# Save the validation dataset --------------------

save(dataprotValidation, file = "./Output files/dataprotValidation.Rdata")
rio::export(dataprotValidation, file = "./Output files/dataprotValidation.xlsx")

rm(absentid, newprotid, olink_mcl)
dataprotValidation <- list(All=dataprotValidation)




# Perform the validation -----------------

# Load the finally selected best model which is treeBagRFE_NoDupProt_Filtered_NF_Size5
load("./Results/Results_Final_Scripts_2024_04_10/Filtered_NF_2024_04_10/treeBagRFE_NoDupProt_Filtered_NF_Size5_station.Rdata")
# The best model is the filtered one
bestfinalmodel <- treeBagRFE_NoDupProt_results_list$All_filtered



## Predict class on the test data using the final model 
pred_class_test <- predict(bestfinalmodel, dataprotValidation$All)

## Create confusion matrix 

(confmatTest <- confusionMatrix(data = relevel(pred_class_test$pred, ref = 'Yes'),
                                reference = relevel(dataprotValidation$All$POD24, ref = 'Yes')))


(Accuracy_test <- paste0(sprintf("%.1f", confmatTest$overall[["Accuracy"]]*100), " (", 
                         sprintf("%.1f", confmatTest$overall[["AccuracyLower"]]*100), "-", 
                         sprintf("%.1f", confmatTest$overall[["AccuracyUpper"]]*100), ")"))

# Get the AUC
# Use my own defined function
modAuc_test <- GetAuc(Actual_class = dataprotValidation$All$POD24, Pred_class = pred_class_test$pred)
modAuc_test

## Export model performances 

perfexport_test <- ModPerfExport(ConfMatrix = confmatTest,
                                 Accuracy = Accuracy_test, 
                                 Model_AUC = modAuc_test)
perfexport_test


