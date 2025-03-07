


# Author: Antoine KD ABALO
# Date of creation: 20 February 2024
# Last modification date: 09 July 2024

#' Modifications:
#' Apply the final code to the whole data of 90 patients
#' Remove the duplicates proteins and missing POD24
#' Apply median value for the duplicated proteins



# Objectives: The main objective is to create the final dataset incoporating 
# the clinical and protein data all together - and secondly to provide an 
# overall description of the dataset to guide model selection decisions.


# Import my functions
source("./Scripts/_00_List_of_functions.R")
rm(list=setdiff(ls(), "install"))

# Install the packages
install(c("rio", "magrittr", "dplyr", "lubridate", "readxl", "tidyr", "corrplot", "caret", "rsample"))


# 1- Load the clinical dataset ----------

cleanMCL <- readxl::read_excel("Datasets/clean_MCL_data.xlsx", 
                               col_types = c("text", "numeric", "text", 
                                             "numeric", "text", "numeric", "numeric", 
                                             "text", "numeric", "text", "numeric", 
                                             "text", "text", "text", "text", "text", 
                                             "text", "text", "date", "date", "numeric", 
                                             "numeric", "date", "text", "text", 
                                             "date", "numeric", "numeric", "text", 
                                             "date", "text", "text"))
# View(cleanMCL)


# 2- Data checking -----------
# # Renaming variables - cleaning

names(cleanMCL)
cleanMCL %<>% 
  dplyr::rename("SampleID" = "prov-id",
         "MIPI_b_cat" = "MIPI-b cat",
         "ki67_morphology_desc" = "morphology_descr",
         "p53_high" = "p53.high.low", # what the levels 1 vs 0 means?
         "last_FU_dt" = "VITALSTATUS_DAT",
         "FU_yrs" = "os.y",
         "Death" = "os",
         "relpase_dt" = "REC1_dat",
         "ttt_dt_after_relapse" = "REC1_beh_dat", 
         "typ_TTT_after_relpase"= "REC1_beh", 
         "result_TTT_after_relapse" = "REC1_beh_resultat",
         "POD24" = "POD24 or PD") %>% 
  dplyr::mutate(PERSNR = gsub(pattern = "-", replacement="", x=PERSNR),
         Death = as.factor(Death), 
         female = ifelse(SEX=="female", 1, 0),
         age_dx_cat = cut(age_dx, breaks=c(-Inf, 60, 70, 80, Inf),
                          right=FALSE, labels=c("<60", "60-69", "70-79", "80+")), # with Patrick
         ki67 = gsub(pattern = "<", replacement = "", x=ki67),
         POD24 = ifelse(POD24 %in% c("na", "missing"), NA, POD24), 
         POD24 = ifelse(POD24 == "no", "No", POD24), 
         POD24 = ifelse(POD24 == "yes", "Yes", POD24), 
         ki67_morphology_desc = ifelse(ki67_morphology_desc %in% c("blastoid", "c i pellet och b i bm") | ki67_morphology_desc == "pleomorf", 
                                       "Blastoid_Pleomorphic", 
                                       ifelse(ki67_morphology_desc == "classic", 
                                              "Classic", ki67_morphology_desc))) %>% 
  dplyr::select(!c("above_65", "rec1.days", "rec1.y")) 


var <- c("MIPI_b_cat", "SR_dx", "LPK_dx", "Hb_dx", "CRP_dx", "Albumin_dx", "ki67", 
         "p53_high", "ki67_morphology_desc", "STADIUM_AA_dx")

for(i in var){
  dat = cleanMCL[[i]]
  dat[dat %in% grep("n/a", dat, ignore.case=T, value=T)] <- NA
  cleanMCL[, i] <- dat
}

# 3- Correct the date variables -----------

grep("_dat|_dt", names(cleanMCL), ignore.case=TRUE, value = TRUE)

var <- grep("_dat|_dt", names(cleanMCL), ignore.case=TRUE, value = TRUE)

for(i in var){
  dat = cleanMCL[[i]]
  print(class(dat))
}


table(cleanMCL$POD24, exclude = NULL) %>% addmargins()

# Import the new clinic data that has a complete POD24 information
## Load the new data 

newprot <- rio::import("./Datasets/New_data/New_POD24_2024_04_10.xlsx") %>% 
  rename(SampleID = 'sample-id')

cleanMCL %>% 
  select(-POD24) %>% 
  left_join(newprot) -> newCleanMCL

table(newCleanMCL$POD24, exclude = NULL) %>% addmargins()
# Missing      NA      No     Yes     Sum 
#      1       6      59      24      90 

# remove the missing
newCleanMCL %<>% 
  filter(!POD24 %in% c("Missing", "NA"))

table(newCleanMCL$POD24, exclude = NULL) %>% addmargins()
# No  Yes  Sum 
# 59  24   83 


# Overwrite cleanMCL data 
cleanMCL <- newCleanMCL; rm(newCleanMCL, newprot)


# 4- Save the patients clinic data -----------
save(cleanMCL, file = "./Datasets/cleanMCL.Rdata")
rio::export(cleanMCL, file = "./Datasets/cleanMCL_Antoine.xlsx")




# 5- Import the dataset that identify specifically each protein ----------------


olink_mcl <- readxl::read_excel("Datasets/olink_mcl_final_EXTRACTION.xlsx")


# Identify proteins measured in all the 4 panels
olink_mcl %>% 
  dplyr::summarise(n = dplyr::n(), .by = c(SampleID, Assay)) |>
  dplyr::filter(n > 1L) |> View()

# Define as a dataframe
olink_mcl %>% 
  dplyr::summarise(n = dplyr::n(), .by = c(SampleID, Assay)) |>
  dplyr::filter(n > 1L) |>
  as.data.frame() -> Dup_Prot_4Panels



# NOW: select variables that identify uniquely the proteins

olink_mcl %<>% 
  filter(!is.na(NPX)) %>% 
  dplyr::select("SampleID", "Assay", "Panel", "NPX") %>% 
  mutate(Assay_panel = paste0(Assay,"__",  substr(Panel, 1, 4))) %>%
  dplyr::select(!c("Assay", "Panel")) 

# and re-Pivot all wider to a rectangular dataset format:
olink_mcl %>% 
  # filter(Assay_panel %in% c("CXCL8__Card", "CXCL8__Infl", "CXCL8__Neur", "CXCL8__Onco")) %>% 
  dplyr::summarise(n = dplyr::n(), .by = c(SampleID, Assay_panel)) |>
  dplyr::filter(n > 1L) |> View()

olink_mcl %>% 
  tidyr::pivot_wider(names_from=Assay_panel, values_from=NPX) -> olink_mcl

# 6- Merge the tabular protein data to the clinical data ----------


# Select important clinical variables:
cat(names(cleanMCL), sep = ", \n")

cleanMCL %>% 
  dplyr::select(PERSNR, 
                SampleID, 
                age_dx,  
                MIPI_cat, 
                STADIUM_AA_dx, 
                LPK_dx, 
                SR_dx, 
                Hb_dx, 
                CRP_dx, 
                Albumin_dx, 
                ki67_morphology_desc, 
                p53_high, 
                POD24,  
                female) %>% 
  left_join(olink_mcl) -> clinicProtData


# 7- Save the final dataset --------

clinicProtData %>% 
  dplyr::select(PERSNR, 
                SampleID, 
                age_dx,  
                MIPI_cat, 
                STADIUM_AA_dx, 
                LPK_dx, 
                SR_dx, 
                Hb_dx, 
                CRP_dx, 
                Albumin_dx, 
                ki67_morphology_desc, 
                p53_high, 
                POD24,  
                female) %>% View

table(clinicProtData$POD24, exclude = NULL) %>% addmargins()

# Save the clinical and protein data combined:
save(clinicProtData, file = "./Output files/clinicProtData.Rdata")


# 8- Description -------------

# Overall description of all proteins
clinicProtData %>% 
  dplyr::select(contains(c("_Card", "_Onco", "_Infl", "_Neur"))) %>% 
  mutate_all(as.numeric) %>% 
  summarise(across(everything(), 
                   list("Min" = ~min(.), 
                        "Median" = ~median(.),
                        "Mean" = ~mean(.), 
                        "Max" = ~max(.), 
                        "Var" = ~var(.), 
                        "SD" = ~sd(.),
                        "ShapiroPvalue" = ~shapiro.test(.)$p.value))) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(stat = rownames(.)) %>% 
  separate(stat, into=c("Proteins", "Panel"), sep="__") %>% 
  separate(Panel, into=c("Panel", "Statistics"), sep="_") %>% 
  pivot_wider(names_from="Statistics", values_from="V1") %>%  
  as.data.frame() -> Description_Overall



# Description_Overall %>% 
#   View


# Overall description of proteins by POD24
clinicProtData %>% 
  dplyr::select(contains(c("_Card", "_Onco", "_Infl", "_Neur", "POD24"))
  ) %>% 
  mutate(across(tidyr::contains(c("_Card", "_Onco", "_Infl", "_Neur")), as.numeric)) %>% 
  mutate(POD24 = as.factor(POD24)) %>% 
  group_by(POD24) %>% 
  summarise(across(tidyr::everything(), 
                   list("Min" = ~min(.), 
                        "Median" = ~median(.),
                        "Mean" = ~mean(.), 
                        "Max" = ~max(.), 
                        "Var" = ~var(.), 
                        "SD" = ~sd(.),
                        "ShapiroPvalue" = ~shapiro.test(.)$p.value))) %>% 
  t() %>% 
  as.data.frame() %>% 
  rename(POD24_Yes = V1, 
         POD24_No = V2) %>% 
  slice(-1) %>% 
  mutate(stat = rownames(.)) %>% 
  separate(stat, into=c("Proteins", "Panel"), sep="__") %>% 
  separate(Panel, into=c("Panel", "Statistics"), sep="_") %>% 
  pivot_wider(names_from="Statistics", values_from=c("POD24_Yes", "POD24_No")) %>%  
  as.data.frame() %>%
  select(contains(c("Proteins", "Panel", "Min", "Mean", "Median", "Max", 
                    "Var", "SD", "Pvalu"))) -> Description_POD24



# Description_POD24 %>% 
#   View



# 9- Compute t-test of wilcoxon rank test given the distribution of the protein ---------------------


# Object to store p-values
p_values <- list()

# Limit the data to necessary variables
clinicProtData %>% 
  dplyr::select(contains(c("_Card", "_Onco", "_Infl", "_Neur", "POD24"))) %>% 
  mutate(across(contains(c("_Card", "_Onco", "_Infl", "_Neur")), as.numeric)) -> Red_clinic_ProtData

# Loop over each variable to perform tests
for (variable in colnames(Red_clinic_ProtData)) {
  if (variable != "POD24" && is.numeric(Red_clinic_ProtData[[variable]])) {
    # Check normality using Shapiro-Wilk test
    shapiro_test <- shapiro.test(Red_clinic_ProtData[[variable]])
    
    if (shapiro_test$p.value > 0.05) {
      # If data is normally distributed, perform t-test
      test_result <- t.test(Red_clinic_ProtData[[variable]] ~ Red_clinic_ProtData$POD24)
    } else {
      # If data is not normally distributed, perform Wilcoxon rank-sum test
      test_result <- wilcox.test(Red_clinic_ProtData[[variable]] ~ Red_clinic_ProtData$POD24)
    }
    
    p_values[[variable]] <- test_result$p.value
  }
}

# Adjust p-values for multiple comparisons using FDR
adjusted_p_values <- p.adjust(unlist(p_values), method = "fdr")

# Combine results into a data frame
Pvalues_ttest_Wilcoxon <- data.frame(
  Variable = names(adjusted_p_values),
  P_Value = unlist(p_values),
  Adjusted_P_Value = adjusted_p_values
)

str(Pvalues_ttest_Wilcoxon)

table(Pvalues_ttest_Wilcoxon$P_Value <= 0.05, exclude = NULL)
table(Pvalues_ttest_Wilcoxon$Adjusted_P_Value <= 0.05, exclude = NULL)

# Print significant results
significant_Ajusted_Pvalues <- Pvalues_ttest_Wilcoxon %>%
  filter(Adjusted_P_Value <= 0.05)

head(significant_Ajusted_Pvalues)



# Export the descriptions:
rio::export(Description_Overall, file = "./Results/Description_Overall.xlsx")
rio::export(Description_POD24, file = "./Results/Description_by_POD24.xlsx")
rio::export(Pvalues_ttest_Wilcoxon, file = "./Results/Pvalues_ttest_Wilcoxon.xlsx")
rio::export(significant_Ajusted_Pvalues, file = "./Results/Significant_Ajusted_Pvalues.xlsx")

rm(significant_Ajusted_Pvalues, p_values, Description_Overall, Description_POD24,
   shapiro_test, test_result, variable, adjusted_p_values, i, dat, var)





# Check for correlation among the significant proteins -------------------------

# Retain the significant proteins at univariate comparison p-value (non-adjusted)
Signft_prot <- Pvalues_ttest_Wilcoxon$Variable[Pvalues_ttest_Wilcoxon$P_Value <= 0.05]

# Reduce the protein set to the significant ones
Red_clinic_ProtData %<>% 
  dplyr::select(POD24, paste(Signft_prot)) 



# Test for correlation:
(co <- cor.test(Red_clinic_ProtData[[Signft_prot[2]]], Red_clinic_ProtData[[Signft_prot[3]]]))
names(co)
co$estimate
plot(Red_clinic_ProtData[[Signft_prot[2]]], Red_clinic_ProtData[[Signft_prot[3]]])

cor_matrix <- cor(Red_clinic_ProtData[,names(Red_clinic_ProtData[, names(Red_clinic_ProtData)!="POD24"])])
# Plot the correlation matrix using corrplot
corrplot(cor_matrix, method = "color", insig = "blank", tl.pos='n')
# Save
# png(file = "./Results/Plotcorrelations.png", units ="cm", height = 15, width = 15, res = 500)
# corrplot(cor_matrix, method = "color", insig = "blank", tl.pos='n'); dev.off()


# # Create a data frame to store results
# cor_test_results <- data.frame(Column1 = character(), Column2 = character(), Estimate = numeric(),
#                                PValue = numeric(), stringsAsFactors = FALSE)
# 
# df <- Red_clinic_ProtData[,names(Red_clinic_ProtData[, names(Red_clinic_ProtData)!="POD24"])]
# 
# for (i in 1:(ncol(df)-1)) {
#   for (j in (i+1):ncol(df)) {
#     result <- cor.test(df[[i]], df[[j]])
#     cor_test_results <- rbind(cor_test_results, 
#                               data.frame(Column1 = names(df)[i], 
#                                          Column2 = names(df)[j], 
#                                          Estimate = result$estimate,
#                                          PValue = result$p.value))
#   }
# }
# 
# # Filter significant correlations 
# significant_correlations <- cor_test_results[cor_test_results$Estimate >0.5 & cor_test_results$PValue < 0.05 |
#                                                cor_test_results$Estimate < -0.5 & cor_test_results$PValue < 0.05, ]
# significant_correlations %>% View
# 
# rio::export(significant_correlations, file = "./Results/significant_correlations.xlsx")
# 
# rm(co, significant_correlations, df, cor_test_results, cor_matrix, Pvalues_ttest_Wilcoxon, result)





# 10- Split the data into test and training set -----------------

set.seed(2024)

# Select the protein data only 
clinicProtData %>% 
  dplyr::select(paste(Signft_prot), contains(c("POD24"))) %>% 
  mutate(across(contains(c("_Card", "_Onco", "_Infl", "_Neur")), as.numeric)) %>%
  mutate(across(contains(c("POD24")), as.factor)) %>% 
  initial_split(strata = "POD24", prop = .75) -> prot_split

# Create initial partition of the dataset
# Create the training data
Prot_train <- training(prot_split) %>% as.data.frame()
# Create the test data
Prot_test  <- testing(prot_split) %>% as.data.frame()


# Check balance in the train and test set
# In the original data
table(clinicProtData$POD24) %>% addmargins()
prop.table(table(clinicProtData$POD24))

# In the training data
table(Prot_train$POD24) %>% addmargins()
prop.table(table(Prot_train$POD24))

# In the test data
table(Prot_test$POD24) %>% addmargins()
prop.table(table(Prot_test$POD24))


# 11- Investigate NA in the final data ------------------------------------------

sum(is.na(clinicProtData)) # 87 In the whole data
sum(is.na(Prot_train)) # 0 In the training data
sum(is.na(Prot_test)) # 0 In the test data


# 11.1- Investigate the distribution of the fifty first prot using a boxplot ------
names(Prot_train)
boxplot(Prot_train[8])

par(mfrow = c(2, 3))
cols <- colors()
plotboxp <- function(start){
  dat <- Prot_train[, names(Prot_train) != "POD24"]
  boxplot(dat[[start]], main = names(dat[start]), col = cols[start], horizontal = F)
}
for(start in 1:(6*5))plotboxp(start)


# Save the plots
# Directory for saving plots
save_dir <- "./Results/Boxplot_panel"

# Check the directory exists, create if not
if (!file.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}



(n <- ncol(Prot_train)-1)

# Create a panel of 9 plots and save them directly 
for (i in seq(1, n, by = 9)) {
  # Set up the layout for multiple plots
  par(mfrow = c(3, 3))
  
  # Loop through each individual plot and add it to the panel
  for (start in i:(i + 8)) {
    if(start <= ncol(Prot_train)-1){
    # Create a new plot
    plotboxp(start)
      }
  }
  
  # Save the combined plot
  filename <- sprintf("combined_plot_%03d.png", i)  # Adjust filename as needed
  filepath <- file.path(save_dir, filename)
  
  # Save the entire panel
  dev.copy(png, filepath, units = "cm", height = 15, width = 15, res = 300)
  dev.off()
  
  # Reset the layout to its original state
  par(mfrow = c(1, 1))
}



# 12- Investigate the distribution of the fifty first prot using an histogram/density plot -------


# Create a function to plot the density

plotdensity <- function(start, Density = FALSE) {
  dat <- Prot_train[, names(Prot_train) != "POD24"]
  if(Density){
    hist(dat[[start]], main = names(dat[start]), xlab = "", ylab="", prob=TRUE )
    lines(density(dat[[start]]), col = cols[start], lwd = 2)
  }else{
    plot(density(dat[[start]]), col = cols[start], xlab = "",
         ylab="", lwd = 2, main = names(dat[start]))
    }
}

plotdensity(6, Density = FALSE)
plotdensity(6, Density = TRUE)
plotdensity(6)

for(i  in 1:24)plotdensity(i, Density = TRUE)

# Save the graphic of density plot only ------- 

# Directory for saving plots
save_dir <- "./Results/Density_panel"

# Ensure the directory exists, create if not
if (!file.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}


(n <- ncol(Prot_train)-1)

# Create a panel of 9 plots and save them directly 
for (i in seq(1, n, by = 9)) {
  # Set up the layout for multiple plots
  par(mfrow = c(3, 3))
  
  # Loop through each individual plot and add it to the panel
  for (start in i:(i + 8)) {
    if(start <= ncol(Prot_train)-1){
      # Create a new plot
      plotdensity(start)
    }
  }
  
  # Save the combined plot
  filename <- sprintf("Density_panel_%03d.png", i)  # Adjust filename as needed
  filepath <- file.path(save_dir, filename)
  
  # Save the entire panel
  dev.copy(png, filepath, units = "cm", height = 15, width = 15, res = 300)
  dev.off()
  
  # Reset the layout to its original state
  par(mfrow = c(1, 1))
}




# Save the graphic of barplot and density plot only ------- 

save_dir <- "./Results/Histo_panel"

# Ensure the directory exists, create if not
if (!file.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}




# Create a panel of 9 plots and save them directly 
for (i in seq(1, n, by = 9)) {
  # Set up the layout for multiple plots
  par(mfrow = c(3, 3))
  
  # Loop through each individual plot and add it to the panel
  for (start in i:(i + 8)) {
    if(start <= ncol(Prot_train)-1){
      # Create a new plot
      plotdensity(start, Density = TRUE)
    }
  }
  
  # Save the combined plot
  filename <- sprintf("Histo_panel_%03d.png", i)  # Adjust filename as needed
  filepath <- file.path(save_dir, filename)
  
  # Save the entire panel
  dev.copy(png, filepath, units = "cm", height = 15, width = 15, res = 300)
  dev.off()
  
  # Reset the layout to its original state
  par(mfrow = c(1, 1))
}





# 13- Investigate near zero and zero variances -------


caret::nearZeroVar(clinicProtData[, names(clinicProtData)!="POD24"], saveMetrics = TRUE) %>%
  tibble::rownames_to_column() %>%
  filter(nzv)

# [1] rowname       freqRatio     percentUnique zeroVar       nzv          
# <0 rows> (or 0-length row.names)


# Save the partiionned data ----------------

save(Prot_test, file = "./Output files/Prot_test.Rdata")
save(Prot_train, file = "./Output files/Prot_train.Rdata")






# 14- Further managements removing age and sex, and creating different datasets -----------------

# 14.1 - For the training data -------------

# Convert all features into numeric

Prot_train %<>%
  mutate(across(tidyr::contains(c("_Card", "_Onco", "_Infl", "_Neur")), as.numeric)) %>% 
  mutate(POD24 = as.factor(POD24)) 


# Define the datasets ------------------

cardio_tr <- Prot_train %>% 
  select(POD24, contains("card"))

neuro_tr <- Prot_train %>% 
  select(POD24, contains("neur"))

onco_tr <- Prot_train %>% 
  select(POD24, contains("onco"))

inflam_tr <- Prot_train %>% 
  select(POD24, contains("infl"))

Onco_Inflam_tr <- Prot_train %>% 
  select(POD24, contains(c("onco", "infl")))


all_tr <- Prot_train

trainlist <- list(Cardiology = cardio_tr, 
                  Neurology = neuro_tr, 
                  Oncology = onco_tr, 
                  Inflamation = inflam_tr,
                  Onco_Inflamation = Onco_Inflam_tr,
                  All = all_tr)
rm(cardio_tr, neuro_tr, inflam_tr, onco_tr, Onco_Inflam_tr, all_tr)







# 14.2 - For the test data -------------


# Convert all features into numeric


Prot_test %<>%
  mutate(across(tidyr::contains(c("_Card", "_Onco", "_Infl", "_Neur")), as.numeric)) %>% 
  mutate(POD24 = as.factor(POD24)) 


# Define the datasets ------------------

cardio_tr <- Prot_test %>% 
  select(POD24, contains("card"))

neuro_tr <- Prot_test %>% 
  select(POD24, contains("neur"))

onco_tr <- Prot_test %>% 
  select(POD24, contains("onco"))

inflam_tr <- Prot_test %>% 
  select(POD24, contains("infl"))

Onco_Inflam_tr <- Prot_test %>% 
  select(POD24, contains(c("onco", "infl")))

all_tr <- Prot_test

testlist <- list(Cardiology = cardio_tr, 
                 Neurology = neuro_tr, 
                 Oncology = onco_tr, 
                 Inflamation = inflam_tr, 
                 'Onco Inflamation' = Onco_Inflam_tr,
                 All = all_tr)
rm(cardio_tr, neuro_tr, inflam_tr, Onco_Inflam_tr, onco_tr, all_tr)





# 14.3- Save Further managements removing age and sex, and creating different datasets -----------------

save(trainlist, file = "./Output files/trainlist.Rdata")
save(testlist, file = "./Output files/testlist.Rdata")



# 15- Further managements removing common proteins, and creating different datasets -----------------

# 15.1 - For the training data -------------

(ProtCard <- gsub("__Card", "", names(trainlist$Cardiology), ignore.case = T))
(ProtOnco <- gsub("__Onco", "", names(trainlist$Oncology), ignore.case = T))
(ProtNeur <- gsub("__Neur", "", names(trainlist$Neurology), ignore.case = T))
(ProtInfl <- gsub("__Infl", "", names(trainlist$Inflamation), ignore.case = T))

CommonProt <- Reduce(intersect, list(ProtCard, ProtOnco, ProtNeur, ProtInfl))
CommonProt

unique(Dup_Prot_4Panels$Assay)


# remove POD24
CommonProt <- CommonProt[CommonProt!="POD24"]
CommonProt

#---- DO NOT RUN IF CommonProt = 0 ------------------


# # Remove these proteins and keep them only if they are from Oncology
# Prottoremove <- c(paste0(CommonProt, "__Neur"),
#                   paste0(CommonProt, "__Card"),
#                   paste0(CommonProt, "__Infl"))
# length(Prottoremove)
# # [1] 3
# ncol(Prot_train)
# # [1] 318
# Prot_train_nodup <- Prot_train[, !(names(Prot_train) %in% Prottoremove)]
# ncol(Prot_train_nodup)
# # [1] 318
# ncol(Prot_train) - ncol(Prot_train_nodup)
# # [1] 0
# 
# names(Prot_train_nodup)[names(Prot_train_nodup) %in% paste0(CommonProt, "__Onco")]
# 
# # Include them into the train list object 
# names(trainlist)
# trainlist$`All no dup Onco` <- Prot_train_nodup
# names(trainlist)
# 
# # 15.2 - For the test data -------------
# 
# Prot_test_nodup <- Prot_test[, !(names(Prot_test) %in% Prottoremove)]
# 
# # Include them into the test list object 
# names(testlist)
# testlist$`All no dup Onco` <- Prot_test_nodup
# names(testlist)
# 
# 
# 
# save(Prot_train_nodup, file = "./Output files/Prot_train_nodup.Rdata")
# save(Prot_test_nodup, file = "./Output files/Prot_test_nodup.Rdata")
# 
# 
# 
# rm(list=setdiff(ls(), c("cleanMCL", "Dup_Prot_4Panels", "Prot_test", "Prot_train", "testlist", "trainlist", "CommonProt")))
# 
# 
# 
# 
# 
# # 16- Further managements take the median value of common proteins, and creating different datasets -----------------
# 
# CommonProt
# # [1] "CXCL8" "IL6"   "TNF" 
# unique(Dup_Prot_4Panels$Assay)
# # [1] "CXCL8" "IL6"   "TNF" 
# 
# # For the training data
# Prot_train$CXCL8_med <- apply(Prot_train[, grep("CXCL8__", names(Prot_train), ignore.case = T, value = T)], 1, median)
# Prot_train$IL6_med <- apply(Prot_train[, grep("IL6__", names(Prot_train), ignore.case = T, value = T)], 1, median)
# Prot_train$TNF_med <- apply(Prot_train[, grep("TNF__", names(Prot_train), ignore.case = T, value = T)], 1, median)
# 
# 
# # Remove the duplicated prot
# Prot_train_noDup_med <- Prot_train %>% 
#   select(!contains(c("CXCL8__", "IL6__", "TNF__")))
# 
# # Include the median values of the duplicated pro into the train list object 
# names(trainlist)
# trainlist$`All no dup Median` <- Prot_train_noDup_med
# names(trainlist)
# 
# 
# 
# 
# # For the test data
# Prot_test$CXCL8_med <- apply(Prot_test[, grep("CXCL8__", names(Prot_test), ignore.case = T, value = T)], 1, median)
# Prot_test$IL6_med <- apply(Prot_test[, grep("IL6__", names(Prot_test), ignore.case = T, value = T)], 1, median)
# Prot_test$TNF_med <- apply(Prot_test[, grep("TNF__", names(Prot_test), ignore.case = T, value = T)], 1, median)
# 
# 
# # Remove the duplicated prot
# Prot_test_noDup_med <- Prot_test %>% 
#   select(!contains(c("CXCL8__", "IL6__", "TNF__")))
# 
# # Include the median values of the duplicated pro into the test list object 
# names(testlist)
# testlist$`All no dup Median` <- Prot_test_noDup_med
# names(testlist)
# 
# 
# 
# 
# # 17- Save Further managements removing common proteins, and creating different datasets -----------------
# 
# save(trainlist, file = "./Output files/trainlist.Rdata")
# save(testlist, file = "./Output files/testlist.Rdata")
# 



#---- END OF: DO NOT RUN IF CommonProt = 0 ------------------





rm(list=setdiff(ls(), c("cleanMCL", "Dup_Prot_4Panels", "Prot_test", "Prot_train", "testlist", "trainlist", "CommonProt")))
