



# Author: Antoine KD ABALO
# Date of creation: 20 August 2024
# Last modification date: 


# Objective: Investigating the distribution of the best feature set in the final model





# Import my functions
source("./Scripts/_00_List_of_functions.R")
rm(list=setdiff(ls(), "install"))


# Install the packages
install(c("magrittr", "dplyr", "ggplot2"))


# load the best predictors from the unfiltered best model

All_Median_pred <- rio::import("./Results/Results_LOOCV_Parallelized_KI_Computer_2024_08_11/Non_Filtered_Only_2024_08_11/All Median for dup_NonFiltered_Best_Predictors.xlsx") %>% 
  mutate(Variables = gsub("__Neur|__Infl|__Onco|__Card", "", Variables)) %>% 
  select(Variables) %>%  
  mutate(Model = "All Median for dup_NonFiltered")


Onco_Infl_pred <- rio::import("./Results/Results_LOOCV_Parallelized_KI_Computer_2024_08_11/Non_Filtered_Only_2024_08_11/Onco Infl_NonFiltered_Best_Predictors.xlsx") %>% 
  mutate(Variables = gsub("__Neur|__Infl|__Onco|__Card", "", Variables)) %>% 
  filter(Variables %in% All_Median_pred$Variables) %>% 
  select(Variables) %>%  
  mutate(Model = "Onco Infl_NonFiltered")


All_Onco_prot_for_dup_pred <- rio::import("./Results/Results_LOOCV_Parallelized_KI_Computer_2024_08_11/Non_Filtered_Only_2024_08_11/All Onco prot for dup_NonFiltered_Best_Predictors.xlsx") %>% 
  mutate(Variables = gsub("__Neur|__Infl|__Onco|__Card", "", Variables)) %>% 
  filter(Variables %in% All_Median_pred$Variables) %>% 
  select(Variables) %>%  
  mutate(Model = "All Onco prot for dup_NonFiltered")


# Load the filtered best predictors 



Filtered_All_Median_pred <- rio::import("./Results/Results_LOOCV_Parallelized_KI_Computer_2024_08_11/Filtered_Only_2024_08_11/All Median for dup_Filtered_Best_Predictors.xlsx") %>% 
  mutate(Variables = gsub("__Neur|__Infl|__Onco|__Card", "", Variables)) %>% 
  filter(Variables %in% All_Median_pred$Variables) %>% 
  select(Variables) %>%  
  mutate(Model = "All Median for dup_Filtered")


Filtered_Onco_Infl_pred <- rio::import("./Results/Results_LOOCV_Parallelized_KI_Computer_2024_08_11/Filtered_Only_2024_08_11/Onco Infl_Filtered_Best_Predictors.xlsx") %>% 
  mutate(Variables = gsub("__Neur|__Infl|__Onco|__Card", "", Variables)) %>% 
  filter(Variables %in% All_Median_pred$Variables) %>% 
  select(Variables) %>%  
  mutate(Model = "Onco Infl_Filtered")


Filtered_All_Onco_prot_for_dup_pred <- rio::import("./Results/Results_LOOCV_Parallelized_KI_Computer_2024_08_11/Filtered_Only_2024_08_11/All Onco prot for dup_Filtered_Best_Predictors.xlsx") %>% 
  mutate(Variables = gsub("__Neur|__Infl|__Onco|__Card", "", Variables)) %>% 
  filter(Variables %in% All_Median_pred$Variables) %>% 
  select(Variables) %>%  
  mutate(Model = "All Onco prot for dup_Filtered")


# Create a dataframe 
ls()

namestocall <- setdiff(ls(), c("install"))
namestocall

Varinothermodels <- data.frame()

for(o in namestocall) Varinothermodels <- rbind(Varinothermodels, get(o, envir = .GlobalEnv))

rm(list = setdiff(ls(), c("install", "Varinothermodels")))

length(unique(Varinothermodels$Model))


Varinothermodels %>% 
  group_by(Variables) %>% 
  summarise(N = n()) %>% 
  ungroup() %>% 
  right_join(Varinothermodels) %>% 
  mutate(Percent = (N/length(unique(Varinothermodels$Model)))*100) %>% 
  arrange(-desc(N)) %>% 
  distinct(Variables, .keep_all = T)-> toplot



# Plot the barplot

barplot(height = toplot$Percent, 
        names.arg = toplot$Variables, 
        col = viridis::viridis(length(unique(toplot$Variables)), option = "rocket"), 
        space = 1, 
        horiz = T, 
        las = 1, 
        cex.names = .57, 
        xlab = "% of resampling", 
        xpd = TRUE, 
        xlim = c(0, max(toplot$Percent)+0.5))

abline(v=50, lty = "dashed", lwd = 2)


# using ggplot2 

p <- ggplot(toplot, aes(x = Percent, y = reorder(Variables, Percent))) +
  geom_bar(stat = "identity", fill = viridis::viridis(length(unique(toplot$Variables)), option = "rocket")) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "black", size = 1) +  # Dashed line at 50%
  scale_x_continuous(limits = c(0, max(toplot$Percent) + 0.5), name = "% of selection in other models") +
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
       file = file.path("./Results/Results_LOOCV_Parallelized_KI_Computer_2024_08_11/Predictors_selected_in_other_models_in_the_final_model.png"), 
       units = "cm", 
       bg = "white",
       height = 25, 
       width = 15)

dev.off()

