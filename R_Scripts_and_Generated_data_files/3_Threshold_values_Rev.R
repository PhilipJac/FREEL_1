#################
# Jacobson et al.
# Navigating the salinity gradient: Individual variation in habitat use and migration of European eel revealed by otolith microchemistry 
#################

## Last update: 2026-10-28

### R-Script: 3_Threshold_values.R

## Load needed packages:
library(tidyverse)
library(segclust2d)
library(stringr)
library(flextable)

########### Load data ############

### Load the data generated in Script 1_Input_data.R ####

# 1. LA-ICP-MS otolith chemistry transect data combined with individual level data from the databases Sötebasen and KUL:
LA_ICP_MS_KUL_Sotebasen = readRDS("R_Scripts_and_Generated_data_files/LA_ICP_MS_KUL_Sotebasen.Rds") # containing 29 columns with 950 477 rows of data 
head(LA_ICP_MS_KUL_Sotebasen)

# 2. FE-EPMA otolith chemistry transect data combined with individual level data from the database Sötebasen:
FE_EPMA_Sotebasen = readRDS("R_Scripts_and_Generated_data_files/FE_EPMA_Sotebasen.Rds") # containing 22 colummns with 58 740 rows of data
head(FE_EPMA_Sotebasen)

# 3. All individual level data combined for both methods, one row = one individual
FE_EPMA_LA_ICP_MS_Individual_data = readRDS("R_Scripts_and_Generated_data_files/FE_EPMA_LA_ICP_MS_Individual_data.Rds") # containing 13 columns with 3648 rows of data 
head(FE_EPMA_LA_ICP_MS_Individual_data)

## FE-EPMA Sr:Ca threshold values ##

# Select reference eel from lake Mälaren and lake Hjälmaren:
## FE-EPMA threshold values:
FE_EPMA_Sotebasen_reference_eel<-FE_EPMA_Sotebasen |> 
    filter(IndividNr %in% c("31645", "31704", "31945", "31950", "31951", "31953", "31955", #Hjälmaren
                            "31957", "31958", "31960", "31961", "31962", # Hjälmaren, 12 inds
                            "31586", "31587", "31602", "31603", "31609", "31611", "31612", #Mälaren
                            "31627", "32074", "32077", "32079", "32086", "32596", "32598", "32601")) #Mälaren, 15 inds)
print(FE_EPMA_Sotebasen_reference_eel, n=Inf)

sort(unique(FE_EPMA_Sotebasen_reference_eel$Site))

FE_EPMA_Sotebasen_reference_eel |>
  group_by(Water_area, Site) |>
   summarise(n=n_distinct(IndividNr))


# Save dataframe for Script 7_Suppl_Info.R
saveRDS(FE_EPMA_Sotebasen_reference_eel, file = "R_scripts_and_Generated_data_files/FE_EPMA_Sotebasen_reference_eel.Rds")

### Make a loop for the segmentation analyses:
result_list <- list()

## Create a folder to save the segmentation and clustering output in the output folder:
output_dir <- "Output/Reference_eel_FE_EPMA_segclust"
dir.create(output_dir, showWarnings = FALSE)

individual <- unique(FE_EPMA_Sotebasen_reference_eel$IndividNr)

# Run the loop for every individual:
for (individ in individual) {
  
  ## Calculate the Sr:Ca ratio for every individual, remove Na-values and run the loop
  individ_data <- FE_EPMA_Sotebasen_reference_eel %>%
    filter(IndividNr == individ) %>%
    select(IndividNr, COMMENT, Sr, Ca, Date, Water_area) %>%
    mutate(Sr_Ca_ratio = (as.numeric(as.character(Sr)) / as.numeric(as.character(Ca))) * 1000) %>%
    filter(!is.na(Sr_Ca_ratio)) %>%
    as.data.frame()
  
  #### Runf the segclust function
  mode_segclust <- tryCatch({
   segclust(individ_data,
             Kmax = 6,
             lmin = 5,
             ncluster = 2:5,
             seg.var = "Sr_Ca_ratio")}, 
             error = function(e) {
    message(paste("ERROR", individ, ":", e$message))
    return(NULL)
   })
  
 ### If the loop failed:
  if (is.null(mode_segclust)) next
  
    # Extract the number of selected clusters according to BIC
  ncluster_BIC <- mode_segclust$ncluster.BIC[[1]]
  
  # Get the mean and sd value for ind and cluster
  state_info <- states(mode_segclust, ncluster = ncluster_BIC)$states
  
  # extract the individualID nr, sampling site and the number of clusters 
  state_info$Individ <- individ
  state_info$Water_area <- unique(individ_data$Water_area)[1]
  state_info$ncluster_BIC <- ncluster_BIC
  
  # Save the results as a list
  result_list[[individ]] <- state_info

  ## Name each plot
  Water_name <- gsub("[^A-Za-z0-9_åäöÅÄÖ]","_", unique(individ_data$Water_area)[1])
  Individual_ID <- gsub("[^A-Za-z0-9_åäöÅÄÖ]", "_", individ)
  
## Plot 1: Segmentation and clustering Sr.Ca transect
 png(filename = file.path(output_dir, paste0(Individual_ID, "_",  Water_name, "_segments.png")),
    width = 1000, height = 600)
  plot(mode_segclust)  # BIC-selected output
  print(plot(mode_segclust))
  dev.off()
  
}
### Combine loop output data:

### Make one data-frame wiht all the results (1 row per individual):
Results_threshold_FE_EPMA <- do.call(rbind, result_list)
rownames(Results_threshold_FE_EPMA) <- NULL

### Results
head(Results_threshold_FE_EPMA) # State = Cluster

## Save dataframe for upcoming scripts
saveRDS(Results_threshold_FE_EPMA, file = "R_Scripts_and_Generated_data_files/Results_threshold_FE_EPMA.Rds")

## Extract the mean and SD per cluster:
Summary_per_Cluster_ref_FE_EPMA<-Results_threshold_FE_EPMA %>% 
  group_by(state) %>% 
  summarise(Mean_Min=min(mu.Sr_Ca_ratio, na.rm = TRUE),
            Mean_Max=max(mu.Sr_Ca_ratio, na.rm = TRUE),
            n= n(),
            sd_Min= min(sd.Sr_Ca_ratio, na.rm = TRUE),
            sd_Max= max(sd.Sr_Ca_ratio, na.rm = TRUE))

Summary_per_Cluster_ref_FE_EPMA

saveRDS(Summary_per_Cluster_ref_FE_EPMA, file = "R_Scripts_and_Generated_data_files/Summary_threshold_values_FE_EPMA.Rds")

# ------ Done with deriving Threshold values for FE-EPMA data ------ #

## LA-ICP-MS threshold values:
# Load data from Script 1_Input_data.R:

LA_ICP_MS_KUL_Sotebasen = readRDS("R_Scripts_and_Generated_data_files/LA_ICP_MS_KUL_Sotebasen.Rds")

LA_ICP_MS_KUL_Sotebasen<-LA_ICP_MS_KUL_Sotebasen |> 
 group_by(IndividNr) |> 
   mutate(Ca_ppm_calc=as.numeric(as.character(Ca43_CPS, na.rm=T))*380000/mean(as.numeric(as.character(Ca43_CPS, na.rm=T))),#mmol mol–1 
          Sr_Ca_ppm=(as.numeric(as.character(Sr_ppm_m88, na.rm=T))*1000)/Ca_ppm_calc) |> 
  filter(Ca_ppm_calc > (mean(Ca_ppm_calc) - 1.25 * sd(Ca_ppm_calc)), 
         Ca_ppm_calc < (mean(Ca_ppm_calc) + 1.25 * sd(Ca_ppm_calc))) |>
  ungroup()

# Subset naturally recruited reference eel:
Reference_eel_LA_ICP_MS <- LA_ICP_MS_KUL_Sotebasen |> 
  filter(IndividNr %in% c("50463", "50473", "50479", "51097", "51103", # Enningsdalälven
                          "51108", "51114", "51116", "51118", "51120", # Enningsdalälven
                          "2019_10", "2015_110", "2019_116", "2019_119", "2021_191", # Stenungssund #2019_110 is wierd!
                          "2019_136", "2019_147", "2019_167", "2019_171", "2019_179", # Stenungssund
                          "41252", "44420", "44425", "44427", "44429", # Göta älv
                          "44441", "44449", "44450", "46705", "50318", # Göta älv
                          "47882", "47892", "48474", "49338", "49339", # Viskan
                          "49342", "49343", "49946", "49947", "49948", # Viskan
                          "48044", "48734", "48736", "48737", "48738", # Kävlingeån
                          "49286", "48742", "49330", "50173", "50174", # Kävlingeån
                          "47995", "47998", "47999", "48002", "48003", # Nybroån
                          "48489", "48490", "48492", "50198", "50200", # Nybroån
                          "51473", "51478", "51482", "51483", "51484", # Karlshman, Blekinge
                          "51494", "51496", "51500", "51502", "51489", # Karlshamn, Blekinge
                          "48767", "51214", "51186", "46596", "46615", # Eight clear natural recruits Mälaren
                          "51245", "50848", "50161"))

head(Reference_eel_LA_ICP_MS)

# Save dataframe for Script 10_Suppl_Info.R
saveRDS(Reference_eel_LA_ICP_MS, file = "R_scripts_and_Generated_data_files/Reference_eel_LA_ICP_MS.Rds")

### Make a loop for the segmentation analyses:
# Sr:Ca ratio loop:

output_dir <- "Output/Reference_eel_LA_ICP_MS_segclust" # Make a folder to save the plots in
dir.create(output_dir, showWarnings = FALSE)

result_list <- list()

individual <- unique(Reference_eel_LA_ICP_MS$IndividNr)

# Run the segmentation for each individual
for (individ in individual) {
  
  ##  Filter out data and make a smaller df to speed it up:
   individ_data_lund <- Reference_eel_LA_ICP_MS %>%
    filter(IndividNr == individ) %>%
     select(IndividNr, Sr_Ca_ppm, Ca_ppm_calc, ElapsedTime_s, Water_area) %>%
     filter(!is.na(Sr_Ca_ppm)) %>% # Remove NAs if present
     mutate(ElapsedTime_s = as.numeric(ElapsedTime_s)) %>%    # convert to numeric
     arrange(ElapsedTime_s) %>%   #Sort
     as.data.frame()
  
 # # skip if to little data
  if (nrow(individ_data_lund) < 10) {
    message(paste("Skipped", individ, "- too few datapoints"))
    next
   }

  #### Run the segclust modell:
  mode_segclust <- tryCatch({
   segclust(individ_data_lund,
             Kmax = 20,
             lmin = 10,
             ncluster = 2:5,
             seg.var = "Sr_Ca_ppm")
   }, error = function(e) {
    message(paste("Error for", individ, ":", e$message))
    return(NULL)
   })
  
 ### If the model failed:
  if (is.null(mode_segclust)) next
  
    # Extract number of BIC cluster per ind:
  ncluster_BIC <- mode_segclust$ncluster.BIC[[1]]
  
  # Get the state for each cluster:
  state_info <- states(mode_segclust, ncluster = ncluster_BIC)$states
  
  # Add IndividNr to the data:
  state_info$Individ <- individ
  state_info$Water_area <- unique(individ_data_lund$Water_area)[1]
  state_info$ncluster_BIC <- ncluster_BIC
  
  # Save as a list
  result_list[[individ]] <- state_info

  ## Name the plots so it is easy to find them in the folder:
  Water_name <- gsub("[^A-Za-z0-9_åäöÅÄÖ]","_", unique(individ_data_lund$Water_area)[1])
  individual_ID <- gsub("[^A-Za-z0-9_åäöÅÄÖ]", "_", individ)
  
## Plot 1: Segclust plot
 png(filename = file.path(output_dir, paste0(individual_ID, "_", Water_name, "_segments.png")),
      width = 1000, height = 600)
 plot(mode_segclust)  # BIC-plot
  print(plot(mode_segclust))
  dev.off() 
}

### Combine all output data:
### Merge into one Df
result_LA_ICP_MS_Threshold <- do.call(rbind, result_list)
rownames(result_LA_ICP_MS_Threshold) <- NULL

### Check the results:
head(result_LA_ICP_MS_Threshold) # state = Cluster

## Save data-file for upcoming scripts
saveRDS(result_LA_ICP_MS_Threshold, file = "R_Scripts_and_Generated_data_files/Results_threshold_LA_ICP_MS.Rds")

## Threshold-values for freshwater caught eel:
result_LA_ICP_MS_Threshold_values_summary<-result_LA_ICP_MS_Threshold |> 
  mutate(Catch_habitat=ifelse(Water_area=="Östersjön (hela) ICES SD 22-32","Coastal_Baltic",
                        ifelse(Water_area=="Västerhavet (hela) ICES SD 20-21","Coastal_Skagerrak", "Freshwater"))) |> 
  filter(ncluster_BIC == 5) |>
     group_by(Catch_habitat, state) %>% 
     summarise(mu_5 = quantile(mu.Sr_Ca_ppm, 0.05, na.rm = TRUE),
               mu_95 = quantile(mu.Sr_Ca_ppm, 0.95, na.rm = TRUE),
               n = n(),
               sd_5 = quantile(sd.Sr_Ca_ppm, 0.05, na.rm = TRUE),
               sd_95 = quantile(sd.Sr_Ca_ppm, 0.95, na.rm = TRUE))|> 
  mutate(Analutical_method="LA_ICP_MS")

result_LA_ICP_MS_Threshold_values_summary

saveRDS(result_LA_ICP_MS_Threshold_values_summary, file = "R_Scripts_and_Generated_data_files/Summary_threshold_values_LA_ICP_MS.Rds")

# ------ Table 2, main text -------- #

## Load the generated .RDS output table if not running the script from start again:

Summary_per_Cluster_ref_FE_EPMA = readRDS("R_Scripts_and_Generated_data_files/Summary_threshold_values_FE_EPMA.Rds")

result_LA_ICP_MS_Threshold_values_summary = readRDS("R_Scripts_and_Generated_data_files/Summary_threshold_values_LA_ICP_MS.Rds")

print(Summary_per_Cluster_ref_FE_EPMA, n=Inf) # Three clusters
print(result_LA_ICP_MS_Threshold_values_summary, n=Inf) # 5 Clusters, select the combination of clusters and known sampling sites for the threshold values

# Pull out the min and max values, merge and create the table for the main text:

## FE-EPMA:

FE_EPMA_Table2 <- Summary_per_Cluster_ref_FE_EPMA %>%
  mutate(Growth_habitat = c("Freshwater",
                             "Coastal/Brackish",
                             "Marine"), # Change the name of each row
         `Sr:Ca threshold values` = c(paste0("0 - ", sprintf("%.2f", Mean_Max[1])), # Upper limit of cluster 1 corresponding to freshwater
                                      paste0(">", sprintf("%.2f", Mean_Max[1], 2), " - ", sprintf("%.2f", Mean_Max[2], 2)), # Upper limit of cluster 1 corresponding to freshwater and upper limit for state 2 = Coastal/Brackish
                                      paste0(">", sprintf("%.2f", Mean_Max[2], 2))), # Upper limit for state 2 (Coastal/Brackish)
         Analytical_method = "FE-EPMA") %>%
  select(Growth_habitat, `Sr:Ca threshold values`, Analytical_method)

FE_EPMA_Table2

## LA-ICP-MS:

Freshwater_upper <- result_LA_ICP_MS_Threshold_values_summary %>%
  filter(Catch_habitat == "Freshwater", state == 1) %>%
  pull(mu_95)
Freshwater_upper # should be 1.67 (1.667376)

Estuary_upper <- result_LA_ICP_MS_Threshold_values_summary %>%
  filter(Catch_habitat == "Freshwater", state == 2) %>%
  pull(mu_95)
Estuary_upper # Should be 2.60 (2.596779)

Marine_lower <- result_LA_ICP_MS_Threshold_values_summary %>%
  filter(Catch_habitat == "Coastal_Baltic", state == 5) %>%
  pull(mu_5)
Marine_lower # Should be 9.13 (9.130356)

# Combine:
LA_ICP_MS_Table2 <- tibble(Growth_habitat = c("Freshwater",
                                              "Estuary",
                                              "Coastal/Brackish",
                                              "Marine"),
                           `Sr:Ca threshold values` = c(paste0("0 - ", sprintf("%.2f", Freshwater_upper)),
                                                       paste0(">", sprintf("%.2f", Freshwater_upper), " - ", sprintf("%.2f", Estuary_upper)),
                                                       paste0(">", sprintf("%.2f", Estuary_upper), " - ", sprintf("%.2f", Marine_lower)),
                                                       paste0(">", sprintf("%.2f", Marine_lower))),
                           Analytical_method = "LA-ICP-MS")

LA_ICP_MS_Table2

## Combine:

Table2 <- bind_rows(FE_EPMA_Table2, LA_ICP_MS_Table2)

names(Table2)

# Create the flextable for the main text:

Table_2_ft <- Table2 |>
  rename(`Growth habitat` = Growth_habitat,
          `Analytical method` = Analytical_method) |> 
  flextable() |>
  theme_booktabs() |>
  flextable::align(j = "Growth habitat", align = "left", part = "all") |>
  flextable::align(j = c("Sr:Ca threshold values", "Analytical method"),
         align = "right", part = "all") |>
  autofit()

Table_2_ft

# Save the table in the Output folder:
# As an image:
save_as_image(Table_2_ft, path = "Output/Table_2.png")

# As a word-file:
save_as_docx(Table_2_ft, path = "Output/Table_2.docx")

# --------- End of Script 3_Threshold_values.R ------------- #