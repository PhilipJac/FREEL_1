#################
# Jacobson et al.
# Navigating the salinity gradient: Individual variation in habitat use and migration of European eel revealed by otolith microchemistry 
#################

## Last update: 2026-06-16

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

#### FE-EPMA:

# Select reference eels from additional sites:

## FE-EPMA threshold values:
# Bolmen
# Roxen
# Vänern
# Ymsen
# Norrtälje
# Kalmarsund
# Höganäs kullen

set.seed(345) # to make the random selection reproducable:

sites_to_sample <- c(
  "Bolmen - Ingen uppgift/ej definierad",
  "Roxen - Ingen uppgift/ej definierad",
  "Vänern - Ingen uppgift/ej definierad",
  "Ymsen - Ingen uppgift/ej definierad",
  "Hjälmaren - Ingen uppgift/ej definierad",
  "Mälaren - Galten",
  "Mälaren - Stallarholmen",
  "Östersjön (hela) ICES SD 22-32 - Kalmarsund",
  "Östersjön (hela) ICES SD 22-32 - Norrtälje-Bylehamn",
  "Västerhavet (hela) ICES SD 20-21 - Höganäs-Kullen"
)

selected_ids <- FE_EPMA_Sotebasen |>
  mutate(Water_site = paste(Water_area, Site, sep = " - ")) |>
  distinct(Water_site, IndividNr) |>
  filter(Water_site %in% sites_to_sample) |>
  group_by(Water_site) |>
  group_modify(~ slice_sample(.x, n = min(10, nrow(.x)))) |>
  ungroup()

FE_EPMA_Sotebasen_sample_ref <- FE_EPMA_Sotebasen |>
  mutate(Water_site = paste(Water_area, Site, sep = " - ")) |>
  semi_join(
    selected_ids,
    by = c("Water_site", "IndividNr")
  )
head(FE_EPMA_Sotebasen_sample_ref) ## There we go!

### Make a loop for the segmentation analyses:
result_list <- list()

## Create a folder to save the segmentation and clustering output in the output folder:
output_dir <- "Output/Reference_eel_FE_EPMA_segclust_sensitivity"
dir.create(output_dir, showWarnings = FALSE)

individual <- unique(FE_EPMA_Sotebasen_sample_ref$IndividNr)

# Run the loop for every individual:
for (individ in individual) {
  
  ## Calculate the Sr:Ca ratio for every individual, remove Na-values and run the loop
  individ_data <- FE_EPMA_Sotebasen_sample_ref %>%
    filter(IndividNr == individ) %>%
    select(IndividNr, COMMENT, Sr, Ca, Date, Water_site) %>%
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
  state_info$Water_site <- unique(individ_data$Water_site)[1]
  state_info$ncluster_BIC <- ncluster_BIC
  
  # Save the results as a list
  result_list[[individ]] <- state_info

  ## Name each plot
  Water_name <- gsub("[^A-Za-z0-9_åäöÅÄÖ]","_", unique(individ_data$Water_site)[1])
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
Results_threshold_FE_EPMA_sensitivity <- do.call(rbind, result_list)
rownames(Results_threshold_FE_EPMA_sensitivity) <- NULL

### Results
head(Results_threshold_FE_EPMA_sensitivity) # State = Cluster

## Extract the mean and SD per cluster:
Summary_per_Cluster_ref_FE_EPMA_sensitivity<-Results_threshold_FE_EPMA_sensitivity %>% 
  filter(ncluster_BIC=="3") %>% 
  group_by(Water_site, state) %>% 
  summarise(Mean_Min=min(mu.Sr_Ca_ratio, na.rm = TRUE),
            mu_5 = quantile(mu.Sr_Ca_ratio, 0.05, na.rm = TRUE),
            mu_95 = quantile(mu.Sr_Ca_ratio, 0.95, na.rm = TRUE),
            Mean_Max=max(mu.Sr_Ca_ratio, na.rm = TRUE),
            n= n(),
            sd_Min= min(sd.Sr_Ca_ratio, na.rm = TRUE),
            sd_5 = quantile(sd.Sr_Ca_ratio, 0.05, na.rm = TRUE),
            sd_95 = quantile(sd.Sr_Ca_ratio, 0.95, na.rm = TRUE),
            sd_Max= max(sd.Sr_Ca_ratio, na.rm = TRUE))

print(Summary_per_Cluster_ref_FE_EPMA_sensitivity, n=Inf)

## Make one estimate for freahwater, one for coastal:

## Extract the mean and SD per cluster:
Summary_per_Cluster_and_hab_ref_FE_EPMA_sensitivity<-Results_threshold_FE_EPMA_sensitivity %>% 
  filter(ncluster_BIC == "3") %>%
  mutate(Sampling_hab = case_when(Water_site %in% c("Östersjön (hela) ICES SD 22-32 - Norrtälje-Bylehamn",
                                                    "Östersjön (hela) ICES SD 22-32 - Kalmarsund") ~ "Baltic_sea",
                                  Water_site == "Västerhavet (hela) ICES SD 20-21 - Höganäs-Kullen" ~ "West_coast",
                                  TRUE ~ "Freshwater")) %>%
  group_by(Sampling_hab, state) %>% 
     summarise(Mean_Min=min(mu.Sr_Ca_ratio, na.rm = TRUE),
               Mean_Max=max(mu.Sr_Ca_ratio, na.rm = TRUE),
               n = n(),
               sd_Min= min(sd.Sr_Ca_ratio, na.rm = TRUE),
               sd_Max= max(sd.Sr_Ca_ratio, na.rm = TRUE))

print(Summary_per_Cluster_and_hab_ref_FE_EPMA_sensitivity, n=Inf)

#### Individual habitat assignment for eel with FE-EPMA otolith chemistry data #####

### Calculate distance from core for each datapoint using the X and Y coordinates for each individual:
FE_EPMA_Sotebasen_Distance <- FE_EPMA_Sotebasen %>%
  group_by(IndividNr) %>%
  arrange(Date, .by_group = TRUE) %>%
  mutate(Sr_Ca_ratio = as.numeric(as.character(Sr)) / as.numeric(as.character(Ca)) * 1000,
         Step_dist = sqrt((X - lag(X))^2 + (Y - lag(Y))^2),  # Distance between consecutive points -> distance between measure-points
         Distance_from_start = cumsum(replace_na(Step_dist, 0))) %>%  # Cumulative distance from first point which is set to 0 (the center of the core), which is the start of the transect
  ungroup()

# Check the distance calculation:
# Distance distriburution
FE_EPMA_Sotebasen_Distance |> 
  group_by(IndividNr) |> 
  summarise(Distance_core_to_edge=max(Distance_from_start)) |> # Get the maximum value for each 
ggplot()+
  geom_histogram(aes(x=Distance_core_to_edge))+ # Looks good
  theme_bw(base_size=22)+
  xlab("Distance from otolith core (µm)")+
  ylab("Amount")

# An example of an Sr:Ca transect:
FE_EPMA_Sotebasen_Distance |> 
  filter(IndividNr %in% c("31643")) |> 
    ggplot(aes(y = Sr_Ca_ratio,
                   x = Distance_from_start)) +
       geom_line()+
       facet_wrap(~IndividNr*Water_area*filename) +
      ylab("(Sr/Ca)*1000") +
      xlab("Distance from core")+
 scale_x_continuous(breaks = seq(0, 2000, by = 200), limits = c(0, 2000)) +
  theme_bw(base_size=22) # Looks good!

## Assign a habitat for each data point in each transect:

# >20 cm filter and transect longer than 513 micrometers and no Sr-mark partially or fully > 513 micrometer from the otolith core:
FE_EPMA_Sotebasen_Distance_Habitat_filter_20cm_513mm_No_Sr_Sensitivity<-FE_EPMA_Sotebasen_Distance |> 
  filter(Total_length_mm>=200, # Only include eel longer than 20 cm
         Distance_from_start>=513) |> # Remove the first section of the Sr:Ca transect, corresponding to the arrival from the Sargasso sea
  filter(!IndividNr %in% c("34151", "34152", "34153", "34154", "34155", 
                           "34156", "34157", "34158", "34159", "34160", 
                           "34161", "34162", "34163", "34164", "34685", "34686")) |>
   filter(!is.na(Sr_Ca_ratio)) |> 
     mutate(Habitat_use = case_when(Sr_Ca_ratio <= 0.94 ~ "Freshwater",
                                   Sr_Ca_ratio > 0.94 & Sr_Ca_ratio <= 5.98 ~ "Coastal/Brackish",
                                   Sr_Ca_ratio > 5.98 ~"Marine/Core"))

n_distinct(FE_EPMA_Sotebasen_Distance_Habitat_filter_20cm_513mm_No_Sr_Sensitivity$IndividNr) # 1877 inds
# Rename the data frame:
FE_EPMA_Sotebasen_Distance_Habitat_Sensitivity<-FE_EPMA_Sotebasen_Distance_Habitat_filter_20cm_513mm_No_Sr_Sensitivity

## Assign life-history to each individual:

FE_EPMA_Sotebasen_Life_history_Sens<-FE_EPMA_Sotebasen_Distance_Habitat_Sensitivity |>
  group_by(Water_area, Site, IndividNr) |>
  summarise(has_fresh = any(Habitat_use == "Freshwater"),
            has_coastal = any(Habitat_use == "Coastal/Brackish"),
            has_marine = any(Habitat_use == "Marine/Core"),
            .groups = "drop") |>
  mutate(Life_history = case_when(has_fresh & !has_coastal & !has_marine ~ "Freshwater", # Only freshwater Sr:Ca values
                                  has_fresh & (has_coastal | has_marine) ~ "Habitat shifter", # Both fresheater and coastal or marine signal
                                  !has_fresh ~ "Coastal", # No freshwater signal
                                   TRUE ~ NA_character_))

n_distinct(FE_EPMA_Sotebasen_Life_history_Sens$IndividNr) # Should be 1877 inds

## Summary table major water areas:
FE_EPMA_Sotebasen_Life_history_Sens<-FE_EPMA_Sotebasen_Life_history_Sens|> 
  mutate(Catch_area_habitat=ifelse(Water_area %in% c("Bolmen", "Vänern", "Ymsen", "Ätran"), "Freshwater_Kattegat", # Freshwater sites connected to Kattegat 
                      ifelse(Water_area %in% c("Hjälmaren", "Mälaren", "Ottsjön", "Ringsjöarna, Skåne", 
                                               "Roxen", "Ångermanälven", "Åssjön"), "Freshwater_Baltic_Sea", # Freshwater sites connected to the Baltic Sea
                      ifelse(Water_area == "Västerhavet (hela) ICES SD 20-21", "Coastal_Kattegat",    # Coastal sites, Kattegat
                      ifelse(Water_area == "Östersjön (hela) ICES SD 22-32", "Coastal_Baltic_Sea",    # Coastal sites, Baltic Sea
                                           "No_classification"))))) |> 
  group_by(Catch_area_habitat, Life_history) |> 
  summarise(Individuals_n=n()) |> 
  group_by(Catch_area_habitat) |> 
  mutate(Total_amount_of_inds=sum(Individuals_n),
         Percent_per_life_history=round((Individuals_n/Total_amount_of_inds)*100, 1),
         Analythical_method="FE-EPMA")

flextable(FE_EPMA_Sotebasen_Life_history_Sens) ## There we go!


## Done with FE-EPMA

## LA-ICP-MS:

### Randomize the individual selection from each sample size for LA-ICP-MS:
LA_ICP_MS_KUL_Sotebasen = readRDS("R_Scripts_and_Generated_data_files/LA_ICP_MS_KUL_Sotebasen.Rds")

head(LA_ICP_MS_KUL_Sotebasen)

LA_ICP_MS_KUL_Sotebasen<-LA_ICP_MS_KUL_Sotebasen |> 
 group_by(IndividNr) |> 
   mutate(Ca_ppm_calc=as.numeric(as.character(Ca43_CPS, na.rm=T))*380000/mean(as.numeric(as.character(Ca43_CPS, na.rm=T))),#mmol mol–1 
          Sr_Ca_ppm=(as.numeric(as.character(Sr_ppm_m88, na.rm=T))*1000)/Ca_ppm_calc) |> 
  filter(Ca_ppm_calc > (mean(Ca_ppm_calc) - 1.25 * sd(Ca_ppm_calc)), 
         Ca_ppm_calc < (mean(Ca_ppm_calc) + 1.25 * sd(Ca_ppm_calc))) |>
  ungroup()

LA_ICP_MS_KUL_Sotebasen |> 
 group_by(Water_area, Site) |>
  summarise(n=n_distinct(IndividNr)) |>
  print(n=Inf)

## Make  anew column combining Water area and site:
LA_ICP_MS_KUL_Sotebasen_sens<-LA_ICP_MS_KUL_Sotebasen |>
mutate(Water_site = paste(Water_area, Site, sep = " - "))

sort(unique(LA_ICP_MS_KUL_Sotebasen_sens$Water_site))

# We want the same sites as in the main text, just randomize the individuals from each site:
sites_to_sample_LA_ICP_MS <- c("Enningdalsälven - Kynne älv-Sundshult",
                              "Västerhavet (hela) ICES SD 20-21 - Stenungsund",
                              "Göta älv - Olidans kraftverk",
                              "Viskan - Kullagård ål",
                              "Nybroån - Övningsfältet",
                              "Östersjön (hela) ICES SD 22-32 - Karlshamn - Veggabadet",
                              "Mälaren - Sotholmen A")

selected_inds_sens_LA_ICP_MS <- LA_ICP_MS_KUL_Sotebasen_sens  |>
  distinct(Water_site, IndividNr) |>
  filter(Water_site %in% sites_to_sample_LA_ICP_MS) |>
  group_by(Water_site) |>
  group_modify(~ slice_sample(.x, n = min(10, nrow(.x)))) |>
  ungroup() # A list of randomly selected inds:

LA_ICP_MS_KUL_Sotebasen_sens_samples <- LA_ICP_MS_KUL_Sotebasen_sens |>
  mutate(Water_site = paste(Water_area, Site, sep = " - ")) |>
  semi_join(
    selected_inds_sens_LA_ICP_MS,
    by = c("Water_site", "IndividNr")
  )
head(LA_ICP_MS_KUL_Sotebasen_sens_samples) ## There we go!

### Make a loop for the segmentation analyses:
result_list <- list()

### Make a loop for the segmentation analyses:
# Sr:Ca ratio loop:

output_dir <- "Output/Reference_eel_LA_ICP_MS_segclust_sensitivty" # Make a folder to save the plots in
dir.create(output_dir, showWarnings = FALSE)

result_list <- list()

individual <- unique(LA_ICP_MS_KUL_Sotebasen_sens_samples$IndividNr)

# Run the segmentation for each individual
for (individ in individual) {
  
  ##  Filter out data and make a smaller df to speed it up:
   individ_data_lund <- LA_ICP_MS_KUL_Sotebasen_sens_samples %>%
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
 plot(mode_segclust)  # automatiskt BIC-val av ncluster/nseg
  print(plot(mode_segclust))
  dev.off() 
}

### Combine all output data:
### Merge into one Df
result_LA_ICP_MS_Threshold_sensitivity <- do.call(rbind, result_list)
rownames(result_LA_ICP_MS_Threshold_sensitivity) <- NULL

### Check the results:
head(result_LA_ICP_MS_Threshold_sensitivity) # state = Cluster

result_LA_ICP_MS_Threshold_sensitivity|>
   filter(ncluster_BIC == 5)|>
  group_by(Water_area)|>
  summarise(n=(n_distinct(Individ)))

#saveRDS(result_LA_ICP_MS_Threshold_sensitivity, file = "R_scripts_and_Generated_data_files/result_LA_ICP_MS_Threshold_sensitivity.Rds")
result_LA_ICP_MS_Threshold_sensitivity = readRDS("R_scripts_and_Generated_data_files/result_LA_ICP_MS_Threshold_sensitivity.Rds")

## Threshold-values for freshwater caught eel:

set.seed(789)
result_LA_ICP_MS_Threshold_values_sensitivity_summary<-result_LA_ICP_MS_Threshold_sensitivity |> 
  mutate(Catch_habitat=ifelse(Water_area=="Östersjön (hela) ICES SD 22-32","Coastal_Baltic",
                        ifelse(Water_area=="Västerhavet (hela) ICES SD 20-21","Coastal_Skagerrak", "Freshwater"))) |> 
  filter(ncluster_BIC == 5) |>
  filter(!Water_area=="Göta älv") |>
  group_by(Water_area, state) |> # group by sampling site/Waterarea:
  slice_sample(n = 10) |> # randomly select 10 inds
  ungroup() |> 
     group_by(Catch_habitat, state) |> 
     summarise(mu_5 = quantile(mu.Sr_Ca_ppm, 0.05, na.rm = TRUE),
               mu_95 = quantile(mu.Sr_Ca_ppm, 0.95, na.rm = TRUE),
               n = n(),
               sd_5 = quantile(sd.Sr_Ca_ppm, 0.05, na.rm = TRUE),
               sd_95 = quantile(sd.Sr_Ca_ppm, 0.95, na.rm = TRUE))|> 
  mutate(Analutical_method="LA_ICP_MS")

result_LA_ICP_MS_Threshold_values_sensitivity_summary

#saveRDS(result_LA_ICP_MS_Threshold_values_sensitivity_summary, file = "R_scripts_and_Generated_data_files/result_LA_ICP_MS_Threshold_values_sensitivity_summary.Rds")

result_LA_ICP_MS_Threshold_values_sensitivity_summary = readRDS("R_scripts_and_Generated_data_files/result_LA_ICP_MS_Threshold_values_sensitivity_summary.Rds")

result_LA_ICP_MS_Threshold_values_sensitivity_summary 

## Habitat assignment:

LA_ICP_MS_KUL_Sotebasen_Distance=readRDS(file="R_Scripts_and_Generated_data_files/LA_ICP_MS_KUL_Sotebasen_Distance.Rds")

#Step 8. Assign habitat for all Sr:Ca value, excluding all points between 0-513 micrometer for all eel above 20 cm in length and without any Sr-mark affecting the habitat classification:

LA_ICP_MS_KUL_Sotebasen_Distance_Habitat_sensitivity<-LA_ICP_MS_KUL_Sotebasen_Distance |>
  filter(Total_length_mm>=200) |> # Only include eel longer than 20 cm
  filter(Distance_from_core>=513) |> # Remove the first section of the Sr:Ca transect, corresponding to the arrival from the Sargasso sea
   filter(!IndividNr %in% c("47849", "47851", "47852", "47853", "47854", "47855", "47856", "47857", "47858"))|> # Remove the first section of the Sr:Ca transect, corresponding to the arrival from the Sargasso sea
 mutate(Habitat_use = case_when(Sr_Ca_ppm >=0 & Sr_Ca_ppm <= 2.53 ~ "Freshwater",
                                   Sr_Ca_ppm > 2.53 & Sr_Ca_ppm <= 2.91  ~ "Estuary", # Lowest of state 1 coastal inds
                                   Sr_Ca_ppm > 2.91 & Sr_Ca_ppm <= 9.10 ~ "Coastal/Brackish",
                                   Sr_Ca_ppm > 9.10 ~"Marine/Core"))


n_distinct(LA_ICP_MS_KUL_Sotebasen_Distance_Habitat_sensitivity$IndividNr) # A total of 1036 inds with a transect longer than 513 micrometers and a body length >20 cm
## Shorten the df name:
LA_ICP_MS_KUL_Sotebasen_Distance_Habitat_sens<-LA_ICP_MS_KUL_Sotebasen_Distance_Habitat_sensitivity
n_distinct(LA_ICP_MS_KUL_Sotebasen_Distance_Habitat_sens$IndividNr) # Should be 1036 inds.


## Assign life-history to each individual:
LA_ICP_MS_KUL_Sotebasen_Life_history_sens<-LA_ICP_MS_KUL_Sotebasen_Distance_Habitat_sens|>
  group_by(Water_area, Site, IndividNr) |>
  summarise(has_fresh = any(Habitat_use == "Freshwater", na.rm=T),
            has_estuary = any(Habitat_use == "Estuary", na.rm=T),
            has_brackish = any(Habitat_use == "Coastal/Brackish", na.rm=T),
            has_marine = any(Habitat_use == "Marine/Core", na.rm=T),
           .groups = "drop") |>
 mutate(Life_history = case_when(has_fresh & !has_estuary & !has_brackish & !has_marine ~ "Freshwater", # Only freshwater
                                 has_fresh & has_estuary & !has_brackish & !has_marine ~ "Freshwater/Estuary",  # Only freshwater + estuary
                                 has_fresh & (has_brackish | has_marine) ~ "Habitat shifter",   # Freshwater + Coastal and/or marine)
                                 has_estuary & !has_fresh & !has_brackish & !has_marine ~ "Estuary", # Only estuary (no freshwater, no coastal)
                                 !has_fresh ~ "Coastal", # No freswhater singal (all estuarine, coastal/brackish and marine combos)
                                 TRUE ~ NA_character_))

LA_ICP_MS_KUL_Sotebasen_Life_history_sens |> 
  group_by(Life_history) |> 
   summarise(n=n()) |> 
   mutate(sum=sum(n)) # 1036 inds, all have an assigned life-hisotry

head(LA_ICP_MS_KUL_Sotebasen_Life_history_sens)

# Make summary table for the four major geographical regions:
LA_ICP_MS_KUL_Sotebasen_Life_history_Major_regions<-LA_ICP_MS_KUL_Sotebasen_Life_history_sens  |> 
  mutate(Catch_area_habitat=ifelse(Water_area %in% c("Aröd å", "Bolmen", "Enningdalsälven", "Göta älv", 
                                                      "Knapebäcken", "Kvarnabäcken", "Löftaån", "Rönne å", 
                                                      "Strömsån", "Viskan", "Ymsen", "Örekilsälven"), "Freshwater_Kattegat", # Freshwater sites connected to Kattegat 
                            ifelse(Water_area %in% c("Dalälven", "Dybäcksån", "Helge å", "Höje å", 
                                                     "Kävlingeån", "Motala Ström", "Mälaren", "Nybroån", 
                                                      "Nyköpingsån", "Skivarpsån", "Ången"), "Freshwater_Baltic_Sea", # Freshwater sites connected to the Baltic Sea
                            ifelse(Water_area == "Västerhavet (hela) ICES SD 20-21", "Coastal_Kattegat",    # Coastal sites, Kattegat
                            ifelse(Water_area == "Östersjön (hela) ICES SD 22-32", "Coastal_Baltic_Sea",    # Coastal sites, Baltic Sea
                                                  "No_classification"))))) |> 
  group_by(Catch_area_habitat, Life_history) |> 
  summarise(Individuals_n=n()) |> 
  group_by(Catch_area_habitat) |> 
  mutate(Total_amount_of_inds=sum(Individuals_n),
         Percent_per_life_history=round((Individuals_n/Total_amount_of_inds)*100, 1),
         Analythical_method="LA-ICP-MS")

head(LA_ICP_MS_KUL_Sotebasen_Life_history_Major_regions)

flextable(LA_ICP_MS_KUL_Sotebasen_Life_history_Major_regions)

sum(LA_ICP_MS_KUL_Sotebasen_Life_history_Major_regions$Individuals_n)

