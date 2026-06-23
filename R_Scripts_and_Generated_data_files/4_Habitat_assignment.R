#################
# Jacobson et al.
# Navigating the salinity gradient: Individual variation in habitat use and migration of European eel revealed by otolith microchemistry 
#################

## Last update: 2026-06-16

#### R-Script: 4_Habitat_assignment.R

## Load packages:

library(tidyverse)
library(zoo)
library(flextable)

# ---------- LA-iCP-MS data ----------------- #

### Load data generated in the script "1_Input_data.R":

# 1. LA-ICP-MS otolith chemistry transect data combined with individual level data from the databases Sötebasen and KUL:
LA_ICP_MS_KUL_Sotebasen = readRDS("R_Scripts_and_Generated_data_files/LA_ICP_MS_KUL_Sotebasen.Rds") # Containing 29 columns with 950477 rows of data 

## Load spot size and scanning speed data for the LA_ICP-MS data:
Scan_speed_spot_size_LA_ICP_MS<-read.delim("Input_data/Scan_speed_Spot_size_LA_ICP_MS.csv", sep=";") # All spot-size and scanning speed data for LA-ICP-MS analyses

#### Habitat assigment LA-ICP-MS data:

## Identify the center of the otolith core and the associated costal signal corresponding 
#  to the birth and arrival of eel from the Sargasso sea to Skagerrak and Kattegat area.

# Start with spot-size and scan-speed data:
head(Scan_speed_spot_size_LA_ICP_MS)

# Clean up the table column names:
Scan_speed_spot_size_LA_ICP_MS <- Scan_speed_spot_size_LA_ICP_MS |> 
  rename(Plate_number=Plate.number,
         IndividNr=Individ..Nr.,
         Scan_speed_micrometer_sec=Scan.speed..µm.s.,
         Spot_size=Spot.size) |> 
  select(Plate_number, IndividNr, Scan_speed_micrometer_sec, Spot_size)

# Merge the Scan speed data with the LA_ICP_MS data:
LA_ICP_MS_KUL_Sotebasen_Scan_speed <- merge(x = LA_ICP_MS_KUL_Sotebasen, 
                                            y = Scan_speed_spot_size_LA_ICP_MS, 
                                            by = "IndividNr", all.x = TRUE)

head(LA_ICP_MS_KUL_Sotebasen_Scan_speed)

## Step 1 and 2. Calculate Sr:Ca values for each individual and make a Sr:Ca smoother:                                 
LA_ICP_MS_KUL_Sotebasen_Scan_speed_core<-LA_ICP_MS_KUL_Sotebasen_Scan_speed |> 
 group_by(IndividNr) |> 
  arrange(ElapsedTime_s, .by_group = TRUE) |> # make sure the data is arrange according to elapsed time for each individual, otherwise, the rolling mean will not work
   mutate(Ca_ppm_calc=as.numeric(as.character(Ca43_CPS, na.rm=T))*380000/mean(as.numeric(as.character(Ca43_CPS, na.rm=T))),#mmol mol–1 
          Sr_Ca_ppm=(as.numeric(as.character(Sr_ppm_m88, na.rm=T))*1000)/Ca_ppm_calc) |> 
  filter(Ca_ppm_calc > (mean(Ca_ppm_calc) - 1.25 * sd(Ca_ppm_calc)), 
         Ca_ppm_calc < (mean(Ca_ppm_calc) + 1.25 * sd(Ca_ppm_calc))) |>
   mutate(SrCa_smooth = zoo::rollmean(Sr_Ca_ppm, 
                                      k = 8, 
                                      fill = "extend", # To avoid NA values in the first and last data-points, extending the first and last value of the window
                                      align = "center"),
         Core_signal_smooth = ifelse(SrCa_smooth > 9.13, "Core", "Outside_core"), # Use the smoothed Sr:Ca ratio to identify the regions of the transect correspondng to the core
         Core_signal_raw=ifelse(Sr_Ca_ppm > 9.13, "Core", "Outside_core")) |> # Use the raw Sr:Ca ratio to identify the regions of the transect corresponding to the core
  ungroup()

# Have a look at the Sr:Ca raw values vs the smoothed values, here is one example:: 
sort(unique(LA_ICP_MS_KUL_Sotebasen_Scan_speed_core$Water_area))

LA_ICP_MS_KUL_Sotebasen_Scan_speed_core |> 
  filter(Water_area=="Göta älv") |>
  filter(IndividNr==44462) |> 
    ggplot(aes(y = SrCa_smooth,
               x = as.numeric(as.character(ElapsedTime_s)))) +
     geom_line(aes(y = Sr_Ca_ppm, 
                   x = as.numeric(as.character(ElapsedTime_s))), colour="black")+
      geom_point(aes(colour=Core_signal_smooth))+
      facet_wrap(~IndividNr*Total_length_mm*Water_area, scales = "free") +
      ylab("(Sr/Ca)*1000") +
      xlab("Elapsed_time") # Smoothing looks good (data-points, raw data is shown as a black line).

# Identify the otolith core region ans its mid-point, used to define the zero-point corresponding to the centre of the otolith core for all individuals:
# Step 3: Find the core section of the transect using cumsum to identify contiguous Sr:Ca regions of the transects above 9.12 of the smoothed Sr:Ca ratio:
LA_ICP_MS_KUL_Sotebasen_Scan_speed_otolith_core <- LA_ICP_MS_KUL_Sotebasen_Scan_speed_core |>
  arrange(IndividNr, ElapsedTime_s) |> 
  group_by(IndividNr) |>
  mutate(Core_run = cumsum(Core_signal_smooth != lag(Core_signal_smooth, default = first(Core_signal_smooth))))

# Look at how many contiguous regions that have been identified:
core_regions <- LA_ICP_MS_KUL_Sotebasen_Scan_speed_otolith_core |>
  filter(Core_signal_smooth == "Core") |> # Select only the regions with a core signal
  group_by(IndividNr, Core_run) |> # group the data per individual and identified regions with core-signals
  summarise(n_points = n(), # number of data-points per core region
            mean_SrCa_smooth = mean(SrCa_smooth, na.rm = TRUE), 
            max_SrCa_smooth = max(SrCa_smooth, na.rm = TRUE),
            min_time = min(ElapsedTime_s),
            max_time = max(ElapsedTime_s),
            center_time = (min_time + max_time)/2, # this is the center-point, works because ElapsedTime_s is in seconds
           .groups = "drop")

# Step 4: Select the longest contiguous region (most consequitive points) with a Sr:Ca value >9.12 of each transect:

core_centers <- core_regions |>
  filter(n_points > 20) |> # Remove spikes and short segments to avoiding estimate for otoliths lacking a clear core Sr:Ca signal.
  arrange(IndividNr, desc(n_points), desc(max_SrCa_smooth)) |> # Sort the estimates per individual, so the first row is the estimated center based on the longest contiguous region. 
                                                         # If there are several estimates with the same number of data points, sort on highest Sr:Ca value and keep the first row of data
  group_by(IndividNr) |> 
  slice(1) |> # Keep only the first row of data per individual
  ungroup() |> 
  mutate(Core_center_time = center_time)

# Chech the data:
core_centers |> 
  ggplot(aes(x=n_points))+
  geom_histogram(fill="white", colour="black", bins=100)+ # Looks good!
  xlim(0,250)+
  theme_classic() # Looks good!

# Add the time corresponding to the Core center time for each individual:
core_centers_merge<-core_centers |> 
  select(IndividNr, Core_center_time)

# Here we have the center-point of the otolith core region in Elapsed_time for each individual
head(core_centers_merge)

# Step 5: Merge the center-point (expressed in time) corresponding to the otolith core to the LA_ICP_MS dataframe:
LA_ICP_MS_KUL_Sotebasen_Scan_speed_zero_point<- merge(x = LA_ICP_MS_KUL_Sotebasen_Scan_speed_core, 
                                                      y = core_centers_merge, 
                                                      by = "IndividNr", all.x = TRUE)

n_distinct(LA_ICP_MS_KUL_Sotebasen_Scan_speed_zero_point$IndividNr)

## Check sample size:
LA_ICP_MS_KUL_Sotebasen_Scan_speed_zero_point_Check<- merge(x = LA_ICP_MS_KUL_Sotebasen_Scan_speed_core, 
                                                      y = core_centers_merge, 
                                                      by = "IndividNr")

head(LA_ICP_MS_KUL_Sotebasen_Scan_speed_zero_point) # Looks good, all inds are there

## Convert ElapsedTime_s to distance (micrometer) from the center-point for each transect:
LA_ICP_MS_KUL_Sotebasen_Distance<-LA_ICP_MS_KUL_Sotebasen_Scan_speed_zero_point |> 
  group_by(IndividNr) |> 
  arrange(IndividNr,ElapsedTime_s) |> 
  mutate(Distance=ElapsedTime_s*Scan_speed_micrometer_sec,
         Zero_point_micrometer=Core_center_time*Scan_speed_micrometer_sec,
         Distance_from_core = Distance - Zero_point_micrometer)
  
# Check data, here are 4 eels from river Kävlingeån with the identified 0-point:
LA_ICP_MS_KUL_Sotebasen_Distance |> 
  filter(IndividNr %in% c("48734", "48736", "48746", "50028")) |> 
     ggplot(aes(y = Sr_Ca_ppm,
                   x = Distance_from_core)) +
       geom_point(aes(colour=Core_signal_raw)) +
       facet_wrap(~IndividNr*Total_length_mm*Water_area, scales = "free") +
       geom_vline(xintercept=0)+
       ylab("(Sr/Ca)*1000") +
      xlab("Distance from core")

# Looks good, we identify the otolith core as zero and then calculated the distance from the core for each data-point of the transect
  
# Step 6: Calculate the mean distance from core to edge of the otolith for eel betwween 10-20 cm in length ###

LA_ICP_MS_KUL_Sotebasen_Distance |> 
  filter(Total_length_mm>=100,
         Total_length_mm<200,
         Distance_from_core>=0) |> 
  group_by(IndividNr) |>
  summarise(Max_distance_from_core = max(Distance_from_core)) |> 
  ungroup() |> 
  summarise(n=n_distinct(IndividNr),
         Mean_max_size=mean(Max_distance_from_core),
         Sd_max_size=sd(Max_distance_from_core)) |> 
  mutate(across(where(is.numeric), ~ round(.x, 0)))

# n=348, mean maximum length = 513 micrometers, sd=206 micrometers.

# Check the derived value for small eel caught in a river flowing to Kattegat 
LA_ICP_MS_KUL_Sotebasen_Distance |> 
 filter(Water_area=="Nybroån") |> 
 filter(Total_length_mm>200) |> 
     ggplot(aes(y = Sr_Ca_ppm,
                   x = Distance_from_core)) +
       geom_point(aes(colour=Core_signal_raw)) +
       facet_wrap(~IndividNr*Total_length_mm, scales = "free") +
       geom_vline(aes(xintercept=513))+
       ylab("(Sr/Ca)*1000") +
      xlab("Distance from core") # Looks good, we have some data in the ned of transect, and the inital core signal will be excluded.

# Save the data frame for script 6_Figures.R
saveRDS(LA_ICP_MS_KUL_Sotebasen_Distance, file="R_Scripts_and_Generated_data_files/LA_ICP_MS_KUL_Sotebasen_Distance.Rds")

#Step 8. Assign habitat for all Sr:Ca value, excluding all points between 0-513 micrometer afor all eel above 20 cm in length

LA_ICP_MS_KUL_Sotebasen_Distance_Habitat_no_filter<-LA_ICP_MS_KUL_Sotebasen_Distance |>
  #filter(Total_length_mm>=200) |> # Only include eel longer than 20 cm
  #filter(Distance_from_core>=513) |> # Remove the first section of the Sr:Ca transect, corresponding to the arrival from the Sargasso sea
 mutate(Habitat_use = case_when(Sr_Ca_ppm >=0 & Sr_Ca_ppm <= 1.67 ~ "Freshwater",
                                   Sr_Ca_ppm > 1.67 & Sr_Ca_ppm <= 2.60  ~ "Estuary",
                                   Sr_Ca_ppm > 2.60 & Sr_Ca_ppm <= 9.13 ~ "Coastal/Brackish",
                                   Sr_Ca_ppm > 9.13 ~"Marine/Core"))

LA_ICP_MS_KUL_Sotebasen_Distance_Habitat_filter20<-LA_ICP_MS_KUL_Sotebasen_Distance |>
  filter(Total_length_mm>=200) |> # Only include eel longer than 20 cm
  #filter(Distance_from_core>=513) |> # Remove the first section of the Sr:Ca transect, corresponding to the arrival from the Sargasso sea
 mutate(Habitat_use = case_when(Sr_Ca_ppm >=0 & Sr_Ca_ppm <= 1.67 ~ "Freshwater",
                                   Sr_Ca_ppm > 1.67 & Sr_Ca_ppm <= 2.60  ~ "Estuary",
                                   Sr_Ca_ppm > 2.60 & Sr_Ca_ppm <= 9.13 ~ "Coastal/Brackish",
                                   Sr_Ca_ppm > 9.13 ~"Marine/Core"))

LA_ICP_MS_KUL_Sotebasen_Distance_Habitat_filter513_20cm<-LA_ICP_MS_KUL_Sotebasen_Distance |>
  filter(Total_length_mm>=200) |> # Only include eel longer than 20 cm
  filter(Distance_from_core>=513) |> # Remove the first section of the Sr:Ca transect, corresponding to the arrival from the Sargasso sea
 mutate(Habitat_use = case_when(Sr_Ca_ppm >=0 & Sr_Ca_ppm <= 1.67 ~ "Freshwater",
                                   Sr_Ca_ppm > 1.67 & Sr_Ca_ppm <= 2.60  ~ "Estuary",
                                   Sr_Ca_ppm > 2.60 & Sr_Ca_ppm <= 9.13 ~ "Coastal/Brackish",
                                   Sr_Ca_ppm > 9.13 ~"Marine/Core"))

n_distinct(LA_ICP_MS_KUL_Sotebasen_Distance_Habitat_no_filter$IndividNr) # total of 1689 when all inds are included

n_distinct(LA_ICP_MS_KUL_Sotebasen_Distance_Habitat_filter20$IndividNr) # A total of 1166 inds with a transect longer than 513 micrometers

n_distinct(LA_ICP_MS_KUL_Sotebasen_Distance_Habitat_filter513_20cm$IndividNr) # A total of 1045 inds with a transect longer than 513 micrometers and a body length >20 cm

sort(unique(LA_ICP_MS_KUL_Sotebasen_Distance_Habitat_filter513_20cm$Habitat_use))

## Shorten the df name:
LA_ICP_MS_KUL_Sotebasen_Distance_Habitat<-LA_ICP_MS_KUL_Sotebasen_Distance_Habitat_filter513_20cm
n_distinct(LA_ICP_MS_KUL_Sotebasen_Distance_Habitat$IndividNr) # Should be 1045 inds.

## Save for script 6_Habitat_shifting.R
saveRDS(LA_ICP_MS_KUL_Sotebasen_Distance_Habitat, file="R_Scripts_and_Generated_data_files/LA_ICP_MS_KUL_Sotebasen_Distance_Habitat.Rds")

## Assign life-history to each individual:
LA_ICP_MS_KUL_Sotebasen_Life_history<-LA_ICP_MS_KUL_Sotebasen_Distance_Habitat|>
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

LA_ICP_MS_KUL_Sotebasen_Life_history |> 
  group_by(Life_history) |> 
   summarise(n=n()) |> 
   mutate(sum=sum(n)) # 1045 inds, all hav an assigned life-hisotry

## Save data-frame to read in to script 6_Figures.R:
saveRDS(LA_ICP_MS_KUL_Sotebasen_Life_history, file="R_Scripts_and_Generated_data_files/LA_ICP_MS_KUL_Sotebasen_Life_history.RDS")

## Combine the assigned life-history with the transect data and look at the classification:
LA_ICP_MS_KUL_Sotebasen_Life_history|> 
  filter(IndividNr %in% c("51262", "51256", "50161",
                          "48767", "51214", "51186", "46596", "46615", "51245", "50848", "50161"))

LA_ICP_MS_KUL_Sotebasen_Distance_Habitat|> 
  filter(IndividNr %in% c("51262", "51256", "50161",
                          "48767", "51214", "51186", "46596", "46615", "51245", "50848", "50161")) |> 
      ggplot(aes(y = Sr_Ca_ppm,
                   x = Distance_from_core)) +
       geom_point(aes(colour=Habitat_use)) +
       facet_wrap(~IndividNr*Total_length_mm*Water_area) +
      ylab("Sr:Ca ratio") +
      xlab("Distance from core")

# Looks good, continute and make a summary table of all inds.
View(LA_ICP_MS_KUL_Sotebasen_Distance_Habitat)

# Make summary table for the four major geographical regions:
sort(unique(LA_ICP_MS_KUL_Sotebasen_Life_history$Water_area))

LA_ICP_MS_KUL_Sotebasen_Life_history_Major_regions<-LA_ICP_MS_KUL_Sotebasen_Life_history |> 
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

flextable(LA_ICP_MS_KUL_Sotebasen_Life_history_Major_regions)

sum(LA_ICP_MS_KUL_Sotebasen_Life_history_Major_regions$Individuals_n)

## Per sampling site:
LA_ICP_MS_KUL_Sotebasen_Life_history_all_sites<-LA_ICP_MS_KUL_Sotebasen_Life_history |> 
  group_by(Water_area, Life_history) |> 
  summarise(Individuals_n=n()) |> 
  group_by(Water_area) |> 
  mutate(Total_amount_of_inds=sum(Individuals_n),
         Percent_per_life_history=round((Individuals_n/Total_amount_of_inds)*100, 1))

LA_ICP_MS_KUL_Sotebasen_Life_history_all_sites

flextable(LA_ICP_MS_KUL_Sotebasen_Life_history_all_sites)

#### Done with LA-ICP-MS Data ####

# -------------------------------------------- #

#### Individual habitat assignment for eel with FE-EPMA otolith chemistry data #####

## Load FE-EPMA otolith chemistry transect data combined with individual level data from the database Sötebasen:
FE_EPMA_Sotebasen = readRDS("R_Scripts_and_Generated_data_files/FE_EPMA_Sotebasen.Rds") # Containing 22 columns with 58740 rows of data 
head(FE_EPMA_Sotebasen)
str(FE_EPMA_Sotebasen)

### Calculate distance from core for each datapoint using the X and Y coordinates for each individual:
FE_EPMA_Sotebasen_Distance <- FE_EPMA_Sotebasen %>%
  group_by(IndividNr) %>%
  arrange(Date, .by_group = TRUE) %>%
  mutate(Sr_Ca_ratio = as.numeric(as.character(Sr)) / as.numeric(as.character(Ca)) * 1000,
         Step_dist = sqrt((X - lag(X))^2 + (Y - lag(Y))^2),  # Distance between consecutive points -> distance between measure-points
         Distance_from_start = cumsum(replace_na(Step_dist, 0))) %>%  # Cumulative distance from first point which is set to 0 (the center of the core), which is the start of the transect
  ungroup()

## Save df for script 6_Figures.R
saveRDS(FE_EPMA_Sotebasen_Distance, file = "R_Scripts_and_Generated_data_files/FE_EPMA_Sotebasen_Distance.Rds")

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

# No filter:
FE_EPMA_Sotebasen_Distance_Habitat_no_filter<-FE_EPMA_Sotebasen_Distance |> 
  #filter(Total_length_mm>=200, # Only include eel longer than 20 cm
    #     Distance_from_start>=513) |> # Remove the first section of the Sr:Ca transect, corresponding to the arrival from the Sargasso sea
   filter(!is.na(Sr_Ca_ratio)) |> 
    mutate(Habitat_use = case_when(Sr_Ca_ratio <= 1.12 ~ "Freshwater",
                                   Sr_Ca_ratio > 1.12 & Sr_Ca_ratio <= 6.97 ~ "Coastal/Brackish",
                                   Sr_Ca_ratio > 6.97 ~"Marine/Core"))
# >20 cm filter:
FE_EPMA_Sotebasen_Distance_Habitat_filter_20cm<-FE_EPMA_Sotebasen_Distance |> 
  filter(Total_length_mm>=200) |> # Only include eel longer than 20 cm
    #     Distance_from_start>=513) |> # Remove the first section of the Sr:Ca transect, corresponding to the arrival from the Sargasso sea
   filter(!is.na(Sr_Ca_ratio)) |> 
    mutate(Habitat_use = case_when(Sr_Ca_ratio <= 1.12 ~ "Freshwater",
                                   Sr_Ca_ratio > 1.12 & Sr_Ca_ratio <= 6.97 ~ "Coastal/Brackish",
                                   Sr_Ca_ratio > 6.97 ~"Marine/Core"))

# >20 cm filter and transect longer than 513 micrometers:
FE_EPMA_Sotebasen_Distance_Habitat_filter_20cm_513mm<-FE_EPMA_Sotebasen_Distance |> 
  filter(Total_length_mm>=200, # Only include eel longer than 20 cm
         Distance_from_start>=513) |> # Remove the first section of the Sr:Ca transect, corresponding to the arrival from the Sargasso sea
   filter(!is.na(Sr_Ca_ratio)) |> 
     mutate(Habitat_use = case_when(Sr_Ca_ratio <= 1.12 ~ "Freshwater",
                                   Sr_Ca_ratio > 1.12 & Sr_Ca_ratio <= 6.97 ~ "Coastal/Brackish",
                                   Sr_Ca_ratio > 6.97 ~"Marine/Core"))

n_distinct(FE_EPMA_Sotebasen_Distance_Habitat_no_filter$IndividNr) # 1985 inds
n_distinct(FE_EPMA_Sotebasen_Distance_Habitat_filter_20cm$IndividNr) # 1893 inds
n_distinct(FE_EPMA_Sotebasen_Distance_Habitat_filter_20cm_513mm$IndividNr) # 1893 inds (no inds longer than 20 cm in length had a shorter transect than 513 micro meter)

# Rename the data frame:
FE_EPMA_Sotebasen_Distance_Habitat<-FE_EPMA_Sotebasen_Distance_Habitat_filter_20cm_513mm

## Save df for script 5_Habitat_shifts and script 6_Habitat_shifting.R
saveRDS(FE_EPMA_Sotebasen_Distance_Habitat, file="R_Scripts_and_Generated_data_files/FE_EPMA_Sotebasen_Distance_Habitat.Rds")

## Assign life-history to each individual:

FE_EPMA_Sotebasen_Life_history<-FE_EPMA_Sotebasen_Distance_Habitat |>
  group_by(Water_area, Site, IndividNr) |>
  summarise(has_fresh = any(Habitat_use == "Freshwater"),
            has_coastal = any(Habitat_use == "Coastal/Brackish"),
            has_marine = any(Habitat_use == "Marine/Core"),
            .groups = "drop") |>
  mutate(Life_history = case_when(has_fresh & !has_coastal & !has_marine ~ "Freshwater", # Only freshwater Sr:Ca values
                                  has_fresh & (has_coastal | has_marine) ~ "Habitat shifter", # Both fresheater and coastal or marine signal
                                  !has_fresh ~ "Coastal", # No freshwater signal
                                   TRUE ~ NA_character_))

n_distinct(FE_EPMA_Sotebasen_Life_history$IndividNr)

# Save data-frame to read in to script 3_Figures.R and 6_Habitat_assignment.R:
saveRDS(FE_EPMA_Sotebasen_Life_history, file="R_Scripts_and_Generated_data_files/FE_EPMA_Sotebasen_Life_history.Rds")

## Summary table major water areas:
sort(unique(FE_EPMA_Sotebasen_Life_history$Water_area))

FE_EPMA_Sotebasen_Life_history_Major_regions<-FE_EPMA_Sotebasen_Life_history|> 
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

flextable(FE_EPMA_Sotebasen_Life_history_Major_regions)

sum(FE_EPMA_Sotebasen_Life_history_Major_regions$Individuals_n)


## Make summary table for all sites:
FE_EPMA_Sotebasen_Life_history_all_sites<-FE_EPMA_Sotebasen_Life_history |>
  filter(!Water_area=="Ingen uppgift/ej definierat") |> 
  group_by(Water_area, Life_history) |> 
  summarise(Individuals_n=n()) |> 
  group_by(Water_area) |> 
  mutate(Total_amount_of_inds=sum(Individuals_n),
         Percent_per_life_history=round((Individuals_n/Total_amount_of_inds)*100, 1))

FE_EPMA_Sotebasen_Life_history_all_sites

flextable(FE_EPMA_Sotebasen_Life_history_all_sites)

# -------- Done with FE-EPMA life-hisotry assignment --------- #

# -------- Table 4, main text --------- #

### Combined life_history categorization results:
flextable(LA_ICP_MS_KUL_Sotebasen_Life_history_Major_regions)
sum(LA_ICP_MS_KUL_Sotebasen_Life_history_Major_regions$Individuals_n) # Should be 1045

flextable(FE_EPMA_Sotebasen_Life_history_Major_regions)
sum(FE_EPMA_Sotebasen_Life_history_Major_regions$Individuals_n) # Should be 1893

# Combine into one table, fix format and save as a flextable for main text in the manuscript:
Complete_Life_history_Major_regions<-rbind(FE_EPMA_Sotebasen_Life_history_Major_regions, LA_ICP_MS_KUL_Sotebasen_Life_history_Major_regions)

Complete_Life_history_Major_regions |> 
  arrange(Catch_area_habitat, Analythical_method, Life_history)
sum(Complete_Life_history_Major_regions$Individuals_n) # Should be 2938 which it is

#Check which life-histories are included in the data:
Complete_Life_history_Major_regions |> 
  group_by(Life_history) |> 
    summarise(n=sum(Individuals_n)) |> 
    mutate(tot_n=sum(n))

# Flextable:
Table_4<- Complete_Life_history_Major_regions %>%
  ungroup() |> 
  mutate(Sampling_habitat = case_when(grepl("^Coastal", Catch_area_habitat) ~ "Coastal",
                                      grepl("^Freshwater", Catch_area_habitat) ~ "Freshwater"),
         Sampling_region = case_when(grepl("Baltic_Sea", Catch_area_habitat) ~ "Baltic Sea",
                                     grepl("Kattegat", Catch_area_habitat) ~ "Skagerrak/Kattegat"),
         Assigned_life_history = recode(Life_history, "Coastal" = "Coastal resident",
                                                       "Freshwater" = "Freshwater resident",
                                                       "Freshwater/Estuary" = "Freshwater/Estuary",
                                                       "Habitat shifter" = "Habitat shifter")) %>%
  select(Sampling_habitat,
         Sampling_region,
         Assigned_life_history,
         Total_N = Total_amount_of_inds,
         Assigned_n = Individuals_n,
         Assigned_percent = Percent_per_life_history,
         Analytical_method = Analythical_method) %>%
  arrange(desc(Sampling_region),
               Sampling_habitat,
               Analytical_method,
               Assigned_life_history)
Table_4

# Make the Flextable:

Table_4_ft <- flextable(Table_4) %>%
  set_header_labels(Sampling_habitat = "Sampling habitat",
                    Sampling_region = "Sampling region",
                    Assigned_life_history = "Assigned life history",
                    Total_N = "Total N",
                    Assigned_n = "Assigned n",
                    Assigned_percent = "Assigned %",
                    Analytical_method = "Analytical method") %>% # Remove "_" and add "%" sign
  colformat_double(j = "Assigned_percent",
                   digits = 1) %>% # Round to 1 digit
   flextable::align(j = "Analytical_method",
         align = "right", part = "all") |>
  autofit()

Table_4_ft

# Save the table in the Output folder:
# As an image:
save_as_image(Table_4_ft, path = "Output/Table_4.png")

# As a word-file:
save_as_docx(Table_4_ft, path = "Output/Table_4.docx")

# --------- End of Script 4_Habitat_assignment.R ------------- #


