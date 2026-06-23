#################
# Jacobson et al.
# Navigating the salinity gradient: Individual variation in habitat use and migration of European eel revealed by otolith microchemistry 
#################

## Last update: 2026-06-16

#### R-Script: 5_Habitat_shifts.R

# Load packages
library(tidyverse)

## Load data generated in Script 4_Habitat_assignment.R

# Otolith transect data with converted habitat values:
FE_EPMA_Sotebasen_Distance_Habitat<-readRDS("R_Scripts_and_Generated_data_files/FE_EPMA_Sotebasen_Distance_Habitat.Rds")
LA_ICP_MS_KUL_Sotebasen_Distance_Habitat<-readRDS("R_Scripts_and_Generated_data_files/LA_ICP_MS_KUL_Sotebasen_Distance_Habitat.Rds")

# Assigned life_history data, one row per individual:
FE_EPMA_Sotebasen_Life_history = readRDS("R_Scripts_and_Generated_data_files/FE_EPMA_Sotebasen_Life_history.Rds")
LA_ICP_MS_KUL_Sotebasen_Life_history = readRDS("R_Scripts_and_Generated_data_files/LA_ICP_MS_KUL_Sotebasen_Life_history.Rds")

# FE-EPMA data ######

# Make a left_join so thath every ind gets its assigned life-history:
FE_EPMA_Sotebasen_Life_history_selected<-FE_EPMA_Sotebasen_Life_history |> 
  select(IndividNr, Life_history)

# Combine with transect data:
FE_EPMA_Sotebasen_Distance_Habitat_Assigned_Life_history<-left_join(FE_EPMA_Sotebasen_Distance_Habitat, FE_EPMA_Sotebasen_Life_history_selected, by="IndividNr")

# Extract all Habitat shifters and count the number of shifts between freshwater and coastal habitats:
FE_EPMA_Sotebasen_Distance_Habitat_shift<-FE_EPMA_Sotebasen_Distance_Habitat_Assigned_Life_history |> 
  filter(Life_history=="Habitat shifter")

sort(unique(FE_EPMA_Sotebasen_Distance_Habitat_shift$Habitat_use))
sort(unique(FE_EPMA_Sotebasen_Distance_Habitat_shift$Water_area))

## Quantify the number of shifts between freshwater and coastal habitats per individual:
Habitat_shifts_FE_EPMA <- FE_EPMA_Sotebasen_Distance_Habitat_shift %>%
  arrange(IndividNr, Distance_from_start) %>%
  mutate(Habitat_simple = case_when(Habitat_use == "Freshwater" ~ "Freshwater",
                                    Habitat_use %in% c("Coastal/Brackish", "Marine/Core") ~ "Coastal/Marine",
                                    TRUE ~ NA_character_)) %>%
  group_by(IndividNr) %>%
  mutate(run_id = cumsum(Habitat_simple != lag(Habitat_simple, default = first(Habitat_simple)))) %>%
  group_by(IndividNr, run_id) %>%
  mutate(run_length = n()) %>%
  ungroup() %>%
  distinct(IndividNr, run_id, Habitat_simple, run_length) %>%
  arrange(IndividNr, run_id) %>%
  group_by(IndividNr) %>% 
  mutate(prev_habitat = lag(Habitat_simple),
         prev_run_length = lag(run_length),
         valid_shift = run_length >= 1 & prev_run_length >= 1,
         shift_from_FW_to_B = valid_shift & Habitat_simple == "Coastal/Marine" & prev_habitat == "Freshwater",
         shift_from_B_to_FW = valid_shift & Habitat_simple == "Freshwater" & prev_habitat == "Coastal/Marine") %>%
  summarise(n_shifts_total = sum(shift_from_FW_to_B | shift_from_B_to_FW, na.rm = TRUE),
            n_from_FW_to_B = sum(shift_from_FW_to_B, na.rm = TRUE),
            n_from_B_to_FW = sum(shift_from_B_to_FW, na.rm = TRUE),
            .groups = "drop")

# Summary table for the results section:
Habitat_shifts_FE_EPMA |> 
  group_by(n_shifts_total) |> 
   summarise(n=n_distinct(IndividNr)) |> 
    mutate(tot_n=sum(n),
           percentage=(n/tot_n)*100)

## Save the datatable for the script 6_Figures.R
saveRDS(Habitat_shifts_FE_EPMA , file="R_Scripts_and_Generated_data_files/Habitat_shifts_FE_EPMA_Mixed.Rds")

# ------ Done with FE-EPMA data ------- #

### LA-ICP-MS #####

# First, select the individuals assigned with a mixed life-history:
LA_ICP_MS_KUL_Sotebasen_Life_history_selected<-LA_ICP_MS_KUL_Sotebasen_Life_history |> 
  select(IndividNr, Life_history)

# Combine with transect data for each individual:
LA_ICP_MS_KUL_Sotebasen_Distance_Habitat_Assigned_Life_history<-left_join(LA_ICP_MS_KUL_Sotebasen_Distance_Habitat, LA_ICP_MS_KUL_Sotebasen_Life_history_selected, by="IndividNr")

# Filter out eel with an assigned mixed life-history:
LA_ICP_MS_KUL_Sotebasen_Distance_Habitat_shifter<-LA_ICP_MS_KUL_Sotebasen_Distance_Habitat_Assigned_Life_history |> 
  filter(Life_history=="Habitat shifter")

names(LA_ICP_MS_KUL_Sotebasen_Distance_Habitat_shifter)

Habitat_shifts_LA_ICP_Mixed <- LA_ICP_MS_KUL_Sotebasen_Distance_Habitat_shifter %>%
  arrange(IndividNr, Distance_from_core) %>% # Make sure the data is sorted correctly
  mutate(Habitat_simple = case_when(Habitat_use %in% c("Freshwater", "Estuary") ~ "Freshwater",
                                    Habitat_use %in% c("Coastal/Brackish", "Marine/Core") ~ "Coastal/Marine",
                                    TRUE ~ NA_character_)) %>%
  group_by(IndividNr) %>%
  mutate(run_id = cumsum(Habitat_simple != lag(Habitat_simple, default = first(Habitat_simple)))) %>%
  group_by(IndividNr, run_id) %>%
  mutate(run_length = n()) %>%
  ungroup() %>%
  distinct(IndividNr, run_id, Habitat_simple, run_length) %>%
  arrange(IndividNr, run_id) %>%
  group_by(IndividNr) %>% 
  mutate(prev_habitat = lag(Habitat_simple),
         prev_run_length = lag(run_length),
         valid_shift = run_length >= 1 & prev_run_length >= 1,
         shift_from_FW_to_B = valid_shift & Habitat_simple == "Coastal/Marine" & prev_habitat == "Freshwater",
         shift_from_B_to_FW = valid_shift & Habitat_simple == "Freshwater" & prev_habitat == "Coastal/Marine") %>%
  summarise(n_shifts_total = sum(shift_from_FW_to_B | shift_from_B_to_FW, na.rm = TRUE),
            n_from_FW_to_B = sum(shift_from_FW_to_B, na.rm = TRUE),
            n_from_B_to_FW = sum(shift_from_B_to_FW, na.rm = TRUE),
            .groups = "drop")

Habitat_shifts_LA_ICP_Mixed|> 
  group_by(n_shifts_total) |> 
   summarise(n=n_distinct(IndividNr)) |> 
    mutate(tot_n=sum(n),
           percentage=(n/tot_n)*100) |> 
  print(n=Inf)

## Save the datatable for the script 6_Figures.R
saveRDS(Habitat_shifts_LA_ICP_Mixed, file="R_Scripts_and_Generated_data_files/Habitat_shifts_LA_ICP_Mixed.Rds")

# -------- End of script 5_Habitat_shifts.R -------- #