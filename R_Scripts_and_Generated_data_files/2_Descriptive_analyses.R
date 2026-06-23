#################
# Jacobson et al.
# Navigating the salinity gradient: Individual variation in habitat use and migration of European eel revealed by otolith microchemistry 
#################

## Last update: 2026-06-16

#### R-Script: 2_Descriptive analyses ####

### This script contains the code for creating table 1 in the main text.

### Load needed packages:

library(tidyverse)
library(flextable)
library(officer)
library(webshot2)  # only needed for saving Table 1 as an image, using the function save_as_image()

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

### Descriptive summary table ###

head(FE_EPMA_LA_ICP_MS_Individual_data)

sort(unique(FE_EPMA_LA_ICP_MS_Individual_data$Water_area))

# ------ Table 1, main text --------- #

FE_EPMA_LA_ICP_MS_Individual_data  |> 
  group_by(Water_area, Site) |> 
   summarise(n=n_distinct(IndividNr)) |> 
  filter(Water_area %in% c("Östersjön (hela) ICES SD 22-32", "Västerhavet (hela) ICES SD 20-21")) |> 
  print(n=35) # these are the coastal sites whihc needs to be re-categorize and translate from Swedish to English 

# Translate, recategorize costal sampling sites and make the summary table:
Table_1 <-FE_EPMA_LA_ICP_MS_Individual_data |> 
  filter(Water_area!="Artificiellt vatten") |> # Remove these in script 1!
  filter(Water_area!="Ingen uppgift/ej definierat") |> # Remove these in script 1
   mutate(Sampling_site=ifelse(Water_area %in% c("Östersjön (hela) ICES SD 22-32", 
                                                 "Västerhavet (hela) ICES SD 20-21"), Site, Water_area)) |> 
   mutate(Sampling_area=ifelse(Sampling_site=="Borgholms kommuns kustvatten", "Brackish, Baltic Sea, Sweden",
                        ifelse(Sampling_site=="Faxe Bugt-Jungshoved (Danmark)", "Brackish, Baltic Sea, Denmark",
                        ifelse(Sampling_site=="Hamina (Finland)", "Brackish, Baltic Sea, Finland",
                        ifelse(Sampling_site=="Höganäs-Kullen", "Brackish, Baltic Sea, Sweden",
                        ifelse(Sampling_site=="Jyllinge (Danmark)", "Brackish, Kattegat, Denmark",
                        ifelse(Sampling_site=="Kalmarsund", "Brackish, Baltic Sea, Sweden",
                        ifelse(Sampling_site=="Karlshamn - Boön", "Brackish, Baltic Sea, Sweden",
                        ifelse(Sampling_site=="Karlshamn - Veggabadet", "Brackish, Baltic Sea, Sweden",
                        ifelse(Sampling_site=="Kattegatt (ICES SD 21)", "Brackish, Kattegat, Denmark", # Caught in Roskildefjorden according to comment in the database
                        ifelse(Sampling_site=="Kemi-Montaja (Finland)", "Brackish, Baltic Sea, Finland",
                        ifelse(Sampling_site=="Klamilanlahti (Finland)", "Brackish, Baltic Sea, Finland",     
                        ifelse(Sampling_site=="Kõiguste (Estland)", "Brackish, Baltic Sea, Estonia",
                        ifelse(Sampling_site=="Kävlinge-Barsebäck", "Brackish, Baltic Sea, Sweden",
                        ifelse(Sampling_site=="Norrtälje-Bylehamn", "Brackish, Baltic Sea, Sweden",
                        ifelse(Sampling_site=="Norrtälje-Väddö", "Brackish, Baltic Sea, Sweden",
                        ifelse(Sampling_site=="Nynäshamns kommuns kustvatten", "Brackish, Baltic Sea, Sweden",
                        ifelse(Sampling_site=="Oskarshamn-Kråkelund", "Brackish, Baltic Sea, Sweden",
                        ifelse(Sampling_site=="Rymättylä (Finland)", "Brackish, Baltic Sea, Finland",
                        ifelse(Sampling_site=="Skurup-Abbekås", "Brackish, Baltic Sea, Sweden",
                        ifelse(Sampling_site=="Stege (Danmark)", "Brackish, Baltic Sea, Denmark",
                        ifelse(Sampling_site=="Stenungsund", "Marine/Brackish, Skagerrak, Sweden",
                        ifelse(Sampling_site=="Söderköping-Aspöja", "Brackish, Baltic Sea, Sweden",
                        ifelse(Sampling_site=="Sölvesborg-Hörvik", "Brackish, Baltic Sea, Sweden",
                        ifelse(Sampling_site=="Sölvesborg-Valjeviken", "Brackish, Baltic Sea, Sweden",
                        ifelse(Sampling_site=="Sölvesborgs kommuns kustvatten", "Brackish, Baltic Sea, Sweden",
                        ifelse(Sampling_site=="Tanum-Havstenssund", "Marine, Skagerrak, Sweden",
                        ifelse(Sampling_site=="Tanum-Veddökilen", "Marine, Skagerrak, Sweden",
                        ifelse(Sampling_site=="Vellinge-Kämpingebukten", "Brackish, Baltic Sea, Sweden",
                        ifelse(Sampling_site=="Vilsandi (Estland)", "Brackish, Baltic Sea, Estonia",
                        ifelse(Sampling_site=="Virolahti (Finland)", "Brackish, Baltic Sea, Finland",
                        ifelse(Sampling_site=="Västervik-Marsö", "Brackish, Baltic Sea, Sweden",
                        ifelse(Sampling_site=="Weichsel lagun (Polen)", "Brackish, Baltic Sea, Poland", 
                                                                        "Freshwater, Inland, Sweden"))))))))))))))))))))))))))))))))) |> 
   separate(Sampling_area, into = c("Sampling_habitat", "Area", "Country"), sep = ",\\s*") |> # Separate at the comma and remove blanks/white-space.
   group_by(Sampling_habitat, Area, Country) |> 
   summarise(Number_of_eel=n_distinct(IndividNr),
            Mean_length_mm=mean(Total_length_mm, na.rm=T),
            Sd_length_mm=sd(Total_length_mm, na.rm=T),
            Mean_weight_gram=mean(Total_weight_g, na.rm=T),
            Sd_weight_gram=sd(Total_weight_g, na.rm=T)) |> # Add sampling years
  mutate(across(where(is.numeric), ~ round(.x, 1))) # Round all numeric values to 1 digit.

data.frame(Table_1) # Looks good.

sum(Table_1$Number_of_eel) # Should be 3648 which it is.

## Create a flextable with wanted layout, sorted from Marine -> Freshwater:
Table_1_ft <- Table_1 |>
  mutate(Sampling_habitat = factor(Sampling_habitat, 
                                   levels = c("Marine",
                                              "Marine/Brackish",
                                              "Brackish",
                                              "Freshwater"))) |>
    arrange(Sampling_habitat, Area, Country) |> # Fix the order
  mutate(`Mean length, mm (±SD)` = sprintf("%.1f (±%.1f)", Mean_length_mm, Sd_length_mm),
         `Mean weight, g (±SD)` = sprintf("%.1f (±%.1f)", Mean_weight_gram, Sd_weight_gram)) |> # Fix the table header format
  select(Sampling_habitat,
         Area,
         Country,
         Number_of_eel,
         `Mean length, mm (±SD)`,
         `Mean weight, g (±SD)`) |>
  rename(`Sampling habitat` = Sampling_habitat,
         `n eel` = Number_of_eel) |>
  flextable() |>
  theme_booktabs() |>
    align(j = c("Sampling habitat", "Area", "Country"),
        align = "left", part = "all") |>
  align(j = c("n eel",
              "Mean length, mm (±SD)",
              "Mean weight, g (±SD)"),
        align = "left", part = "header") |>
    align(j = c("n eel",
              "Mean length, mm (±SD)",
              "Mean weight, g (±SD)"),
              align = "right", part = "body") |>
  autofit()

Table_1_ft

# Save the table in the output folder:

# As an image:
save_as_image(Table_1_ft, path = "Output/Table_1.png")

# As a word-file:
save_as_docx(Table_1_ft, path = "Output/Table_1.docx")

## Number of individuals per analytical method, for the M&M section ###

names(FE_EPMA_LA_ICP_MS_Individual_data)
sort(unique(FE_EPMA_LA_ICP_MS_Individual_data$Water_area))

FE_EPMA_LA_ICP_MS_Individual_data |> 
  group_by(Analytical_method) |> 
   summarise(Number_of_inds=n()) |> 
   ungroup() |> 
   mutate(Tot_inds=sum(Number_of_inds))

## Should be 1958 for FE-EPMA and 1690 for LA-ICP-MS (total = 3648)

# --------- End of Script 2_Descriptive_analyses.R ------------- #
