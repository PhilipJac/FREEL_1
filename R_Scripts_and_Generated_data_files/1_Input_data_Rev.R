
#################
# Jacobson et al.
# Navigating the salinity gradient: Individual variation in habitat use and migration of European eel revealed by otolith microchemistry 
#################

## Last update: 2026-10-28

#### R-Script: 1_Input_data ####

### This script loads and merge otolith microchemistry transect data (LA-ICP-MS and FE-EPMA) with 
### individual level eel data from the databases Sötebasen and KUL and saves these combined data-files as .RDS files into the "R_Scripts_and_Generated_data_files" folder
### for later usage

rm(list = ls())

### Load packages (Check which once are needed for this script later on!):
 
library(tidyverse)
library(stringr)
library(DBI)
library(odbc)

##### LA-ICP-MS data #####

# Path to the folder that contains all the Output_1_TimeSeries folders with transect data, stored as individual .txt files
main_path <- "Input_data/LA_ICP_MS/Sotebasen"

# 1. Import a list of all .txt files
files <- list.files(main_path,
                    pattern = "\\.txt$",
                    full.names = TRUE,
                    recursive = TRUE)

## Check for duplicates in file-name, should not be any duplicates because there should be one file for each individual:

file_names <- basename(files)
file_counts <- table(file_names)
duplicate_files <- file_counts[file_counts > 1]

duplicate_files  # should be 0, no duplicates of raw-data transects in the folder

# 3. Import all files
df_list <- list()

# Define the columns you are interested in (e.g., the first 6 columns)
# Select the elements we want to speed up the import loop
columns_of_interest <- c("ElapsedTime_s", "Ca43_CPS", "IntStdWv", "Sr_ppm_m88") 

## A loop to identify, read in and combine all LA-ICP-MS data:

for (file in files) {

  df <- read.delim(file, sep = "\t", header = TRUE)

  ### Select wanted columns:

  df <- df[, columns_of_interest]

  # Add the filename (e.g., "50025.txt"), whihc is the Id-nr of each individual, used to combine with ind. level data
  df$filename <- basename(file)

  # Add the folder name (Output_1_TimeSeries)
  df$folder <- basename(dirname(file))

  df_list[[length(df_list) + 1]] <- df
}

# 4. Combine all text files into one data frame
combined_df_LA_ICP_MS <- do.call(rbind, df_list)

combined_df_LA_ICP_MS$IndividNr <- as.integer(gsub("[^0-9]", "", combined_df_LA_ICP_MS$filename)) # Extract the Individual ID nr fto enable merge with ind. level data in the database Sötebasen

head(combined_df_LA_ICP_MS, 30) # Looks good

#### Load Individual level data from Sötebasen ####

## Sötebasen:
Sotebasen_inds <- read.delim(file = "Input_data/Ind_level_data/ÅlindividerSötebasen.csv",
                                     sep = ";",
                                     dec = ",",
                                     fileEncoding = "Latin1")

## Check for dublicates in IndividNr, correspdonding to the .txt file name containing the otolith microchemistry data:

table(table(Sotebasen_inds$IndividNr) > 1) # We have 20 duplciates in IndividNr

# Remove duplicates:
Sotebasen_inds_unique <- Sotebasen_inds %>%
  group_by(IndividNr) %>%
  slice(1) %>%  # Keep the first row of the duplicated row
  ungroup()

table(table(Sotebasen_inds_unique$IndividNr) > 1) # There we go, all duplicated rows removed, 22575 individuals in total

## Combine Individual level data and the LA-ICP-MS data:

LA_ICP_MS_combined_Sotebasen <- merge(x = combined_df_LA_ICP_MS, 
                                      y = Sotebasen_inds_unique, 
                                      by = "IndividNr", all.x = TRUE)

head(LA_ICP_MS_combined_Sotebasen)

LA_ICP_MS_combined_Sotebasen |> 
  distinct(IndividNr, .keep_all = TRUE) |> 
    group_by(Behandling) |> 
     summarise(n=n())

## Filter out the eel that have been sampled from the restocking quarantine facility: 

LA_ICP_MS_combined_Sotebasen<-LA_ICP_MS_combined_Sotebasen %>% 
  filter(!VattenNamn=="Artificiellt vatten")

# Check the number of individuals with LA-ICP-MS data:

n_distinct(LA_ICP_MS_combined_Sotebasen$IndividNr) # 1481 individuals.

## Done loading and combining LA-ICP-MS data with individual level data stored in Sötebasen.

### Load LA-ICP-MS data for eel with data in the database KUL

 # Define path where the data is stored:
main_path_2 <- "Input_data/LA_ICP_MS/KUL"

# 1. Import a list of all .txt files
files <- list.files(main_path_2,
                    pattern = "\\.txt$",
                    full.names = TRUE,
                    recursive = TRUE)

## Check for duplicates in file-name, should not be any duplicates because there should be one file for each individual:

file_names <- basename(files)
file_counts <- table(file_names)
duplicate_files <- file_counts[file_counts > 1]

duplicate_files  # should be 0, no duplicates of raw-data transects in the folder

# 3. Import all files
df_list <- list()

# Define the columns you are interested in (e.g., the first 6 columns)
# Select the elements we want to speed up the import loop
columns_of_interest <- c("ElapsedTime_s", "Ca43_CPS", "IntStdWv", "Sr_ppm_m88") 
 
 ## For Loop to find, read in and combine all LA-ICP-MS data:
 
 for (file in files) {
 
   df <- read.delim(file, sep = "\t", header = TRUE)
 
   ### Select wanted columns:
   df <- df[, columns_of_interest]
 
   # Add the filename (e.g., "50025.txt")
   df$filename <- basename(file)
 
   # Add the folder name (Output_1_TimeSeries)
   df$folder <- basename(dirname(file))
 
   df_list[[length(df_list) + 1]] <- df
 }

 # 4. Combine all text files from the LOOP above into one data frame
 LA_ICP_MS_KUL <- do.call(rbind, df_list)

sort(unique(LA_ICP_MS_KUL$filename)) # These .txt have another individualID-nr, where every year is a letter (except for 2019, when there is a mix) and then the ind.nr starts at 1.

head(LA_ICP_MS_KUL)

 # Conver the filename to its correpsonding ID-nr for merging with ind.level data in the database KUL: 
LA_ICP_MS_KUL <- LA_ICP_MS_KUL %>%
  mutate(filename_clean = sub("\\.txt$", "", filename),
         IndividNr = case_when(grepl("-19$", filename_clean) ~ as.integer(sub("-19$", "", filename_clean)),
                     TRUE ~ as.integer(gsub("[^0-9]", "", filename_clean))),
         Letter = gsub("[^A-Za-z]", "", filename_clean),
         Year = case_when(grepl("-19$", filename_clean) ~ 2019,
                                Letter == "A" ~ 2013,
                                Letter == "B" ~ 2014,
                                Letter == "C" ~ 2015,
                                Letter == "D" ~ 2016,
                                Letter == "E" ~ 2017,
                                Letter == "F" ~ 2018,
                                Letter == "G" ~ 2019,
                                Letter == "H" ~ 2021,
                                TRUE ~ NA_integer_),
          Comb_Year_IndividNr = paste(Year, IndividNr, sep = "_"))

sort(unique(LA_ICP_MS_KUL$Comb_Year_IndividNr)) # Looks good.

LA_ICP_MS_KUL |>
    summarise(number_of_inds=n_distinct(Comb_Year_IndividNr)) # 209 inds.

### Import individual level data for Stenungssund, extracted from KUL #####

## This section is commented out since you need to have an internal SLU log-in account and a permit to connect to the the internal server:

#AQ_con <- dbConnect(odbc(), 
#                   Driver = "SQL Server", 
#                  Server = "aquadw_prod.db.slu.se", 
#                    USER = "phja0001",
#                    Database = "AquaDW_prod", 
#                    Trusted_Connection = "True",
#                    encoding = "utf-8")

#All_eels_caught_after_2012_KUL<-dbGetQuery(AQ_con,"select * from AquaDW_prod.mart.KUL_INDIVID_NEW where Art = 'Ål' AND Fiskeår > 2012") # extracted 2026-06-15, extract all indivudual-level data for all eel caught after 2012 in the whole database

#head(All_eels_caught_after_2012_KUL)
## Select the wanted columns and save as and RDS-file:

#All_eels_caught_after_2012_KUL_SS <- All_eels_caught_after_2012_KUL |> 
 # filter(Lokal=="Stenungsund", 
 #       Artbestämning=="Gulål") |> # Filter out individual level data of yeelow eel from Stenungsund, the coastal site we have LA-ICP-MS data from. We have only analysed yellow eel.
 # mutate(Comb_Year_IndividNr=paste(Fiskeår, Löpnr, sep = "_")) |> # Create a column combining Year and IndividNr which start at 1 for each year, usde to combine with LA-ICP-MS data
 # select(Artbestämning, Lokal, Fiskeår, Löpnr, Slutlängd_mm, Totalvikt_gram, Kön, Ålder, Bäst_longitud_decimal_grader, Bäst_latitud_decimal_grader, Comb_Year_IndividNr) |> # Selected wanted columns
 #   rename(Life_stage = Artbestämning, # Translate columns to english
 #          Site = Lokal,
 #          Year = Fiskeår, # Sampling year
 #          IndividNr = Löpnr,
 #          Total_length_mm = Slutlängd_mm,
 #          Total_weight_g = Totalvikt_gram,
 #          Sex = Kön,
 #          Age = Ålder,
 #          Longitude = Bäst_longitud_decimal_grader, # Representative Longitude
 #          Latitude = Bäst_latitud_decimal_grader) # Representative Lat

# Save the export from KUL as not everyone can access it using the databse connection above:
# saveRDS(All_eels_caught_after_2012_KUL_SS, file = "Input_data/Ind_level_data/Ind_level_data_KUL.Rds")

# Read in the data-file extracted from KUL using the code above:
Ind_level_data_KUL <- readRDS("Input_data/Ind_level_data/Ind_level_data_KUL.Rds")

head(Ind_level_data_KUL) # Looks good.

## Merge with LA-ICP-MS data:
LA_ICP_MS_combined_KUL <- merge(x = LA_ICP_MS_KUL, 
                                y = Ind_level_data_KUL,
                                by = "Comb_Year_IndividNr", all.x = TRUE)

## select wanted columns and translate them:
LA_ICP_MS_combined_KUL<- LA_ICP_MS_combined_KUL |> 
  select(Comb_Year_IndividNr, ElapsedTime_s, Ca43_CPS, IntStdWv, Sr_ppm_m88, filename, 
         Year.x, Life_stage, Site, Total_length_mm, Total_weight_g, Sex, Age, Longitude, Latitude) |> 
  rename(Year=Year.x,
         IndividNr=Comb_Year_IndividNr)

## Check import:
n_distinct(LA_ICP_MS_combined_KUL$IndividNr) # n=209, all individuals have been combined

## Done importing and combining LA-ICP-MS data with individual level data from the database KUL ##

### Combine the two dataframes with LA-ICP-MS and individual level data from KUL and Sotebasen ######

str(LA_ICP_MS_combined_KUL) # This is the column-format we want
str(LA_ICP_MS_combined_Sotebasen) # From Sötebasen, select the columns we want, translate from Swedish to English and convert the data formats:

LA_ICP_MS_combined_Sotebasen<-LA_ICP_MS_combined_Sotebasen |> 
  select(IndividNr, ElapsedTime_s, Ca43_CPS, IntStdWv, Sr_ppm_m88, 
         filename, Årtal1, Stadium, VattenNamn, Plats, Längd1, Vikt1, Kön, Ålder, 
         WGS84_N_Vatten, WGS84_E_Vatten, WGS84_N_Plats, WGS84_E_Plats) |> 
  rename(Year=Årtal1,
         Life_stage=Stadium,
         Water_area=VattenNamn,
         Site=Plats,
         Total_length_mm=Längd1,
         Total_weight_g=Vikt1,
         Sex=Kön,
         Age=Ålder,
         Latitude_water_area=WGS84_N_Vatten,
         Longitude_water_area=WGS84_E_Vatten,
         Latitude_site=WGS84_N_Plats,
         Longitude_site=WGS84_E_Plats)
  
# Fix data structure:
LA_ICP_MS_combined_Sotebasen$ElapsedTime_s<-as.numeric(as.character(LA_ICP_MS_combined_Sotebasen$ElapsedTime_s))
LA_ICP_MS_combined_Sotebasen$Ca43_CPS<-as.numeric(as.character(LA_ICP_MS_combined_Sotebasen$Ca43_CPS))
LA_ICP_MS_combined_Sotebasen$Sr_ppm_m88<-as.numeric(as.character(LA_ICP_MS_combined_Sotebasen$Sr_ppm_m88))
LA_ICP_MS_combined_Sotebasen$Latitude_water_area<-as.numeric(as.character(LA_ICP_MS_combined_Sotebasen$Latitude_water_area))
LA_ICP_MS_combined_Sotebasen$Longitude_water_area<-as.numeric(as.character(LA_ICP_MS_combined_Sotebasen$Longitude_water_area))
LA_ICP_MS_combined_Sotebasen$Latitude_site<-as.numeric(as.character(LA_ICP_MS_combined_Sotebasen$Latitude_site))
LA_ICP_MS_combined_Sotebasen$Longitude_site<-as.numeric(as.character(LA_ICP_MS_combined_Sotebasen$Longitude_site))
LA_ICP_MS_combined_Sotebasen$IndividNr<-as.character(LA_ICP_MS_combined_Sotebasen$IndividNr)

# Fix KUL data structure:

## Add the lacking column Water_area, and fix the longitude and latitude columns for the dataframe containing data from KUL so it matches with Sotebasen:
LA_ICP_MS_combined_KUL<-LA_ICP_MS_combined_KUL |> 
  mutate(Water_area="Västerhavet (hela) ICES SD 20-21") |> # Same label as used in Sotebasen: sort(unique(LA_ICP_MS_combined_Sotebasen$Water_area)), and we only have KUl data from this site (STenungsund, Swedish west coast)
  rename(Latitude_site=Latitude,
         Longitude_site=Longitude) |> 
  mutate(Latitude_water_area=NA,  # Add this column that exist in Sotebasen but not in KUL
         Longitude_water_area=NA) |> # Add this column that exist in Sotebasen but not in KUL
  select(IndividNr, ElapsedTime_s, Ca43_CPS, IntStdWv, Sr_ppm_m88, filename, Year, Life_stage, Water_area, Site,
         Total_length_mm, Total_weight_g, Sex, Age, Latitude_water_area, Longitude_water_area, Latitude_site, Longitude_site) # Change the order of columns so they match

# Check so that the columnname matches in the two dataframes:
identical(names(LA_ICP_MS_combined_KUL), names(LA_ICP_MS_combined_Sotebasen)) # TRUE = the columns match

# Combine into one data frame:
LA_ICP_MS_KUL_Sotebasen <- bind_rows(LA_ICP_MS_combined_KUL, LA_ICP_MS_combined_Sotebasen)

# Check combined data:
str(LA_ICP_MS_KUL_Sotebasen) # Structure as we want it to be.
n_distinct(LA_ICP_MS_KUL_Sotebasen$IndividNr) # Should be 1481 (Sotebasen) + 209 (KUL) = 1690 individuals.

# Add a column with the otolith chemistry analytical method:
LA_ICP_MS_KUL_Sotebasen<-LA_ICP_MS_KUL_Sotebasen |> 
  mutate(Analytical_method="LA-ICP-MS")

head(LA_ICP_MS_KUL_Sotebasen) # Looks good.

# A final check that each IndividNr correpsonds to one indvidual by checking for conflicts in individual level data:
LA_ICP_MS_check <- LA_ICP_MS_KUL_Sotebasen %>%
  group_by(IndividNr) %>%
  summarise(across(c(Life_stage, Water_area, Site, Total_length_mm, Total_weight_g,
                     Sex, Age, Latitude_water_area, Longitude_water_area,
                     Latitude_site, Longitude_site), ~ n_distinct(.) > 1,   # TRUE if there is more than one unique value
                     .names = "conflict_{col}"),
                     .groups = "drop") %>%
 filter(if_any(starts_with("conflict_"), ~ . == TRUE))  # keep only rows with conflicts

# If we have any individual with conflicting information, it will be shown here:
LA_ICP_MS_check # No conflicts, the data frame is 0

LA_ICP_MS_check <- LA_ICP_MS_KUL_Sotebasen %>%
  group_by(IndividNr) %>%
  summarise(across(c(Life_stage, Water_area, Site, Total_length_mm, Total_weight_g,
                     Sex, Age, Latitude_water_area, Longitude_water_area,
                     Latitude_site, Longitude_site), ~ n_distinct(.) > 1,   # TRUE if there is more than one unique value
                     .names = "conflict_{col}"),
                     .groups = "drop") %>%
 filter(if_any(starts_with("conflict_"), ~ . == FALSE)) 

LA_ICP_MS_check # Contains 1690 rows, so the loop works and no conflicts, all are FALSE.

# Save the dataframe as a .RDS datafile for the upcoming scripts
saveRDS(LA_ICP_MS_KUL_Sotebasen, file = "R_scripts_and_Generated_data_files/LA_ICP_MS_KUL_Sotebasen.Rds")

#LA_ICP_MS_KUL_Sotebasen = readRDS("LA_ICP_MS_KUL_Sotebasen.Rds")

#------------------ Done with LA-ICP-MS and individual level data import ------------------#

####### Load FE-EMPA otolith chemistry data and combine with Sotebasen ###########

## Path to folder containing the FE-EPMA data
main_path_3 <- "Input_data/FE_EPMA"

# List all .txt files in the directory
files <- list.files(main_path_3, pattern = "\\.txt$", full.names = TRUE)

# Initialize an empty list to store data frames
df_list <- list()

# Loop through each file and import the data
for (file in files) {
  # Read the file, skipping the first two rows
  df <- read.delim(file, header = FALSE, sep = "\t", skip = 1)
  
  # Keep only the first 6 columns and the 17th column containing date and time, needed for sorting the data per transect
  df <- df[, c(1:6, 14:17)]
  
  # Remove rows with header-like entries
  df <- df[!df$V1 %in% c("WT%", "AT%", "COORDINATES", "DATE", "COMMENT"), ]
  
  # Add a column with the filename
  df$filename <- basename(file)    # Add just the file name
  
  # Assign column names
  colnames(df) <- c("COMMENT", "Sr", "Ca", "C", "O", "SUM", "X", "Y", "Z","Date", "filename")
  
  # Append the data frame to the list, using a numeric index
  df_list[[length(df_list) + 1]] <- df   # Use index instead of file name
}

# Combine the data frames generated in the LOOP into one:
combined_df_Fe_EPMA <- do.call(rbind, df_list)
head(combined_df_Fe_EPMA)
# Make a letter lower case and remove hidden white-space etc in the filename:
combined_df_Fe_EPMA$COMMENT <- tolower(trimws(combined_df_Fe_EPMA$COMMENT))

# mMke a column with the IndividID-nr without any letters, used for grouping and merging with individual level data from Sötebasen:
combined_df_Fe_EPMA$IndividNr <- as.integer(gsub("[^0-9]", "", combined_df_Fe_EPMA$COMMENT)) # remove all letters from the IndividNR

## For some individuals, there are more than one transect, denoted as _a (first) and _b (second). 
# A second transect was done if there was something wrong with the first transect. We want to use the second transect for these individuals:

transect_summary <- combined_df_Fe_EPMA %>%
  mutate(suffix = str_extract(COMMENT, "_[ab]$")) %>%
  distinct(IndividNr, suffix) %>%     # one row per individual per suffix
  group_by(IndividNr) %>%
  summarise(has_a = any(suffix == "_a", na.rm = TRUE),
            has_b = any(suffix == "_b", na.rm = TRUE),
            no_suffix = all(is.na(suffix)))

# Summarise it:
transect_summary %>%
  summarise(total_individuals = n(),
            both_ab = sum(has_a & has_b),
            only_a  = sum(has_a & !has_b),
            only_b  = sum(!has_a & has_b),
            no_suffix = sum(no_suffix)) # 2085 out of 2109 individuals have one transect, 5 have two, for which the second is denoted using _b. 
                                        # No ind have only _a (which is what we expect).

# Filter so we only keep the second transect (denoted as "_b") for the 5 individuals with two transects:
combined_df_Fe_EPMA_cleaned <- combined_df_Fe_EPMA %>%
  mutate(suffix = str_extract(COMMENT, "_[ab]$")) %>%
  group_by(IndividNr) %>%
  filter(if (any(suffix == "_b", na.rm = TRUE)) {
    suffix == "_b"} 
    else 
    {TRUE}) %>%
  ungroup() %>%
  select(-suffix)

combined_df_Fe_EPMA_cleaned %>%
  filter(grepl("_a$", COMMENT)) %>%
  distinct(IndividNr) # Should be 0 because there should not be any _a transects in the df.

# Check number of individuals:
combined_df_Fe_EPMA_cleaned %>%
  distinct(IndividNr) %>%
  nrow() #Still 2109 individuals which is correct, we have only removed duplicate transects, not individuals.

# Check number of transect per individual:
combined_df_Fe_EPMA_cleaned %>%
  distinct(IndividNr, COMMENT) %>%
  count(IndividNr) %>%
  filter(n > 1) # Duplicate remain for NA, which should be removed as we cannot link the transect files to any individual:

combined_df_Fe_EPMA_cleaned %>%
  filter(is.na(IndividNr)) %>%
  distinct(COMMENT, filename) # These COMMENTS are hard to convert into individuals, remove them.

#  Remove two the NA transect files:
combined_df_Fe_EPMA_cleaned<-combined_df_Fe_EPMA_cleaned|> 
  filter(!is.na(IndividNr))

# Dubble-check so that we have one transect per individual:
combined_df_Fe_EPMA_cleaned %>%
  distinct(COMMENT) %>%
  nrow() # 2108 unique transects

combined_df_Fe_EPMA_cleaned %>%
  distinct(IndividNr) %>%
  nrow() # 2108 unique individuals

# Check data structure:
str(combined_df_Fe_EPMA_cleaned)

# Fix the data structure:
combined_df_Fe_EPMA_cleaned$Sr<-as.numeric(combined_df_Fe_EPMA_cleaned$Sr)
combined_df_Fe_EPMA_cleaned$Ca<-as.numeric(combined_df_Fe_EPMA_cleaned$Ca)
combined_df_Fe_EPMA_cleaned$C <-as.numeric(combined_df_Fe_EPMA_cleaned$C)
combined_df_Fe_EPMA_cleaned$O <- as.numeric(combined_df_Fe_EPMA_cleaned$O)
combined_df_Fe_EPMA_cleaned$SUM<- as.numeric(combined_df_Fe_EPMA_cleaned$SUM)
combined_df_Fe_EPMA_cleaned$X <-as.numeric(combined_df_Fe_EPMA_cleaned$X)
combined_df_Fe_EPMA_cleaned$Y <- as.numeric(combined_df_Fe_EPMA_cleaned$Y)
combined_df_Fe_EPMA_cleaned$Z<- as.numeric(combined_df_Fe_EPMA_cleaned$Z)
combined_df_Fe_EPMA_cleaned$Date<- as.POSIXct(combined_df_Fe_EPMA_cleaned$Date, format = "%Y-%m-%d %H:%M:%S")
combined_df_Fe_EPMA_cleaned$IndividNr<- as.character(combined_df_Fe_EPMA_cleaned$IndividNr)

str(combined_df_Fe_EPMA_cleaned)
### Combine Fe-EPMA and individual level data in Sötebasen using IndividNr:

FE_EPMA_combined_Sotebasen <- merge(x = combined_df_Fe_EPMA_cleaned, 
                                      y = Sotebasen_inds_unique, 
                                      by = "IndividNr", all.x = TRUE)

## Filter out all rows without plate-nr to ensure that we only use individual with complete data
FE_EPMA_combined_Sotebasen<-FE_EPMA_combined_Sotebasen %>% 
  filter(PlattaNr>=0)

n_distinct(FE_EPMA_combined_Sotebasen$IndividNr) # 1973 individuals with otolith chemistry and indvidual level data

head(FE_EPMA_combined_Sotebasen) # Looks good
tail(FE_EPMA_combined_Sotebasen) # Looks good
str(FE_EPMA_combined_Sotebasen)
### Sort the data according to date and time for each transect:

FE_EPMA_combined_Sotebasen<-FE_EPMA_combined_Sotebasen %>% 
  group_by(IndividNr) %>% 
  arrange(Date) %>% # this is the order fo the transect, expressed in date and time.
  ungroup()

# Select, translate, add analytical method, fix data structure:

FE_EPMA_Sotebasen<-FE_EPMA_combined_Sotebasen |> 
  select(IndividNr, Date, COMMENT, Ca, Sr, X, Y, Z,
         filename, Årtal1, Stadium, VattenNamn, Plats, Längd1, Vikt1, Kön, Ålder, 
         WGS84_N_Vatten, WGS84_E_Vatten, WGS84_N_Plats, WGS84_E_Plats) |> 
  rename(Year=Årtal1,
         Life_stage=Stadium,
         Water_area=VattenNamn,
         Site=Plats,
         Total_length_mm=Längd1,
         Total_weight_g=Vikt1,
         Sex=Kön,
         Age=Ålder,
         Latitude_water_area=WGS84_N_Vatten,
         Longitude_water_area=WGS84_E_Vatten,
         Latitude_site=WGS84_N_Plats,
         Longitude_site=WGS84_E_Plats) |> 
  mutate(Analytical_method="FE-EPMA")

# Fix data classification
FE_EPMA_Sotebasen$Latitude_water_area<-as.numeric(as.character(FE_EPMA_Sotebasen$Latitude_water_area))
FE_EPMA_Sotebasen$Longitude_water_area<-as.numeric(as.character(FE_EPMA_Sotebasen$Longitude_water_area))
FE_EPMA_Sotebasen$Latitude_site<-as.numeric(as.character(FE_EPMA_Sotebasen$Latitude_site))
FE_EPMA_Sotebasen$Longitude_site<-as.numeric(as.character(FE_EPMA_Sotebasen$Longitude_site))

str(FE_EPMA_Sotebasen) # Looks good.

## Additional test to see that import and merge of FE-EPMA and individual-level data from Sötebasen worked:

# Check that each IndividNr correpsonds to one indvidual by checking for conflicts in individual level data:
FE_EPMA_check <- FE_EPMA_Sotebasen %>%
  group_by(IndividNr) %>%
  summarise(across(c(Life_stage, Water_area, Site, Total_length_mm, Total_weight_g,
                     Sex, Age, Latitude_water_area, Longitude_water_area,
                     Latitude_site, Longitude_site), ~ n_distinct(.) > 1,   # TRUE if there is more than one unique value
                     .names = "conflict_{col}"),
                     .groups = "drop") %>%
  filter(if_any(starts_with("conflict_"), ~ . == TRUE))  # keep only rows with conflicts

# See which individuals have conflicting info
FE_EPMA_check # No conflicts.

# We should now have n=30 for every individual, as one transect from the FE-EPMA always should be 30 data-points
FE_EPMA_Sotebasen|>
  count(IndividNr, name = "n_points") |>
  arrange(n_points) |> 
  filter(!n_points==30) # A few that have more or less than n=30, which all should have.

# Have a look at the complete data for these in Sötebasen
FE_EPMA_Sotebasen |> 
  filter(IndividNr %in% c("32595", "31880", "35040", "31633", "34924", "31604", "31607")) |> 
  group_by(filename) |> 
   count(IndividNr, name = "n_points")

## Plot the Sr-concentration of the transects, from the 1st to the last data-point:
FE_EPMA_Sotebasen |> 
  filter(IndividNr %in% c("32595", "31880", "35040", "31633", "34924", "31604", "31607")) |> 
    ggplot(aes(x=Date, y=Sr, group=filename))+
    geom_point()+
    geom_line()+
    facet_wrap(~IndividNr*filename, scales="free")

# For ind 31604 and 31607, both otoliths have been analyzed and data should be in two separate files (PL33 and PL36) which is why we have n=60. The proiles look very similar.
# select the profiles in plate 36, if there was something wrong with the first analysis (data in PL 33).
# For ind 31633 and 34924, the last data-points (the 31st), seems to be off in time. Filter out this data-point for these two individuals:

# Here we remove the two profiles from plate 33 for Individuals 31604 and 31607 and then we remove tha last point for the two individuals with 31 data-points.
FE_EPMA_Sotebasen <- FE_EPMA_Sotebasen |>
  filter(!(IndividNr %in% c("31604", "31607") &
           filename == "HW1_PL33.txt")) |>
  arrange(IndividNr, Date) |>   # Just to ensure that data is arranged properly
  group_by(IndividNr) |>
  filter(!(n() == 31 & row_number() == n())) |> # remove the last data-point for the two inds with 31 data-points
  filter(n() == 30) |> # Remove the three individuals with less than 30 data-points, as all transects derived using FE-EPMA should have 30 data-points
  ungroup()

# We should now have n=30 for every individual
FE_EPMA_Sotebasen|>
  count(IndividNr, name = "n_points") |>
  arrange(n_points) |> 
  filter(!n_points==30) # =0, which it should be.

FE_EPMA_Sotebasen|>
  count(IndividNr, name = "n_points") |>
  arrange(n_points) |> 
  filter(n_points==30) # SHould be 1970 rows, which it is.

n_distinct(FE_EPMA_Sotebasen$IndividNr) # Should be 1970

## Check sampling sites:
sort(unique(FE_EPMA_Sotebasen$Water_area)) # Filter out this site "Ingen uppgift/ej definierat", as it means "No data/not defined"

FE_EPMA_Sotebasen<-FE_EPMA_Sotebasen |> 
  filter(!Water_area=="Ingen uppgift/ej definierat")

# Check updated sample size
n_distinct(FE_EPMA_Sotebasen$IndividNr) # 1958 individuals.

# There we go, now the data table contains one transect (n=30) per individual for 1958 individuals with FE-EPMA and ind. level data with known sampling site.

# Save the dataframe as a .RDS datafile for usage in upcoming scripts
saveRDS(FE_EPMA_Sotebasen, file = "R_scripts_and_Generated_data_files/FE_EPMA_Sotebasen.Rds")

###### Combine individual level data for FE-EPMA and LA-ICP-MS data into one data frame:
names(FE_EPMA_Sotebasen)
names(LA_ICP_MS_KUL_Sotebasen)

# Remove chemistry data columns and combine the two tables:
FE_EPMA_ind_data <- FE_EPMA_Sotebasen |> 
  select(IndividNr, 
         Life_stage,
         Water_area,
         Site,
         Total_length_mm,
         Total_weight_g,
         Sex,
         Age,
         Latitude_water_area,
         Longitude_water_area,
         Latitude_site,
         Longitude_site,
         Analytical_method) |> 
  group_by(IndividNr) |>
  summarise(across(everything(), first), .groups = "drop") # Gives the value in the first row of data for each column per individual

head(FE_EPMA_ind_data) # One row = one individual -> the wanted format.

#LA-ICP-MS data:
LA_ICP_MS_ind_data <- LA_ICP_MS_KUL_Sotebasen |> 
  select(IndividNr, 
         Life_stage,
         Water_area,
         Site,
         Total_length_mm,
         Total_weight_g,
         Sex,
         Age,
         Latitude_water_area,
         Longitude_water_area,
         Latitude_site,
         Longitude_site,
         Analytical_method) |> 
  group_by(IndividNr) |>
  summarise(across(everything(), first), .groups = "drop") # Gives the value in the first row of data for each column per individual

head(LA_ICP_MS_ind_data) # One row = one individual

# Combine the two dfs wiht individual level data:

FE_EPMA_LA_ICP_MS_Individual_data <- bind_rows(LA_ICP_MS_ind_data, FE_EPMA_ind_data)

# Save it as as RDS file.
saveRDS(FE_EPMA_LA_ICP_MS_Individual_data, file = "R_scripts_and_Generated_data_files/FE_EPMA_LA_ICP_MS_Individual_data.Rds")

## This script generates the following three dataframes:

# 1. LA-ICP-MS otolith chemistry transect data combined with individual level data from the databases Sötebasen and KUL:
LA_ICP_MS_KUL_Sotebasen = readRDS("R_Scripts_and_Generated_data_files/LA_ICP_MS_KUL_Sotebasen.Rds") # containing 29 columns with 950 477 rows of data 
head(LA_ICP_MS_KUL_Sotebasen)

# 2. FE-EPMA otolith chemistry transect data combined with individual level data from the database Sötebasen:
FE_EPMA_Sotebasen = readRDS("R_Scripts_and_Generated_data_files/FE_EPMA_Sotebasen.Rds") # containing 22 colummns with 58 740 rows of data
head(FE_EPMA_Sotebasen)

# 3. All individual level data combined for both methods, one row = one individual
FE_EPMA_LA_ICP_MS_Individual_data = readRDS("R_Scripts_and_Generated_data_files/FE_EPMA_LA_ICP_MS_Individual_data.Rds") # containing 13 columns with 3648 rows of data 
head(FE_EPMA_LA_ICP_MS_Individual_data)

#-------------- End of script 1_Input_data.R------------#

