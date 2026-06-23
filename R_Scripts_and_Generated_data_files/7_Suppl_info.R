#################
# Jacobson et al.
# Navigating the salinity gradient: Individual variation in habitat use and migration of European eel revealed by otolith microchemistry 
#################

## Last update: 2026-06-16

#### R-Script: 7_Suppl_info.R

## Load needed packages:
library(tidyverse)
library(flextable)
library(ggmagnify)
library(patchwork)
library(viridis)
library(terra)
library(sf)
library(rnaturalearth)
#library(rnaturalearthdata)
library(ggrepel)
library(ggnewscale)
library(ggspatial)
library(tidytext)  
library(forcats)
library(segclust2d)
library(ggtext)

##########

#### Supplementary material ###

# Create an output folder for the suppl.material
output_dir <- "Output/Suppl_Inf"
dir.create(output_dir, showWarnings = FALSE)

# Fig. S1. Map of reference sites

## Salinity data:

file <- "Input_data/cmems_mod_bal_phy_my_P1Y-m_1773566080862.nc" #sob 2019 Baltic Sea

# Terra automatically finds variables inside the file
r <- rast(file)

sal <- r$sob   # The data we want to plot
sal <- sal[[1]]   # first time step
xmin <- 8; xmax <- 31
ymin <- 53; ymax <- 70

sal_baltic <- crop(sal, ext(xmin, xmax, ymin, ymax))
df_sal <- as.data.frame(sal_baltic, xy = TRUE)
names(df_sal) <- c("lon", "lat", "salinity")

df_sal <- subset(df_sal, salinity < 50)

# Coastlines:
coast <- ne_countries(scale = "medium", returnclass = "sf")
coast_baltic <- st_crop(coast,
                        xmin = xmin, xmax = xmax,
                        ymin = ymin, ymax = ymax)

## Lake data:
## Read in lake shape-file:
lakes_sweden <- read_sf("Input_data/Topography/Swedish rivers and lakes, shapefiles/Lakes2009.shp") # Lakes
print(st_crs(lakes_sweden)$epsg)

## Convert to WGS 84
lakes_sweden_wgs84 <- st_transform(lakes_sweden, crs = 4326) # Lakes
print(st_crs(lakes_sweden_wgs84)$epsg)

# Sample sites of reference eel::
# Load the individual level data from the script 1_Input_data
FE_EMPA_LA_ICP_MS_Individual_data = readRDS("R_Scripts_and_Generated_data_files/FE_EPMA_LA_ICP_MS_Individual_data.Rds") # containing 13 columns with 3567 rows of data 
head(FE_EMPA_LA_ICP_MS_Individual_data)

# Extract the latitude and longitude for each sampling site for coastal sites, each water area for freshwater sites, add missing and fix coordinates:
Reference_Sample_sites_map_figure1<-FE_EMPA_LA_ICP_MS_Individual_data |> 
  mutate(lat = case_when(Water_area %in% c("Östersjön (hela) ICES SD 22-32","Västerhavet (hela) ICES SD 20-21") ~ Latitude_site,
                         TRUE ~ Latitude_water_area),
         long = case_when(Water_area %in% c("Östersjön (hela) ICES SD 22-32", "Västerhavet (hela) ICES SD 20-21") ~ Longitude_site,
                           TRUE ~ Longitude_water_area),
         label=ifelse(Water_area=="Östersjön (hela) ICES SD 22-32", Site,
               ifelse(Water_area=="Västerhavet (hela) ICES SD 20-21", Site, Water_area)),   
         Sampling_site_category=ifelse(Water_area %in% c("Östersjön (hela) ICES SD 22-32", 
                                                         "Västerhavet (hela) ICES SD 20-21"), "Coastal", "Freshwater")) |> 
         filter(label %in% c("Mälaren", "Hjälmaren", 
                             "Stenungsund", "Göta älv", "Viskan", 
                            "Kävlingeån", "Enningdalsälven", 
                             "Nybroån", "Karlshamn - Boön")) |> 
  group_by(label) |> 
    summarise(lat=mean(lat), # Get the latitude per site, for Stenungsund there are several latitudes close to each other, why we use the mean value
              long=mean(long), # get the longitude per site, for Stenungsund there are several longitudes close to each other, why we use the mean value
              Sampling_site_category=unique(Sampling_site_category), # Get the sampling site category 
              .groups = "drop")

sort(unique(Reference_Sample_sites_map_figure1$label)) #

Reference_Sample_sites_map_figure1<-Reference_Sample_sites_map_figure1 |> 
  mutate(Analytical_method=ifelse(label=="Hjälmaren", "FE-EPMA",
                            ifelse(label=="Mälaren", "FE-EPMA and LA-ICP-MS", "LA-ICP-MS")),
        label=ifelse(label=="Karlshamn - Boön", "Karlshamn", label))

## Get colours fo the analytical method points:
plasma_cols <- viridis(7, option = "plasma")

# Create the map:
Fig_S1<-ggplot() +
  geom_raster(data = df_sal, 
              aes(x = lon, 
                  y = lat, 
                  fill = salinity)) + # Salinity data
  scale_fill_viridis(option = "plasma",
                     name = "Salinity (PSU)",
                     guide = guide_colourbar(title.position = "top")) +
   # Reset fill scale
  ggnewscale::new_scale_fill() +
  ## Land on top of raster
  geom_sf(data = coast_baltic,
          fill = "grey80", color = "black") +
   ## Lakes on top of land
   geom_sf(data = lakes_sweden_wgs84, fill="#0d0887", colour=NA)+ # 0 color of the plasma palette
  
  ## Sample-sites and number of eel per site
  geom_point(data = Reference_Sample_sites_map_figure1,
    aes(x = long,
        y = lat,
       shape = Sampling_site_category,
       fill=Analytical_method),
       colour = "white",
       stroke = 1,
       size = 8,
       alpha=0.9) +
## Type of symbol
  scale_shape_manual(values=c(Coastal = 21, 
                              Freshwater = 24),
                              name = "Sampling habitat") +
 scale_fill_manual(values= c("FE-EPMA" = plasma_cols[3],
                            "FE-EPMA and LA-ICP-MS" = plasma_cols[6],
                             "LA-ICP-MS" = plasma_cols[1]),
                    name="Analytical method")+
  ## Labels
  geom_label_repel(data = Reference_Sample_sites_map_figure1,
                    aes(x = long, y = lat, label = label), 
                        colour="black",
                        size = 4,
                        max.overlaps = Inf,
                        box.padding = 0.5,
                        point.padding=0.5,
                        force=4,
                        nudge_x=2,
                        nudge_y=0.6,
                        fill = alpha(c("white"),0.8))+
  ## Coordinate system
  coord_sf(xlim = c(xmin, xmax),
           ylim = c(ymin, ymax),
           expand = FALSE) +
  ## Add a north arrow:
   annotation_north_arrow(location = "br",   # bottom right
                          which_north = "true", # geographic north
                          style = north_arrow_fancy_orienteering(text_size=15),
                          height = unit(3, "cm"),
                          width = unit(2, "cm")) +
  ## Scale_bar
   annotation_scale(location = "br",      # bottom right
                    width_hint = 0.5,     # relative size
                    text_cex = 1) +
  ## Add text labels showing Skagerrak, Kattegat and the Baltic Sea:
geom_text(aes(x = 9.5, y = 58, label = "Skagerrak"),
              size = 4.5,
              fontface = "italic",
              color = "Black",
              angle = 35)+
geom_text(aes(x = 11.25, y = 56.9, label = "Kattegat"),
              size = 4.5,
              fontface = "italic",
              color = "Black",
              angle = -70)+
geom_text(aes(x = 18.5, y = 56.3, label = "Baltic Sea"),
              size = 4.5,
              fontface = "italic",
              color = "Black",
              angle = 50)+
theme_bw(base_size = 16) +
labs(x = "Longitude (°E)",
     y = "Latitude (°N)") +
theme(legend.position = c(0.5, 0.98),   # inside plot, top-center
      legend.justification = c(0.5, 1),
      legend.box = "horizontal",        # legends next to each other
      legend.direction = "horizontal",
      legend.background = element_rect(fill = alpha("white", 0.9), colour = "black"))+
guides(fill = guide_legend(title.position = "top", nrow=3, override.aes = list(shape = 21,
                                                                               colour = "black",
                                                                              size = 5)),
       shape = guide_legend(title.position = "top",
                             nrow = 2,
                             override.aes = list(size = 5, colour="black"),
                             keyheight = unit(1.35, "cm")))
 
Fig_S1

## Save it: 
ggsave("Output/Suppl_Inf/Fig_S1.tiff",
  plot = Fig_S1,
  width = 8,
  height = 12,
  units = "in",
  dpi = 300,
  compression = "lzw")


## Table S2, Reference eel:

## Import data generated in Script 3_Threshold_Values.R

FE_EPMA_Sotebasen_reference_eel = readRDS("R_scripts_and_Generated_data_files/FE_EPMA_Sotebasen_reference_eel.Rds")

FE_EPMA_Sotebasen_reference_eel_Suppl_table<-FE_EPMA_Sotebasen_reference_eel |> 
  group_by(Water_area) |> 
    summarise(Inds=n_distinct(IndividNr),
               Mean_length=mean(Total_length_mm, na.rm=T),
               SD_length=sd(Total_length_mm, na.rm=T),
               Mean_weight=mean(Total_weight_g, na.rm=T),
               SD_weight=sd(Total_weight_g, na.rm=T),
               mean_age=mean(Age, na.rm=T),
               SD_age=sd(Age, na.rm=T)) |> 
  mutate(Analytical_method="FE-EPMA")

FE_EPMA_Sotebasen_reference_eel_Suppl_table

## LA-ICP-MS:
LA_ICP_MS_reference_eel = readRDS("R_scripts_and_Generated_data_files/Reference_eel_LA_ICP_MS.Rds")

sort(unique(Reference_eel_LA_ICP_MS$Water_area))
sort(unique(Reference_eel_LA_ICP_MS$IndividNr))

### Descriptive data table for Supplementary material:

LA_ICP_MS_reference_eel_Suppl_table<-Reference_eel_LA_ICP_MS |> 
  group_by(Water_area) |> 
    summarise(Inds=n_distinct(IndividNr),
               Mean_length=mean(Total_length_mm, na.rm=T),
               SD_length=sd(Total_length_mm, na.rm=T),
               Mean_weight=mean(Total_weight_g, na.rm=T),
               SD_weight=sd(Total_weight_g, na.rm=T),
               mean_age=mean(Age, na.rm=T),
               SD_age=sd(Age, na.rm=T)) |> 
  mutate(across(where(is.numeric), ~ round(.x, 1))) |> # Round all numeric values to 1 digit.
  mutate(Analytical_method="LA-ICP-MS")

LA_ICP_MS_reference_eel_Suppl_table

## Combine the two data-frames:

TableS2<-rbind(FE_EPMA_Sotebasen_reference_eel_Suppl_table, LA_ICP_MS_reference_eel_Suppl_table)

# Fix layout and save as a flextable:

TableS2_ft_data <- TableS2 |>
  mutate(Site = case_when(Water_area == "Västerhavet (hela) ICES SD 20-21" ~ "Stenungsund",
                                   Water_area == "Östersjön (hela) ICES SD 22-32" ~ "Karlshamn",
                                   TRUE ~ Water_area),
         Habitat = case_when(Water_area == "Västerhavet (hela) ICES SD 20-21" ~ "Coastal, Skagerrak",
                             Water_area == "Östersjön (hela) ICES SD 22-32" ~ "Coastal, southern Baltic Sea",
                             TRUE ~ "Freshwater"),
        `Mean length, mm (± 1 SD)` = sprintf("%.1f (± %.1f)", Mean_length, SD_length),
        `Mean weight, g (± 1 SD)` = sprintf("%.1f (± %.1f)", Mean_weight, SD_weight),
        `Mean age (± 1 SD)` = sprintf("%.1f (± %.1f)", mean_age, SD_age)) |>
  transmute(Site,
            n = Inds,
            `Mean length, mm (± 1 SD)`,
            `Mean weight, g (± 1 SD)`,
            `Mean age (± 1 SD)`,
             Habitat,
            `Analytical method` = Analytical_method)

TableS2_ft_data

## Fix formatting:

TableS2_ft <- TableS2_ft_data |>
  flextable() |>
  theme_booktabs() |>
  autofit()

TableS2_ft

# Save it in the correct Output folder:
# As an image:
save_as_image(TableS2_ft, path = "Output/Suppl_Inf/Table_S2.png")

# As a word-file:
save_as_docx(TableS2_ft, path = "Output/Suppl_Inf/Table_S2.docx")

#### Comparing LA-ICP-MS and FE-EPMA data, same individuals #####

### Load the newly aqcuired LA-ICP-MS data for which we have FE-EPMA data since before:

 # Define path where the data is stored:
 Path_Compare_data <- "Input_data/Suppl_LA_ICP_MS_FE_EPMA_Comp"
 # 1. Find all subfolders named "Output_1_TimeSeries"
 Comp_folders <- list.dirs(Path_Compare_data, full.names = TRUE, recursive = TRUE)
 
 # 2. List all .txt files inside those folders
 files <- unlist(lapply
                (Comp_folders, 
                function(x) { list.files(x, pattern = "\\.txt$", full.names = TRUE)
 }))
 
 # 3. Import all files
 df_list <- list()
 
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
 LA_ICP_MS_Compare <- do.call(rbind, df_list)

 # Extract the IndividNr and Year from the filename to make it mergable: 
 LA_ICP_MS_Compare <- LA_ICP_MS_Compare %>%
   mutate(filename_clean = sub("\\.txt$", "", filename), # Remove the .txt suffix of the file-name
          IndividNr = as.integer(gsub("[^0-9]", "", filename_clean))) # Save only the number from the filename
      
 head(LA_ICP_MS_Compare)

LA_ICP_MS_Compare |>
  group_by(IndividNr) |> 
     summarise(Data_points=n()) |> 
  arrange(Data_points) |> 
  print(n=Inf)

LA_ICP_MS_Compare<-LA_ICP_MS_Compare |> 
  mutate(Analythical_method="LA-ICP-MS")

## Calculate the Sr:Ca ratio for each transect, combine data and plot the transects:

# First, identify and specify the zero-point corresponding to the centre of the otolith core:

Scan_speed_LA_ICP_MS<-6 # 6 micrometer/s was the scan-speed for all of these inds.

LA_ICP_MS_Compare_core<-LA_ICP_MS_Compare|> 
 group_by(IndividNr) |> 
   mutate(Ca_ppm_calc=as.numeric(as.character(Ca43_CPS, na.rm=T))*380000/mean(as.numeric(as.character(Ca43_CPS, na.rm=T))),
          Sr_Ca_ppm=(as.numeric(as.character(Sr_ppm_m88, na.rm=T))*1000)/Ca_ppm_calc,
          Distance=ElapsedTime_s*Scan_speed_LA_ICP_MS) |> 
  filter(Ca_ppm_calc > (mean(Ca_ppm_calc) - 1.25 * sd(Ca_ppm_calc)), 
         Ca_ppm_calc < (mean(Ca_ppm_calc) + 1.25 * sd(Ca_ppm_calc))) |>
  arrange(ElapsedTime_s, .by_group = TRUE) |> # make sure the data is arrange according to elapsed time for each individual, otherwise, the rolling mean will not work
   mutate(SrCa_smooth = zoo::rollmean(Sr_Ca_ppm, 
                                      k = 8, 
                                      fill = "extend", # To avoid NA values in the first and last data-points, extending the first and last value of the window
                                      align = "center"),
         Core_signal_smooth = ifelse(SrCa_smooth > 9.13, "Core", "Outside_core"), # Use the smoothed Sr:Ca ratio to identify the regions of the transect correspondng to the core
         Core_signal_raw=ifelse(Sr_Ca_ppm > 9.13, "Core", "Outside_core")) |> # Use the raw Sr:Ca ratio to identify the regions of the transect correspondng to the core
  ungroup()

# Step 2: Find the core section of the transect using cumsum to identify contiguous Sr:Ca regions of the transects above 9.12 of the smoothed Sr:Ca ratio:
LA_ICP_MS_Compare_core <- LA_ICP_MS_Compare_core|>
  arrange(IndividNr, ElapsedTime_s) |> 
  group_by(IndividNr) |>
  mutate(Core_run = cumsum(Core_signal_smooth != lag(Core_signal_smooth, default = first(Core_signal_smooth)))) 

# Step 2: Look at how many contiguous regions that have been identified:
LA_ICP_MS_core_regions <- LA_ICP_MS_Compare_core |>
  filter(Core_signal_smooth == "Core") |> # Select only the regions with a core signal
  group_by(IndividNr, Core_run) |> # group the data per individual and identified egions with core-signals
  summarise(
    n_points = n(), # number of data-points per core region
    mean_SrCa_smooth = mean(SrCa_smooth, na.rm = TRUE), 
    max_SrCa_smooth = max(SrCa_smooth, na.rm = TRUE),
    min_time = min(ElapsedTime_s),
    max_time = max(ElapsedTime_s),
    center_time = (min_time + max_time)/2, # this is the center-point, works because ElapsedTime_s is in seconds
    .groups = "drop")

# Filter out the longest contiguous region (most consequitive points) with a Sr:Ca value >9.13 of each transect:
LA_ICP_MS_core_centers <- LA_ICP_MS_core_regions |>
  filter(n_points > 20) |> # Remove spikes and short segments to avoiding estimate for otoliths lacking a clear core Sr:Ca signal.
  arrange(IndividNr, desc(n_points), desc(max_SrCa_smooth)) |> # Sort the estimates per individual, so the first row is the estimated center based on the longest contiguous region. 
                                                         # If there are several estimates with the sanme number of data points, sort on highest Sr:Ca value
  group_by(IndividNr) |> 
  slice(1) |> # Keep only the first row of data per individual
  ungroup() |> 
  mutate(Core_center_time = center_time)

# Add the time corresponding to the Core center time for each individual:
core_centers_merge<-LA_ICP_MS_core_centers  |> 
  select(IndividNr, Core_center_time)

LA_ICP_MS_Compare_zero_point<- merge(x = LA_ICP_MS_Compare_core, 
                                     y = core_centers_merge, 
                                     by = "IndividNr", all.x = TRUE)

# Check the derived zero-point:
LA_ICP_MS_Compare_zero_point|> 
    ggplot(aes(y = Sr_Ca_ppm,
                  x = as.numeric(as.character(ElapsedTime_s)))) +
      geom_point() +
      facet_wrap(~IndividNr)+
      geom_vline(aes(xintercept=Core_center_time))+
      ylab("Sr:Ca Ratio") +
      xlab("Elapsed_time") # Looks good!

# Caluclate distance from core center, set to be 0.
LA_ICP_MS_Compare_zero_point_Distance<-LA_ICP_MS_Compare_zero_point |> 
  group_by(IndividNr) |> 
  arrange(IndividNr,ElapsedTime_s) |> 
  mutate(Distance=ElapsedTime_s*Scan_speed_LA_ICP_MS,
         Zero_point_micrometer=Core_center_time*Scan_speed_LA_ICP_MS,
         Distance_from_core = Distance - Zero_point_micrometer)

# Plot it to check that it worked:
LA_ICP_MS_Compare_zero_point_Distance|> 
    ggplot(aes(y = Sr_Ca_ppm,
                  x = Distance_from_core)) +
      geom_point() +
      facet_wrap(~IndividNr)+
      ylab("(Sr/Ca)*1000") +
      xlab("Distance_from_core") # Looks good!

#### Assign habitat for each data-point:
LA_ICP_MS_Compare_Life_history <- LA_ICP_MS_Compare_zero_point_Distance|> 
          mutate(Habitat_use = case_when(Distance_from_core < 513 ~ "Excluded",
                                         Sr_Ca_ppm >=0 & Sr_Ca_ppm <= 1.67 ~ "Freshwater",
                                         Sr_Ca_ppm > 1.67 & Sr_Ca_ppm <= 2.60  ~ "Estuary",
                                         Sr_Ca_ppm > 2.60 & Sr_Ca_ppm <= 9.13 ~ "Coastal/Brackish",
                                         Sr_Ca_ppm > 9.13 ~"Marine/Core"))

## Assign life history:
LA_ICP_MS_Compare_Assinged_life_history_summary<-LA_ICP_MS_Compare_Life_history |>
  group_by(Analythical_method, IndividNr) |>
  summarise(has_fresh = any(Habitat_use == "Freshwater", na.rm = TRUE),
            has_estuary = any(Habitat_use == "Estuary", na.rm = TRUE),
            has_brackish = any(Habitat_use == "Coastal/Brackish", na.rm = TRUE),
            has_marine = any(Habitat_use == "Marine/Core", na.rm = TRUE),
           .groups = "drop") |>
 mutate(Life_history = case_when(has_fresh & !has_estuary & !has_brackish & !has_marine ~ "Freshwater", # Only freshwater
                                 has_fresh & has_estuary & !has_brackish & !has_marine ~ "Freshwater/Estuary",  # Only freshwater + estuary
                                 has_fresh & (has_brackish | has_marine) ~ "Habitat shifter",   # Freshwater + Coastal and/or marine
                                 has_estuary & !has_fresh & !has_brackish & !has_marine ~ "Estuary", # Only estuary (no freshwater, no coastal)
                                 !has_fresh ~ "Coastal", # No freswhater singal (all estuarine, coastal/brackish and marine combos)
                                 TRUE ~ NA_character_))

print(LA_ICP_MS_Compare_Assinged_life_history_summary, n=Inf) # All classified as habitat shifters.

## Plot the data:

LA_ICP_MS_Compare_Life_history |> 
  mutate(Habitat_simple = case_when(Habitat_use %in% c("Freshwater", "Estuary") ~ "Freshwater",
                                    Habitat_use %in% c("Coastal/Brackish", "Marine/Core") ~ "Coastal",
                                    Habitat_use == "Excluded" ~"Excluded",
         TRUE ~ NA_character_)) |> 
  ggplot(aes(x=Distance_from_core, y=Sr_Ca_ppm))+
   geom_line()+
  geom_point(aes(colour=Habitat_simple))+
  facet_wrap(~IndividNr+Analythical_method)+
  geom_hline(yintercept=2.60)+ # Threshold value for Freshwater+Estuary
  theme_bw(base_size=18)

## Habitat assignement Fe-EPMA:
# Extract this data from the FE-EPMA dataframe generated in script 1.Input_data.R
FE_EPMA_Sotebasen = readRDS("R_Scripts_and_Generated_data_files/FE_EPMA_Sotebasen.Rds") 

FE_EPMA_Sotebasen_Compare<-FE_EPMA_Sotebasen |> 
  filter(IndividNr %in% c("31584", "31585","31586", "31587", "31588", "31591", "31596", "31599", 
                          "31600", "31601", "31602", "31606", "31616", "31617", "31618", "31620", 
                          "31621", "31623", "31625", "31627", "33916", "33917", "33919", "33920", "33924"))

FE_EPMA_Sotebasen_Compare |> 
    group_by(IndividNr) |> 
     summarise(Data_points=n()) |> 
  print(n=Inf) # There we go, 30 points per ind.

FE_EPMA_Sotebasen_Compare<-FE_EPMA_Sotebasen_Compare |> 
  group_by(IndividNr) %>%
  arrange(Date, .by_group = TRUE) %>%
  mutate(Sr_Ca_ratio = as.numeric(as.character(Sr)) / as.numeric(as.character(Ca)) * 1000,
         Step_dist = sqrt((X - lag(X))^2 + (Y - lag(Y))^2),  # Distance between consecutive points -> distance between measure-points
         Distance_from_start = cumsum(replace_na(Step_dist, 0))) %>%  # Cumulative distance from first point which is set to 0 (the center of the core), which is the start of the transect
  ungroup()

## Assign habitat use for each data-point:
FE_EPMA_Sotebasen_Compare<-FE_EPMA_Sotebasen_Compare |> 
   filter(!is.na(Sr_Ca_ratio)) %>%
    mutate(Habitat_use = case_when(Distance_from_start < 515 ~ "Excluded",
                                   Sr_Ca_ratio <= 1.12 ~ "Freshwater",
                                   Sr_Ca_ratio > 1.12 & Sr_Ca_ratio <= 6.97 ~ "Coastal/Brackish",
                                   Sr_Ca_ratio > 6.97 ~"Marine/Core"))

## Assign life-history:

FE_EPMA_Sotebasen_Compare_Assigned_Life_history<-FE_EPMA_Sotebasen_Compare |>
  group_by(Analytical_method, IndividNr) |>
summarise(has_fresh = any(Habitat_use == "Freshwater", na.rm = TRUE),
          has_brackish = any(Habitat_use == "Coastal/Brackish", na.rm = TRUE),
          has_marine = any(Habitat_use == "Marine/Core", na.rm = TRUE),
         .groups = "drop") |>
  mutate(Life_history = case_when(has_fresh & !has_brackish ~ "Freshwater",
                                  has_fresh & (has_brackish | has_marine) ~ "Habitat shifter",
                                  !has_fresh ~ "Coastal"))

print(FE_EPMA_Sotebasen_Compare_Assigned_Life_history, n=Inf) # All classified as habitat shifters

# Plot the data:
FE_EPMA_Sotebasen_Compare|>
  mutate(Habitat_simple = case_when(Habitat_use == "Freshwater" ~ "Freshwater",
                                    Habitat_use %in% c("Coastal/Brackish", "Marine/Core") ~ "Coastal",
                                    Habitat_use=="Excluded" ~ "Excluded",
                                    TRUE ~ NA_character_)) %>%
  ggplot(aes(x=Distance_from_start, y=Sr_Ca_ratio))+
   geom_line()+
  geom_point(aes(colour=Habitat_simple))+
  facet_wrap(~IndividNr+Analytical_method)+
  geom_hline(yintercept=1.12)+
  theme_bw(base_size=18)

### All eel classified as mixed, irrespective of method.

#### ---------- 

## Quantify and add the number of shifts between freshwater and coastal habitats:

## FE-EPMA:

Habitat_shifts_FE_EPMA_Compare_summary <- FE_EPMA_Sotebasen_Compare |> 
  filter(Distance_from_start>=513) |> 
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

Habitat_shifts_FE_EPMA_Compare_summary

# Summary table:
Habitat_shifts_FE_EPMA_Compare_summary |> 
  group_by(n_shifts_total) |> 
   summarise(n=n_distinct(IndividNr)) |> 
    mutate(tot_n=sum(n))

### Quantify habitat shifts for LA-ICP-MS data:

  LA_ICP_MS_Compare_Habitat_shift_summary <- LA_ICP_MS_Compare_Life_history |> 
  arrange(IndividNr, Distance_from_core) |> # Make sure the data is sorted correctly
  filter(Distance_from_core>=513) |> 
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

  print(LA_ICP_MS_Compare_Habitat_shift_summary, n=Inf)

  LA_ICP_MS_Compare_Habitat_shift_summary |> 
  group_by(n_shifts_total) |> 
   summarise(n=n_distinct(IndividNr)) |> 
    mutate(tot_n=sum(n))


 LA_ICP_MS_Compare_Life_history|> 
  mutate(Habitat_simple = case_when(Habitat_use %in% c("Freshwater", "Estuary") ~ "Freshwater/Estuary",
         Habitat_use %in% c("Coastal/Brackish","Marine/Core") ~ "Coastal",
         TRUE ~ "Excluded")) |> 
  filter(IndividNr=="31585") |> 
  ggplot(aes(x=Distance_from_core, y=Sr_Ca_ppm))+
   geom_line()+
  geom_point(aes(colour=Habitat_simple))+
  facet_wrap(~IndividNr+Analythical_method)+
  geom_hline(yintercept=2.60)+ # Threshold value for freshwater+Estuary
  theme_bw(base_size=18)

### Make the results table:
# FE-EPMA:
head(FE_EPMA_Sotebasen_Compare_Assigned_Life_history)

FE_EPMA_Assigned_life_history<-FE_EPMA_Sotebasen_Compare_Assigned_Life_history |> 
  mutate(Assigned_Life_history_FE_EPMA=Life_history) |> 
   select(IndividNr, Assigned_Life_history_FE_EPMA)

FE_EPMA_Assigned_life_history

# LA-ICP-MS:
head(LA_ICP_MS_Compare_Assinged_life_history_summary)

LA_ICP_MS_Assigned_life_history<-LA_ICP_MS_Compare_Assinged_life_history_summary |> 
  mutate(Assigned_Life_history_LA_ICP_MS=Life_history) |> 
   select(IndividNr, Assigned_Life_history_LA_ICP_MS)

LA_ICP_MS_Assigned_life_history

# Merge:
FE_EPMA_LA_ICP_MS_Assigned_life_history<-merge(FE_EPMA_Assigned_life_history, LA_ICP_MS_Assigned_life_history, by="IndividNr")

FE_EPMA_LA_ICP_MS_Assigned_life_history

## Add number of shifts between freshwater and coastal areas:

# LA_ICP_MS:

LA_ICP_MS_Compare_Habitat_shift_summary_Only_tot<-LA_ICP_MS_Compare_Habitat_shift_summary |> 
  mutate(Number_of_habitat_shifts_LA_ICP_MS=n_shifts_total) |> 
   select(IndividNr, Number_of_habitat_shifts_LA_ICP_MS)

LA_ICP_MS_Compare_Habitat_shift_summary_Only_tot

# FE-EPMA:

Habitat_shifts_FE_EPMA_Compare_summary_Only_tot<-Habitat_shifts_FE_EPMA_Compare_summary |> 
   mutate(Number_of_habitat_shifts_FE_EPMA=n_shifts_total) |> 
   select(IndividNr, Number_of_habitat_shifts_FE_EPMA)

Habitat_shifts_FE_EPMA_LA_ICP_MS<-merge(Habitat_shifts_FE_EPMA_Compare_summary_Only_tot, 
                                        LA_ICP_MS_Compare_Habitat_shift_summary_Only_tot,
                                        by="IndividNr")
Habitat_shifts_FE_EPMA_LA_ICP_MS

## Merge

Table_Suppl<-merge(FE_EPMA_LA_ICP_MS_Assigned_life_history, Habitat_shifts_FE_EPMA_LA_ICP_MS, by="IndividNr")

# Add ind level data:

names(FE_EPMA_Sotebasen)

FE_EPMA_Sotebasen_Compare_Ind_data<-FE_EPMA_Sotebasen |> 
  filter(IndividNr %in% c("31584", "31585","31586", "31587", "31588", "31591", "31596", "31599", 
                          "31600", "31601", "31602", "31606", "31616", "31617", "31618", "31620", 
                          "31621", "31623", "31625", "31627", "33916", "33917", "33919", "33920", "33924")) |> 
  select(IndividNr, Water_area, Total_length_mm, Total_weight_g) |> 
  group_by(IndividNr) |> 
  slice(1) # Take the first row for each ind.

   FE_EPMA_Sotebasen_Compare_Ind_data

Table_Suppl_Complete<-merge(FE_EPMA_Sotebasen_Compare_Ind_data, Table_Suppl, by="IndividNr")

Table_Suppl_Complete

Table_S3_ft<-flextable(Table_Suppl_Complete) # There we go.

Table_S3_ft

## Save the table:
# As an image:
save_as_image(Table_S3_ft, path = "Output/Suppl_Inf/Table_S3.png")

# As a word-file:
save_as_docx(Table_S3_ft, path = "Output/Suppl_Inf/Table_S3.docx")

### Make a figure for appendix:

## Combine and show in one plot both profiles:

FE_EPMA_Sotebasen_Compare_33916<-FE_EPMA_Sotebasen_Compare|>
  mutate(Habitat_simple = case_when(Habitat_use == "Freshwater" ~ "Freshwater",
                                    Habitat_use %in% c("Coastal/Brackish", "Marine/Core") ~ "Coastal",
                                    Habitat_use=="Excluded" ~ "Not included",
                                    TRUE ~ NA_character_)) %>%
  filter(IndividNr=="33916") |> 
  select(IndividNr, Distance_from_start, Sr_Ca_ratio, Analytical_method) |> 
  rename(Distance=Distance_from_start,
         Sr_Ca_Ratio=Sr_Ca_ratio,
         Method=Analytical_method) |> 
    mutate(IndividNr=as.character(IndividNr))


LA_ICP_MS_Sotebasen_Compare_33916<-LA_ICP_MS_Compare_Life_history|> 
  mutate(Habitat_simple = case_when(Habitat_use %in% c("Freshwater", "Estuary") ~ "Freshwater/Estuary",
         Habitat_use %in% c("Coastal/Brackish","Marine/Core") ~ "Coastal",
         TRUE ~ "Not included")) |> 
  filter(IndividNr=="33916") |>
  select(IndividNr, Distance_from_core, Sr_Ca_ppm, Analythical_method) |> 
  rename(Distance=Distance_from_core,
         Sr_Ca_Ratio=Sr_Ca_ppm,
         Method=Analythical_method) |> 
  mutate(IndividNr=as.character(IndividNr))

Compare_33916_plot<-rbind(LA_ICP_MS_Sotebasen_Compare_33916, FE_EPMA_Sotebasen_Compare_33916)
unique(Compare_33916_plot$Method)

Compare_33916_plot<-Compare_33916_plot |> 
  mutate(Method=ifelse(Method=="Fe-EPMA", "FE-EPMA", Method))

Compare_33916_plot_save<-Compare_33916_plot |> 
ggplot(aes(x=Distance, y=Sr_Ca_Ratio))+
  geom_line(data = subset(Compare_33916_plot, Method == "LA-ICP-MS"), aes(colour = Method)) +
  geom_line(data = subset(Compare_33916_plot, Method == "FE-EPMA"), aes(colour = Method)) +
  ylim(0,30)+
  xlim(-400, 2000)+
  theme_classic(base_size=18)+
  labs(colour="Analytical method")+
   ylab("Sr:Ca ratio")+
   xlab("Distance from core (µm)")+
  scale_colour_manual(values = c("LA-ICP-MS"="Orange",
                                 "FE-EPMA"="Dodgerblue"),
                      breaks = c("FE-EPMA", "LA-ICP-MS"))+
    annotate("text", x = -Inf, y = Inf, label = "C", 
           hjust = -0.1, vjust = 1.1, size = 8, fontface="bold")+
  theme(legend.position="top")

Compare_33916_plot_save

### Habiats shift:

 LA_ICP_33916<-LA_ICP_MS_Compare_Life_history|> 
  mutate(Habitat_simple = case_when(Habitat_use %in% c("Freshwater", "Estuary") ~ "Freshwater/Estuary",
         Habitat_use %in% c("Coastal/Brackish","Marine/Core") ~ "Coastal",
         TRUE ~ "Excluded")) |> 
  mutate(Habitat_simple = factor(Habitat_simple, levels = c("Excluded", "Coastal", "Freshwater/Estuary"))) |> # Set the order of apperance
  filter(IndividNr=="33916")

Habitat_simple<- levels(LA_ICP_33916$Habitat_simple) 
cols_LA_ICP_MS <- c("Excluded" = "grey",
        setNames(rev(viridis(length(Habitat_simple) - 1, option = "plasma", begin = 0, end = 0.85)),
                                    Habitat_simple[Habitat_simple != "Excluded"]))


 LA_ICP_33916_Plot<- LA_ICP_33916 |> 
  ggplot(aes(x=Distance_from_core, y=Sr_Ca_ppm))+
   geom_line()+
  geom_point(aes(colour=Habitat_simple))+
  geom_hline(yintercept=1.67)+
   geom_hline(yintercept=2.60)+
  theme_classic(base_size=18)+
  scale_colour_manual(values = cols_LA_ICP_MS, name="Assigned\nhabitat") +
  geom_magnify(from = c(xmin=1160, xmax=1210, ymin=1, ymax=4), # 930 -10000
               to=c(1050, 2000, 14, 28), axes="y")+
   ylab("Sr:Ca ratio")+
   xlab("Distance from core (µm)")+
   labs(colour = NULL)+
   annotate("text", x = 0.1, y = 1.67+0.5, label ="Sr:Ca = 1.67, Freshwater", hjust = 0, size=4)+
    annotate("text", x = 0.1, y = 2.60+0.5, label ="Sr:Ca = 2.60, Estuary", hjust = 0, size=4)+
   ylim(0,30)+
   xlim(0,2000)+
      annotate("text", x = -Inf, y = Inf, label = "E", 
           hjust = -0.1, vjust = 1.1, size = 8, fontface="bold")+
   theme(legend.position="top")

 LA_ICP_33916_Plot

# FE-EPMA:
FE_EPMA_33916_Plot<-FE_EPMA_Sotebasen_Compare|>
  mutate(Habitat_simple = case_when(Habitat_use == "Freshwater" ~ "Freshwater",
                                    Habitat_use %in% c("Coastal/Brackish", "Marine/Core") ~ "Coastal",
                                    Habitat_use=="Excluded" ~ "Excluded",
                                    TRUE ~ NA_character_)) %>%
  mutate(Habitat_simple = factor(Habitat_simple,
                                    levels = c("Excluded", "Coastal", "Freshwater"))) |> # Set the order of apperance
  filter(IndividNr=="33916")

Habitat_simple_FE_EPMA<- levels(FE_EPMA_33916_Plot$Habitat_simple) 
cols_FE_EPMA <- c("Excluded" = "grey",
        setNames(rev(viridis(length(Habitat_simple_FE_EPMA) - 1, option = "plasma", begin = 0, end = 0.85)),
                                    Habitat_simple_FE_EPMA[Habitat_simple_FE_EPMA != "Excluded"]))

  FE_EPMA_33916_Plot<-FE_EPMA_33916_Plot |> 
  ggplot(aes(x=Distance_from_start, y=Sr_Ca_ratio))+
   geom_line()+
  geom_point(aes(colour=Habitat_simple))+
 # facet_wrap(~IndividNr+Analytical_method)+
  geom_hline(yintercept=1.12)+
  theme_classic(base_size=18)+
  scale_colour_manual(values = cols_FE_EPMA, name="Assigned\nhabitat") +
   ylab("Sr:Ca ratio")+
   xlab("Distance from core (µm)")+
   labs(colour = NULL)+
   annotate("text", x = 0.1, y = 1.12+0.5, label = "Sr:Ca = 1.12, Freshwater", hjust = 0, size=4)+
   ylim(0,30)+
   xlim(0,2000)+
  annotate("text", x = -Inf, y = Inf, label = "D", 
           hjust = -0.1, vjust = 1.1, size = 8, fontface="bold")+
  theme(legend.position="top")

FE_EPMA_33916_Plot

## Cobine the three plots:
Fig_S2 <- Compare_33916_plot_save / (FE_EPMA_33916_Plot+ LA_ICP_33916_Plot)

Fig_S2

# Save the plot:
ggsave("Output/Suppl_Inf/Fig_S2.tiff",
  plot = Fig_S2,
  width = 14,
  height = 14,
  units = "in",
  dpi = 300,
  compression = "lzw")

### Make a plot showing all Sr:Ca ratio transects, one facet per individual:

# FE-EPMA data:
FE_EPMA_Sotebasen_Compare_All_inds<-FE_EPMA_Sotebasen_Compare|>
  mutate(Habitat_simple = case_when(Habitat_use == "Freshwater" ~ "Freshwater",
                                    Habitat_use %in% c("Coastal/Brackish", "Marine/Core") ~ "Coastal",
                                    Habitat_use=="Excluded" ~ "Not included",
                                    TRUE ~ NA_character_)) %>%
  select(IndividNr, Distance_from_start, Sr_Ca_ratio, Analytical_method) |> 
  rename(Distance=Distance_from_start,
         Sr_Ca_Ratio=Sr_Ca_ratio,
         Method=Analytical_method) |> 
    mutate(IndividNr=as.character(IndividNr))

# LA-ICP-MS data:
LA_ICP_MS_Sotebasen_Compare_All_inds<-LA_ICP_MS_Compare_Life_history|> 
  mutate(Habitat_simple = case_when(Habitat_use %in% c("Freshwater", "Estuary") ~ "Freshwater/Estuary",
         Habitat_use %in% c("Coastal/Brackish","Marine/Core") ~ "Coastal",
         TRUE ~ "Not included")) |> 
  select(IndividNr, Distance_from_core, Sr_Ca_ppm, Analythical_method) |> 
  rename(Distance=Distance_from_core,
         Sr_Ca_Ratio=Sr_Ca_ppm,
         Method=Analythical_method) |> 
  mutate(IndividNr=as.character(IndividNr))

# Combine the data-tables:
Compare_All_inds_plot<-rbind(LA_ICP_MS_Sotebasen_Compare_All_inds, FE_EPMA_Sotebasen_Compare_All_inds)
unique(Compare_All_inds_plot$Method)


Compare_All_inds_plot<-Compare_All_inds_plot |> 
  mutate(Method=ifelse(Method=="Fe-EPMA", "FE-EPMA", Method))

# Fix spelling of analythical method
Fig_S3<-Compare_All_inds_plot |> 
ggplot(aes(x=Distance, y=Sr_Ca_Ratio))+
  geom_line(data = subset(Compare_All_inds_plot, Method == "LA-ICP-MS"),
            aes(colour = Method)) +
  geom_line(data = subset(Compare_All_inds_plot, Method == "FE-EPMA"),
            aes(colour = Method)) +
  ylim(0,30)+
  xlim(-500, 2000)+
  theme_classic(base_size=36)+
  labs(colour="Analytical method")+
  ylab("Sr:Ca ratio")+
  xlab("Distance from core (µm)")+
  scale_colour_manual(values = c("LA-ICP-MS"="Orange",
                                 "FE-EPMA"="Dodgerblue"),
                      breaks = c("FE-EPMA", "LA-ICP-MS"))+
  theme(legend.position="top")+
  facet_wrap(~IndividNr, ncol=3)

Fig_S3

# Save the plot:
ggsave("Output/Suppl_Inf/Fig_S3.tiff",
  plot = Fig_S3,
  width = 20,
  height = 32,
  units = "in",
  dpi = 600,
  compression = "lzw")
