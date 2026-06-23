#################
# Jacobson et al.
# Navigating the salinity gradient: Individual variation in habitat use and migration of European eel revealed by otolith microchemistry 
#################

## Last update: 2026-06-16

#### R-Script: 6_Figures.R

## Load needed packages:

library(terra)
library(ggplot2) # In tidyverse
library(sf)
library(viridis)
library(rnaturalearth)
#library(rnaturalearthdata)
library(ggrepel)
library(dplyr) # In tidyverse
library(ggnewscale)
library(ggspatial)
library(tidytext)  
library(forcats) # In tidyverse
library(patchwork)
library(segclust2d)
library(ggtext)

######

## Figure 1 - Map of sampling sites, number of eel per site and the salinity gradient of the Baltic Sea region ####

# Sea floor salinity data from CMEMS: 
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

# Sample sites and number of eel per site:
# Load the individual level data from the script 1_Input_data
FE_EMPA_LA_ICP_MS_Individual_data = readRDS("R_Scripts_and_Generated_data_files/FE_EPMA_LA_ICP_MS_Individual_data.Rds") # containing 13 columns with 3648 rows of data 

# Extract the latitude and longitude for each sampling site for coastal sites, each water area for freshwater sites, add missing and fix coordinates:
Sample_sites_map_figure1<-FE_EMPA_LA_ICP_MS_Individual_data |> 
  mutate(lat = case_when(Site == "Kattegatt (ICES SD 21)" ~ 55.872206, # Kattegatt (ICES SD 21) should be in Roskildefjorden in Denmark according to free-text comment in Sötebasen
                         Water_area == "Ringsjöarna, Skåne" ~ 55.873999, # Ringsjöarna lack coordinates in Sötebasen
                         Water_area %in% c("Östersjön (hela) ICES SD 22-32","Västerhavet (hela) ICES SD 20-21") ~ Latitude_site,
                         TRUE ~ Latitude_water_area),
          long = case_when(Site == "Kattegatt (ICES SD 21)" ~ 12.031596, # Kattegatt (ICES SD 21) should be in Roskildefjorden in Denmark according to free-text comment in Sötebasen
                           Water_area == "Ringsjöarna, Skåne" ~ 13.506483, # Ringsjöarna lack coordinates in Sötebasen
                           Water_area %in% c("Östersjön (hela) ICES SD 22-32", "Västerhavet (hela) ICES SD 20-21") ~ Longitude_site,
                           TRUE ~ Longitude_water_area),
         label=ifelse(Water_area=="Östersjön (hela) ICES SD 22-32", Site,
               ifelse(Water_area=="Västerhavet (hela) ICES SD 20-21", Site, Water_area)),   
         Sampling_site_category=ifelse(Water_area %in% c("Östersjön (hela) ICES SD 22-32", 
                                                         "Västerhavet (hela) ICES SD 20-21"), "Coastal", "Freshwater")) |> 
  group_by(label) |> 
    summarise(n=n_distinct(IndividNr), # number of individuals
              lat=mean(lat), # Get the latitude per site, for Stenungsund there are several latitudes close to each other, why we use the mean value
              long=mean(long), # get the longitude per site, for Stenungsund there are several longitudes close to each other, why we use the mean value
              Sampling_site_category=unique(Sampling_site_category), # Get the sampling site category 
              .groups = "drop")

sum(Sample_sites_map_figure1$n) # Should be 3648 which it is.

# Create the map:
Figure_1 <- ggplot() +
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
  geom_point(data = Sample_sites_map_figure1,
             aes(x = long,
                 y = lat,
                 fill = n,
                 shape = Sampling_site_category),
                 colour = "black",
                 stroke = 1,
                 size = 6,
                 alpha=0.9) +
## Type of symbol
  scale_shape_manual(values=c(Coastal = 21, 
                              Freshwater = 24),
                              name = "Sampling habitat") +
  # Gradient fill colour 
 scale_fill_gradient(low="#d7b5d8", 
                     high="forestgreen",
                     name = "Number of eel",
                     breaks=c(100, 300, 500),
                     guide = guide_colourbar(title.position = "top"))+
  
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
  ## Add scale bar
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
      guides(fill = guide_colourbar(title.position = "top"),
             shape = guide_legend(title.position = "top", nrow = 1,
             override.aes = list(size = 5),
             keyheight = unit(1.35, "cm")))
 
Figure_1

## Save it: 
ggsave("Output/Figure_1.tiff",
  plot = Figure_1,
  width = 8,
  height = 12,
  units = "in",
  dpi = 300,
  compression = "lzw")


## Figure 2 ####

## Figure 2 are photos of analyzed otoliths, not generated using this script and already in the output folder.

## Figure 3 ####

## Segclust output plots, one for FE_EPMA and one for LA-ICP-MS:

## FE-EPMA data:
# Load Sr:Ca and individual leve data from Script: 1_Input_data.R
FE_EPMA_Sotebasen = readRDS("R_Scripts_and_Generated_data_files/FE_EPMA_Sotebasen.Rds")

# Filter out one of the reference eels:
FE_EPMA_ref_Ind_32079<-data.frame(FE_EPMA_Sotebasen |> 
  filter(IndividNr == "32079") |> 
  select("IndividNr", "Sr", "Ca", "COMMENT", "Date", "X", "Y") %>% 
  mutate(Sr_Ca_ratio=(as.numeric(as.character(Sr)) / as.numeric(as.character(Ca))) * 1000))

# Perform segmentation using the segclust function
mode_segclust_FE_EPMA_ref_Ind_32079<- segclust(FE_EPMA_ref_Ind_32079,
                                               Kmax = 6, 
                                               lmin = 5, 
                                               ncluster = 2:6,
                                               seg.var = c("Sr_Ca_ratio"))

plot_BIC(mode_segclust_FE_EPMA_ref_Ind_32079) # 3 clusters, 3 segments according to BIC
plot(mode_segclust_FE_EPMA_ref_Ind_32079) # BIC-selected output plot we want to make publishable

## Extract values from the segclust function for the 3 cluster, 3 segment config:
segclust_output_FE_EPMA <- mode_segclust_FE_EPMA_ref_Ind_32079$outputs$`3 class - 3 segments`
result <- merge(segclust_output_FE_EPMA$segments[, c("begin", "end", "state")], # Start, end and range of each segment
                 segclust_output_FE_EPMA$states[, c("state", "mu.Sr_Ca_ratio", "sd.Sr_Ca_ratio")], by = "state") # Mean and SD in Sr:Ca per segment
result <- result[order(result$begin), ] # Order the results
result # There we go, now combine with the Sr:Ca transect data in a plot:

## Make A-F panels:

# Subplot A:
Figure3_A<-FE_EPMA_ref_Ind_32079 |> 
  mutate(Row_number=row_number()) |> 
  ggplot(aes(x=Row_number, y=Sr_Ca_ratio))+
  geom_line()+ # Sr:Ca ratio
  geom_rect(data = result,    # SD per segment
            aes(xmin = begin, 
                xmax = end,
                ymin = mu.Sr_Ca_ratio - sd.Sr_Ca_ratio, 
                ymax = mu.Sr_Ca_ratio + sd.Sr_Ca_ratio,
                fill = factor(state)),
                alpha = 0.5,
                inherit.aes = FALSE) +
  geom_segment(data = result,   # mean value as a line (per segment)
              aes(x = begin, 
                   xend = end,
                   y = mu.Sr_Ca_ratio, 
                   yend = mu.Sr_Ca_ratio,
                   color = factor(state)),
                   linewidth = 1.2,
                   inherit.aes = FALSE) +
    scale_fill_viridis_d(option = "plasma")+
    scale_colour_viridis_d(option = "plasma")+
   labs(color = "Cluster", fill = "Cluster")+
  theme_classic(base_size=22) +
  ggtitle("FE-EPMA")+
  ylab("Sr:Ca ratio")+
  xlab("Transect data-point")+
  annotate("text", x = -Inf, y = Inf, label = " A", 
           hjust = -0.1, vjust = 1.1, size = 8, fontface="bold")+
  theme(plot.title = element_text(hjust = 0.5))

Figure3_A

# Subplot B

## LA-ICP-MS
# Load transect data from script 1_Input_data.R
LA_ICP_MS_KUL_Sotebasen = readRDS("R_Scripts_and_Generated_data_files/LA_ICP_MS_KUL_Sotebasen.Rds")

# Filter out one of the reference eels:
LA_ICP_MS_ref_Ind_50161 <- LA_ICP_MS_KUL_Sotebasen |> 
  filter(IndividNr == "50161") |> 
   mutate(Ca_ppm_calc=as.numeric(as.character(Ca43_CPS, na.rm=T))*380000/mean(as.numeric(as.character(Ca43_CPS, na.rm=T))),#mmol mol–1 
          Sr_Ca_ppm=(as.numeric(as.character(Sr_ppm_m88, na.rm=T))*1000)/Ca_ppm_calc) |> 
  filter(Ca_ppm_calc > (mean(Ca_ppm_calc) - 1.25 * sd(Ca_ppm_calc)), 
         Ca_ppm_calc < (mean(Ca_ppm_calc) + 1.25 * sd(Ca_ppm_calc))) |>
  ungroup() |> 
     select(IndividNr, Sr_Ca_ppm, Ca_ppm_calc, ElapsedTime_s, Water_area) %>%
     filter(!is.na(Sr_Ca_ppm)) %>% # Remove NAs if present
     mutate(ElapsedTime_s = as.numeric(ElapsedTime_s)) %>%    # convert to numeric
     arrange(ElapsedTime_s) %>%   #Sort
     as.data.frame()
  
  #### Run the segclust modell:
  mode_segclust_LA_ICP_MS_ref_Ind_50161 <- segclust(LA_ICP_MS_ref_Ind_50161,
                                                    Kmax = 20,
                                                    lmin = 10,
                                                    ncluster = 2:5,
                                                    seg.var = "Sr_Ca_ppm")

plot_BIC(mode_segclust_LA_ICP_MS_ref_Ind_50161)
plot(mode_segclust_LA_ICP_MS_ref_Ind_50161) # 5 clusters, 11 segments in total.

# Extract the segclust model output and plot it in a ggplot:

segclust_output_LA_ICP_MS<-mode_segclust_LA_ICP_MS_ref_Ind_50161$outputs$`5 class - 11 segments`
result_LA_ICP_MS <- merge(segclust_output_LA_ICP_MS$segments[, c("begin", "end", "state")], # Start, end and range of each segment
                segclust_output_LA_ICP_MS$states[, c("state", "mu.Sr_Ca_ppm", "sd.Sr_Ca_ppm")], 
                by = "state") # Mean and SD in Sr:Ca per segment
result_LA_ICP_MS <- result_LA_ICP_MS[order(result_LA_ICP_MS$begin), ] # Order the results
result_LA_ICP_MS

# Create the plot:

head(segclust_output_LA_ICP_MS)

Figure3_B<-LA_ICP_MS_ref_Ind_50161|> 
  mutate(Row_number=row_number()) |> 
  ggplot(aes(x=Row_number, y=Sr_Ca_ppm))+
  geom_line(linewidth=0.5)+ # Sr:Ca ratio
  geom_rect(data = result_LA_ICP_MS,    # SD per segment
            aes(xmin = begin, 
                xmax = end,
                ymin = mu.Sr_Ca_ppm - sd.Sr_Ca_ppm, 
                ymax = mu.Sr_Ca_ppm + sd.Sr_Ca_ppm,
                fill = factor(state)),
                alpha = 0.5,
                inherit.aes = FALSE) +
  geom_segment(data = result_LA_ICP_MS,   # mean value as a line (per segment)
              aes(x = begin, 
                   xend = end,
                   y = mu.Sr_Ca_ppm, 
                   yend = mu.Sr_Ca_ppm,
                   color = factor(state)),
                   linewidth = 1.2,
                   inherit.aes = FALSE) +
  labs(color = "Cluster", fill = "Cluster")+
      scale_fill_viridis_d(option = "plasma")+
    scale_colour_viridis_d(option = "plasma")+
  theme_classic(base_size=22) +
  ggtitle("LA-ICP-MS")+
  ylab("Sr:Ca ratio")+
  xlab("Transect data-point")+
 annotate("text", x = -Inf, y = Inf, label = " B", 
          hjust = -0.1, vjust = 1.1, size = 8, fontface="bold")+
  theme(plot.title = element_text(hjust = 0.5))

Figure3_B

###

# Subplot C

Results_threshold_FE_EPMA = readRDS("R_Scripts_and_Generated_data_files/Results_threshold_FE_EPMA.Rds") # Generated in Script 3_Threshold_values.R

Summary_per_segment_ref_FE_EPMA<-Results_threshold_FE_EPMA %>% 
  group_by(state) %>% 
  summarise(mu_Min=min(mu.Sr_Ca_ratio, na.rm = TRUE),
            mu_Max=max(mu.Sr_Ca_ratio, na.rm = TRUE),
            Mean_n=n(),
            sd_Min=min(sd.Sr_Ca_ratio, na.rm = TRUE),
            sd_Max=max(sd.Sr_Ca_ratio, na.rm = TRUE),
            sd_n=n())

Summary_per_segment_ref_FE_EPMA

## Make the two subplots with the Sr:Ca ratios per segment:

Figure3_C<-Results_threshold_FE_EPMA |> 
  ggplot(aes(x=mu.Sr_Ca_ratio, y=sd.Sr_Ca_ratio))+
    geom_point(aes(fill=as.factor(state)), 
                   shape=21, 
                   alpha=0.8, 
                   size=5,
                   show.legend=FALSE) +
    geom_rect(data = Summary_per_segment_ref_FE_EPMA,
              aes(xmin = mu_Min,
                  xmax = mu_Max,
                  ymin = sd_Min,
                  ymax = sd_Max,
                  color = as.factor(state)),
                  fill = NA,
                  linewidth = 1.25,
                  inherit.aes = FALSE) +
   labs(fill="Cluster number",
        colour="Min-max cluster\nboundaries")+
   theme_classic(base_size=22)+
   ylab("1 S.D. of Sr:Ca ratio")+
   xlab("Mean Sr:Ca ratio")+
   scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, 5))+
     annotate("text", x = -Inf, y = Inf, label = " C", 
           hjust = -0.1, vjust = 1.1, size = 8, fontface="bold")+
      scale_fill_viridis_d(option = "plasma")+
      scale_colour_viridis_d(option = "plasma")

Figure3_C

# Subplot D

# LA-ICP-MS data:
## 5th and 95th percentile mean Sr:Ca ~ SD Sr:Ca:

Results_threshold_LA_ICP_MS=readRDS("R_Scripts_and_Generated_data_files/Results_threshold_LA_ICP_MS.Rds") ## Generated in Script 3_Threshold_values.R

head(Results_threshold_LA_ICP_MS)
sort(unique(Results_threshold_LA_ICP_MS$Water_area))

Results_threshold_LA_ICP_MS_Habitat<-Results_threshold_LA_ICP_MS |> 
  mutate(Catch_habitat=ifelse(Water_area=="Östersjön (hela) ICES SD 22-32","Coastal_Baltic",
                        ifelse(Water_area=="Västerhavet (hela) ICES SD 20-21","Coastal_Skagerrak", "Freshwater"))) |> 
  filter(Catch_habitat=="Freshwater") |>
  filter(ncluster_BIC == 5) |>
     group_by(state) %>% 
  summarise(mu_5 = quantile(mu.Sr_Ca_ppm, 0.05, na.rm = TRUE),
            mu_95 = quantile(mu.Sr_Ca_ppm, 0.95, na.rm = TRUE),
            sd_5 = quantile(sd.Sr_Ca_ppm, 0.05, na.rm = TRUE),
            sd_95 = quantile(sd.Sr_Ca_ppm, 0.95, na.rm = TRUE))

# Threshold values (5th and 95th percentile))
Results_threshold_LA_ICP_MS_Habitat

## Subset coastal eel:
Results_threshold_LA_ICP_MS_habitat_coast<-Results_threshold_LA_ICP_MS|> 
    mutate(Catch_habitat=ifelse(Water_area=="Östersjön (hela) ICES SD 22-32","Coastal, southern Baltic Sea",
                        ifelse(Water_area=="Västerhavet (hela) ICES SD 20-21","Coastal, Skagerrak", "Freshwater"))) |> 
 # mutate(VattenNamn_ordered=recode(VattenNamn_ordered, "Stenungssund"="Stenungsund")) |> 
  filter(Catch_habitat %in% c("Coastal, southern Baltic Sea", "Coastal, Skagerrak")) |>
  filter(ncluster_BIC == 5) |>
     group_by(Catch_habitat, state) %>% 
  summarise(mu_5 = quantile(mu.Sr_Ca_ppm, 0.05, na.rm = TRUE),
           mu_95 = quantile(mu.Sr_Ca_ppm, 0.95, na.rm = TRUE),
           sd_5 = quantile(sd.Sr_Ca_ppm, 0.05, na.rm = TRUE),
           sd_95 = quantile(sd.Sr_Ca_ppm, 0.95, na.rm = TRUE))


Figure3_D<-Results_threshold_LA_ICP_MS |> 
  mutate(Catch_habitat=ifelse(Water_area=="Östersjön (hela) ICES SD 22-32","Coastal, southern Baltic Sea",
                       ifelse(Water_area=="Västerhavet (hela) ICES SD 20-21","Coastal, Skagerrak", "Freshwater"))) |> 
  filter(Catch_habitat=="Freshwater") |>
  filter(ncluster_BIC == 5) |>
  ggplot(aes(x=mu.Sr_Ca_ppm, y=sd.Sr_Ca_ppm))+
    geom_point(aes(fill=as.factor(state)), 
                   shape=21, 
                   alpha=0.8, 
                   size=5,
                   show.legend=FALSE) +
    geom_rect(data = Results_threshold_LA_ICP_MS_Habitat,
              aes(xmin = mu_5,
                  xmax = mu_95,
                  ymin = sd_5,
                  ymax = sd_95,
                  color = as.factor(state)),
                  fill = NA,
                  linewidth = 1.25,
                  inherit.aes = FALSE) +
   geom_rect(data = Results_threshold_LA_ICP_MS_habitat_coast, # Add the coastal data
             aes(xmin = mu_5,
                 xmax = mu_95,
                 ymin = sd_5,
                 ymax = sd_95,
                 color = as.factor(state),
                 linetype=Catch_habitat),
                 fill = NA,
                 linewidth = 1.5,
                 inherit.aes = FALSE)+
  scale_linetype_manual(values = c("Coastal, southern Baltic Sea"="dashed",
                                   "Coastal, Skagerrak"="dotted")) +
  guides(colour = guide_legend(order = 1),
         linetype = guide_legend(order = 2,
                                 title = "Coastal sampling site",
                                 override.aes = list(color = "black",
                                                    linewidth = 0.6,
                                                     shape = NA)))+
   labs(fill="Cluster of Sr:Ca\nratio transect",
        colour = "5<sup>th</sup> – 95<sup>th</sup> cluster<br>percentile boundaries")+ # ggtext needed for this
   theme_classic(base_size=22)+
   ylab("1 S.D. of Sr:Ca ratio")+
   xlab("Mean Sr:Ca ratio")+
   scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, 5))+
     annotate("text", x = -Inf, y = Inf, label = " D", 
               hjust = -0.1, vjust = 1.1, size = 8, fontface="bold")+
        scale_fill_viridis_d(option = "plasma")+
        scale_colour_viridis_d(option = "plasma")+
 theme(legend.title = ggtext::element_markdown(hjust = 0),
       legend.text.align = 0)

Figure3_D

## Subplot E

head(Results_threshold_FE_EPMA)

Figure3_E<-Results_threshold_FE_EPMA  |> 
  mutate(Catch_habitat="Freshwater") |> 
  group_by(Catch_habitat, state_ordered) |> 
  summarise(min = min(mu.Sr_Ca_ratio),
            median = quantile(mu.Sr_Ca_ratio, 0.5),
            max = max(mu.Sr_Ca_ratio)) |> 
  ggplot(aes(x = factor(state_ordered),
             y = median, 
             fill = Catch_habitat, 
             group = Catch_habitat,
             colour=factor(state_ordered))) +
  geom_ribbon(aes(ymin = min, 
                  ymax = max), 
                  alpha=0.5, colour="black") +
  geom_errorbar(aes(ymin = min, 
                    ymax = max),
                    width = 0.05,
                    linewidth = 2.4,
                    colour="white") +
  geom_errorbar(aes(ymin = min, 
                    ymax = max),
                    width = 0.05,
                    linewidth = 1.6,
                    colour="black") +
  labs(x = "Cluster of Sr:Ca transect", 
       y = "Range of mean Sr:Ca ratio\nvalues per transect cluster",
       fill = "Sampling habitat")+
      # colour = "Cluster")+
  theme_classic(base_size=22)+
  geom_hline(yintercept=1.12)+ # Freshwater upper limit
  geom_hline(yintercept=7.27)+ # Marine
  ## Add text to each threshold area:
  annotate("text", x = 3.5, y = 0.3, label = "Freshwater",  fontface = "italic", hjust = 1, size=7)+
  annotate("text", x = 3.5, y = 3, label = "Coastal/Brackish",  fontface = "italic", hjust = 1, size=7)+
  annotate("text", x = 1, y = 10, label = "Marine", fontface = "italic", hjust = 0, size=7)+
  ## Add the threshold values in Sr:Ca
  annotate("text", x = 0.5, y = 1.12+0.4, label = "Sr:Ca = 1.12", hjust = 0, size=5)+
  annotate("text", x = 0.5, y = 7.27+0.4, label = "Sr:Ca = 6.97", hjust = 0, size=5)+
  # Add panel label
  annotate("text", x = -Inf, y = Inf, label = " E", 
           hjust = -0.1, vjust = 1.1, size = 8, fontface="bold")+
   scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 5))+
   scale_fill_viridis_d(option = "plasma", name = "Sampling habitat") +
   #scale_colour_viridis_d(option = "plasma", name = "Cluster")+
  guides(fill = guide_legend(order = 1))

Figure3_E


### Show the distribtion of data points per segment using range-bars:

Figure3_F <- Results_threshold_LA_ICP_MS |>  
  mutate(Catch_habitat=ifelse(Water_area=="Östersjön (hela) ICES SD 22-32","Coastal, southern Baltic Sea",
                       ifelse(Water_area=="Västerhavet (hela) ICES SD 20-21","Coastal, Skagerrak", "Freshwater"))) |> 
  mutate(Catch_habitat = factor(Catch_habitat, 
                                levels = c("Freshwater",
                                           "Coastal, southern Baltic Sea",
                                           "Coastal, Skagerrak"))) |> 
  filter(ncluster_BIC == 5) |>
  group_by(Catch_habitat, state) |> 
  summarise(mu_5 = quantile(mu.Sr_Ca_ppm, 0.05, na.rm = TRUE),
            median = quantile(mu.Sr_Ca_ppm, 0.5, na.rm = TRUE),
            mu_95 = quantile(mu.Sr_Ca_ppm, 0.95, na.rm = TRUE),
            .groups = "drop") |> 
  ggplot(aes(x = factor(state),
             y = median,
             group = Catch_habitat,
             fill = Catch_habitat,
             colour = Catch_habitat)) +
  geom_ribbon(aes(ymin = mu_5, ymax = mu_95, fill=Catch_habitat),
              alpha = 0.5, colour="black") +
    geom_errorbar(aes(ymin = mu_5, 
                    ymax = mu_95),
                    position = position_dodge(width = 0.1),
                    width = 0.3,
                    linewidth = 2.4,
                    colour="white") +
  geom_errorbar(aes(ymin = mu_5, ymax = mu_95),
                position = position_dodge(width = 0.1),
                width = 0.3,
                linewidth = 1.6,
                colour="black") +
  # Axes and labels:
  labs(x="Cluster of Sr:Ca transect",
       y = "Range of mean Sr:Ca ratio values per<br>transect cluster (5<sup>th</sup>–95<sup>th</sup> percentile)",
       fill = "Sampling habitat",
       colour = "Cluster") +
  scale_fill_viridis_d(option = "plasma", name = "Sampling habitat") +
 # scale_colour_viridis_d(option = "plasma", name = "Cluster") +
  theme_classic(base_size = 22) +
  # Horizontal thresholds;
  geom_hline(yintercept = 1.67) +
  geom_hline(yintercept = 2.60) +
  #geom_hline(yintercept = 3.90) +
  #geom_hline(yintercept = 4.42) +
  geom_hline(yintercept = 9.13) +
  # Annotations:
  annotate("text", x = 3.5, y = 0.83, label = "Freshwater", fontface = "italic", hjust = 0, size = 7) +
  annotate("text", x = 3.5, y = 2.1, label = "Estuary", fontface = "italic", hjust = 0, size = 7) +
  annotate("text", x = 0.5, y = 8, label = "Coastal/Brackish", fontface = "italic", hjust = 0, size = 7) +
  annotate("text", x = 1, y = 12, label = "Marine", fontface = "italic", hjust = 0, size = 7) +
  # Sr:Ca threshold labels:
  annotate("text", x = 0.45, y = 1.67 + 0.4, label = "Sr:Ca = 1.67", hjust = 0, size = 5) +
  annotate("text", x = 0.45, y = 2.60 + 0.4, label = "Sr:Ca = 2.60", hjust = 0, size = 5) +
  annotate("text", x = 0.45, y = 9.13 + 0.4, label = "Sr:Ca = 9.13", hjust = 0, size = 5) +
  # Panel label:
  annotate("text", x = -Inf, y = Inf, label = " F", hjust = -0.1, vjust = 1.1, size = 8, fontface = "bold")+
  guides(fill = guide_legend(order = 1))+
  theme(axis.title.y = ggtext::element_markdown())
  
Figure3_F

Figure_3<-(Figure3_A + Figure3_B) / (Figure3_C + Figure3_D) / (Figure3_E + Figure3_F)

Figure_3

ggsave("Output/Figure_3.tiff",
  plot = Figure_3,
  width = 20,
  height = 20,
  units = "in",
  dpi = 300,
  compression = "lzw")

#### Figure 4 ####

## FE-EPMA:
FE_EPMA_Sotebasen_Distance = readRDS("R_Scripts_and_Generated_data_files/FE_EPMA_Sotebasen_Distance.Rds")

# Conver Sr:Ca to habitat use
FE_EPMA_Sotebasen_Distance_Habitat<-FE_EPMA_Sotebasen_Distance |> 
    mutate(Habitat_use = case_when(Distance_from_start < 513 ~ "Excluded",
                                   Sr_Ca_ratio <= 1.12 ~ "Freshwater",
                                   Sr_Ca_ratio > 1.12 & Sr_Ca_ratio <= 6.97 ~ "Coastal/Brackish",
                                   Sr_Ca_ratio > 6.97 ~"Marine/Core"))

# Select representative profiles for assigned life-histories:
Selected_examples_FE_EPMA<-FE_EPMA_Sotebasen_Distance_Habitat|> 
  filter(IndividNr %in% c("33645","35324","34144")) |> 
  mutate(Analytical_method="FE-EPMA",
        Site_plot=ifelse(Water_area=="Östersjön (hela) ICES SD 22-32", "Coastal Baltic Sea,",
                  ifelse(Water_area=="Mälaren", "Freshwater lake,",
                  ifelse(Water_area=="Ätran", "Freshwater river,", "ERROR"))),
        facet_group = factor(paste(Site_plot, Analytical_method)))

sort(unique(Selected_examples_FE_EPMA$facet_group))

## Sort sampling site according to assigned life-history 
Selected_examples_FE_EPMA$facet_group<- factor(Selected_examples_FE_EPMA$facet_group, # Sort according to salinity gradient, high -> low
                   levels = c("Coastal Baltic Sea, FE-EPMA", "Freshwater river, FE-EPMA", "Freshwater lake, FE-EPMA"))

## Fix order of assigned habitats, high to low salinity:
Selected_examples_FE_EPMA$Habitat_use<- factor(Selected_examples_FE_EPMA$Habitat_use,
                   levels = c("Excluded", 
                              "Marine/Core", 
                              "Coastal/Brackish", 
                              "Freshwater"))

# Save factor-order for colour scheme
levels_hab_FE_EPMA <- levels(Selected_examples_FE_EPMA$Habitat_use) 

## Creat a color vector:
habitat_cols <- c("Excluded" = "grey",
                  "Marine/Core" = rev(viridis(4, option = "plasma"))[1],
                  "Coastal/Brackish" = rev(viridis(4, option = "plasma"))[2],
                  "Estuary" = rev(viridis(4, option = "plasma"))[3],
                  "Freshwater" = rev(viridis(4, option = "plasma"))[4])

## Create facet-text, one for each assigned life-history:
facet_labels_FE_EPMA <- data.frame(facet_group = c("Coastal Baltic Sea, FE-EPMA", 
                                                   "Freshwater river, FE-EPMA",
                                                   "Freshwater lake, FE-EPMA"),
                           label = c("italic('Coastal resident')",
                                     "italic('Freshwater resident')",
                                     "italic('Habitat shifter')"),
                           x = Inf,
                           y = Inf)

facet_labels_FE_EPMA$facet_group <- factor(facet_labels_FE_EPMA$facet_group,
                          levels = levels(Selected_examples_FE_EPMA$facet_group))

Figure_4_FE_EPMA <- Selected_examples_FE_EPMA |>
  ggplot(aes(y = Sr_Ca_ratio,
             x = Distance_from_start)) +
  geom_line()+
  geom_point(aes(colour = Habitat_use), size=2) +
  ylab("Sr:Ca ratio") +
  xlab("Distance from otolith core (µm)") +
  theme_classic(base_size = 22) +
  facet_wrap(~facet_group) +
  scale_colour_manual(values = habitat_cols, name="Assigned habitat") +
  scale_x_continuous(breaks = c(0, 500, 1000, 1500, 2000),
                   labels = c("0", "500", "1000", "1500", "2000"),
                   limits = c(0, 2000))+
  theme(legend.position = "top")+
  geom_text(data = facet_labels_FE_EPMA,
            aes(x = 100, y = 15, label = label),
            inherit.aes = FALSE,
            hjust = -0.1, vjust = 1.5,
            parse = TRUE,
            size=6)+
  geom_hline(yintercept=1.12, linetype="dashed")+
  ylim(0, 20)

Figure_4_FE_EPMA

## LA-ICP-MS:

LA_ICP_MS_KUL_Sotebasen_Distance = readRDS("R_Scripts_and_Generated_data_files/LA_ICP_MS_KUL_Sotebasen_Distance.Rds")

# Convert Sr:Ca ratio to habitat use:
LA_ICP_MS_KUL_Sotebasen_Life_hist <- LA_ICP_MS_KUL_Sotebasen_Distance|> 
          mutate(Habitat_use = case_when(Distance_from_core < 513 ~ "Excluded",
                                         Sr_Ca_ppm >=0 & Sr_Ca_ppm <= 1.67 ~ "Freshwater",
                                         Sr_Ca_ppm > 1.67 & Sr_Ca_ppm <= 2.60  ~ "Estuary",
                                         Sr_Ca_ppm > 2.60 & Sr_Ca_ppm <= 9.13 ~ "Coastal/Brackish",
                                         Sr_Ca_ppm > 9.13 ~"Marine/Core"))

# Check classification:
Selected_examples<-LA_ICP_MS_KUL_Sotebasen_Life_hist|> 
  filter(IndividNr %in% c("48474","50161","2018_33")) # Select three characteristic plots for assigned life-histories
sort(unique(Selected_examples$Water_area))

Selected_examples$Habitat_use<- factor(Selected_examples$Habitat_use,
                   levels = c("Excluded", "Marine/Core", "Coastal/Brackish", "Estuary", "Freshwater"))

Selected_examples_LA_ICP_MS<-Selected_examples|> 
 mutate(Site_plot=ifelse(Water_area=="Västerhavet (hela) ICES SD 20-21", "Coastal Kattegat,", 
                  ifelse(Water_area=="Viskan", "Freshwater river,",
                  ifelse(Water_area=="Mälaren", "Freshwater lake,"))), # Fix sampling site label
        Analytical_method="LA-ICP-MS",
        facet_group = factor(paste(Site_plot, Analytical_method)))

sort(unique(Selected_examples_LA_ICP_MS$facet_group))

## Sort sampling site according to assigned life-history 
Selected_examples_LA_ICP_MS$facet_group<- factor(Selected_examples_LA_ICP_MS$facet_group, # Sort according to salinity gradient, high -> low
                   levels = c("Coastal Kattegat, LA-ICP-MS", "Freshwater river, LA-ICP-MS", "Freshwater lake, LA-ICP-MS"))

levels_hab <- levels(Selected_examples_LA_ICP_MS$Habitat_use) # To ensure thath colour scheme follow the factor order

# Create color vector
cols <- c("Excluded" = "grey",
          setNames(rev(viridis(length(levels_hab) - 1, option = "plasma")),
                                      levels_hab[levels_hab != "Excluded"]))

## Create facet-text, one for each assigned life-history:
facet_labels <- data.frame(facet_group = c("Coastal Kattegat, LA-ICP-MS", 
                                            "Freshwater river, LA-ICP-MS", 
                                            "Freshwater lake, LA-ICP-MS"),
                           label = c("italic('Coastal resident')",
                                     "italic('Freshwater resident')",
                                     "italic('Habitat shifter')"),
                           x = Inf,
                           y = Inf)

facet_labels$facet_group <- factor(facet_labels$facet_group,
                          levels = levels(Selected_examples_LA_ICP_MS$facet_group))

Figure_4_LA_ICP_MS <- Selected_examples_LA_ICP_MS |>
  ggplot(aes(y = Sr_Ca_ppm,
             x = Distance_from_core)) +
  geom_line()+
  geom_point(aes(colour = Habitat_use), size=2) +
  ylab("Sr:Ca ratio") +
  xlab("Distance from otolith core (µm)") +
  theme_classic(base_size = 22) +
  facet_wrap(~facet_group) +
  scale_colour_manual(values = habitat_cols, name="Assigned habitat") +
  scale_x_continuous(breaks = c(0, 500, 1000, 1500),
                   labels = c("0", "500", "1000", "1500"),
                   limits = c(0, 1500))+
  theme(legend.position = "top")+
  geom_text(data = facet_labels,
            aes(x = 100, y = 15, label = label),
            inherit.aes = FALSE,
            hjust = -0.1, vjust = 1.5,
            parse = TRUE,
            size=6)+
  geom_hline(yintercept=1.67, linetype="dashed")+
  ylim(0,25)

Figure_4_LA_ICP_MS

## Combine the two plots:

Figure_4 <-Figure_4_FE_EPMA / Figure_4_LA_ICP_MS

Figure_4

## Save it:
ggsave("Output/Figure_4.tiff",
       plot = Figure_4,
       width = 14,
       height = 12,
       units = "in",
       dpi = 300,
       compression = "lzw")

### Figure 5 ####

## The proportion of life-histories per site, per major geographical region, sorted from north to south:

# Read in assigned life-history data from the R-script: 5_Habitat_assignment.R
LA_ICP_MS_KUL_Sotebasen_Life_history = readRDS("R_Scripts_and_Generated_data_files/LA_ICP_MS_KUL_Sotebasen_Life_history.Rds")
FE_EPMA_Sotebasen_Life_history = readRDS("R_Scripts_and_Generated_data_files/FE_EPMA_Sotebasen_Life_history.Rds")

## Select Water_area, site, IndividNr and assigned-lifehistory and combine the LA_ICP_MS and the FE_EPMA data, 
# then merge with Indivudual-level data:
LA_ICP_MS_KUL_Selected_Life_history<-LA_ICP_MS_KUL_Sotebasen_Life_history |> 
  select(Water_area, Site, IndividNr, Life_history) |> 
  mutate(Analythical_method="LA-ICP-MS") # Add method

FE_EPMA_Sotebasen_Selected_Life_history<-FE_EPMA_Sotebasen_Life_history |>
  select(Water_area, Site, IndividNr, Life_history) |> 
  mutate(Analythical_method="FE-EPMA") # Add method

Selected_life_history_complete<- rbind(LA_ICP_MS_KUL_Selected_Life_history, FE_EPMA_Sotebasen_Selected_Life_history)

head(Selected_life_history_complete)
# Combine with individual level data:

# Load the individual level data from the R-script: 1_Input_data.R
FE_EPMA_LA_ICP_MS_Individual_data = readRDS("R_Scripts_and_Generated_data_files/FE_EPMA_LA_ICP_MS_Individual_data.Rds") 
head(FE_EPMA_LA_ICP_MS_Individual_data)

Selected_life_history_complete_Ind_data<-merge(x = Selected_life_history_complete, 
                                            y = FE_EPMA_LA_ICP_MS_Individual_data, 
                                            by = "IndividNr", all.x = TRUE)

# Check merge:
all(Selected_life_history_complete_Ind_data$Water_area.x == Selected_life_history_complete_Ind_data$Water_area.y) # Complete match
all(Selected_life_history_complete_Ind_data$Site.x == Selected_life_history_complete_Ind_data$Site.y) # Complete match

# Group the data based on major geographical region, sites and sort in accordance to the general salinity gradient:
# (North -> South for Kattegat and Skagerrak); South -> North fro the Baltic Sea, using the same method as for the Map of all sample sizes (Figure 3, the map).

names(Selected_life_history_complete_Ind_data)

# Extract the latitude and longitude for each sampling site for coastal sites, each water area for freshwater sites, add missing and fix coordinates:
Life_history_proportions_per_site<-Selected_life_history_complete_Ind_data |> 
  rename(Site=Site.x,
        Water_area=Water_area.x) |> 
  mutate(lat = case_when(Site == "Kattegatt (ICES SD 21)" ~ 55.872206, # Kattegatt (ICES SD 21) should be in Roskildefjorden in Denmark according to free-text comment in Sötebasen
                         Water_area == "Ringsjöarna, Skåne" ~ 55.873999, # Ringsjöarna lack coordinates in Sötebasen
                         Water_area %in% c("Östersjön (hela) ICES SD 22-32","Västerhavet (hela) ICES SD 20-21") ~ Latitude_site,
                         TRUE ~ Latitude_water_area),
         long = case_when(Site == "Kattegatt (ICES SD 21)" ~ 12.031596, # Kattegatt (ICES SD 21) should be in Roskildefjorden in Denmark according to free-text comment in Sötebasen
                           Water_area == "Ringsjöarna, Skåne" ~ 13.506483, # Ringsjöarna lack coordinates in Sötebasen
                           Water_area %in% c("Östersjön (hela) ICES SD 22-32", "Västerhavet (hela) ICES SD 20-21") ~ Longitude_site,
                           TRUE ~ Longitude_water_area),
         Site = ifelse(Site == "Kattegatt (ICES SD 21)", "Roskildefjorden", Site), 
         label = ifelse(Water_area=="Östersjön (hela) ICES SD 22-32", Site,
                  ifelse(Water_area=="Västerhavet (hela) ICES SD 20-21", Site, Water_area)),
         Catch_area_habitat=ifelse(Water_area %in% c("Aröd å", "Bolmen", "Enningdalsälven", "Göta älv", 
                                               "Knapebäcken", "Kvarnabäcken", "Löftaån", "Rönne å", 
                                               "Strömsån", "Viskan", "Ymsen", "Örekilsälven", "Vänern", "Ätran"), "Freshwater Skagerrak/Kattegat", # Freshwater sites connected to Kattegat 
                      ifelse(Water_area %in% c("Dalälven", "Dybäcksån", "Helge å", "Höje å", 
                                               "Kävlingeån", "Motala Ström", "Mälaren", "Nybroån", 
                                               "Nyköpingsån", "Skivarpsån", "Ången", "Hjälmaren", "Ottsjön", 
                                               "Ringsjöarna, Skåne", "Roxen", "Ångermanälven", "Åssjön"), "Freshwater Baltic Sea", # Freshwater sites connected to the Baltic Sea
                      ifelse(Water_area == "Västerhavet (hela) ICES SD 20-21", "Coastal Skagerrak/Kattegat",    # Coastal sites, Kattegat
                      ifelse(Water_area == "Östersjön (hela) ICES SD 22-32", "Coastal Baltic Sea",    # Coastal sites, Baltic Sea
                                           "No_classification")))),
          Life_history = ifelse(Life_history=="Mixed", "Habitat shifter", Life_history)) |> 
   group_by(Catch_area_habitat, label, Analythical_method, Life_history) |> 
            summarise(Individuals_n=n(),
                      lat=mean((as.numeric(as.character(lat)))), # Get the latitude per site, for Stenungsund there are several latitudes close to each other, why we use the mean value
                      long=mean((as.numeric(as.character(long)))),
                      .groups = "drop") |> 
  group_by(Catch_area_habitat, label, Analythical_method) |> 
  mutate(Total_amount_of_inds=sum(Individuals_n), # total number of individuals
         Percent_per_life_history=round((Individuals_n/Total_amount_of_inds)*100, 1), # Proportion per life history 
         lat=mean((as.numeric(as.character(lat)))), # Get the latitude per site, for Stenungsund there are several latitudes close to each other, why we use the mean value
         long=mean((as.numeric(as.character(long))))) # get the longitude per site, for Stenungsund there are several longitudes close to each other, why we use the mean value
      
## The plot:

# Fix the data so we can order the facet panels and format and translate site names:
plot_data <- Life_history_proportions_per_site |>
  mutate(facet_group = factor(paste(Catch_area_habitat, Analythical_method),
                              levels = c("Coastal Skagerrak/Kattegat FE-EPMA",
                                         "Coastal Skagerrak/Kattegat LA-ICP-MS",
                                         "Freshwater Skagerrak/Kattegat FE-EPMA",
                                         "Freshwater Skagerrak/Kattegat LA-ICP-MS",
                                         "Coastal Baltic Sea FE-EPMA",
                                         "Coastal Baltic Sea LA-ICP-MS",
                                         "Freshwater Baltic Sea FE-EPMA",
                                         "Freshwater Baltic Sea LA-ICP-MS"))) |> 
  mutate(label=ifelse(label=="Borgholms kommuns kustvatten", "Borgholm",
               ifelse(label=="Faxe Bugt-Jungshoved (Danmark)", "Faxe Bugt-Jungshoved",
               ifelse(label=="Hamina (Finland)", "Hamina",
               ifelse(label=="Karlshamn - Boön", "Karlshamn",
               ifelse(label=="Karlshamn - Veggabadet", "Karlshamn",
               ifelse(label=="Kemi-Montaja (Finland)", "Kemi-Montaja",
               ifelse(label=="Klamilanlahti (Finland)", "Klamilanlahti",
               ifelse(label=="Kõiguste (Estland)", "Kõiguste",
               ifelse(label=="Köge Bugt (Danmark)", "Köge Bugt",
               ifelse(label=="Köge Bugt (Danmark)", "Köge Bugt",
               ifelse(label=="Nynäshamns kommuns kustvatten", "Nynäshamn",
               ifelse(label=="Rymättylä (Finland)", "Rymättylä",
               ifelse(label=="Stege (Danmark)", "Stege",
               ifelse(label=="Sölvesborgs kommuns kustvatten", "Sölvesborg",
               ifelse(label=="Vilsandi (Estland)", "Vilsandi",
               ifelse(label=="Virolahti (Finland)", "Virolahti",
               ifelse(label=="Weichsel lagun (Polen)", "Vistula Lagoon",
               ifelse(label=="Ringsjöarna, Skåne", "Ringsjöarna",
               ifelse(label=="Jyllinge (Danmark)", "Jyllinge", label)))))))))))))))))))) |> 
  mutate(label = fct_drop(label))

sort(unique(plot_data$label))

## Make the plot, with the sites sorted from north to south within each facet:

# Define the colours:
 cols_life <- setNames(viridis(5, option = "plasma"),
                                c("Freshwater", "Freshwater/Estuary", "Estuary", "Habitat shifter", "Coastal"))

# Make the plot:
Figure_5 <- ggplot(plot_data, aes(x = reorder_within(label, desc(lat), facet_group), 
                              y = Percent_per_life_history, 
                              fill = Life_history)) +
  geom_col(position = "stack", colour="black") + 
  geom_text(data = plot_data |> distinct(facet_group, label, lat, Total_amount_of_inds),
            aes(x = reorder_within(label, desc(lat), facet_group), y = 109,
                                   label = Total_amount_of_inds),
                                   inherit.aes = FALSE, size=4) + # To get the sampling sites sorthed north to south within each panel
  facet_wrap(~facet_group, scales = "free_x", ncol = 2) +
  scale_x_reordered(drop=TRUE) +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 20)) +
  labs(y = "Assigned life-history (%)", x = "Site", fill = "Assigned life history") +
  theme_classic(base_size=20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        legend.key.size = unit(1.2, "cm"))+
  scale_fill_manual(values = cols_life,
                    breaks = c("Coastal", "Habitat shifter", "Estuary", "Freshwater/Estuary", "Freshwater"))

Figure_5

# Save it:
ggsave("Output/Figure_5.png",
  plot = Figure_5,
  width = 14,
  height = 16,
  units = "in",
  dpi = 300)

  
### Figure 6 #####

# Number of habitat shifts between frehwater and any other habitat (estuary and coastal)

# Load habitat shift data from script 6_Habitat_shifting.R:
Habitat_shifts_FE_EPMA_Mixed=readRDS("R_Scripts_and_Generated_data_files/Habitat_shifts_FE_EPMA_Mixed.Rds")
Habitat_shifts_LA_ICP_Mixed=readRDS("R_Scripts_and_Generated_data_files/Habitat_shifts_LA_ICP_Mixed.Rds")

Habitat_shifts_FE_EPMA_Mixed<-Habitat_shifts_FE_EPMA_Mixed |> 
  mutate(Analytical_method="FE-EPMA")

Habitat_shifts_LA_ICP_Mixed<-Habitat_shifts_LA_ICP_Mixed|> 
  mutate(Analytical_method="LA-ICP-MS")

Habitat_shifts_FE_EPMA_Mixed_combined<-rbind(Habitat_shifts_FE_EPMA_Mixed, 
                                             Habitat_shifts_LA_ICP_Mixed)

# Make a histrogram:
Habitat_shifts_FE_EPMA_Mixed_combined |>  
  ggplot(aes(x=n_shifts_total)) +
   geom_histogram(bins=59)+
    facet_wrap(~Analytical_method)


## Make one histogram per analytical method to make the bins=1:
p1<-Habitat_shifts_FE_EPMA_Mixed |> 
    ggplot(aes(x=n_shifts_total)) +
   geom_histogram(bins=14, fill="grey85", colour="black")+
    theme_classic(base_size=22)+
    xlab("Number of shifts between freshwater and coastal habitats")+
    ylab("Number of eel")+
    scale_y_continuous(breaks = seq(0, 250, by = 25))+
    scale_x_continuous(breaks = seq(0, 14, by = 1))+
    annotate("text", x = -Inf, y = Inf, label = " A", 
           hjust = -0.1, vjust = 1.1, size = 8, fontface="bold")

p2<-Habitat_shifts_LA_ICP_Mixed |> 
    filter(n_shifts_total<=100) |> 
    ggplot(aes(x=n_shifts_total)) +
   geom_histogram(bins=59, fill="grey85", colour="black")+
    theme_classic(base_size=22)+
    xlab("Number of shifts between freshwater/estuary and coastal habitats")+
    ylab("Number of eel")+
    scale_y_continuous(breaks = seq(0, 100, by = 10), limits=c(0,100))+
    scale_x_continuous(breaks = seq(0, 60, by = 5))+
    annotate("text", x = -Inf, y = Inf, label = " B", 
           hjust = -0.1, vjust = 1.1, size = 8, fontface="bold")
  
# Combine the two plots into one panel:
Figure_6 <- p1 / p2

Figure_6

# Save the plot:
ggsave("Output/Figure_6.tiff",
  plot = Figure_6,
  width = 12,
  height = 10,
  units = "in",
  dpi = 300,
  compression = "lzw")

# ------- End of Script 6_Figures.R --------- #



