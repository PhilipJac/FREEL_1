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

# Plot males on a map:

## Sötebasen:
Sotebasen_inds <- read.delim(file = "Input_data/Ind_level_data/ÅlindividerSötebasen.csv",
                                     sep = ";",
                                     dec = ",",
                                     fileEncoding = "Latin1")
head(Sotebasen_inds) 
sort(unique(Sotebasen_inds$Kön)) # f = female, o = unclassified, iu = To be checked, m = male

  Sotebasen_inds_males <- Sotebasen_inds  |> 
  filter(Kön=="m")
names(Sotebasen_inds_males)

Sotebasen_inds_males <- Sotebasen_inds_males |>
  mutate(across(c(WGS84_N_Plats, WGS84_N_Vatten,
                  WGS84_E_Plats, WGS84_E_Vatten), ~ readr::parse_number(as.character(.x))))

Sotebasen_inds_males_map <-Sotebasen_inds_males |> 
   mutate(lat = case_when(VattenNamn %in% c("Östersjön (hela) ICES SD 22-32","Västerhavet (hela) ICES SD 20-21") ~ WGS84_N_Plats,
                         TRUE ~ WGS84_N_Vatten),
          lon = case_when(VattenNamn %in% c("Östersjön (hela) ICES SD 22-32","Västerhavet (hela) ICES SD 20-21") ~ WGS84_E_Plats,
                         TRUE ~ WGS84_E_Vatten),
          label=ifelse(VattenNamn=="Östersjön (hela) ICES SD 22-32", Plats,
                ifelse(VattenNamn=="Västerhavet (hela) ICES SD 20-21", Plats, VattenNamn)),   
          Sampling_site_category=ifelse(VattenNamn %in% c("Östersjön (hela) ICES SD 22-32", 
                                                         "Västerhavet (hela) ICES SD 20-21"), "Coastal", "Freshwater")) |>
   group_by(label) |> 
    summarise(n=n_distinct(IndividNr), # number of individuals
              lat=mean(lat, na.rm = TRUE), # Get the latitude per site, for Stenungsund there are several latitudes close to each other, why we use the mean value
              lon=mean(lon, na.rm = TRUE), # get the longitude per site, for Stenungsund there are several longitudes close to each other, why we use the mean value
              Sampling_site_category=unique(Sampling_site_category), # Get the sampling site category 
              .groups = "drop") |>
  filter(lat>45)|> # Remove NA and Italian males
  mutate(Databse="Sötebasen")

Sotebasen_inds_males_map # This is the strucutre we want

## KUL:
# To do
  




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

  geom_point(data = Sotebasen_inds_males_map,
             aes(x = lon,
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
   geom_text(aes(x = 20.25, y = 59, label = "Baltic Sea"),
              size = 5,
              fontface = "italic",
              color = "Black",
              angle = 10)+
    geom_text(aes(x = 18.3, y = 56, label = "Baltic Proper"),
              size = 4.5,
              fontface = "italic",
              color = "grey",
              angle = 40)+
  geom_text(aes(x = 19.55, y = 61.75, label = "Gulf of Bothnia"),
              size = 4.5,
              fontface = "italic",
              color = "grey",
              angle = 75)+
   geom_text(aes(x = 23, y = 64.7, label = "Bothnian Bay"),
              size = 4.5,
              fontface = "italic",
              color = "grey",
              angle = 50)+
  geom_text(aes(x = 25, y = 59.75, label = "Gulf of Finland"),
              size = 4.5,
              fontface = "italic",
              color = "grey",
              angle = 10)+
  
  ## Add countries abbreviatons:

  # Denmark
geom_text(aes(x = 9.25, y = 56.2, label = "DNK"),
              size = 5, 
              color = "Black", angle = 0) +
# Estonia
geom_text(aes(x = 25.5, y = 58.7, label = "EST"),
               size = 5, 
               color = "Black", 
               angle = 0) +
# Finland
geom_text(aes(x = 24.75, y = 62.0, label = "FIN"),
              size = 5, 
              color = "Black", 
              angle = 0) +
# Germany
geom_text(aes(x = 11.5, y = 53.5, label = "DEU"),
              size = 5,
              color = "Black", 
              angle = 0) +
# Latvia
geom_text(aes(x = 25.25, y = 56.9, label = "LVA"),
              size = 5,
              color = "Black", 
              angle = 0) +
# Lithuania
geom_text(aes(x = 24.0, y = 55.5, label = "LTU"),
               size = 5, 
               color = "Black", 
               angle = 0) +
# Poland
geom_text(aes(x = 19.0, y = 53.75, label = "POL"),
              size = 5,
              color = "Black", 
              angle = 0) +
# Russia
geom_text(aes(x = 29.5, y = 58.5, label = "RUS"),
          size = 5,
          color = "Black", 
          angle = 0) +
# Sweden
geom_text(aes(x = 14, y = 61.5, label = "SWE"),
              size = 5,
              color = "Black", 
              angle = 0) +
# Norway
geom_text(aes(x = 10.25, y = 61.0, label = "NOR"),
          size = 5,
          color = "Black", 
          angle = 0) +
  
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
ggsave("Output/Figure_1_rev.tiff",
  plot = Figure_1,
  width = 8,
  height = 12,
  units = "in",
  dpi = 300,
  compression = "lzw")
