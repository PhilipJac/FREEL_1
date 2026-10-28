#################
# Jacobson et al.
# Navigating the salinity gradient: Individual variation in habitat use and migration of European eel revealed by otolith microchemistry 
#################

## Last update: 2026-10-21

#### R-Script: 6_Statistical_analyses_life_hist_variation.R 

# Load needed packages:
library(tidyverse)
library(brms) # Might need to install Rtools43 to make it work https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html
library(purrr)

## Load input data:
FE_EPMA_Sotebasen_Life_history = readRDS("R_Scripts_and_Generated_data_files/FE_EPMA_Sotebasen_Life_history.Rds")
LA_ICP_MS_KUL_Sotebasen_Life_history = readRDS("R_Scripts_and_Generated_data_files/LA_ICP_MS_KUL_Sotebasen_Life_history.Rds")

# Prepaire the data in the same way as for figure 5:

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


## Set up the data for model-fitting, same way as figure 5:
model_data <- Selected_life_history_complete_Ind_data |>
  rename(Site = Site.x,
         Water_area = Water_area.x) |>
  mutate(Life_history = if_else(Life_history == "Mixed", "Habitat shifter", Life_history),
         Site = if_else(Site == "Kattegatt (ICES SD 21)", "Roskildefjorden", Site),
         label = if_else(Water_area %in% c("Östersjön (hela) ICES SD 22-32", "Västerhavet (hela) ICES SD 20-21"), Site, Water_area),
         label = if_else(label %in% c("Karlshamn - Boön", "Karlshamn - Veggabadet"), "Karlshamn", label),
         Catch_area_habitat = case_when(Water_area %in% c("Aröd å", "Bolmen", "Enningdalsälven", "Göta älv",
                                                          "Knapebäcken", "Kvarnabäcken", "Löftaån", "Rönne å",
                                                          "Strömsån", "Viskan", "Ymsen", "Örekilsälven",
                                                          "Vänern", "Ätran") ~ "Freshwater Skagerrak/Kattegat",
                                        Water_area %in% c("Dalälven", "Dybäcksån", "Helge å", "Höje å",
                                                          "Kävlingeån", "Motala Ström", "Mälaren", "Nybroån",
                                                          "Nyköpingsån", "Skivarpsån", "Ången", "Hjälmaren",
                                                          "Ottsjön", "Ringsjöarna, Skåne", "Roxen",
                                                          "Ångermanälven", "Åssjön") ~ "Freshwater Baltic Sea",

      Water_area == "Västerhavet (hela) ICES SD 20-21" ~
        "Coastal Skagerrak/Kattegat",

      Water_area == "Östersjön (hela) ICES SD 22-32" ~
        "Coastal Baltic Sea",

      TRUE ~ NA_character_
    ),
    Region = factor(Catch_area_habitat),
    Site_ID = interaction(Catch_area_habitat, label, drop = TRUE),
    Life_history = factor(Life_history),
    Analythical_method = factor(Analythical_method)
  ) |>
  filter(
    !is.na(Region),
    !is.na(Site_ID),
    !is.na(Life_history),
    !is.na(Analythical_method)
  ) |>
  droplevels()

head(model_data)

# Check pipe function:
names(Selected_life_history_complete_Ind_data)

Selected_life_history_complete_Ind_data |>
  distinct(Water_area.x) |>
  mutate(
    Catch_area_habitat = case_when(
      Water_area.x %in% c(
        "Aröd å", "Bolmen", "Enningdalsälven", "Göta älv",
        "Knapebäcken", "Kvarnabäcken", "Löftaån", "Rönne å",
        "Strömsån", "Viskan", "Ymsen", "Örekilsälven",
        "Vänern", "Ätran"
      ) ~ "Freshwater Skagerrak/Kattegat",
      Water_area.x %in% c(
        "Dalälven", "Dybäcksån", "Helge å", "Höje å",
        "Kävlingeån", "Motala Ström", "Mälaren", "Nybroån",
        "Nyköpingsån", "Skivarpsån", "Ången", "Hjälmaren",
        "Ottsjön", "Ringsjöarna, Skåne", "Roxen",
        "Ångermanälven", "Åssjön"
      ) ~ "Freshwater Baltic Sea",
      Water_area.x == "Västerhavet (hela) ICES SD 20-21" ~
        "Coastal Skagerrak/Kattegat",
      Water_area.x == "Östersjön (hela) ICES SD 22-32" ~
        "Coastal Baltic Sea",
      TRUE ~ NA_character_
    )
  ) |>
  filter(is.na(Catch_area_habitat)) # Should be 0

head(model_data) # looks good!

model_data |>
  distinct(Analythical_method, Region, Site_ID) |>
  count(Analythical_method, Region) # good!

model_data |>
  count(Analythical_method, Region, Life_history) |>
  print() # looks good!

## Fit the models, one for each analythical model:

data_FE_EPMA <- model_data |>
  filter(Analythical_method == "FE-EPMA") |>
  droplevels()

head(data_FE_EPMA)

data_LA_ICP_MS<- model_data |>
  filter(Analythical_method == "LA-ICP-MS") |>
  droplevels()

# Test these to see if the Rtools43 is correctly installed and works:

#Sys.which("make")
#Sys.which("g++")
#Sys.getenv("PATH")
#file.exists("C:/rtools43/usr/bin/make.exe")
#file.exists("C:/rtools43/x86_64-w64-mingw32.static.posix/bin/g++.exe")
#file.exists("C:/rtools43/x86_64-w64-mingw32.static.posix/bin/gcc.exe")
#system2("make", "--version")

# FE-EPMA model:
#fit_epma <- brm(Life_history ~ Region + (1 | Site_ID),
 #               data = data_FE_EPMA,
  #              family = categorical(),
   #             chains = 2,
    #            iter = 1000,
     #           cores = 1,
      #          seed = 123,
       #         refresh=100)

#summary(fit_epma)
#plot(fit_epma)

#nuts_params(fit_epma) |>
#  dplyr::filter(Parameter == "divergent__") |>
#  dplyr::summarise(
#    divergences = sum(Value))

#data_FE_EPMA |>
#  count(Region, Life_history) |>
#  tidyr::pivot_wider(
#    names_from = Life_history,
#    values_from = n,
#    values_fill = 0)

## Try with priors: 

get_prior(Life_history ~ Region + (1 | Site_ID),
  data = data_FE_EPMA,
  family = categorical())

priors_epma <- c(prior(normal(0, 1.5), class = "b", dpar = "muFreshwater"),
                 prior(normal(0, 1.5), class = "b", dpar = "muHabitatshifter"))

validate_prior(priors_epma,
               Life_history ~ Region + (1 | Site_ID),
               data = data_FE_EPMA,
               family = categorical())

fit_epma_regularised_099 <- brm(Life_history ~ Region + (1 | Site_ID),
                            data = data_FE_EPMA,
                            family = categorical(),
                            prior = priors_epma,
                            chains = 2,
                            iter = 2000,
                            warmup = 1000,
                            cores = 1,
                            seed = 123,
                            refresh = 100,
                            control = list(adapt_delta = 0.99))

summary(fit_epma_regularised_099)

# Save model output to avoid the need of running it again:
saveRDS(fit_epma_regularised_099, "fit_epma_regularised_099.rds")

fit_epma_regularised_099<-readRDS("fit_epma_regularised_099.rds")

summary(fit_epma_regularised_099)

## Looks quite good, check model diagnostics:
# Count divergent transitions
brms::nuts_params(fit_epma_regularised_099) |>
  dplyr::filter(Parameter == "divergent__") |>
  dplyr::summarise(
    divergences = sum(Value)
  ) # looks good

# Check maximum tree-depth warnings
brms::nuts_params(fit_epma_regularised_099) |>
  dplyr::filter(Parameter == "treedepth__") |>
  dplyr::summarise(
    max_treedepth_observed = max(Value)
  ) # looks good

## Check predicted life-history probabilities:

epma_region_predictions <- fitted(
  fit_epma_regularised_099,
  newdata = data.frame(
    Region = factor(
      levels(data_FE_EPMA$Region),
      levels = levels(data_FE_EPMA$Region)
    ),
    Site_ID = NA
  ),
  re_formula = NA,
  summary = TRUE
)

epma_region_table <- tibble(
  Region = levels(data_FE_EPMA$Region),
  Coastal = epma_region_predictions[, , "P(Y = Coastal)"][, "Estimate"],
  Coastal_low = epma_region_predictions[, , "P(Y = Coastal)"][, "Q2.5"],
  Coastal_high = epma_region_predictions[, , "P(Y = Coastal)"][, "Q97.5"],
  Freshwater = epma_region_predictions[, , "P(Y = Freshwater)"][, "Estimate"],
  Freshwater_low = epma_region_predictions[, , "P(Y = Freshwater)"][, "Q2.5"],
  Freshwater_high = epma_region_predictions[, , "P(Y = Freshwater)"][, "Q97.5"],
  Habitat_shifter = epma_region_predictions[, , "P(Y = Habitat shifter)"][, "Estimate"],
  Habitat_shifter_low = epma_region_predictions[, , "P(Y = Habitat shifter)"][, "Q2.5"],
  Habitat_shifter_high = epma_region_predictions[, , "P(Y = Habitat shifter)"][, "Q97.5"]
)

epma_region_table

### Next step:

# Create one prediction row for each geographical region
new_regions <- data.frame(Region = factor(levels(data_FE_EPMA$Region),
                          levels = levels(data_FE_EPMA$Region)))

# Generate posterior predicted probabilities
epma_epred <- posterior_epred(fit_epma_regularised_099,
                              newdata = new_regions,
                              re_formula = NA)
# Check the dimensions
dim(epma_epred)
# Next step:

# Extract region and life-history names
region_names <- as.character(new_regions$Region)

life_history_names <- dimnames(epma_epred)[[3]]

# Fallback in case category names are not attached to the array
if (is.null(life_history_names)) {
  life_history_names <- levels(data_FE_EPMA$Life_history)
}

# Generate all unique pairs of regions
region_pairs <- combn(
  seq_along(region_names),
  2,
  simplify = FALSE
)

# Calculate posterior contrasts
epma_contrasts <- do.call(
  rbind,
  lapply(seq_along(life_history_names), function(k) {

    do.call(rbind, lapply(region_pairs, function(pair) {

      # Difference in predicted proportions:
      # Region 1 minus Region 2
      difference <- epma_epred[, pair[1], k] -
                    epma_epred[, pair[2], k]

      data.frame(
        Life_history = life_history_names[k],
        Region_1 = region_names[pair[1]],
        Region_2 = region_names[pair[2]],

        # Differences expressed in percentage points
        Difference_pp = 100 * median(difference),

        CI_low = 100 * unname(quantile(difference, 0.025)),
        CI_high = 100 * unname(quantile(difference, 0.975)),

        # Posterior probabilities for the direction of the difference
        Prob_region1_higher = mean(difference > 0),
        Prob_region2_higher = mean(difference < 0)
      )
    }))
  })
)

# Round and inspect the results
epma_contrasts <- epma_contrasts |>
  mutate(
    across(
      c(Difference_pp, CI_low, CI_high),
      ~ round(.x, 1)
    ),
    across(
      c(Prob_region1_higher, Prob_region2_higher),
      ~ round(.x, 3)
    )
  )

print(epma_contrasts) # This is per region

## Variation among sites within regions:

# One row per sampling site
site_data <- data_FE_EPMA |>
  distinct(Region, Site_ID) |>
  mutate(
    Region = factor(
      Region,
      levels = levels(data_FE_EPMA$Region)
    ),
    Site_ID = factor(
      Site_ID,
      levels = levels(data_FE_EPMA$Site_ID)
    )
  )

# Check how many sites are represented in each region
site_data |>
  count(Region, name = "Number_of_sites")

# Posterior predicted probabilities for each site
# Include the estimated site-level random effects
epma_site_epred <- posterior_epred(
  fit_epma_regularised_099,
  newdata = site_data,
  re_formula = NULL
)

# Check dimensions: posterior draws x sites x life-history categories
dim(epma_site_epred)

# Extract category names from the prediction array
category_names <- dimnames(epma_site_epred)[[3]]

# Calculate the among-site SD for every posterior draw
site_variation <- list()
result_index <- 1

for (region_name in levels(site_data$Region)) {

  site_indices <- which(
    as.character(site_data$Region) == region_name
  )

  for (k in seq_along(category_names)) {

    # Draws x sites matrix for this category and region
    site_probabilities <- epma_site_epred[
      , site_indices, k, drop = FALSE
    ]

    dim(site_probabilities) <- c(
      dim(epma_site_epred)[1],
      length(site_indices)
    )

    # SD among sites for each posterior draw
    sd_draws <- apply(
      site_probabilities,
      1,
      sd
    )

    site_variation[[result_index]] <- tibble(
      Region = region_name,
      Life_history = category_names[k],
      Median_SD_pp = 100 * median(sd_draws),
      CI_low = 100 * unname(
        quantile(sd_draws, 0.025)
      ),
      CI_high = 100 * unname(
        quantile(sd_draws, 0.975)
      )
    )

    result_index <- result_index + 1
  }
}

site_variation_results <- bind_rows(site_variation) |>
  mutate(
    across(
      c(Median_SD_pp, CI_low, CI_high),
      ~ round(.x, 1)
    )
  )

print(site_variation_results)

saveRDS(site_variation_results, "site_variation_results_FE_EPMA.rds")

### LA-ICP-MS data #######

# This is without priors:

# fit_laicpms <- brm(Life_history ~ Region + (1 | Site_ID),
#                   data = data_LA_ICP_MS,
#                   family = categorical(),
#                   chains = 2,
#                   iter = 2000,
#                   warmup =1000,
#                   cores = 1,
#                   seed = 123,
#                   refresh=100,
#                   control = list(adapt_delta = 0.99))


#summary(fit_laicpms)

#brms::nuts_params(fit_laicpms) |>
#  dplyr::filter(Parameter == "divergent__") |>
#  dplyr::summarise(
#    divergences = sum(Value)
#  )

#brms::nuts_params(fit_laicpms) |>
#  dplyr::filter(Parameter == "treedepth__") |>
#  dplyr::summarise(
#    max_treedepth_observed = max(Value),
#    transitions_at_10 = sum(Value >= 10)
#  )
                     
#data_LA_ICP_MS |>
#  dplyr::count(Region, Life_history) |>
#  tidyr::pivot_wider(
#    names_from = Life_history,
#    values_from = n,
#    values_fill = 0
#  )

#data_LA_ICP_MS |>
#  dplyr::distinct(Region, Site_ID) |>
#  dplyr::count(Region)

## Try with priors:
priors_laicpms <- c(
  # Regional effects
  prior(normal(0, 1.5),
        class = "b", dpar = "muFreshwater"),

  prior(normal(0, 1.5),
        class = "b", dpar = "muFreshwaterEstuary"),

  prior(normal(0, 1.5),
        class = "b", dpar = "muHabitatshifter"),

  # Intercepts
  prior(normal(0, 2.5),
        class = "Intercept", dpar = "muFreshwater"),

  prior(normal(0, 2.5),
        class = "Intercept", dpar = "muFreshwaterEstuary"),

  prior(normal(0, 2.5),
        class = "Intercept", dpar = "muHabitatshifter"),

  # Site-level variation
  prior(student_t(3, 0, 1.5),
        class = "sd", dpar = "muFreshwater"),

  prior(student_t(3, 0, 1.5),
        class = "sd", dpar = "muFreshwaterEstuary"),

  prior(student_t(3, 0, 1.5),
        class = "sd", dpar = "muHabitatshifter")
)


# Check them: 
get_prior(
  Life_history ~ Region + (1 | Site_ID),
  data = data_LA_ICP_MS,
  family = categorical()
)

 #Fit he model using priors and regularization:

fit_laicpms_regularised <- brm(
  Life_history ~ Region + (1 | Site_ID),
  data = data_LA_ICP_MS,
  family = categorical(),
  prior = priors_laicpms,
  chains = 2,
  iter = 2000,
  warmup = 1000,
  cores = 1,
  seed = 123,
  refresh = 100,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 12
  )
)

summary(fit_laicpms_regularised)

# Save model output to avoid the need of running it again:
saveRDS(fit_laicpms_regularised, "fit_laicpms_regularised.rds")

fit_laicpms_regularised<-readRDS("fit_laicpms_regularised.rds")

# Check divergent transitions
brms::nuts_params(fit_laicpms_regularised) |>
  dplyr::filter(Parameter == "divergent__") |>
  dplyr::summarise(divergences = sum(Value))

# Check maximum treedepth
brms::nuts_params(fit_laicpms_regularised) |>
  dplyr::filter(Parameter == "treedepth__") |>
  dplyr::summarise(
    max_treedepth_observed = max(Value),
    transitions_at_12 = sum(Value >= 12)
  )

summary(fit_laicpms_regularised)

## LA-ICP-MS: regional predicted life-history probabilities

# Get region and life-history levels
laicpms_regions <- levels(data_LA_ICP_MS$Region)
laicpms_life_histories <- levels(data_LA_ICP_MS$Life_history)

# One prediction row per region
laicpms_new_regions <- data.frame(
  Region = factor(
    laicpms_regions,
    levels = levels(data_LA_ICP_MS$Region)
  )
)

# Posterior predicted probabilities, excluding site-specific effects
laicpms_region_predictions <- posterior_epred(
  fit_laicpms_regularised,
  newdata = laicpms_new_regions,
  re_formula = NA
)

# Create a table of posterior medians and 95% credible intervals
# Use posterior draws directly for consistent summaries

laicpms_region_probabilities <- purrr::map_dfr(
  seq_along(laicpms_regions),
  function(r) {
    purrr::map_dfr(
      seq_along(laicpms_life_histories),
      function(h) {

        p <- laicpms_region_predictions[, r, h]

        tibble(
          Region = laicpms_regions[r],
          Life_history = laicpms_life_histories[h],
          Estimate = median(p),
          CI_low = unname(quantile(p, 0.025)),
          CI_high = unname(quantile(p, 0.975))
        )
      }
    )
  }
) |>
  mutate(
    across(
      c(Estimate, CI_low, CI_high),
      ~ round(100 * .x, 1)
    )
  )

print(laicpms_region_probabilities, n = Inf)

## Extract posterior estimates:

regions <- levels(data_LA_ICP_MS$Region)
life_histories <- levels(data_LA_ICP_MS$Life_history)

newdata_regions <- data.frame(
  Region = factor(
    regions,
    levels = regions
  )
)

# Posterior predictions, excluding site-specific random effects
epred_laicpms <- posterior_epred(
  fit_laicpms_regularised,
  newdata = newdata_regions,
  re_formula = NA
)

# Check the dimensions of the prediction array
dim(epred_laicpms)

## LA-ICP-MS: pairwise regional contrasts

laicpms_region_pairs <- combn(
  seq_along(laicpms_regions),
  2,
  simplify = FALSE
)

laicpms_region_contrasts <- purrr::map_dfr(
  seq_along(laicpms_life_histories),
  function(h) {
    purrr::map_dfr(
      laicpms_region_pairs,
      function(pair) {

        r1 <- pair[1]
        r2 <- pair[2]

        # Posterior differences in predicted probabilities
        # Region 1 minus Region 2, in percentage points
        differences <- 100 * (
          laicpms_region_predictions[, r1, h] -
          laicpms_region_predictions[, r2, h]
        )

        tibble(
          Life_history = laicpms_life_histories[h],
          Region_1 = laicpms_regions[r1],
          Region_2 = laicpms_regions[r2],
          Difference_pp = median(differences),
          CI_low = unname(
            quantile(differences, 0.025)
          ),
          CI_high = unname(
            quantile(differences, 0.975)
          ),
          Prob_region1_higher = mean(differences > 0),
          Prob_region2_higher = mean(differences < 0)
        )
      }
    )
  }
) |>
  mutate(
    across(
      c(Difference_pp, CI_low, CI_high),
      ~ round(.x, 1)
    ),
    across(
      c(Prob_region1_higher, Prob_region2_higher),
      ~ round(.x, 3)
    )
  )

print(laicpms_region_contrasts, n = Inf)

## LA-ICP-MS: variation among sites within regions

# One row per sampling site
laicpms_site_data <- data_LA_ICP_MS |>
  distinct(Region, Site_ID) |>
  mutate(
    Region = factor(
      Region,
      levels = levels(data_LA_ICP_MS$Region)
    ),
    Site_ID = factor(
      Site_ID,
      levels = levels(data_LA_ICP_MS$Site_ID)
    )
  )

# Number of sites per region
laicpms_site_data |>
  count(Region, name = "Number_of_sites")

# Posterior predicted probabilities for each site,
# including estimated site-level random effects
laicpms_site_predictions <- posterior_epred(
  fit_laicpms_regularised,
  newdata = laicpms_site_data,
  re_formula = NULL
)

# Dimensions: draws x sites x life-history categories
dim(laicpms_site_predictions)

# Extract category names
laicpms_site_categories <- dimnames(
  laicpms_site_predictions
)[[3]]

if (is.null(laicpms_site_categories)) {
  laicpms_site_categories <- levels(data_LA_ICP_MS$Life_history)
}

## Calculate among-site variation for LA-ICP-MS
## Handle regions with fewer than two sites

laicpms_site_variation_list <- list()
result_index <- 1

for (region_name in levels(laicpms_site_data$Region)) {

  site_indices <- which(
    as.character(laicpms_site_data$Region) == region_name
  )

  n_sites <- length(site_indices)

  for (k in seq_along(laicpms_site_categories)) {

    # If fewer than two sites, SD cannot be estimated
    if (n_sites < 2) {

      laicpms_site_variation_list[[result_index]] <- tibble(
        Region = region_name,
        Life_history = laicpms_site_categories[k],
        Number_of_sites = n_sites,
        Median_SD_pp = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_
      )

    } else {

      # Posterior draws x sites for this life-history category
      site_probabilities <- laicpms_site_predictions[
        , site_indices, k, drop = FALSE
      ]

      dim(site_probabilities) <- c(
        dim(laicpms_site_predictions)[1],
        n_sites
      )

      # SD among sites for each posterior draw
      sd_draws <- apply(
        site_probabilities,
        1,
        sd
      )

      # Remove any non-finite values before summarising
      sd_draws <- sd_draws[is.finite(sd_draws)]

      laicpms_site_variation_list[[result_index]] <- tibble(
        Region = region_name,
        Life_history = laicpms_site_categories[k],
        Number_of_sites = n_sites,
        Median_SD_pp = 100 * median(sd_draws),
        CI_low = 100 * unname(
          quantile(sd_draws, 0.025)
        ),
        CI_high = 100 * unname(
          quantile(sd_draws, 0.975)
        )
      )
    }

    result_index <- result_index + 1
  }
}

# Combine and format results
laicpms_site_variation <- bind_rows(
  laicpms_site_variation_list
) |>
  mutate(
    across(
      c(Median_SD_pp, CI_low, CI_high),
      ~ round(.x, 1)
    )
  )

print(laicpms_site_variation, n = Inf)

saveRDS(laicpms_site_variation, "site_variation_results_LA_ICP_MS.rds")

## Main results of both models:

# Model summaries
summary(fit_epma_regularised_099)
summary(fit_laicpms_regularised)

# FE-EPMA results
print(epma_region_table, n = Inf)
print(epma_contrasts)
print(site_variation_results, n = Inf)

# LA-ICP-MS results
print(laicpms_region_probabilities, n = Inf)
print(laicpms_region_contrasts, n = Inf)
print(laicpms_site_variation, n = Inf)

## Save results table for plotting script:

# Save model-predicted regional probabilities for plotting
saveRDS(epma_region_table,"R_Scripts_and_Generated_data_files/epma_region_table.rds")
saveRDS(laicpms_region_probabilities,"R_Scripts_and_Generated_data_files/laicpms_region_probabilities.rds")
