#################
# Jacobson et al.
# Navigating the salinity gradient: Individual variation in habitat use and migration of European eel revealed by otolith microchemistry 
#################

# Last update: 2026-10-21

#### R-Script: 8_Life_history_ind_level_comparison.R

# Load packages
library(tidyverse)
library(emmeans)

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

head(FE_EPMA_Sotebasen_Distance_Habitat_Assigned_Life_history)

# Keep the first row of each individual and then compare length-weight-age and gender data:
FE_EPMA_Ind_Assinged_life_history <- FE_EPMA_Sotebasen_Distance_Habitat_Assigned_Life_history %>%
  group_by(IndividNr) %>%
  slice(1) %>%
  ungroup() 

n_distinct(FE_EPMA_Ind_Assinged_life_history$IndividNr) # Should be 1877 which it is

## Number of males:
FE_EPMA_Ind_Assinged_life_history %>%
  group_by(Sex) %>%
   summarise(n=n_distinct(IndividNr))

FE_EPMA_Ind_Assinged_life_history %>%
  filter(Sex=="m") %>%
  group_by(Water_area, Site) %>%
  summarise(n=n_distinct(IndividNr)) # 57 males from coastal sampling sites.

# Have a look at the data:

FE_EPMA_Ind_Assinged_life_history |>
     ggplot(aes(y=Total_weight_g, x=Total_length_mm, colour=as.factor(Life_history), group=as.factor(Life_history)))+
       geom_point()

FE_EPMA_Ind_Assinged_life_history |>
     ggplot(aes(y=log(Total_weight_g), x=log(Total_length_mm), colour=as.factor(Life_history), group=as.factor(Life_history)))+
       geom_point(alpha=0.5)+
       stat_smooth(method="lm")

# Save the data for the figure script:
saveRDS(FE_EPMA_Ind_Assinged_life_history, file="R_Scripts_and_Generated_data_files/FE_EPMA_Ind_Assinged_life_history_weight_length.Rds")

## Stat test if the length-weight relationships differs between the life-histories:

Condition_test_FE_EPMA_interaction <- lm(log(Total_weight_g) ~ log(Total_length_mm) * factor(Life_history), data = FE_EPMA_Ind_Assinged_life_history)
Condition_test_FE_EPMA_no_interaction  <- lm(log(Total_weight_g) ~ log(Total_length_mm) + factor(Life_history), data = FE_EPMA_Ind_Assinged_life_history)
Condition_test_FE_EPMA_no_life_history  <- lm(log(Total_weight_g) ~ log(Total_length_mm), data = FE_EPMA_Ind_Assinged_life_history)

anova(Condition_test_FE_EPMA_no_life_history, Condition_test_FE_EPMA_no_interaction, Condition_test_FE_EPMA_interaction)
AIC(Condition_test_FE_EPMA_no_life_history, Condition_test_FE_EPMA_no_interaction, Condition_test_FE_EPMA_interaction)
summary(Condition_test_FE_EPMA_interaction) # Selected model according to AIC

par(mfrow=c(2,2))
plot(Condition_test_FE_EPMA_interaction)

#### Extract model estimates:

df_FE <- FE_EPMA_Ind_Assinged_life_history %>%
  mutate(log_weight = log(Total_weight_g),
         log_length = log(Total_length_mm),
         Life_history = factor(Life_history))

model_FE <- lm(log_weight ~ log_length * Life_history, data = df_FE)

slopes_FE <- emtrends(model_FE,
                      pairwise ~ Life_history,
                      var = "log_length",
                      adjust = "tukey")

slopes_FE$emtrends # Estimated slopes and CI
slopes_FE$contrasts # SLope comparisons
## Predict the weight of a 800mm coastal resident and freshwater resident eel:
# Predict weight for 800 mm eels in each life-history group

pred_800 <- emmeans(model_FE,
                    ~ Life_history,
                    at = list(log_length = log(800)))

# Predictions are on the log-weight scale
pred_800_df <- as.data.frame(pred_800)

# Back-transform to grams
pred_800_df$predicted_weight_g <- exp(pred_800_df$emmean)

pred_800_df[, c("Life_history", "predicted_weight_g")] # There we go:
# a 800 mm coastal eel is predicted to weigh 1098.631, a freshwater resident eel is predcited to weigh 1014.604 grams

### LA-ICP-MS #####

# First, select the individuals assigned with a mixed life-history:
LA_ICP_MS_KUL_Sotebasen_Life_history_selected<-LA_ICP_MS_KUL_Sotebasen_Life_history |> 
  select(IndividNr, Life_history)

# Combine with transect data for each individual:
LA_ICP_MS_KUL_Sotebasen_Distance_Habitat_Assigned_Life_history<-left_join(LA_ICP_MS_KUL_Sotebasen_Distance_Habitat, LA_ICP_MS_KUL_Sotebasen_Life_history_selected, by="IndividNr")

# Extract one row per individual:
LA_ICP_MS_Ind_Assinged_life_history <- LA_ICP_MS_KUL_Sotebasen_Distance_Habitat_Assigned_Life_history %>%
  group_by(IndividNr) %>%
  slice(1) %>%
  ungroup() # 1045, which it should be.

## Number of males:
LA_ICP_MS_Ind_Assinged_life_history  %>%
  filter(Sex %in% c("Hane", "m")) %>%
 group_by(Water_area, Site) %>%
  summarise(n=n_distinct(IndividNr)) %>% 
            ungroup()%>%
  mutate(Tot_n=sum(n)) # 8 males in total, 1 is from a coastal sampling site, the rest from freshwater lakes (Mälaren and Bolmen).

# Plot the data:

LA_ICP_MS_Ind_Assinged_life_history  |>
     ggplot(aes(y=Total_weight_g, x=Total_length_mm, colour=as.factor(Life_history), group=as.factor(Life_history)))+
       geom_point()

LA_ICP_MS_Ind_Assinged_life_history  |>
     ggplot(aes(y=log(Total_weight_g), x=log(Total_length_mm), colour=as.factor(Life_history), group=as.factor(Life_history)))+
       geom_point()+
       stat_smooth(method="lm")

# Save the data for the figure script:
saveRDS(LA_ICP_MS_Ind_Assinged_life_history , file="R_Scripts_and_Generated_data_files/LA_ICP_MS_Ind_Assinged_life_history_weight_length.Rds")

## Stat-testing if the length-weight relationship differs between Life-history categories:

Condition_test_LA_ICP_MS_interaction <- lm(log(Total_weight_g) ~ log(Total_length_mm) * factor(Life_history), data = LA_ICP_MS_Ind_Assinged_life_history)
Condition_test_LA_ICP_MS_no_interaction  <- lm(log(Total_weight_g) ~ log(Total_length_mm) + factor(Life_history), data = LA_ICP_MS_Ind_Assinged_life_history)
Condition_test_LA_ICP_MS_no_life_history  <- lm(log(Total_weight_g) ~ log(Total_length_mm), data = LA_ICP_MS_Ind_Assinged_life_history)

anova(Condition_test_LA_ICP_MS_no_life_history, Condition_test_LA_ICP_MS_no_interaction, Condition_test_LA_ICP_MS_interaction)
AIC(Condition_test_LA_ICP_MS_no_life_history, Condition_test_LA_ICP_MS_no_interaction, Condition_test_LA_ICP_MS_interaction)
summary(Condition_test_LA_ICP_MS_interaction) # Selected model according to AIC

par(mfrow=c(2,2))
plot(Condition_test_LA_ICP_MS_interaction)

## Em-means to extract model parameters:

df_LA_ICP_MS <- LA_ICP_MS_Ind_Assinged_life_history %>%
  mutate(log_weight = log(Total_weight_g),
         log_length = log(Total_length_mm),
         Life_history = factor(Life_history))

model_LA_ICP_MS <- lm(log_weight ~ log_length * Life_history, data = df_LA_ICP_MS)

slopes_LA_ICP_MS <- emtrends(model_LA_ICP_MS,
                             pairwise ~ Life_history,
                             var = "log_length",
                             adjust = "tukey")

slopes_LA_ICP_MS$emtrends
slopes_LA_ICP_MS$contrasts

pred_800_LA <- emmeans(model_LA_ICP_MS,
                    ~ Life_history,
                    at = list(log_length = log(800)))

# Predictions are on the log-weight scale
pred_800_LA_df <- as.data.frame(pred_800_LA)

# Back-transform to grams
pred_800_LA_df$predicted_weight_g <- exp(pred_800_LA_df$emmean)

pred_800_LA_df[, c("Life_history", "predicted_weight_g")] # There we go:
# a 800 mm coastal eel is predicted to weigh 960.1588 grams, a freshwater resident eel is predcited to weigh 895.9660 grams, 
# Freshwater/Estuary to weigh 930.1634 grams and Habitat shifting eel to weigh 974.9967
