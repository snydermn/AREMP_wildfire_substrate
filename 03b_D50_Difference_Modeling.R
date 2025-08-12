
####******** MODELING DIFFERENCE for D50 ********************######
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

#### Broadest model categories will be: Klamath and Cascades

# Load Libraries####
library(randomForest)
library(tidyverse)
library(dplyr)
library(caret)
library(ranger)
library(pdp)
library(patchwork)
library(pdp)
library(lattice)

# Load Data Sets ####

# bankfull processing

bankfull <- read.csv("all_site_survey_covars.csv") %>% 
  group_by(SITE_ID) %>% 
  arrange(SITE_ID) %>% 
  dplyr::select(c(SITE_ID, bankfull_mean)) %>% 
  summarise(bankfull_all_mean = mean(bankfull_mean))

# all post-fire surveys
pf_survey_covars <- read.csv("post_fire_surveys_covar.csv") %>% 
  filter(ppc_burned_perc >= 1) %>% 
  left_join(bankfull, by = c("SITE_ID" = "SITE_ID")) %>% 
  mutate(relative_bankfull = bankfull_mean / bankfull_all_mean) %>% 
  filter(QC_visit.x == "N")

# site id linkages that will be used for the difference models
site_id_link <- read.csv("aremp_site_info/site_id_linkage.csv") %>% 
  filter(QC_visit == "N")

# klamath pred_obs data for testing on all post-fire sites-- produced from D50_Unburned_Modeling script
pred_obs.pf.k <- read_csv("RF_models/FINAL_RF_MODELS/D50/post_fire_testing/predicted_observed/klamath_D50_phi_post_fire_pred_obs.csv")

pred_obs.pf.c <- read_csv("RF_models/FINAL_RF_MODELS/D50/post_fire_testing/predicted_observed/cascades_D50_phi_post_fire_pred_obs.csv")

####### KLAMATH DIFFERENCE MODELING: EXAMINING TEST RUNS ON POST-FIRE SURVEYS ####
#### cREATE POST-FIRE SURVEY DATA FRAME FOR first burned cats ####


pf_covars_filtered <- pf_survey_covars %>% 
  group_by(site_survey_id) %>% 
  arrange(site_survey_id, ppc_burned_sqkm) %>% 
  mutate(area_burned_rank = row_number()) %>% 
  ungroup() %>% 
  relocate(area_burned_rank, .after= site_survey_id) %>% 
  filter(area_burned_rank == "1") %>% 
  dplyr::select(c(ppc_burned_sqkm,rank.x,
                  ppc_himod_sqkm, rip_burned_sqkm, 
                  rip_himod_sqkm,
                  reach_burned_sqkm,
                  reach_himod_sqkm, YPF,
                  CREEK_CODE,BFIWS, elev_mean,
                  SITE_ID, survey_year,
                  slope23d_perc_greater,                    
                  area_sqkm.x, stream_slope, 
                  bankfull_mean,relative_bankfull,
                  Confinement,COMPSTRGTHWS, HYDRLCONDWS, dom_geology,                                    
                  KFFACTWS,PERMWS,RCKDEPWS,OMWS, PRECIP8110WS,
                  TMEAN8110WS, high_flow_count,                                    
                  mean_ad, mean_cc, frequency_km,                                    
                  RDCRSWS,RDDENSWS, AQUA_PROV, site_survey_id)) %>% 
  drop_na() %>% 
  filter(AQUA_PROV == "Klamath/Siskiyou") %>% 
  filter(SITE_ID != "CAEMI0015") %>% #remove catchment area outlier
  dplyr::select(-c(AQUA_PROV)) %>% 
  filter(Confinement != "None")

difference.fines.k <- pred_obs.pf.k %>% 
  mutate(difference = predicted - observed) %>% 
  left_join(site_id_link, by = c("site_id" = "site_id", 
                                 "year" = "survey_year")) %>% 
  select(-c(X, OBJECTID, QC_visit)) %>% 
  left_join(pf_covars_filtered, by = c("site_survey_id" = "site_survey_id"),
            relationship = "many-to-many")

difference_train.fines.k<- difference.fines.k %>% 
  dplyr::select(c(Model, difference, ppc_burned_sqkm,
                  ppc_himod_sqkm, rip_burned_sqkm, 
                  rip_himod_sqkm,rank.x,
                  reach_burned_sqkm,
                  reach_himod_sqkm, YPF,
                  CREEK_CODE,BFIWS, elev_mean,
                  SITE_ID, survey_year,
                  slope23d_perc_greater,                    
                  area_sqkm.x, stream_slope, 
                  bankfull_mean,relative_bankfull,
                  Confinement,COMPSTRGTHWS, HYDRLCONDWS, dom_geology,                                    
                  KFFACTWS,PERMWS,RCKDEPWS,OMWS, PRECIP8110WS,
                  TMEAN8110WS, high_flow_count,                                    
                  mean_ad, mean_cc, frequency_km,                                    
                  RDCRSWS,RDDENSWS)) %>% 
  drop_na() 




##### Difference Modeling: Klamath ######

# Step 1: Split the dataframe by "Model" and store each subset in a list
model_data_list <- split(difference_train.fines.k, difference_train.fines.k$Model)

# Step 2: Initialize lists to store the outputs
accuracy_list.kd <- list()
error_list.kd <- list()
model_results.kd <- list()  # List to store dataframes with observed and predicted together
varimp_list.kd <- list()

# Step 3: Loop through the list of dataframes to fit models and store outputs
for (model_run in names(model_data_list)) {
  # Access the current dataframe
  model_data <- model_data_list[[model_run]]
  
  # Calculate the difference between observed and predicted values
  model_data <- model_data %>%
    select(-c(Model))
  
  # Fit a regression random forest model to predict the difference
  rf_model <- ranger(difference~ . , data= model_data , 
                     num.trees =500, mtry = 2, min.node.size = 12,
                     importance = 'permutation', 
                     scale.permutation.importance = TRUE) 
  
  # Store model outputs
  # Accuracy (e.g., R-squared)
  accuracy <- rf_model$r.squared[length(rf_model$r.squared)]
  accuracy_list.kd[[model_run]] <- accuracy
  
  # variable importance
  var_imp <- rf_model$variable.importance
  varimp_list.kd[[model_run]] <- var_imp
  
  # Error (e.g., MSE)
  error <- rf_model$prediction.error
  error_list.kd[[model_run]] <- error
  
  # Combine observed and predicted into a single dataframe for this model run
  model_output <- data.frame(
    observed = model_data$difference,
    predicted = rf_model$predictions,
    site_id = model_data$SITE_ID,
    survey_year = model_data$survey_year
  )
  
  # Add this dataframe to the results list
  model_results.kd[[model_run]] <- model_output
}

# create variable importance data frame
importance_all.kd <- as.data.frame(do.call(rbind, varimp_list.kd)) %>% 
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "importance")
# create model stats dataframes 
rsq_all.kd <- as.data.frame(do.call(rbind, accuracy_list.kd))
error_all.kd <- as.data.frame(do.call(rbind, error_list.kd)) %>% 
  mutate(RMSE = sqrt(V1)) %>% select(-c(V1))

# create predicted and observed data frame
pred_obs.kd <- as.data.frame(do.call(rbind, model_results.kd)) 

# Calculate mean %IncMSE for each variable
perc_mse.kd <- importance_all.kd %>%
  group_by(variable) %>%
  summarise(importance.m = mean(importance))


## VARIABLE IMPORTANCE PLOTS ####

final.df.kd <- importance_all.kd %>% 
  filter(importance >= 0)

# Reorder Variables based on mean IncNodePurity (largest to smallest)
final.df.kd$variable <- factor(final.df.kd$variable, 
                               levels = perc_mse.kd$variable[order(perc_mse.kd$importance.m, 
                                                                   decreasing = FALSE)]) 
#classify each covariate into categories for plotting
importance_data.kd <- final.df.kd %>% 
  mutate(category = case_when(
    variable == "CREEK_CODE" ~ "Location",
    variable == "SITE_ID" ~"Location",
    variable == "survey_year" ~ "Temporal",
    variable == "slope23d_perc_greater" ~ "Topography",
    variable == "area_sqkm.x" ~ "Topography",
    variable == "stream_slope" ~ "Stream Morphology",
    variable == "relative_bankfull" ~ "Temporal",
    variable == "bankfull_mean" ~ "Stream Morphology",
    variable == "Confinement"~ "Stream Morphology",
    variable == "BFIWS" ~ "Climate",
    variable == "elev_mean" ~ "Topography",
    variable == "COMPSTRGTHWS" ~ "Lithology/Soils",
    variable == "HYDRLCONDWS"  ~ "Lithology/Soils",
    variable == "dom_geology" ~ "Lithology/Soils",
    variable == "KFFACTWS" ~ "Lithology/Soils",
    variable == "PERMWS"~ "Lithology/Soils",
    variable == "RCKDEPWS" ~ "Lithology/Soils",
    variable == "OMWS" ~ "Lithology/Soils",
    variable == "PRECIP8110WS" ~ "Climate",
    variable == "TMEAN8110WS" ~ "Climate",
    variable == "high_flow_count" ~ "Climate",
    variable == "mean_ad" ~ "Land Cover",
    variable == "mean_cc" ~ "Land Cover",
    variable == "frequency_km" ~ "Stream Morphology",
    variable == "RDCRSWS" ~ "Land Cover",
    variable == "RDDENSWS" ~ "Land Cover",
    TRUE ~ "Other"
  )) 

# create dataframe for barplots for plotting
bar_data.kd <- importance_data.kd %>% 
  select(variable, category, importance ) %>% 
  group_by(variable, category) %>% 
  summarize(max.imp = max(importance)) %>% 
  ungroup()

## Labels for variable importance plot
variable_labels <- c(
  CREEK_CODE = "AREMP Watershed ID",
  SITE_ID = "AREMP Site ID",
  survey_year = "Survey Year",
  slope23d_perc_greater = "%subwatershed slopes >23 degrees",
  area_sqkm.x = "subwatershed area",
  stream_slope = "stream slope",
  relative_bankfull = "relative bankfull width",
  bankfull_mean = "mean bankfull width",
  Confinement = "valley confinement",
  BFIWS = "baseflow index",
  elev_mean = "mean subwatershed elevation",
  COMPSTRGTHWS = "rock compressive strength",
  HYDRLCONDWS = "rock hydraulic cond.",
  dom_geology = "dominant lithology",
  KFFACTWS = "soil erodibility factor",
  PERMWS = "soil permeability",
  RCKDEPWS = "soil depth",
  OMWS = "soil organic matter",
  PRECIP8110WS = "mean annual precip.",
  TMEAN8110WS = "mean annual temp.",
  high_flow_count = "high flow count",
  mean_ad = "mean stand age",
  mean_cc = "mean canopy cover",
  frequency_km = "large wood frequency",
  RDCRSWS = "road crossings",
  RDDENSWS = "road density",
  ppc_burned_sqkm = "subwatershed area burned",
  ppc_himod_sqkm = "subwatershed area burned himod.",
  rip_burned_sqkm = "riparian area burned",
  rip_himod_sqkm = "riparian area burned himod.",
  rank.x = "fire rank",
  reach_burned_sqkm = "reach area burned",
  reach_himod_sqkm = "reach area burned himod.",
  YPF = "years post-fire"
)




# FINAL VARIABLE IMPORTANCE PLOT
klam.d.VarImp <- ggplot(importance_data.kd, aes(x= importance, y = variable, fill = category)) +  
  geom_bar(data =bar_data.kd, aes(x= max.imp, y = variable, fill = category),
           stat = "identity", alpha = 0.6)+
  geom_boxplot() + 
  geom_segment(aes(x = importance, xend = importance, 
                   y = variable, yend = variable), linetype = "dashed") +
  scale_fill_manual(breaks = c("Climate", "Land Cover", "Lithology/Soils", 
                               "Location", "Stream Morphology", "Temporal",
                               "Topography", "Other"),
                    label = c("Climate", "Land Cover", "Lithology/Soils", 
                              "Location", "Stream Morphology", "Temporal",
                              "Topography", "Fire"),
                    values = c("lemonchiffon2", "seagreen", "thistle4", "goldenrod2",
                               "steelblue", "pink", "burlywood4", "red3"))+
  xlab(expression(paste("perc_MSE"))) +
  ylab("") + 
  theme_bw() + 
  scale_y_discrete(labels = variable_labels)+
  scale_x_continuous(limits = c(0, 36),
                     breaks = seq(0,36, by = 10))+
  theme( axis.text.y = element_text(angle = 0, hjust = 1), 
          panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         text = element_text(size = 30),
         legend.position = "none")+
  labs(title = "Klamath Difference",
       subtitle = "(pred - obs)")


####### Cascades DIFFERENCE MODELING: EXAMINING TEST RUNS ON POST-FIRE SURVEYS ####
#### cREATE POST-FIRE SURVEY DATA FRAME FOR first burned cats ####


pf_covars_filtered.c <- pf_survey_covars %>% 
  group_by(site_survey_id) %>% 
  arrange(site_survey_id, ppc_burned_sqkm) %>% 
  mutate(area_burned_rank = row_number()) %>% 
  ungroup() %>% 
  relocate(area_burned_rank, .after= site_survey_id) %>% 
  filter(area_burned_rank == "1") %>% 
  dplyr::select(c(ppc_burned_sqkm,
                  ppc_himod_sqkm, rip_burned_sqkm, 
                  rip_himod_sqkm,rank.x,
                  reach_burned_sqkm,
                  reach_himod_sqkm, YPF,
                  CREEK_CODE,BFIWS, elev_mean,
                  SITE_ID, survey_year,
                  slope23d_perc_greater,                    
                  area_sqkm.x, stream_slope, 
                  bankfull_mean,relative_bankfull,
                  Confinement,COMPSTRGTHWS, HYDRLCONDWS, dom_geology,                                    
                  KFFACTWS,PERMWS,RCKDEPWS,OMWS, PRECIP8110WS,
                  TMEAN8110WS, high_flow_count,                                    
                  mean_ad, mean_cc, frequency_km,                                    
                  RDCRSWS,RDDENSWS, AQUA_PROV, site_survey_id)) %>% 
  drop_na() %>% 
  filter(AQUA_PROV == "Western Cascades"|
           AQUA_PROV == "High Cascades" |
           AQUA_PROV == "North Cascades") %>% 
  filter(SITE_ID != "CAEMI0015") %>% #remove catchment area outlier
  dplyr::select(-c(AQUA_PROV)) %>% 
  filter(Confinement != "None")

difference.fines.c <- pred_obs.pf.c %>% 
  mutate(difference = predicted - observed) %>% 
  left_join(site_id_link, by = c("site_id" = "site_id", 
                                 "year" = "survey_year")) %>% 
  select(-c(X, OBJECTID, QC_visit)) %>% 
  left_join(pf_covars_filtered.c, by = c("site_survey_id" = "site_survey_id"),
            relationship = "many-to-many")

difference_train.fines.c <- difference.fines.c %>% 
  dplyr::select(c(Model, difference, ppc_burned_sqkm,
                  ppc_himod_sqkm, rip_burned_sqkm, 
                  rip_himod_sqkm,rank.x,
                  reach_burned_sqkm,
                  reach_himod_sqkm, YPF,
                  CREEK_CODE,BFIWS, elev_mean,
                  SITE_ID, survey_year,
                  slope23d_perc_greater,                    
                  area_sqkm.x, stream_slope, 
                  bankfull_mean,relative_bankfull,
                  Confinement,COMPSTRGTHWS, HYDRLCONDWS, dom_geology,                                    
                  KFFACTWS,PERMWS,RCKDEPWS,OMWS, PRECIP8110WS,
                  TMEAN8110WS, high_flow_count,                                    
                  mean_ad, mean_cc, frequency_km,                                    
                  RDCRSWS,RDDENSWS)) %>% 
  drop_na() 




##### Difference Modeling: Cascades ######

# Step 1: Split the dataframe by "Model" and store each subset in a list
model_data_list3 <- split(difference_train.fines.c, difference_train.fines.c$Model)

# Step 2: Initialize lists to store the outputs
accuracy_list.cd <- list()
error_list.cd <- list()
model_results.cd <- list()  # List to store dataframes with observed and predicted together
varimp_list.cd <- list()

# Step 3: Loop through the list of dataframes to fit models and store outputs
for (model_run in names(model_data_list3)) {
  # Access the current dataframe
  model_data <- model_data_list3[[model_run]]
  
  # Calculate the difference between observed and predicted values
  model_data <- model_data %>%
    select(-c(Model))
  
  # Fit a regression random forest model to predict the difference
  rf_model <- ranger(difference~ . , data= model_data , 
                     num.trees =812, mtry = 10, min.node.size = 3,
                     importance = 'permutation', 
                     scale.permutation.importance = TRUE) 
  
  # Store model outputs
  # Accuracy (e.g., R-squared)
  accuracy <- rf_model$r.squared[length(rf_model$r.squared)]
  accuracy_list.cd[[model_run]] <- accuracy
  
  # variable importance
  var_imp <- rf_model$variable.importance
  varimp_list.cd[[model_run]] <- var_imp
  
  
  # Error (e.g., RMSE)
  error <- rf_model$prediction.error
  error_list.cd[[model_run]] <- error
  
  # Combine observed and predicted into a single dataframe for this model run
  model_output <- data.frame(
    observed = model_data$difference,
    predicted = rf_model$predictions,
    site_id = model_data$SITE_ID,
    survey_year = model_data$survey_year
  )
  
  # Add this dataframe to the results list
  model_results.cd[[model_run]] <- model_output
}

# create variable importance data frame
importance_all.cd <- as.data.frame(do.call(rbind, varimp_list.cd)) %>% 
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "importance")
# create model stats dataframes 
rsq_all.cd <- as.data.frame(do.call(rbind, accuracy_list.cd))
error_all.cd <- as.data.frame(do.call(rbind, error_list.cd)) %>% 
  mutate(RMSE = sqrt(V1)) %>% select(-c(V1))

# create predicted and observed data frame
pred_obs.cd <- as.data.frame(do.call(rbind, model_results.cd)) 

# Calculate mean %IncMSE for each variable
perc_mse.cd <- importance_all.cd %>%
  group_by(variable) %>%
  summarise(importance.m = mean(importance))


## VARIABLE IMPORTANCE PLOTS ####

final.df.cd <- importance_all.cd %>% 
  filter(importance >= 0)

# Reorder Variables based on mean IncNodePurity (largest to smallest)
final.df.cd$variable <- factor(final.df.cd$variable, 
                               levels = perc_mse.cd$variable[order(perc_mse.cd$importance.m, 
                                                                   decreasing = FALSE)]) 
#classify each covariate into categories for plotting
importance_data.cd <- final.df.cd %>% 
  mutate(category = case_when(
    variable == "CREEK_CODE" ~ "Location",
    variable == "SITE_ID" ~"Location",
    variable == "survey_year" ~ "Temporal",
    variable == "slope23d_perc_greater" ~ "Topography",
    variable == "area_sqkm.x" ~ "Topography",
    variable == "stream_slope" ~ "Stream Morphology",
    variable == "relative_bankfull" ~ "Temporal",
    variable == "bankfull_mean" ~ "Stream Morphology",
    variable == "Confinement"~ "Stream Morphology",
    variable == "BFIWS" ~ "Climate",
    variable == "elev_mean" ~ "Topography",
    variable == "COMPSTRGTHWS" ~ "Lithology/Soils",
    variable == "HYDRLCONDWS"  ~ "Lithology/Soils",
    variable == "dom_geology" ~ "Lithology/Soils",
    variable == "KFFACTWS" ~ "Lithology/Soils",
    variable == "PERMWS"~ "Lithology/Soils",
    variable == "RCKDEPWS" ~ "Lithology/Soils",
    variable == "OMWS" ~ "Lithology/Soils",
    variable == "PRECIP8110WS" ~ "Climate",
    variable == "TMEAN8110WS" ~ "Climate",
    variable == "high_flow_count" ~ "Climate",
    variable == "mean_ad" ~ "Land Cover",
    variable == "mean_cc" ~ "Land Cover",
    variable == "frequency_km" ~ "Stream Morphology",
    variable == "RDCRSWS" ~ "Land Cover",
    variable == "RDDENSWS" ~ "Land Cover",
    TRUE ~ "Other"
  )) %>% 
  mutate(AP = "Cascades")

# create dataframe for barplots for plotting
bar_data.cd <- importance_data.cd %>% 
  select(variable, category, importance ) %>% 
  group_by(variable, category) %>% 
  summarize(max.imp = max(importance)) %>% 
  ungroup()

# FINAL VARIABLE IMPORTANCE PLOT
casc.cd <- ggplot(importance_data.cd, aes(x= importance, y = variable, fill = category)) +  
  geom_bar(data =bar_data.cd, aes(x= max.imp, y = variable, fill = category),
           stat = "identity", alpha = 0.6)+
  geom_boxplot() + 
  geom_segment(aes(x = importance, xend = importance, 
                   y = variable, yend = variable), linetype = "dashed") +
  scale_fill_manual(breaks = c("Climate", "Land Cover", "Lithology/Soils", 
                               "Location", "Stream Morphology", "Temporal",
                               "Topography", "Other"),
                    label = c("Climate", "Land Cover", "Lithology/Soils", 
                              "Location", "Stream Morphology", "Temporal",
                              "Topography", "Fire"),
                    values = c("lemonchiffon2", "seagreen", "thistle4", "goldenrod2",
                               "steelblue", "pink", "burlywood4", "red3"))+
  xlab(expression(paste("perc_MSE"))) +
  ylab("Variables") + 
  theme_bw() + 
  scale_y_discrete(labels = variable_labels)+
  scale_x_continuous(limits = c(0, 36),
                     breaks = seq(0,36, by = 10))+
  theme( axis.text.y = element_text(angle = 0, hjust = 1), 
           panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         text = element_text(size = 30), legend.position = "none")+
  labs(title = "Cascades Difference",
       subtitle = "(pred - obs)")


#########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
#### Partial Dependence Plot for Klamath ####
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~#########

# Choose variable of interest
target_variable <- "COMPSTRGTHWS"

# Initialize list to store PDPs
pdp_list.k <- list()

# Loop through each model to compute partial dependence
for (i in seq_along(model_data_list)) {
  model_data <- model_data_list[[i]]
  rf_model <- ranger(difference ~ ., data = model_data, 
                     num.trees = 500, mtry = 2, min.node.size = 12,
                     importance = 'permutation', scale.permutation.importance = TRUE)
  
  # Refit each model with ranger to work with partial() properly
  pd <- partial(rf_model, pred.var = target_variable,
                train = model_data, grid.resolution = 20, progress = "none") %>%
    mutate(model_run = names(model_data_list)[i])
  
  pdp_list.k[[i]] <- pd
}

# Combine all PDPs
pdp_all.k <- bind_rows(pdp_list.k) %>% 
  mutate(LineType = "Individual")


# Compute average PDP across all models
pdp_summary.k <- pdp_all.k %>%
  group_by(!!sym(target_variable)) %>%
  summarise(mean_yhat = mean(yhat), .groups = "drop") %>% 
  mutate(LineType = "Average")

pdp_combined <- bind_rows(
  pdp_all.k %>% rename(y = yhat),
  pdp_summary.k %>% rename(y = mean_yhat)
)


scale <- c("gray80",  "black")

ggplot(pdp_combined, aes_string(x = target_variable, y = "y", group = "interaction(LineType, model_run)", linetype = "LineType")) +
  geom_line(aes(color = LineType), size = 1, alpha = ifelse(pdp_combined$LineType == "Individual", 0.5, 1)) +
  scale_linetype_manual(values = c("Individual" = "solid", "Average" = "dotdash")) +
  scale_color_manual(values = c("Individual" = "gray70", "Average" = "black")) +
  labs(x = target_variable,
       y = "Partial Dependence (difference)",
       #title = paste("Partial Dependence of", target_variable),
       linetype = NULL,
       color = NULL) +
  theme_classic() +
  scale_fill_manual(name="", values= scale)+
  theme(legend.position = c(0.15, 0.85),
        text = element_text(size = 25),   # Move legend box to upper right
                         legend.background = element_rect(fill = "white", color = "black"))






#########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####
#### Partial Dependence Plot for Cascades ####
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~#########


# Choose variable of interest
target_variable.c <- "frequency_km"

# Initialize list to store PDPs
pdp_list.c <- list()

# Loop through each model to compute partial dependence
for (i in seq_along(model_data_list3)) {
  model_data <- model_data_list3[[i]]
  rf_model <- ranger(difference ~ ., data = model_data, 
                     num.trees = 500, mtry = 2, min.node.size = 12,
                     importance = 'permutation', scale.permutation.importance = TRUE)
  
  # Refit each model with ranger to work with partial() properly
  pd <- partial(rf_model, pred.var = target_variable.c,
                train = model_data, grid.resolution = 20, progress = "none") %>%
    mutate(model_run = names(model_data_list3)[i])
  
  pdp_list.c[[i]] <- pd
}

# Combine all PDPs
pdp_all.c <- bind_rows(pdp_list.c) %>% 
  mutate(LineType = "Individual")


# Compute average PDP across all models
pdp_summary.c <- pdp_all.c %>%
  group_by(!!sym(target_variable.c)) %>%
  summarise(mean_yhat = mean(yhat), .groups = "drop") %>% 
  mutate(LineType = "Average")

pdp_combined.c <- bind_rows(
  pdp_all.c %>% rename(y = yhat),
  pdp_summary.c %>% rename(y = mean_yhat)
)


scale <- c("gray80",  "black")

ggplot(pdp_combined.c, aes_string(x = target_variable.c, y = "y", group = "interaction(LineType, model_run)", linetype = "LineType")) +
  geom_line(aes(color = LineType), size = 1, alpha = ifelse(pdp_combined.c$LineType == "Individual", 0.5, 1)) +
  scale_linetype_manual(values = c("Individual" = "solid", "Average" = "dotdash")) +
  scale_color_manual(values = c("Individual" = "gray70", "Average" = "black")) +
  labs(x = target_variable.c,
       y = "Partial Dependence (difference)",
    #   title = paste("Partial Dependence of", target_variable.c),
       linetype = NULL,
       color = NULL) +
  theme_classic() +
  scale_fill_manual(name="", values= scale)+
  theme(legend.position = c(0.15, 0.85),  
        text = element_text(size = 25),   # Move legend box to upper right
                         legend.background = element_rect(fill = "white", color = "black"))

