
##-----------------------------------------------------------##
#### *********** D50 (phi) MODELS ************####
##-----------------------------------------------------------##

## CASCADES AND KLAMATH MODELS ##

## Load Libraries####
library(randomForest)
library(tidyverse)
library(dplyr)
library(caret)
library(ranger)
library(pdp)
library(patchwork)
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
  filter(QC_visit.x == "N")%>% 
  mutate(D50_phi = log2(D50))

# site id linkages that will be used for the difference models
site_id_link <- read.csv("aremp_site_info/site_id_linkage.csv") %>% 
  filter(QC_visit == "N")



##### KLAMATH POSTFIRE MODELING ####
pf_data.k <- pf_survey_covars %>% 
  filter(ppc_burned_perc > 1) %>% 
  group_by(site_survey_id) %>% 
  arrange(site_survey_id, ppc_burned_sqkm) %>% 
  mutate(area_burned_rank = row_number()) %>% 
  ungroup() %>% 
  relocate(area_burned_rank, .after= site_survey_id) %>% 
  filter(area_burned_rank == "1") %>% 
  dplyr::select(c(D50_phi, ppc_burned_sqkm,
                  ppc_himod_sqkm, rip_burned_sqkm, 
                  rip_himod_sqkm,
                  reach_burned_sqkm,rank.x,
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
                  RDCRSWS,RDDENSWS, AQUA_PROV)) %>% 
  filter(AQUA_PROV == "Klamath/Siskiyou") %>% 
  filter(SITE_ID != "CAEMI0015") %>% #remove catchment area outlier
  dplyr::select(-c(AQUA_PROV)) %>% 
  drop_na() %>% 
  filter(Confinement != "None")



#### KLAMATH: Run RF Model 50 times and extract results ####
# variable importance with error bars##
importance_list.k <- list()
r_sq.k <- list()
error.k <- list()
model_list.k <- list()
pred_obs <- list()

## create 50 models and summarize ###
for(n in 1:50){
  set.seed(n)
  # mtry = 10, trees = 975, min_n = 13
  rf <- ranger(D50_phi~ . , data=pf_data.k, 
               num.trees =975, mtry = 10, min.node.size = 13,
               importance = 'permutation', 
               scale.permutation.importance = TRUE) 
  
  importance<- rf$variable.importance
  
  df<-  data.frame(variable = names(importance), importance = importance, iteration = n )
  
  accuracy <- data.frame(r.sq = rf$r.squared, iteration = n)
  
  error <- data.frame(error = rf$prediction.error, iteration = n)
  
  predicted <- data.frame(pred = rf$predictions, obs = pf_data.k$D50_phi, 
                          site_id = pf_data.k$SITE_ID,
                          year = pf_data.k$survey_year)
  
  
  importance_list.k[[n]] <- df
  r_sq.k[[n]]<- accuracy
  error.k[[n]]<- error
  model_list.k[[n]]<- rf
  pred_obs[[n]] <- predicted
} 

# create variable importance data frame
importance_all.k <- as.data.frame(do.call(rbind, importance_list.k))
# create model stats dataframes 
rsq_all.k <- as.data.frame(do.call(rbind, r_sq.k))
error_all.k <- as.data.frame(do.call(rbind, error.k))
results.k <- rsq_all.k %>% 
  left_join(error_all.k, by = c("iteration" = "iteration")) %>% 
  mutate(RMSE = sqrt(error)) %>% 
  select(-c(iteration, error))

# create predicted and observed data frame
pred_obs.k <- as.data.frame(do.call(rbind, pred_obs)) 

# Calculate mean %IncMSE for each variable
perc_mse.k <- importance_all.k %>%
  group_by(variable) %>%
  summarise(importance.m = mean(importance))%>% 
  arrange(desc(importance.m)) %>% 
  mutate(rank = row_number()) %>% 
  filter(rank < 11) %>%  dplyr::select(-rank)

#### CREATE PARTIAL DEPENDENCE PLOTS for Variable Importance Figure ####
# build one RF model to use for the PDPs
rf.k <- ranger(D50_phi~ . , data=pf_data.k, 
                     num.trees =1325, mtry = 5, min.node.size = 1,
                     importance = 'permutation', 
                     scale.permutation.importance = TRUE) 

#identify all covariates to loop through
top_vars.k <- perc_mse.k$variable[1:10]
#create list to store plots
pdp_list.k <- list()
# for loop to cycle through all covariates and create partial dependence dataframe
for (var in top_vars.k) {
  # extract partial dependence info
  pd <- partial(
    object = rf.k,
    pred.var = var,
    train = pf_data.k
  )
  #store partial dependence info for each covariate (covariate values, yhat values, and covariate name)
  pd$var <- var 
  pdp_list.k[[var]] <- pd
}
# standardize the column titles for all dataframes in the pdp_list
pdp_list2.k <- lapply(pdp_list.k, function(df){
  colnames(df)[1]<- "value"
  df})
# combine all dataframes in list into one dataframe
pdp_data.k <- do.call(rbind, pdp_list2.k)

#rescale partial dependence values for plotting
pdp_data2.k <- pdp_data.k %>%
  mutate(value = as.numeric(value)) %>% 
  rename(variable = var) %>% 
  group_by(variable) %>% 
  mutate(x_rescaled = (value - min(value)) /(max(value)-min(value))) %>% 
  ungroup() %>% 
  drop_na(yhat) 

## create pdp skeleton plot
pdp_plot.k <- ggplot(pdp_data2.k, aes(x = x_rescaled, y = yhat))+
  geom_line()+
  geom_line(linewidth=0.5)+
  facet_grid(scales = "free", 
             rows = vars(factor(variable,
                                perc_mse.k$variable[order(perc_mse.k$importance.m,
                                                          decreasing = TRUE)])))+
  labs(x="", y="") +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_line(color="transparent"),
    axis.line = element_line(color="transparent"),
    rect = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
  )

## VARIABLE IMPORTANCE PLOTS ####
final.df.k <- importance_all.k 

# Reorder Variables based on mean IncNodePurity (largest to smallest)
final.df.k$variable <- factor(final.df.k$variable, 
                              levels = perc_mse.k$variable[order(perc_mse.k$importance.m, 
                                                                 decreasing = FALSE)]) 
#classify each covariate into categories for plotting
importance_data.k <- final.df.k %>% 
  mutate(category = case_when(
    variable == "CREEK_CODE" ~ "Location",
    variable == "SITE_ID" ~ "Location",
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
  mutate(AP = "Klamath")%>% 
  filter(variable != "<NA>")

# create dataframe for barplots for plotting
bar_data <- importance_data.k %>% 
  select(variable, category, importance ) %>% 
  group_by(variable, category) %>% 
  summarize(max.imp = max(importance)) %>% 
  ungroup()

# FINAL VARIABLE IMPORTANCE PLOT
klam <- ggplot(importance_data.k, aes(x= importance, y = variable, fill = category)) +  
  geom_bar(data =bar_data, aes(x= max.imp, y = variable, fill = category),
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
  theme( axis.text.y = element_text(angle = 0, hjust = 1), 
           panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         text = element_text(size = 25), 
         legend.position = "none")+
  scale_x_continuous(limits = c(0, 44),
                     breaks = seq(0,44, by = 10))+
#  ggtitle("Klamath Postfire D50 Model")+
  inset_element(pdp_plot.k, left = -0.0005, bottom = -0.03, right = 0.15, top = 1)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
##### CASCADES POSTFIRE MODELING ####
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~#####

pf_data.c <- pf_survey_covars %>% 
  filter(ppc_burned_perc > 1) %>% 
  group_by(site_survey_id) %>% 
  arrange(site_survey_id, ppc_burned_sqkm) %>% 
  mutate(area_burned_rank = row_number()) %>% 
  ungroup() %>% 
  relocate(area_burned_rank, .after= site_survey_id) %>% 
  filter(area_burned_rank == "1") %>% 
  filter(AQUA_PROV == "Western Cascades"|
           AQUA_PROV == "High Cascades" |
           AQUA_PROV == "North Cascades") %>% 
  dplyr::select(c(D50_phi, ppc_burned_sqkm,
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
                  RDCRSWS,RDDENSWS, AQUA_PROV)) %>% 
  filter(SITE_ID != "CAEMI0015") %>% #remove catchment area outlier
  dplyr::select(-c(AQUA_PROV)) %>% 
  drop_na() %>% 
  filter(Confinement != "None")


#### CASCADES: Run RF Model 50 times and extract results ####
# variable importance with error bars##
importance_list.c <- list()
r_sq.c <- list()
error.c <- list()
model_list.c <- list()
pred_obs.c <- list()

## create 50 models and summarize ###
for(n in 1:50){
  set.seed(n)
  # mtry = 5, trees = 1325, min_n = 1
  rf <- ranger(D50_phi~ . , data=pf_data.c, 
               num.trees =1325, mtry = 5, min.node.size = 1,
               importance = 'permutation', 
               scale.permutation.importance = TRUE) 
  
  importance<- rf$variable.importance
  
  df<-  data.frame(variable = names(importance), importance = importance, iteration = n )
  
  accuracy <- data.frame(r.sq = rf$r.squared, iteration = n)
  
  error <- data.frame(error = rf$prediction.error, iteration = n)
  
  predicted <- data.frame(pred = rf$predictions, obs = pf_data.c$D50_phi,
                          site_id = pf_data.c$SITE_ID,
                          year = pf_data.c$survey_year)
  
  importance_list.c[[n]] <- df
  r_sq.c[[n]]<- accuracy
  error.c[[n]]<- error
  model_list.c[[n]]<- rf
  pred_obs.c[[n]] <- predicted
} 

# create variable importance data frame
importance_all.c <- as.data.frame(do.call(rbind, importance_list.c))
# create model stats dataframes 
rsq_all.c <- as.data.frame(do.call(rbind, r_sq.c))
error_all.c <- as.data.frame(do.call(rbind, error.c))
results.c <- rsq_all.c %>% 
  left_join(error_all.c, by = c("iteration" = "iteration")) %>% 
  mutate(RMSE = sqrt(error)) %>% 
  select(-c(iteration, error))
# create predicted and observed data frame
pred_obs.c <- as.data.frame(do.call(rbind, pred_obs.c))

# Calculate mean %IncMSE for each variable
perc_mse.c <- importance_all.c %>%
  group_by(variable) %>%
  summarise(importance.m = mean(importance))%>% 
  arrange(desc(importance.m)) %>% 
  mutate(rank = row_number()) %>% 
  filter(rank < 11) %>%  dplyr::select(-rank)

#### CREATE PARTIAL DEPENDENCE PLOTS for Variable Importance Figure ####
# build one RF model to use for the PDPs
rf.c <- ranger(D50_phi~ . , data=pf_data.c, 
               num.trees =1325, mtry = 5, min.node.size = 1,
               importance = 'permutation', 
               scale.permutation.importance = TRUE) 

#identify all covariates to loop through
top_vars.c <- perc_mse.c$variable[1:10]
#create list to store plots
pdp_list.c <- list()
# for loop to cycle through all covariates and create partial dependence dataframe
for (var in top_vars.c) {
  # extract partial dependence info
  pd <- partial(
    object = rf.c,
    pred.var = var,
    train = pf_data.c
  )
  #store partial dependence info for each covariate (covariate values, yhat values, and covariate name)
  pd$var <- var 
  pdp_list.c[[var]] <- pd
}
# standardize the column titles for all dataframes in the pdp_list
pdp_list2.c <- lapply(pdp_list.c, function(df){
  colnames(df)[1]<- "value"
  df})
# combine all dataframes in list into one dataframe
pdp_data.c <- do.call(rbind, pdp_list2.c)

#rescale partial dependence values for plotting
pdp_data2.c <- pdp_data.c %>%
  mutate(value = as.numeric(value)) %>% 
  rename(variable = var) %>% 
  group_by(variable) %>% 
  mutate(x_rescaled = (value - min(value)) /(max(value)-min(value))) %>% 
  ungroup() %>% 
  drop_na(yhat) 

## create pdp skeleton plot
pdp_plot.c <- ggplot(pdp_data2.c, aes(x = x_rescaled, y = yhat))+
  geom_line()+
  geom_line(linewidth=0.5)+
  facet_grid(scales = "free", 
             rows = vars(factor(variable,
                                perc_mse.c$variable[order(perc_mse.c$importance.m,
                                                          decreasing = TRUE)])))+
  labs(x="", y="") +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_line(color="transparent"),
    axis.line = element_line(color="transparent"),
    rect = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
  )

## VARIABLE IMPORTANCE PLOTS ####
final.df.c <- importance_all.c 

# Reorder Variables based on mean IncNodePurity (largest to smallest)
final.df.c$variable <- factor(final.df.c$variable, 
                              levels = perc_mse.c$variable[order(perc_mse.c$importance.m, 
                                                                 decreasing = FALSE)])
#classify each covariate into categories for plotting
importance_data.c <- final.df.c %>% 
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
  mutate(AP = "Cascades") %>% 
  filter(variable != "<NA>")

# create dataframe for barplots for plotting
bar_data.c <- importance_data.c %>% 
  select(variable, category, importance ) %>% 
  group_by(variable, category) %>% 
  summarize(max.imp = max(importance)) %>% 
  ungroup()

# FINAL VARIABLE IMPORTANCE PLOT
casc <- ggplot(importance_data.c, aes(x= importance, y = variable, fill = category)) +  
  geom_bar(data =bar_data.c, aes(x= max.imp, y = variable, fill = category),
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
  scale_x_continuous(limits = c(0, 44),
                     breaks = seq(0,44, by = 10))+
  theme( axis.text.y = element_text(angle = 0, hjust = 1), 
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         text = element_text(size = 25),
         legend.position = "none")+
#  ggtitle("Cascades postfire D50 Model")+
  inset_element(pdp_plot.c, left = -0.0005, bottom = -0.032, right = 0.15, top = 1)





###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#######
#### PARTIAL DEPENDENCE PLOTS ####
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#######

## CASCADES ####

# Create a directory to save the plots
output_dir <- "FIGURES/supplement/PDP/postfire_casc_D50"
dir.create(output_dir, showWarnings = FALSE)

# Define top ten target variables
target_variables <- c("bankfull_mean", "area_sqkm.x", "PERMWS", 
                      "RCKDEPWS", "COMPSTRGTHWS","ppc_himod_sqkm", "elev_mean", 
                      "survey_year", "KFFACTWS", "SITE_ID")

# Loop over variables
for (var in target_variables) {
  
  pdp_list <- list()
  
  for (i in seq_along(model_list.c)) {
    rf_model <- ranger(D50_phi ~ ., data = pf_data.c,
                       num.trees = 875, mtry = 7, min.node.size = 2,
                       importance = 'permutation', scale.permutation.importance = TRUE)
    
    pd <- partial(rf_model, pred.var = var,
                  train = pf_data.c, grid.resolution = 20, progress = "none") %>%
      mutate(model_run = names(model_list.c)[i])
    
    pdp_list[[i]] <- pd
  }
  
  # Combine all partials
  pdp_all <- bind_rows(pdp_list)
  
  # Summarize across models
  pdp_summary <- pdp_all %>%
    group_by_at(var) %>%
    summarise(mean_yhat = mean(yhat), .groups = "drop") %>%
    mutate(LineType = "Average")
  
  # Plot
  p <- ggplot(pdp_summary, aes_string(x = var, y = "mean_yhat")) +
    geom_line(color = "black", size = 1) +
    labs(x = var,
         y = "Partial Dependence",
         title = paste("Partial Dependence of", var)) +
    theme_classic() +
    theme(
      text = element_text(size = 16),
      legend.position = "none"
    )
  
  # Save plot
  ggsave(filename = file.path(output_dir, paste0("PDP_", var, ".png")),
         plot = p, width = 6, height =4.5)
}



## Klamath ####

# Create a directory to save the plots
output_dir.k <- "FIGURES/supplement/PDP/postfire_klam_D50"
dir.create(output_dir.k, showWarnings = FALSE)

# Define your ten target variables
target_variables.k <- c("bankfull_mean", "OMWS", "COMPSTRGTHWS", 
                        "survey_year", "PRECIP8110WS","slope23d_perc_greater", "ppc_burned_sqkm", 
                        "area_sqkm.x", "RDDENSWS", "ppc_himod_sqkm")

# Loop over variables
for (var in target_variables.k) {
  
  pdp_list <- list()
  
  for (i in seq_along(model_list.k)) {
    rf_model <- ranger(D50_phi ~ ., data = pf_data.k,
                       num.trees = 875, mtry = 7, min.node.size = 2,
                       importance = 'permutation', scale.permutation.importance = TRUE)
    
    pd <- partial(rf_model, pred.var = var,
                  train = pf_data.k, grid.resolution = 20, progress = "none") %>%
      mutate(model_run = names(model_list.k)[i])
    
    pdp_list[[i]] <- pd
  }
  
  # Combine all partials
  pdp_all <- bind_rows(pdp_list)
  
  # Summarize across models
  pdp_summary <- pdp_all %>%
    group_by_at(var) %>%
    summarise(mean_yhat = mean(yhat), .groups = "drop") %>%
    mutate(LineType = "Average")
  
  # Plot
  p <- ggplot(pdp_summary, aes_string(x = var, y = "mean_yhat")) +
    geom_line(color = "black", size = 1) +
    labs(x = var,
         y = "Partial Dependence",
         title = paste("Partial Dependence of", var)) +
    theme_classic() +
    theme(
      text = element_text(size = 16),
      legend.position = "none"
    )
  
  # Save plot
  ggsave(filename = file.path(output_dir.k, paste0("PDP_", var, ".png")),
         plot = p, width = 6, height =4.5)
}
