####******** FINAL D50 phi RF MODEL RUNS ********************######


#### Broadest model categories will be: Klamath and Cascades
##### We will test each category on all post-fire sites as well as fire 
##### thresholds: 0-25%, 25-50%, 50-75%, 75-100%
##### In the case that mutliple fires have burned a pour point catchment,
##### Fire x Survey interactions will be selected for the fire that burned 
##### the most area in the catchment. No repeat site_survey_ids.

# Load Libraries####
library(randomForest)
library(tidyverse)
library(dplyr)
library(caret)
library(ranger)
library(pdp)
library(patchwork)
library(lattice)
library(cowplot)

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

# all unburned and pre-fire surveys
site_survey_covars <- read.csv("all_site_survey_covars.csv") %>% 
  anti_join(pf_survey_covars, by = c("site_survey_id" = "site_survey_id")) %>% 
  left_join(bankfull, by = c("SITE_ID" = "SITE_ID")) %>% 
  mutate(relative_bankfull = bankfull_mean / bankfull_all_mean) %>% 
  filter(QC_visit.x == "N") %>% 
  mutate(D50_phi = log2(D50))

# site id linkages that will be used for the difference models
site_id_link <- read.csv("aremp_site_info/site_id_linkage.csv") %>% 
  filter(QC_visit == "N") 


#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#### ***** MODEL TESTING FUNCTION ***** ####
#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

test_models <- function(models_list, new_df, response){
  # covariates
  x_new <- new_df %>%  dplyr::select(-all_of(response))
  # response variable
  y_new <- new_df[[response]]
  
  #store results
  results <- data.frame(
    Model = integer(),
    R2 = numeric(),
    MSE = numeric()
  )
  # create predictions list
  predictions_list <- list()
  
  #For loop to cycle through stored models and test on new postfire datasets
  for (i in 1:length(models_list)){
    #stored model list
    model <- models_list[[i]]
    #predict response variable using stored model and new dataset (x_new)
    pred_test <- predict(model, x_new)
    #store Rsq
    variance_explained <- cor(pred_test$predictions, y_new)^2*100
    #store rmse
    rmse <- sqrt(mean((pred_test$predictions - y_new)^2))
    #create results dataframe
    results <- rbind(results, data.frame(Model = i, R2 = variance_explained,
                                         RMSE = rmse))
    #store in predictions list along with other information
    predictions_list[[i]] <- data.frame(
      Model = i,
      observed = y_new,
      predicted = pred_test$predictions,
      site_id = x_new$SITE_ID,
      year = x_new$survey_year.x
    )
  }
  
  return(list(
    summary_stats = results,
    predictions = predictions_list
  ))
}

#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#### ***** KLAMATH MODELING ***** ####
#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## KLAMATH UNBURNED TRAINING Dataset ####
D50_phi_data.k=subset(site_survey_covars,select=c(D50_phi, CREEK_CODE,
                                                     SITE_ID, survey_year.x,
                                                     slope23d_perc_greater,                    
                                                     area_sqkm.x, stream_slope, 
                                                     relative_bankfull, bankfull_mean,
                                                     AQUA_PROV,
                                                     Confinement,BFIWS, elev_mean,
                                                     
                                                     # lithology                  
                                                     COMPSTRGTHWS, HYDRLCONDWS, dom_geology,                                    
                                                     
                                                     # soils                  
                                                     KFFACTWS,PERMWS,RCKDEPWS,OMWS,                                    
                                                     
                                                     # climate                  
                                                     PRECIP8110WS,TMEAN8110WS, high_flow_count,                                    
                                                     
                                                     # vegetation                 
                                                     mean_ad, mean_cc, frequency_km,                                    
                                                     
                                                     # roads                  
                                                     RDCRSWS,RDDENSWS)) %>% 
  drop_na() %>% 
  filter(AQUA_PROV == "Klamath/Siskiyou") %>% 
  filter(SITE_ID != "CAEMI0015") %>% #remove catchment area outlier
  mutate(SITE_ID = as.factor(SITE_ID),
         CREEK_CODE = as.factor(CREEK_CODE)) %>% 
  # mutate(D50_phi2 = log10(D50_phi+1)) %>%
  dplyr::select(-c(AQUA_PROV)) %>% 
  filter(Confinement != "None")


#### KLAMATH: Run Unburned RF Model 50 times and extract results ####

#create lists 
importance_list.k <- list()
r_sq.k <- list()
error.k <- list()
model_list.k <- list()
pred_obs <- list()

## create 50 models and summarize ##
for(n in 1:50){
  
  # run RF model 
  rf <- ranger(D50_phi~ . , data=D50_phi_data.k, 
               num.trees =750, mtry = 2, min.node.size = 1,
               importance = 'permutation', 
               scale.permutation.importance = TRUE) 
  #get variable importance info and store in list
  importance<- rf$variable.importance
  df<-  data.frame(variable = names(importance), importance = importance, iteration = n )
  importance_list.k[[n]] <- df
  
  #get model accuracy and error info and store in list
  accuracy <- data.frame(r.sq = rf$r.squared, iteration = n)
  error <- data.frame(error = rf$prediction.error, iteration = n)
  r_sq.k[[n]]<- accuracy
  error.k[[n]]<- error
  
  #get predicted values from model runs and store along with observed values and site info
  predicted <- data.frame(pred = rf$predictions, obs = D50_phi_data.k$D50_phi, 
                          site_id = D50_phi_data.k$SITE_ID,
                          year = D50_phi_data.k$survey_year.x)
  pred_obs[[n]] <- predicted
  
  
  #store all 50 models into a list
  model_list.k[[n]]<- rf
  
} 

#create variable importance dataframe
importance_all.k <- as.data.frame(do.call(rbind, importance_list.k))
#create model stats dataframes: Rsq and error (model reports MSE as error)
rsq_all.k <- as.data.frame(do.call(rbind, r_sq.k))
error_all.k <- as.data.frame(do.call(rbind, error.k))
#create dataframe for the predicted and observed values
pred_obs.k <- as.data.frame(do.call(rbind, pred_obs)) 

#create full results dataframe for model performance (Rsq. and RMSE)
results.k <- rsq_all.k %>% 
  left_join(error_all.k, by = c("iteration" = "iteration")) %>% 
  mutate(RMSE = sqrt(error)) %>% 
  select(-c(iteration, error)) 

# Calculate mean %IncMSE for each variable
perc_mse.k <- importance_all.k %>%
  group_by(variable) %>%
  summarise(importance.m = mean(importance))%>% 
  arrange(desc(importance.m)) %>% 
  mutate(rank = row_number()) %>% 
  filter(rank < 11) %>%  dplyr::select(-rank)



#### CREATE PARTIAL DEPENDENCE PLOTS for Variable Importance Figure ####
# build one RF model to use for the PDPs
rf.k <- ranger(D50_phi~ . , data=D50_phi_data.k, 
               num.trees =750, mtry = 2, min.node.size = 1,
               importance = 'permutation', 
               scale.permutation.importance = TRUE)

#identify all covariates to loop through
top_vars <- perc_mse.k$variable[1:10]
#create list to store plots
pdp_list <- list()
# for loop to cycle through all covariates and create partial dependence dataframe
for (var in top_vars) {
  # extract partial dependence info
  pd <- partial(
    object = rf.k,
    pred.var = var,
    train = D50_phi_data.k
  )
  #store partial dependence info for each covariate (covariate values, yhat values, and covariate name)
  pd$var <- var 
  pdp_list[[var]] <- pd
}
# standardize the column titles for all dataframes in the pdp_list
pdp_list2 <- lapply(pdp_list, function(df){
  colnames(df)[1]<- "value"
  df})
# combine all dataframes in list into one dataframe
pdp_data <- do.call(rbind, pdp_list2)

#rescale partial dependence values for plotting
pdp_data2 <- pdp_data %>%
  mutate(value = as.numeric(value)) %>% 
  rename(variable = var) %>% 
  group_by(variable) %>% 
  mutate(x_rescaled = (value - min(value)) /(max(value)-min(value))) %>% 
  ungroup() %>% 
  drop_na(yhat) 

## create pdp skeleton plot
pdp_plot <- ggplot(pdp_data2, aes(x = x_rescaled, y = yhat))+
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
# create dataframe of importance values to analyze
final.df.k <- importance_all.k %>% 
  filter(importance >=0)
# Reorder Variables based on mean perc_mse (largest to smallest)
final.df.k$variable <- factor(final.df.k$variable, 
                              levels = perc_mse.k$variable[order(perc_mse.k$importance.m, 
                                                                 decreasing = FALSE)]) 
# classify each covariate into categories for plotting
importance_data.k <- final.df.k %>% 
  mutate(category = case_when(
    variable == "CREEK_CODE" ~ "Location",
    variable == "SITE_ID" ~"Location",
    variable == "survey_year.x" ~ "Temporal",
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
    variable == "RDDENSWS" ~ "Land Cover"
  )) %>% 
  mutate(AP = "Klamath")%>% 
  filter(variable != "<NA>")

# create dataframe for barplots for plotting
bar_data <- importance_data.k %>% 
  select(variable, category, importance ) %>% 
  group_by(variable, category) %>% 
  summarize(max.imp = max(importance)) %>% 
  ungroup()

# FINAL VARIABLE IMPORTANCE PLOTTING
klam.VarImp <- ggplot(importance_data.k, aes(x= importance, y = variable, fill = category)) +  
  geom_bar(data =bar_data, aes(x= max.imp, y = variable, fill = category),
           stat = "identity", alpha = 0.6)+
  geom_boxplot() + 
  geom_segment(aes(x = importance, xend = importance, 
                   y = variable, yend = variable), linetype = "dashed") +
  scale_fill_manual(breaks = c("Climate", "Land Cover", "Lithology/Soils", 
                               "Location", "Stream Morphology", "Temporal",
                               "Topography"),
                    label = c("Climate", "Land Cover", "Lithology/Soils", 
                              "Location", "Stream Morphology", "Temporal",
                              "Topography"),
                    values = c("lemonchiffon2", "seagreen", "thistle4", "goldenrod2",
                               "steelblue", "pink", "burlywood4"))+
  xlab(expression(paste("perc_MSE"))) +
  ylab("Variables") + 
  theme_bw() + 
  theme( axis.text.y = element_text(angle = 0, hjust = 1), 
           panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         text = element_text(size = 25),
         legend.position = "none")+
  scale_x_continuous(limits = c(0, 34),
                     breaks = seq(0,34, by = 10))+
 # ggtitle("Klamath Unburned D50 Model")+
  inset_element(pdp_plot, left = -0.0005, bottom = -0.03, right = 0.15, top = 1)

#-----------------------------------------------------------------------#
#### **** KLAMATH UNBURNED MODEL TESTING ON POSTFIRE SURVEYS ******####
#-----------------------------------------------------------------------#

# POST-FIRE SURVEY DATA FRAMEs SELECT Burned Cats for which fire burned the MOST PPC AREA

#### TEST MODELS ON ALL POST-FIRE SURVEY DATA ####
D50_phi_test.pf.k <- pf_survey_covars %>% 
  group_by(site_survey_id) %>% 
  arrange(site_survey_id, ppc_burned_sqkm) %>% 
  mutate(area_burned_rank = row_number()) %>% 
  ungroup() %>% 
  relocate(area_burned_rank, .after= site_survey_id) %>% 
  filter(area_burned_rank == "1") %>% 
  dplyr::select(c(D50_phi, CREEK_CODE,
                  SITE_ID, survey_year.x,
                  slope23d_perc_greater,                    
                  area_sqkm.x, stream_slope, 
                  relative_bankfull, bankfull_mean,
                  AQUA_PROV,
                  Confinement,BFIWS, elev_mean,
                  
                  # lithology                  
                  COMPSTRGTHWS, HYDRLCONDWS, dom_geology,                                    
                  
                  # soils                  
                  KFFACTWS,PERMWS,RCKDEPWS,OMWS,                                    
                  
                  # climate                  
                  PRECIP8110WS,TMEAN8110WS, high_flow_count,                                    
                  
                  # vegetation                 
                  mean_ad, mean_cc, frequency_km,                                    
                  
                  # roads                  
                  RDCRSWS,RDDENSWS)) %>% 
  drop_na() %>% 
  filter(AQUA_PROV == "Klamath/Siskiyou") %>% 
  filter(SITE_ID != "CAEMI0015") %>% #remove catchment area outlier
  mutate(SITE_ID = as.factor(SITE_ID),
         CREEK_CODE = as.factor(CREEK_CODE)) %>% 
  dplyr::select(-c(AQUA_PROV)) %>% 
  filter(Confinement != "None")

### Test on all Klamath post-fire sites 
# Run function test_models to test unburned model on postfire sites
results.pf.k <- test_models(model_list.k, D50_phi_test.pf.k, 'D50_phi')
#compile summary stats for model testing
summary_stats.pf.k <- results.pf.k$summary_stats
#compile dataframe of all predicted and observed values
predictions.pf.k <- do.call(rbind, results.pf.k$predictions)
#write_csv(summary_stats.pf.k, "RF_models/FINAL_RF_MODELS/D50/post_fire_testing/klamath_D50_phi_post_fire_testing.csv")
#write_csv(predictions.pf.k, "RF_models/FINAL_RF_MODELS/D50/post_fire_testing/predicted_observed/klamath_D50_phi_post_fire_pred_obs.csv")



#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
#### ***** CASCADES MODELING ***** ####
#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

## CASCADES UNBURNED TRAINING DATA ####
D50_phi_data.c=subset(site_survey_covars,select=c(D50_phi, CREEK_CODE,
                                                     SITE_ID, survey_year.x,
                                                     slope23d_perc_greater,                    
                                                     area_sqkm.x, stream_slope, 
                                                     relative_bankfull, bankfull_mean,
                                                     AQUA_PROV,
                                                     Confinement,BFIWS, elev_mean,
                                                     
                                                     # lithology                  
                                                     COMPSTRGTHWS, HYDRLCONDWS, dom_geology,                                    
                                                     
                                                     # soils                  
                                                     KFFACTWS,PERMWS,RCKDEPWS,OMWS,                                    
                                                     
                                                     # climate                  
                                                     PRECIP8110WS,TMEAN8110WS, high_flow_count,                                    
                                                     
                                                     # vegetation                 
                                                     mean_ad, mean_cc, frequency_km,                                    
                                                     
                                                     # roads                  
                                                     RDCRSWS,RDDENSWS)) %>% 
  drop_na() %>% 
  filter(AQUA_PROV == "Western Cascades"|
           AQUA_PROV == "High Cascades" |
           AQUA_PROV == "North Cascades") %>% 
  filter(SITE_ID != "CAEMI0015") %>% #remove catchment area outlier
  mutate(SITE_ID = as.factor(SITE_ID),
         CREEK_CODE = as.factor(CREEK_CODE)) %>% 
  dplyr::select(-c(AQUA_PROV)) %>% 
  filter(Confinement != "None")


#### CASCADES: Run RF UNBURNED Model 50 times and extract results ####

#create lists
importance_list.c <- list()
r_sq.c <- list()
error.c <- list()
model_list.c <- list()
pred_obs.c <- list()

## create 50 models and summarize ###
for(n in 1:50){
  
  # Run RF Model
  rf <- ranger(D50_phi~ . , data=D50_phi_data.c, 
               num.trees =750, mtry = 10, min.node.size = 2,
               importance = 'permutation', 
               scale.permutation.importance = TRUE) 
  
  #get variable importance info and store in list
  importance<- rf$variable.importance
  df<-  data.frame(variable = names(importance), importance = importance, iteration = n )
  importance_list.c[[n]] <- df
  
  #get model accuracy and error info and store in list
  accuracy <- data.frame(r.sq = rf$r.squared, iteration = n)
  error <- data.frame(error = rf$prediction.error, iteration = n)
  r_sq.c[[n]]<- accuracy
  error.c[[n]]<- error
  
  #get predicted values from model runs and store along with observed values and site info
  predicted <- data.frame(pred = rf$predictions, obs = D50_phi_data.c$D50_phi, 
                          site_id = D50_phi_data.c$SITE_ID,
                          year = D50_phi_data.c$survey_year.x)
  pred_obs.c[[n]] <- predicted
  
  
  #store all 50 models into a list
  model_list.c[[n]]<- rf
} 

#create variable importance dataframe
importance_all.c <- as.data.frame(do.call(rbind, importance_list.c))
#create model stats dataframes: Rsq and error (model reports MSE as error)
rsq_all.c <- as.data.frame(do.call(rbind, r_sq.c))
error_all.c <- as.data.frame(do.call(rbind, error.c))
#create dataframe for the predicted and observed values
pred_obs.c <- as.data.frame(do.call(rbind, pred_obs.c))

#create full results dataframe for model performance (Rsq. and RMSE)
results.c <- rsq_all.c %>% 
  left_join(error_all.c, by = c("iteration" = "iteration")) %>% 
  mutate(RMSE = sqrt(error)) %>% 
  select(-c(iteration, error))

# Calculate mean %IncMSE for each variable
perc_mse.c <- importance_all.c %>%
  group_by(variable) %>%
  summarise(importance.m = mean(importance)) %>% 
  arrange(desc(importance.m)) %>% 
  mutate(rank = row_number()) %>% 
  filter(rank < 11) %>%  dplyr::select(-rank)




#### CREATE PARTIAL DEPENDENCE PLOTS for Variable Importance Figure ####
# build one RF model to use for the PDPs
rf.c <- ranger(D50_phi~ . , data=D50_phi_data.c, 
               num.trees =750, mtry = 10, min.node.size = 2,
               importance = 'permutation', 
               scale.permutation.importance = TRUE) 
#identify all covariates to loop through
top_vars.c <- perc_mse.c$variable[1:10]
#create list to store plot info
pdp_list.c <- list()
# for loop to cycle through all covariates and create partial dependence dataframe
for (var in top_vars.c) {
  # extract partial dependence info
  pd <- partial(
    object = rf.c,
    pred.var = var,
    train = D50_phi_data.c
  )
  
  #store partial dependence info for each covariate (covariate values, yhat values, and covariate name)
  pd$var <- var 
  pdp_list.c[[var]] <- pd
}
# standardize the column titles for all dataframes in the pdp_list
pdp_list2.c <- lapply(pdp_list.c, function(df){
  colnames(df)[1]<- "value"
  df
})
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
# create dataframe of importance values to analyze
final.df.c <- importance_all.c 

# Reorder Variables based on mean IncNodePurity (largest to smallest)
final.df.c$variable <- factor(final.df.c$variable, 
                              levels = perc_mse.c$variable[order(perc_mse.c$importance.m, 
                                                                 decreasing = FALSE)])
# classify each covariate into categories for plotting
importance_data.c <- final.df.c %>% 
  mutate(category = case_when(
    variable == "CREEK_CODE" ~ "Location",
    variable == "SITE_ID" ~"Location",
    variable == "survey_year.x" ~ "Temporal",
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
    variable == "RDDENSWS" ~ "Land Cover"
  )) %>% 
  mutate(AP = "Cascades") %>% 
  filter(variable != "<NA>")

# create dataframe for barplots for plotting
bar_data.c <- importance_data.c %>% 
  select(variable, category, importance ) %>% 
  group_by(variable, category) %>% 
  summarize(max.imp = max(importance)) %>% 
  ungroup()

# FINAL VARIABLE IMPORTANCE PLOTTING
casc.VarImp <- ggplot(importance_data.c, aes(x= importance, y = variable, fill = category)) +  
  geom_bar(data =bar_data.c, aes(x= max.imp, y = variable, fill = category),
           stat = "identity", alpha = 0.6)+
  geom_boxplot() + 
  geom_segment(aes(x = importance, xend = importance, 
                   y = variable, yend = variable), linetype = "dashed") +
  scale_fill_manual(breaks = c("Climate", "Land Cover", "Lithology/Soils", 
                               "Location", "Stream Morphology", "Temporal",
                               "Topography"),
                    label = c("Climate", "Land Cover", "Lithology/Soils", 
                              "Location", "Stream Morphology", "Temporal",
                              "Topography"),
                    values = c("lemonchiffon2", "seagreen", "thistle4", "goldenrod2",
                               "steelblue", "pink", "burlywood4"))+
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
 # ggtitle("Cascades Unburned D50 Model")+
  inset_element(pdp_plot.c, left = -.0005, bottom = -0.027, right = 0.15, top = 1)



#-----------------------------------------------------------------------#
#### **** CASCADES UNBURNED MODEL TESTING ON POSTFIRE SURVEYS ******####
#-----------------------------------------------------------------------#

# POST-FIRE SURVEY DATA FRAMEs SELECT Burned Cats for which fire burned the MOST PPC AREA

#### TEST MODELS ON ALL POST-FIRE SURVEY DATA ####

D50_phi_test.pf.c <- pf_survey_covars %>% 
  group_by(site_survey_id) %>% 
  arrange(site_survey_id, ppc_burned_sqkm) %>% 
  mutate(area_burned_rank = row_number()) %>% 
  ungroup() %>% 
  relocate(area_burned_rank, .after= site_survey_id) %>% 
  filter(area_burned_rank == "1") %>% 
  filter(AQUA_PROV == "Western Cascades"|
           AQUA_PROV == "High Cascades" |
           AQUA_PROV == "North Cascades") %>%
  dplyr::select(c(D50_phi, CREEK_CODE,
                  SITE_ID, survey_year.x,
                  slope23d_perc_greater,                    
                  area_sqkm.x, stream_slope, 
                  relative_bankfull, bankfull_mean,
                  AQUA_PROV,
                  Confinement,BFIWS, elev_mean,
                  
                  # lithology                  
                  COMPSTRGTHWS, HYDRLCONDWS, dom_geology,                                    
                  
                  # soils                  
                  KFFACTWS,PERMWS,RCKDEPWS,OMWS,                                    
                  
                  # climate                  
                  PRECIP8110WS,TMEAN8110WS, high_flow_count,                                    
                  
                  # vegetation                 
                  mean_ad, mean_cc, frequency_km,                                    
                  
                  # roads                  
                  RDCRSWS,RDDENSWS)) %>% 
  drop_na() %>% 
  filter(SITE_ID != "CAEMI0015") %>% #remove catchment area outlier
  dplyr::select(-c(AQUA_PROV)) %>% 
  filter(Confinement != "None")


### Test on all Cascades post-fire sites 
# Run function test_models to test unburned model on postfire sites
results.pf.c <- test_models(model_list.c, D50_phi_test.pf.c, 'D50_phi')
#compile summary stats for model testing
summary_stats.pf.c <- results.pf.c$summary_stats
#compile dataframe of all predicted and observed values
predictions.pf.c <- do.call(rbind, results.pf.c$predictions)
#write_csv(summary_stats.pf.c, "RF_models/FINAL_RF_MODELS/D50/post_fire_testing/cascades_D50_phi_post_fire_testing.csv")
#write_csv(predictions.pf.c, "RF_models/FINAL_RF_MODELS/D50/post_fire_testing/predicted_observed/cascades_D50_phi_post_fire_pred_obs.csv")

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#######
#### PARTIAL DEPENDENCE PLOTS ####
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#######

## CASCADES ####

# Create a directory to save the plots
output_dir <- "FIGURES/supplement/PDP/unb_casc_D50"
dir.create(output_dir, showWarnings = FALSE)

# Define top ten target variables
target_variables <- c("bankfull_mean", "Confinement", "stream_slope", 
                      "slope23d_perc_greater", "COMPSTRGTHWS","mean_cc", "area_sqkm.x", 
                      "PRECIP8110WS", "survey_year.x", "OMWS")

# Loop over variables
for (var in target_variables) {
  
  pdp_list <- list()
  
  for (i in seq_along(model_list.c)) {
    rf_model <- ranger(D50_phi ~ ., data = D50_phi_data.c,
                       num.trees = 875, mtry = 7, min.node.size = 2,
                       importance = 'permutation', scale.permutation.importance = TRUE)
    
    pd <- partial(rf_model, pred.var = var,
                  train = D50_phi_data.c, grid.resolution = 20, progress = "none") %>%
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
output_dir.k <- "FIGURES/supplement/PDP/unb_klam_D50"
dir.create(output_dir.k, showWarnings = FALSE)

# Define your ten target variables
target_variables.k <- c("bankfull_mean", "area_sqkm.x","PRECIP8110WS","SITE_ID",
                        "CREEK_CODE", "BFIWS", "HYDRLCONDWS", 
                        "RDDENSWS", "RDCRSWS", "OMWS")

# Loop over variables
for (var in target_variables.k) {
  
  pdp_list <- list()
  
  for (i in seq_along(model_list.k)) {
    rf_model <- ranger(D50_phi ~ ., data = D50_phi_data.k,
                       num.trees = 875, mtry = 7, min.node.size = 2,
                       importance = 'permutation', scale.permutation.importance = TRUE)
    
    pd <- partial(rf_model, pred.var = var,
                  train = D50_phi_data.k, grid.resolution = 20, progress = "none") %>%
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
