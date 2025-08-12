## RF MODEL TUNING: D50 Difference####

#### CASCADES ####

# 0. Load  and install packages ----
# List of all packages needed
package_list <- c('tidyverse', 'tidymodels', 'ranger', 'lubridate', 'vip', 'ggpubr', 'DALEXtra','doParallel' , 'ggtext', 'patchwork')

# Check if there are any packacges missing
packages_missing <- setdiff(package_list, rownames(installed.packages()))

# If we find a package missing, install them
if(length(packages_missing) >= 1) install.packages(packages_missing) 

# Now load all the packages
lapply(package_list, require, character.only = TRUE)


# Load Data Sets ####

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

# Cascades pred_obs data for testing on all post-fire sites
pred_obs.pf.c <- read_csv("RF_models/FINAL_RF_MODELS/D50/post_fire_testing/predicted_observed/cascades_D50_phi_post_fire_pred_obs.csv")


## Create Modeling Dataset for the Cascades ####

pf_covars_filtered.C <- pf_survey_covars %>% 
  group_by(site_survey_id) %>% 
  arrange(site_survey_id, fire_year) %>% 
  mutate(fire_rank = row_number()) %>% 
  ungroup() %>% 
  relocate(fire_rank, .after= site_survey_id) %>% 
  filter(fire_rank == "1") %>% 
  dplyr::select(c(ppc_burned_sqkm,
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
  filter(AQUA_PROV == "Western Cascades"|
           AQUA_PROV == "High Cascades" |
           AQUA_PROV == "North Cascades") %>% 
  filter(SITE_ID != "CAEMI0015") %>% #remove catchment area outlier
  dplyr::select(-c(AQUA_PROV)) %>% 
  filter(Confinement != "None")

difference.fines.c <- pred_obs.pf.c %>% 
  mutate(difference = observed - predicted) %>% 
  left_join(site_id_link, by = c("site_id" = "site_id", 
                                 "year" = "survey_year")) %>% 
  select(-c(X, OBJECTID, QC_visit)) %>% 
  left_join(pf_covars_filtered.C, by = c("site_survey_id" = "site_survey_id"),
            relationship = "many-to-many")

difference_train.fines.c<- difference.fines.c %>% 
  dplyr::select(c(Model, difference, ppc_burned_sqkm,
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
                  RDCRSWS,RDDENSWS, site_survey_id)) %>% 
  drop_na() %>% 
  filter(Model == "10") %>% 
  mutate(CREEK_CODE = as.factor(CREEK_CODE),
         SITE_ID = as.factor(SITE_ID),
         dom_geology = as.factor(dom_geology)) %>% 
  select(-c(Model))



set.seed(123)

#prep the dataset into a training and testing one
trees_split <- initial_split(difference_train.fines.c)
trees_train <- training(trees_split)
trees_test <- testing(trees_split)

#Model recipe, predict percent fines
tree_rec <- recipe(difference ~ ., data = trees_train) %>%
  update_role(site_survey_id, new_role = "ID") 

#prepare the testing df by passing on the recipe
tree_prep <- prep(tree_rec)
#prepare the training dataset
juiced <- juice(tree_prep)

# we will tune the hyperparameters of the RF, so we prepare the model for that and put it in a workflow
tune_spec <- rand_forest(
  mtry = tune(),
  trees = tune(),
  min_n = tune()) %>%
  set_mode("regression") %>%
  set_engine("ranger")

tune_wf <- workflow() %>%
  add_recipe(tree_rec) %>%
  add_model(tune_spec)

#prepare the folds for the tuning
set.seed(234)
trees_folds <- vfold_cv(trees_train, v = 5)

# Now run a bunch of models in parallel with different combinations of parameters to find what works best
doParallel::registerDoParallel()

set.seed(345)
tune_res <- tune_grid(
  tune_wf,
  resamples = trees_folds,
  grid = 20
)

tune_res

# lets have a look...
tune_res %>%
  collect_metrics() %>%
  filter(.metric == "rmse") %>%
  select(mean, min_n, mtry, trees) %>%
  pivot_longer(min_n:trees,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") 

#Do a more targeted tuning with a new grid
rf_grid <- grid_regular(
  mtry(range = c(6, 12)),
  min_n(range = c(2, 10)),
  trees(range = c(500, 1700)),
  levels = 5
)
# check grid
rf_grid

#tune the grid again
set.seed(456)
regular_res <- tune_grid(
  tune_wf,
  resamples = trees_folds,
  grid = rf_grid
)

regular_res

# Let's check the results
regular_res %>%
  collect_metrics() %>%
  filter(.metric == "rmse") %>%
  select(mean, min_n, mtry, trees) %>%
  pivot_longer(min_n:trees,
               values_to = "value",
               names_to = "parameter") %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "rmse")

#select the best one and finish the model 
best_model <- select_best(regular_res, metric = "rmse")

final_rf <- finalize_model(
  tune_spec,
  best_model
)

# look at final model parameters: In this case, the best parameters are:
# mtry = 9, trees = 500, min_n = 6
final_rf

# run final model on training and testing datasets
final_wf <- workflow() %>%
  add_recipe(tree_rec) %>%
  add_model(final_rf)

final_res <- final_wf %>%
  last_fit(trees_split)

# see final model stats
final_res %>%
  collect_metrics()


