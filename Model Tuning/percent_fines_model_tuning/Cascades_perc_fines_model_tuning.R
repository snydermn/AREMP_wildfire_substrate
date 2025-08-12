## RF MODEL TUNING: PERCENT FINES ####

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

# all unburned and pre-fire surveys
site_survey_covars <- read.csv("all_site_survey_covars.csv") %>% 
  anti_join(pf_survey_covars, by = c("site_survey_id" = "site_survey_id")) %>% 
  left_join(bankfull, by = c("SITE_ID" = "SITE_ID")) %>% 
  mutate(relative_bankfull = bankfull_mean / bankfull_all_mean) %>% 
  filter(QC_visit.x == "N") 

# site id linkages that will be used for the difference models
site_id_link <- read.csv("aremp_site_info/site_id_linkage.csv") %>% 
  filter(QC_visit == "N")


## Create Modeling Dataset for the Cascades ####


## Cascades UNBURNED TRAINING DATA ##
## CASCADES UNBURNED TRAINING DATA ##
perc_fines_data.c=subset(site_survey_covars,select=c(perc_fines, CREEK_CODE,
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
                                                     RDCRSWS,RDDENSWS, site_survey_id)) %>% 
  drop_na() %>% 
  filter(AQUA_PROV == "Western Cascades"|
           AQUA_PROV == "High Cascades" |
           AQUA_PROV == "North Cascades") %>% 
  filter(SITE_ID != "CAEMI0015") %>% #remove catchment area outlier %>% 
  mutate(SITE_ID = as.factor(SITE_ID),
         CREEK_CODE = as.factor(CREEK_CODE)) %>% 
  dplyr::select(-c(AQUA_PROV)) %>% 
  filter(Confinement != "None")



set.seed(123)

#prep the dataset into a training and testing one
trees_split <- initial_split(perc_fines_data.c)
trees_train <- training(trees_split)
trees_test <- testing(trees_split)

#Model recipe, predict percent fines
tree_rec <- recipe(perc_fines ~ ., data = trees_train) %>%
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
  mtry(range = c(5, 10)),
  min_n(range = c(2, 10)),
  trees(range = c(500, 1000)),
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
# mtry = 7, trees = 875, min_n = 2
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


