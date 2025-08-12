# AREMP_wildfire_substrate
Read Me File for Wall et al., “Effects of wildfire on streambed sediment in the Cascades and Klamath provinces of the Pacific Northwest”
08/06/2025
CSV Files:
**Keys for all abbreviations and units for all variables can be found in Table 1 below. 
All_site_survey_covars.csv
•	Dataframe of all response and explanatory variables for all AREMP subwatersheds
o	Contains more sites and covariates than used in analysis. Dataframe is filtered and cleaned in scripts. 
Post_fire_surveys_covar.csv
•	Dataframe of all response and explanatory variables for all AREMP sites that have had fire in their subwatersheds. 
•	Includes all covariates in “all_site_survey_covars.csv” and includes wildfire metrics.
Site_id_linkage.csv
•	Site id key to connect specific survey ids with specific sites.
R Scripts
01a_Percent_Fines_Unburned_Modeling
•	Constructs Cascades and Klamath unburned percent fines models and tests the unburned models on the post-fire datasets
•	Creates variable importance plots for the unburned percent fines models
***Output data frames from this script are used in “03a_Percent_Fines_Difference_Modeling”
01b_D50_Unburned_Modeling
•	Constructs Cascades and Klamath unburned D50 models and tests the unburned models on the post-fire datasets
•	Creates variable importance plots and partial dependence plots for the unburned D50 models
***Output data frames from this script are used in “03b_D50_Difference_Modeling”
02a_Percent_Fines_Postfire_Modeling
•	Constructs Cascades and Klamath percent fines models using the post-fire dataset
•	Creates variable importance plots and partial dependence plots for the post-fire percent fines models
02b_D50_Postfire_Modeling
•	Constructs Cascades and Klamath D50 models using the post-fire dataset
•	Creates variable importance plots and partial dependence plots for the post-fire D50 models
03a_Percent_Fines_Difference_Modeling
•	Uses dataframes produced in step 01a
•	Constructs Cascades and Klamath Percent Fines Difference Models
•	Creates Variable importance and partial dependence plots
03b_D50_Difference_Modeling
•	Uses dataframes produced in step 01b
•	Constructs Cascades and Klamath D50 Difference Models
•	Creates Variable importance and partial dependence plots
Model Tuning Folder:
Within the folder there are two subfolders: One for the Percent Fines Models and one for the D50 models. These scripts are for tuning each model in the paper to find the best model parameters.
Each script within each folder is for each individual model: 
Percent_fines_model_tuning Folder:
Model:	Script name
Cascades Unburned	Cascades_perc_fines_model_tuning
Klamath Unburned	Klamath_perc_fines_model_tuning
Cascades Post-fire	Cascades_perc_fines_postfire_model_tuning
Klamath Post-fire	Klamath_postfire_perc_fines_model_tuning
Cascades Difference	Cascades_difference_model_tuning
Klamath Difference	Klamath_difference_model_tuning

D50_model_tuning Folder:
Model:	Script name
Cascades Unburned	Cascades_d50_model_tuning
Klamath Unburned	Klamath_d50_model_tuning
Cascades Post-fire	Cascades_D50_postfire_model_tuning
Klamath Post-fire	Klamath_d50_postfire_model_tuning
Cascades Difference	Cascades_D50_difference_model_tuning
Klamath Difference	Klamath_D50_difference_model_tuning



Table 1. Covariates Summary Table
Variables	Abbreviation	Unit	Spatial Extent	Category
SITE ID	SITE_ID	Categorical	Reach	Location
Watershed ID	CREEK_CODE	Categorical	HUC12 Watershed	Location
Survey Year	Survey_year	years	Reach	Temporal
Mean Bankfull Width	Bankfull_mean	m	Reach	Stream Morphology
Relative bankfull	Relative_bankfull	-	Reach	Temporal
Streambed Slope	Stream_slope	degrees	Reach	Stream Morphology
Valley Confinement	Confinement	Categorical	Reach	Stream Morphology
In-Stream Large Wood	Frequency_km	Units/km	Reach	Stream Morphology
High Flow Count	High_flow_count	-	Nearest USGS Gage	Climate
Subwatershed Area	Area_sqkm.x	Km2	Subwatershed	Topography
% of subwatershed with slopes 
> 23 degrees	Slope23d_perc_greater	%	Subwatershed	Topography
Mean Elevation	Elev_mean	m	Subwatershed	Topography
Baseflow Index	BFIWS	%	Subwatershed	Climate
30-year Normal Precipitation	PRECIP8110WS	mm	Subwatershed	Climate
30-year Mean Temperature	TMEAN8110WS	celsius	Subwatershed	Climate
Mean Rock Compressive Strength	COMPSTRGTHWS	megaPascals	Subwatershed	Lithology/Soils
Mean Lithological Hydraulic Conductivity	HYDRLCONDWS	micrometers/s	Subwatershed	Lithology/Soils
Dominant Geology	Dom_geology	Categorical	Subwatershed	Lithology/Soils
Mean Soil Erodibility Factor	KFFACTWS	-	Subwatershed	Lithology/Soils
Soil Permeability	PERMWS	cm/hr	Subwatershed	Lithology/Soils
Mean Soil Depth to Bedrock	RCKDEPWS	cm	Subwatershed	Lithology/Soils
Soil Organic Matter Content	OMWS	%	Subwatershed	Lithology/Soils
Mean Tree Age	Mean_ad	years	Subwatershed	Land Cover
Mean Canopy Cover	Mean_cc	%	Subwatershed	Land Cover
Road Crossings	RDCRSWS	crossings/sqkm	Subwatershed	Land Cover
Upslope Road Density	RDDENSWS	km/sqkm	Subwatershed	Land Cover
Years Between Survey and Fire	YPF	Years	Subwatershed	Fire
Area of Subwatershed  Burned	Ppc_burned_sqkm	Sqkm	Subwatershed	Fire
Area of Subwatershed Burned at Himod	Ppc_himod_sqkm	Sqkm	Subwatershed	Fire
Area of Riparian Area Burned	Rip_burned_sqkm	Sqkm	Riparian	Fire
Area of Riparian Area Burned at Himod	Rip_himod_sqkm	Sqkm	Riparian	Fire
Area of Reach Burned	Reach_burned_sqkm	Sqkm	Reach	Fire
Area of Reach Burned at Himod	Reach_himod_sqkm	Sqkm	Reach	Fire
Fire Rank	Rank.x	Count	Subwatershed	Fire


