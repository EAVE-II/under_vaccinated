######################################################################
## Title: [Insert full title of paper]
## Code author: Steven Kerr steven.kerr@ed.ac.uk
## Description: Analysis of suboptimal vaccination
##              01_data_cleaning and 02_data_prep need to be run first.
##              names_map from  03_cohort_description is also needed.
######################################################################

library(naivebayes)
library(survival)
library(lubridate)
library(broom)
library(survey)

output_dir = "./output/"

# Pre-requisites
#source("./code/01_data_cleaning.R")
#source("./code/02_data_prep.R")
# The only thing needed from 02_cohort_description is names_map
#source('./output/02_cohort_description.R')

# Read in
#df_cohort = readRDS('./data/df_cohort.rds')

# Add more entries to names_map to automatically generate display names
names_map["last_positive_test_group"] = "Last positive test"
names_map["num_tests_6m_group"] = "Number of tests in last 6 months"
names_map["shielding"] = "Ever shielding"
names_map["other_household_shielding"] = "Other household member ever shielding"
names_map["num_doses_time_2"] = "Number of doses"
names_map["num_doses_time_3"] = "Number of doses"
names_map["under_vaccinated_time"] = "Vaccination status"
names_map["age_17cat"] = "Age group"
names_map["covid_mcoa_hosp_ever"] = "Previous COVID-19 hospitalisation"
names_map["urban_rural_2cat"] = "Urban/rural classification"
names_map['competing_hosp'] = "Competing non-COVID-19 hospitalisation"
names_map["num_pos_tests_6m_group"] = "Number of positive tests in last 6 months"
names_map["Q_BMI_imputed"] = "BMI"
names_map["bmi_imputed_cat"] = "BMI"
names_map["health_board"] = "Health board"


#### Functions

# Fit logistic model and saves results table, coefficients and covariance matrix
lr_analysis = function(df_cohort, age_group, dep_var, ind_vars, exclude_non_unit_weight, name) {
  
  ## Logistic regression
  dir = paste0(output_dir, "lr_undervaccination_", age_group, '_', name)
  
  if (exclude_non_unit_weight){
    dir = paste0(dir, "_weight1")
    
    df_cohort = filter(df_cohort, eave_weight == 1)
  } 
  
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  
  df_cohort = filter(
    df_cohort,
    # flag_incon == 0,
    age_4cat == age_group,
    is.na(death_date) | death_date > study_start) %>%
    droplevels()
  
  formula = as.formula(paste(c(paste0(dep_var, " ~ "), ind_vars), collapse = " + "))
  
  # Weighted logistic regression
  # 
  # survey_design = svydesign(id = ~1,
  #                           weights = ~ eave_weight,
  #                           data = df_cohort)
  # 
  # model = svyglm(formula, design = survey_design, data = df_cohort, family = binomial
  
  model <- glm(formula, data = df_cohort, family = "binomial")
  
  # Coefficients
  coef = data.frame("variable" = names(coef(model)), "coef" = coef(model), row.names = 1:length(coef(model)))
  write.csv(coef, paste0(dir, "/coef.csv"), row.names = FALSE)
  
  # Covariance matrix
  cov = vcov(model)
  write.csv(cov, paste0(dir, "/cov.csv"))
  
  # Results table
  write.csv(create_logistic_results_table(df_cohort, model, dep_var, ind_vars), paste0(dir, "/results.csv"), row.names = FALSE)
  
  # Model performance metrics
  write.csv(broom::glance(model), paste0(dir, "/perf_metrics.csv"))
  
  # Variance inflation factor, for testing collinearity
  write.csv(vif(model), paste0(dir, "/vif.csv"))
}


# Does what it says
create_logistic_results_table = function(df_cohort, model, dep_var, ind_vars) {
  df_vars = data.frame()
  
  for (var in ind_vars) {
    
    Levels = as.character(levels( as.factor(pull(df_cohort, !!sym(var)))))
    
    new_df = data.frame(Variable = var, Levels) %>%
      mutate(
        OR = case_when( Levels == as.character(levels( as.factor(pull(df_cohort, !!sym(var))))[1])  ~ 'Ref',
                        TRUE ~ NA_character_) ) 
    
    new_df = df_cohort %>%
      group_by(!!sym(var)) %>%
      # Level 0 gets mapped to 1, and 1 gets mapped to 2 by as.numeric
      summarise(`Under-vaccinated` = mean(as.numeric(under_vaccinated) -1)) %>%
      mutate(`Under-vaccinated` = round(100 * `Under-vaccinated`, 1) ) %>%
      mutate(`Under-vaccinated` = paste0(`Under-vaccinated`, '%')) %>%
      rename(Levels = !!sym(var)) %>%
      right_join(new_df, by = 'Levels')
    
    df_vars = bind_rows(df_vars, new_df)
  }
  
  df_vars = mutate(
    df_vars, 
    term = paste0(Variable, Levels),
    Variable = as.character(Variable)) %>%
    filter(!is.na(Levels))
  
  # Get ORs & 95% CIs
  results = round(exp(cbind("OR" = coef(model), confint.default(model, level = 0.95), digits = 2)), digits = 2) %>%
    as.data.frame() %>%
    select(-digits)
  
  names(results) = c("OR", "lower", "upper")
  
  results = mutate(results, term = rownames(results))
  
  results = mutate(results, OR = paste0(OR, " (", lower, ", ", upper, ")")) %>%
    select(term, OR)
  
  results = left_join(df_vars, results, by = 'term') %>%
    mutate(
      OR = coalesce(OR.x, OR.y),
      Variable = names_map[Variable] ,
      Variable = case_when(
        !duplicated(Variable) ~ Variable,
        TRUE ~ ""
      )
    ) %>%
    select(-term, -OR.x, -OR.y) %>%
    select(Variable, Levels, `Under-vaccinated`, OR) %>%
    rename("Odds ratio" = OR)
}


# Fits naive bayes model, saves graphs of fitted marginal distributions
# ENGLAND NORTHERN IRELAND AND WALES DON'T HAVE TO DO THE NAIVE BAYES ESIMATION
nb_analysis = function(df_cohort, age_group, dep_var, ind_vars) {
  dir = paste0(output_dir, "nb_undervaccination_", age_group)
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  
  df_cohort = filter(df_cohort, age_4cat == age_group) %>%
    droplevels()
  
  formula = as.formula(paste(c(paste0(dep_var, " ~ "), ind_vars), collapse = " + "))
  nb_model = naive_bayes(formula, data = df_cohort, usekernel = T)
  
  for (var in ind_vars) {
    png(paste0(dir, "/", var, ".png"), width = 350, height = 350)
    
    plot(nb_model, which = var)
    
    dev.off()
  }
}


#### Undervaccination analysis

## Logistic regression
lr_dep_var = "under_vaccinated"
lr_ind_vars_minimal = c(
  "sex", 
  "age_17cat", 
  "ethnicity_5cat", 
  "urban_rural_2cat", 
  "simd2020_sc_quintile", 
  "n_risk_gps_6cat",
  "n_risk_gps_3cat")

lr_ind_vars_extended = c(
  "sex", 
  "age_17cat", 
  "ethnicity_5cat",
  "num_ppl_household",
  "health_board", 
  "urban_rural_2cat", 
  "simd2020_sc_quintile",
  #"bmi_imputed_cat",
  "n_risk_gps_6cat",
  "n_risk_gps_3cat",
  "last_positive_test_group")

var_list = list(
  "minimal" =  lr_ind_vars_minimal,
  "extended" = lr_ind_vars_extended)


for (age_group in c("5-11", "12-15", "16-74", "75+")) {
  print(age_group)
  
  for (var_set in 1:length(var_list)){
    
    ind_vars = var_list[[var_set]]
    name = names(var_list)[var_set]
    
    print(name)
    
    if (age_group %in%  c('5-11', '12-15')){
      ind_vars = setdiff(ind_vars, c('age_17cat', 'n_risk_gps_6cat', 'bmi_imputed_cat'))
    } else {
      ind_vars = setdiff(ind_vars, 'n_risk_gps_3cat')
    }
    
    lr_analysis(
      df_cohort, 
      age_group = age_group, 
      dep_var = lr_dep_var, 
      ind_vars = ind_vars,
      exclude_non_unit_weight = TRUE,
      name = name
    )
  }
}


# age_17cat only has one category for 5-11 year olds and 12-15 year olds
# - drop it as a predictor
# lr_analysis(
#   df_cohort, 
#   age_group = "5-11", 
#   dep_var = lr_dep_var, 
#   ind_vars = setdiff(lr_ind_vars, c('age_17cat', "n_risk_gps_6cat")),
#   exclude_non_unit_weight = TRUE
#   )
# 
# lr_analysis(
#   df_cohort, 
#   age_group = "12-15", 
#   dep_var = lr_dep_var, 
#   ind_vars = setdiff(lr_ind_vars, c('age_17cat', "n_risk_gps_6cat")),
#   exclude_non_unit_weight = TRUE
# )
# 
# lr_analysis(
#   df_cohort, 
#   age_group = "16-74", 
#   dep_var = lr_dep_var, 
#   ind_vars = setdiff(lr_ind_vars, "n_risk_gps_3cat"),
#   exclude_non_unit_weight = TRUE
# )
# 
# lr_analysis(
#   df_cohort, 
#   age_group = "75+", 
#   dep_var = lr_dep_var, 
#   ind_vars = setdiff(lr_ind_vars, "n_risk_gps_3cat"),
#   exclude_non_unit_weight = TRUE
# )



## Naive bayes model
# ENGLAND NORTHERN IRELAND AND WALES DON'T HAVE TO DO THE NAIVE BAYES ESIMATION
# nb_dep_var = "under_vaccinated"
# nb_ind_vars = c("sex", "ethnicity_5cat", "urban_rural_2cat", "simd2020_sc_quintile", "n_risk_gps", "age")
# 
# nb_analysis(df_cohort, "5-11", nb_dep_var, nb_ind_vars)
# nb_analysis(df_cohort, "12-15", nb_dep_var, nb_ind_vars)
# nb_analysis(df_cohort, "16-74", nb_dep_var, nb_ind_vars)
# nb_analysis(df_cohort, "75+", nb_dep_var, nb_ind_vars)


## Comparison of logistic regression and naive bayes
# Add predicted probabilities and classes
# df_cohort = mutate(df_cohort,
#                    lr_prob = predict(lr_model, df_cohort, type = "response"),
#                    nb_prob = predict(nb_model, type = 'prob')[, 2]) %>%
#   mutate(lr_class = ifelse(lr_prob  > 0.5, 1, 0),
#          nb_class = ifelse(nb_prob > 0.5, 1, 0))
#
#
# # Confusion matrices
# lr_confusion = table(df_cohort$lr_class, df_cohort$under_vaccinated)
# nb_confusion = table(df_cohort$nb_class, df_cohort$under_vaccinated)
#
# # MSE for class predictions
# mean( ( as.numeric(df_cohort$lr_prob) - (as.numeric(df_cohort$under_vaccinated) -1) )^2, na.rm = TRUE)
# mean( ( as.numeric(df_cohort$nb_prob) - (as.numeric(df_cohort$under_vaccinated) -1))^2, na.rm = TRUE)

