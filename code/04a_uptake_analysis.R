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

# Add more entries to names_map to automatically generate display names
names_map["last_positive_test_group"] = "Last positive test"
names_map["num_tests_6m_group"] = "Number of tests in last 6 months"
names_map["shielding"] = "Ever shielding"
names_map["num_doses_time_2"] = "Number of doses"
names_map["age_17cat"] = "Age group"
names_map["covid_hosp_ever"] = "Previous COVID-19 hospitalisation"
names_map["urban_rural_2cat"] = "Urban/rural classification"

#### Functions

# Fit logistic model and saves results table, coefficients and covariance matrix
lr_analysis = function(df_cohort, age_group, dep_var, ind_vars) {

  ## Logistic regression
  dir = paste0(output_dir, "lr_undervaccination_", age_group)
  if (!dir.exists(dir)) {
    dir.create(dir)
  }

  df_cohort = filter(
    df_cohort, 
      age_4cat == age_group,
      is.na(death_date) | death_date > study_start) %>%
    droplevels()

  # Weighted logistic regression
  formula = as.formula(paste(c(paste0(dep_var, " ~ "), ind_vars), collapse = " + "))
  
  survey_design = svydesign(id = ~1,
                            weights = ~ eave_weight,
                            data = df_cohort)
  
  model = svyglm(formula, design = survey_design, family = binomial, data = df_cohort)
  
  # This takes up too much space
  #saveRDS(model, paste0(dir, "/model.rds"))

  # Coefficients
  coef = data.frame("variable" = names(coef(model)), "coef" = coef(model), row.names = 1:length(coef(model)))
  write.csv(coef, paste0(dir, "/coef.csv"), row.names = FALSE)

  # Covariance matrix
  cov = vcov(model)
  write.csv(cov, paste0(dir, "/cov.csv"))

  # Results table
  write.csv(create_logistic_results_table(df_cohort, model, dep_var, ind_vars), paste0(dir, "/results.csv"), row.names = FALSE)
}


# Does what it says
create_logistic_results_table = function(df_cohort, model, dep_var, ind_vars) {
  df_vars = data.frame()

  for (var in ind_vars) {
    
    Levels = as.character(levels( as.factor(pull(df_cohort, !!sym(var)))))
    
    #Levels = unique(pull(df_cohort, !!sym(var)))
    
    new_df = data.frame(Variable = var, Levels) %>%
      mutate(
        OR = case_when( Levels == as.character(levels( as.factor(pull(df_cohort, !!sym(var))))[1])  ~ 'Ref',
          TRUE ~ NA_character_) ) 
 
    new_df = df_cohort %>%
      group_by(!!sym(var)) %>%
      summarise(Uptake = mean(as.numeric(fully_vaccinated) - 1)) %>%
      mutate(Uptake = round(100 * Uptake, 1) ) %>%
      mutate(Uptake = paste0(Uptake, '%')) %>%
      rename(Levels = !!sym(var)) %>%
      right_join(new_df, by = 'Levels')

    df_vars = bind_rows(df_vars, new_df)
    
    # df_vars = mutate(df_vars,
    #   OR = case_when( Levels == as.character(levels( as.factor(pull(df_cohort, !!sym(var))))[1])  ~ 'Ref',
    #                 TRUE ~ NA_character_) )
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
      #OR = replace_na(OR, "Ref"),
      Variable = names_map[Variable] ,
      Variable = case_when(
        !duplicated(Variable) ~ Variable,
        TRUE ~ ""
      )
    ) %>%
  select(-term, -OR.x, -OR.y) %>%
  select(Variable, Levels, Uptake, OR) %>%
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
lr_dep_var = "fully_vaccinated"
lr_ind_vars = c("sex", "age_17cat", "urban_rural_2cat", "simd2020_sc_quintile", "n_risk_gps")

# age_17cat only has one category for 5-11 year olds and 12-15 year olds
# - drop it as a predictor
lr_analysis(df_cohort, "5-11", lr_dep_var, setdiff(lr_ind_vars, 'age_17cat'))
lr_analysis(df_cohort, "12-15", lr_dep_var, setdiff(lr_ind_vars, 'age_17cat'))
lr_analysis(df_cohort, "16-74", lr_dep_var, lr_ind_vars)
lr_analysis(df_cohort, "75+", lr_dep_var, lr_ind_vars)


## Naive bayes model
# ENGLAND NORTHERN IRELAND AND WALES DON'T HAVE TO DO THE NAIVE BAYES ESIMATION
nb_dep_var = "fully_vaccinated"
nb_ind_vars = c("sex", "urban_rural_2cat", "simd2020_sc_quintile", "n_risk_gps", "age")

nb_analysis(df_cohort, "5-11", nb_dep_var, nb_ind_vars)
nb_analysis(df_cohort, "12-15", nb_dep_var, nb_ind_vars)
nb_analysis(df_cohort, "16-74", nb_dep_var, nb_ind_vars)
nb_analysis(df_cohort, "75+", nb_dep_var, nb_ind_vars)


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
# lr_confusion = table(df_cohort$lr_class, df_cohort$fully_vaccinated)
# nb_confusion = table(df_cohort$nb_class, df_cohort$fully_vaccinated)
#
# # MSE for class predictions
# mean( ( as.numeric(df_cohort$lr_prob) - (as.numeric(df_cohort$fully_vaccinated) -1) )^2, na.rm = TRUE)
# mean( ( as.numeric(df_cohort$nb_prob) - (as.numeric(df_cohort$fully_vaccinated) -1))^2, na.rm = TRUE)

