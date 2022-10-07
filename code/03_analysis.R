######################################################################
## Title: [Insert full title of paper]
## Code author: Steven Kerr steven.kerr@ed.ac.uk
## Description: Prediction of who is undervaccinated,
##              and Cox model of severe outcomes
######################################################################

library(naivebayes)
library(survival)
library(lubridate)
library(broom)

output_dir = "./output/"

#### Functions

# Creates dataframe that is ready for survival model to be fitted.
# df should have a column event_col that is binary and says whether person had an event
# It should also have a column event_col_date that is date of event
create_df_survival = function(df, event_col, calendar_days, study_start, study_end){
  
  event_date_col = paste0(event_col, '_date')
  
  # Add in endpoints
  cohort_endpoints = select(endpoints, EAVE_LINKNO, !!sym(event_col), !!sym(event_date_col)) %>%
    # Only use events that happened in the study period
    filter( !!sym(event_date_col) >= study_start & !!sym(event_date_col) <= study_end) %>%
    # Get date of first event
    group_by(EAVE_LINKNO) %>%
    mutate(!!sym(event_date_col) := min(!!sym(event_date_col), na.rm = TRUE)) %>%
    ungroup() %>%
    distinct()
  
  df = left_join(df, cohort_endpoints) %>%
    mutate(surv_date = pmin(study_end, !!sym(event_date_col), na.rm = TRUE)) %>%
    mutate(surv_date = pmin(surv_date, NRS.Date.Death, na.rm = TRUE)) %>%
    # Censor them if they get a fifth dose
    mutate(surv_date = pmin(surv_date, date_vacc_5, na.rm = TRUE)) %>%
    mutate(event = case_when(surv_date ==  !!sym(event_date_col) ~ 1,
                             TRUE ~ 0)) %>%
    mutate(surv_time = as.numeric(surv_date - (study_start)) ) %>%
    # Some people died before study start and will therefore have negative survival time. Remove them.
    filter(surv_time >= 0) %>%
    # If they had an event on day of vaccination, assume that the event happened shortly after vaccination
    # The exact amount of time you add on doesn't matter when fitting Cox model, only the ordering.
    mutate(surv_time = ifelse(surv_time == 0, 0.001, surv_time)) 
  
  df = tmerge(df, df, id = EAVE_LINKNO, event = event(surv_time, event)) 
  
  # dataframe of start times for different vaccination status
  df_vs = select(df, EAVE_LINKNO, date_vacc_1, date_vacc_2, date_vacc_3, date_vacc_4) %>%
    mutate(uv = -Inf,
           vacc_1 = as.numeric(date_vacc_1 + 14 - study_start),
           vacc_2 = as.numeric(date_vacc_2 + 14 - study_start),
           vacc_3 = as.numeric(date_vacc_3 + 14 - study_start),
           vacc_4 = as.numeric(date_vacc_4 + 14 - study_start)) %>%
    select(-date_vacc_1, -date_vacc_2, -date_vacc_3, -date_vacc_4)
  
  df_vs = pivot_longer(df_vs, cols = c('uv', 'vacc_1', 'vacc_2', 'vacc_3', 'vacc_4'),
                       names_to = 'vs', values_to = 'time') %>%
    filter(!is.na(time))
  
  # dataframe of start times for different calendar periods
  # df_calendar = select(df, EAVE_LINKNO) 
  # 
  # calendar_start = seq(study_start, max(df$surv_date), make_difftime(calendar_days * 24 * 60 * 60, 'days'))
  # 
  # for (i in 1:length(calendar_start)){
  #   df_calendar[as.character(calendar_start[i])] = as.numeric( calendar_start[i] - study_start)
  # }
  # 
  # df_calendar = pivot_longer(df_calendar, cols = as.character(calendar_start),
  #                            names_to = 'calendar', values_to = 'time')
  
  # Add in vaccination status as a time dependent variable
  df = tmerge(df, df_vs, id = EAVE_LINKNO, exposure = tdc(time, vs)) %>%
    rename(vs_time = exposure)
  
  # Add in calendar period as a time dependent variable
  # df = tmerge(df, df_calendar, id = EAVE_LINKNO, exposure = tdc(time, calendar)) %>% 
  #   rename(calendar = exposure)
  
  df = mutate(df, 
              duration = tstop - tstart,
              vaccine_deficit = case_when(age_gp == '5-15' ~ pmax(0, 2 - num_doses),
                                           age_gp == '16-75' ~ pmax(0, 3 - num_doses),
                                           age_gp == '75+' ~ pmax(0, 4 - num_doses) )) %>%
    mutate(vaccine_deficit = factor(vaccine_deficit, levels = c(0, 1, 2, 3, 4)))

}


# Saves table of person years, cumulative incidence curve, fits cox model and 
# saves coefficients and covariance matrix
cox_analysis = function(df_cohort, outcome, age_group, calendar_days, study_start, study_end){
  
  # Hospitalisation or death
  dir = paste0(output_dir, 'cox_', outcome, '_', age_group)
  if (!dir.exists(dir)){dir.create(dir)}
  
  # For testing
  df_cohort = filter(df_cohort, age_group == age_group)
  df_survival = create_df_survival(df_cohort, outcome, calendar_days = calendar_days, study_start = study_start, study_end = study_end)
  
  # Person years by vaccine deficit
  person_years = pyears( Surv(duration, event) ~ vaccine_deficit,
                         data = df_survival, scale=1,  
                         data.frame=TRUE)$data %>%
    # Our unit of time is days, convert to years
    mutate(pyears = pyears/365.25) 
  
  write.csv(person_years, paste0(dir, '/person_years.csv'), row.names = FALSE)
  
  
  # Cumulative incidence curve
  survfit( Surv(tstart, tstop, event) ~ vaccine_deficit,
           # Reset the origin to take account of tim
           data = df_survival %>% 
             group_by(EAVE_LINKNO, vs_time) %>%
             mutate(offset = min(tstart),
                    tstart = tstart - offset,
                    tstop = tstop - offset) %>%
             select(-offset) %>%
             ungroup() ) %>%
    tidy() %>%
    mutate(strata = gsub("vaccine_deficit=0", "0", strata),
           strata = gsub("vaccine_deficit=1", "1", strata),
           strata = gsub("vaccine_deficit=2", "2", strata),
           strata = gsub("vaccine_deficit=3", "3", strata),
           strata = gsub("vaccine_deficit=4", "4", strata)) %>%
    ggplot(aes(x = time, y = 1-estimate, colour = strata, fill = strata)) +
    geom_line() +
    geom_ribbon(aes(ymin = 1 - conf.high, ymax = 1 - conf.low), linetype = 2, alpha = 0.1) +
    xlab('Time (days)') +
    ylab('Cumulative incidence') +
    labs(color = 'Vaccine deficit') +
    labs(fill = 'Vaccine deficit') +
    theme_bw()
  
  ggsave(paste0(dir, '/cumulative_incidence.png'))
  
  
  # Cox model
  cox_model = coxph( Surv(tstart, tstop, event) ~ Sex + age_gp + ur6_2016_name + simd2020_sc_quintile + n_risk_gps  +
                       last_positive_test_group + shielding + num_tests_6m_group + vaccine_deficit,
                     data = df_survival)
  
  # Results table
  write.csv(create_cox_results_table(cox_model), paste0(dir, '/results_table.csv'), row.names = FALSE)
  
  # Coefficients
  cox_coef = data.frame('variable' = names(coef(cox_model)), 'coef' = coef(cox_model), row.names = 1:length(coef(cox_model)))
  write.csv(cox_coef, paste0(dir, '/coef.csv'), row.names = FALSE)
  
  # Covariance matrix
  cox_cov = vcov(cox_model)
  write.csv(cox_cov, paste0(dir, '/cov.csv'), row.names = FALSE)
  
}

# Does what it says
create_logistic_results_table = function(model){
  #get ORs & 95%CIs
  results = round(exp(cbind("OR" = coef(model), confint.default(model, level=0.95), digits=3)), digits = 3) %>%
    as.data.frame() %>%
    select(-digits)
  
  names(results) = c('OR', 'lower', 'upper')
  
  results = mutate(results, variable = rownames(results))
  
  results = mutate( results, OR = paste0(OR, ' (', lower, '-', upper, ')') )  %>%
    select(variable, OR)
  
  names(results) = c('variable', 'OR')
  
  rownames(results) = 1:nrow(results)
  
  results
}


# Creates a summary table of Cox model results
create_cox_results_table = function(model) {
  
  tidy(model) %>%
    filter(!grepl('calendar', term)) %>%
    mutate(ucl = estimate + 1.96 * std.error,
           lcl = estimate - 1.96 * std.error) %>%
    select(term, estimate, lcl, ucl) %>%
    mutate_if(is.numeric, ~formatC(round(exp(.), 2), format = "f", big.mark = ",", drop0trailing = TRUE)) %>%
    mutate(estimate = paste0(estimate, ' (', lcl, ' - ', ucl, ')' )) %>%
    select(-lcl, -ucl) %>%
    rename(HR = estimate)
  
}




#### Undervaccination analyis

## Logistic regression
dir = paste0(output_dir, 'lr_undervaccination')
if (!dir.exists(dir)){dir.create(dir)}

lr_model = glm(fully_vaccinated ~ Sex + age_gp + ur6_2016_name + simd2020_sc_quintile + n_risk_gps,
                family = binomial, data = df_cohort)

# Coefficients
lr_coef = data.frame('variable' = names(coef(lr_model)), 'coef' = coef(lr_model), row.names = 1:length(coef(lr_model)))
write.csv(lr_coef, paste0(dir, '/coef.csv'), row.names = FALSE)

# Covariance matrix
lr_cov = vcov(lr_model)
write.csv(lr_cov, paste0(dir, '/cov.csv'), row.names = FALSE)

# Results table
write.csv(create_logistic_results_table(lr_model),  paste0(dir, '/results.csv'), row.names = FALSE)



## Naive bayes model
dir = paste0(output_dir, 'nb_undervaccination')
if (!dir.exists(dir)){dir.create(dir)}

nb_model = naive_bayes(fully_vaccinated ~ Sex + ur6_2016_name + simd2020_sc_quintile + n_risk_gps + ageYear, 
                       data = df_cohort, usekernel = T) 

vars = c('Sex', 'ur6_2016_name', 'simd2020_sc_quintile', 'n_risk_gps', 'ageYear')

for (var in vars){

  png(paste0(dir, '/', var, '.png'), width = 350, height = 350)
  
  plot(nb_model, which = var)
  
  dev.off()
}




## Comparison of logistic regression and naive bayes
# Add predicted probabilities and classes
df_cohort = mutate(df_cohort,
                   lr_prob = predict(lr_model, df_cohort, type = "response"),
                   nb_prob = predict(nb_model, type = 'prob')[, 2]) %>%
  mutate(lr_class = ifelse(lr_prob  > 0.5, 1, 0),
         nb_class = ifelse(nb_prob > 0.5, 1, 0))


# Confusion matrices
lr_confusion = table(df_cohort$lr_class, df_cohort$fully_vaccinated)
nb_confusion = table(df_cohort$nb_class, df_cohort$fully_vaccinated)

# MSE for class predictions
mean( ( as.numeric(df_cohort$lr_prob) - (as.numeric(df_cohort$fully_vaccinated) -1) )^2, na.rm = TRUE)
mean( ( as.numeric(df_cohort$nb_prob) - (as.numeric(df_cohort$fully_vaccinated) -1))^2, na.rm = TRUE)







#### Severe outcomes analysis

# End dates
study_end_hosp = max(  max(smr$ADMISSION_DATE, na.rm = TRUE),
                       max(smr$DISCHARGE_DATE, na.rm = TRUE), na.rm = TRUE)
study_end_hosp = min(study_end, study_end_hosp)

study_end_death = max(all_deaths$NRS.Date.Death, na.rm = TRUE)
study_end_death = min(study_end, study_end_death)

study_end_hosp_death = min(study_end_hosp, study_end_death, na.rm = TRUE)
study_end_hosp_death = min(study_end, study_end_hosp_death)



cox_analysis(df_cohort, age_group = '5-15', outcome = 'covid_hosp_death',
             calendar_days = NA, study_start = study_start, study_end = study_end_hosp_death)

cox_analysis(df_cohort, age_group = '16-75', outcome = 'covid_hosp_death',
             calendar_days = NA, study_start = study_start, study_end = study_end_hosp_death)

cox_analysis(df_cohort, age_group = '75+', outcome = 'covid_hosp_death',
             calendar_days = NA, study_start = study_start, study_end = study_end_hosp_death)





########### Model experimentation. This is not used in the final results!

# # For testing
# #df_survival = create_df_survival( df = sample_n(df_cohort, 1000), event_col = 'covid_hosp_death', calendar_days = 14, study_start = study_start, study_end = study_end_hosp_death)
# df_survival = create_df_survival( df_cohort, 'covid_hosp_death', calendar_days = 28, study_start = study_start, study_end = study_end_hosp_death)
# 
# # Cox model
# # This worked
# cox_model = coxph( Surv(tstart, tstop, event) ~ Sex + age_gp + ur6_2016_name + simd2020_sc_quintile + n_risk_gps  +
#                  last_positive_test_group + shielding + num_tests_6m_group + vaccine_deficit,
#                   data = df_survival)
# 
# summary(cox_model)
# 
# # This works, but the calendar terms are all NA
# cox_model = coxph( Surv(tstart, tstop, event) ~ Sex + age_gp + ur6_2016_name + simd2020_sc_quintile + n_risk_gps  +
#                      last_positive_test_group + shielding + num_tests_6m_group + vaccine_deficit + calendar,
#                    data = df_survival)
# 
# summary(cox_model)
# 
# # This does not work. The addition of covid_hosp_ever seems to ruin it
# cox_model = coxph( Surv(tstart, tstop, event) ~ Sex + age_gp + ur6_2016_name + simd2020_sc_quintile + n_risk_gps  +
#                      last_positive_test_group + covid_hosp_ever + shielding + num_tests_6m_group + vaccine_deficit + calendar,
#                    data = df_survival)
# 
# summary(cox_model)
# 
# # This also didnt work. Calendar has been removed, and covid_hosp_ever included
# cox_model = coxph( Surv(tstart, tstop, event) ~ Sex + age_gp + ur6_2016_name + simd2020_sc_quintile + n_risk_gps  +
#                      last_positive_test_group + covid_hosp_ever + shielding + num_tests_6m_group + vaccine_deficit,
#                    data = df_survival)
# 
# summary(cox_model)
# 



















