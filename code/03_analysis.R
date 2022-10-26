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

# Add more entries to names_map to automatically generate display names
names_map["last_positive_test_group"] = "Last positive test"
names_map["num_tests_6m_group"] = "Number of tests in last 6 months"
names_map["shielding"] = "Ever shielding"
names_map["num_doses_time"] = "Number of doses"
names_map["age_gp_2"] = "Age group"
names_map["covid_hosp_ever"] = "Previous COVID-19 hospitalisation"

#### Functions

# Fit logistic model and saves results table, coefficients and covariance matrix
lr_analysis = function(df_cohort, age_group, dep_var, ind_vars) {

  ## Logistic regression
  dir = paste0(output_dir, "lr_undervaccination_", age_group)
  if (!dir.exists(dir)) {
    dir.create(dir)
  }

  df_cohort = filter(df_cohort, age_gp == age_group) %>%
    droplevels()

  formula = as.formula(paste(c(paste0(dep_var, " ~ "), ind_vars), collapse = " + "))
  model = glm(formula, family = binomial, data = df_cohort)

  # Coefficients
  coef = data.frame("variable" = names(coef(model)), "coef" = coef(model), row.names = 1:length(coef(model)))
  write.csv(coef, paste0(dir, "/coef.csv"), row.names = FALSE)

  # Covariance matrix
  cov = vcov(model)
  write.csv(cov, paste0(dir, "/cov.csv"), row.names = FALSE)

  # Results table
  write.csv(create_logistic_results_table(df_cohort, model, dep_var, ind_vars), paste0(dir, "/results.csv"), row.names = FALSE)
}


# Does what it says
create_logistic_results_table = function(df_cohort, model, dep_var, ind_vars) {
  df_vars = data.frame()

  for (var in ind_vars) {
    Levels = unique(pull(df_cohort, !!sym(var)))

    df_vars = bind_rows(df_vars, data.frame(Variable = var, Levels))
  }

  df_vars = mutate_all(df_vars, as.character) %>%
    mutate(term = paste0(Variable, Levels)) %>%
    filter(!is.na(Levels))

  # Get ORs & 95% CIs
  results = round(exp(cbind("OR" = coef(model), confint.default(model, level = 0.95), digits = 2)), digits = 2) %>%
    as.data.frame() %>%
    select(-digits)

  names(results) = c("OR", "lower", "upper")

  results = mutate(results, term = rownames(results))

  results = mutate(results, OR = paste0(OR, " (", lower, ", ", upper, ")")) %>%
    select(term, OR)

  results = left_join(df_vars, results) %>%
    select(-term) %>%
    mutate(
      OR = replace_na(OR, "Ref"),
      Variable = names_map[Variable],
      Variable = case_when(
        !duplicated(Variable) ~ Variable,
        TRUE ~ ""
      )
    ) %>%
    rename("Odds ratio" = OR)
}


# Fits naive bayes model, saves graphs of fitted marginal distributions
nb_analysis = function(df_cohort, age_group, dep_var, ind_vars) {
  dir = paste0(output_dir, "nb_undervaccination_", age_group)
  if (!dir.exists(dir)) {
    dir.create(dir)
  }

  df_cohort = filter(df_cohort, age_gp == age_group) %>%
    droplevels()

  formula = as.formula(paste(c(paste0(dep_var, " ~ "), ind_vars), collapse = " + "))
  nb_model = naive_bayes(formula, data = df_cohort, usekernel = T)


  for (var in ind_vars) {
    png(paste0(dir, "/", var, ".png"), width = 350, height = 350)

    plot(nb_model, which = var)

    dev.off()
  }
}


# Saves table of person years, cumulative incidence curve, fits cox model and
# saves results table, coefficients and covariance matrix
cox_analysis = function(df_cohort, age_group, dep_var, ind_vars, calendar_days, study_start, study_end) {

  # Hospitalisation or death
  dir = paste0(output_dir, "cox_", dep_var, "_", age_group)
  if (!dir.exists(dir)) {
    dir.create(dir)
  }

  df_cohort = filter(df_cohort, age_gp == age_group) %>%
    droplevels()
  
  df_survival = create_df_survival(df_cohort, dep_var, calendar_days = calendar_days, study_start = study_start, study_end = study_end)

  # Cumulative incidence curve
  plot = survfit(Surv(tstart, tstop, event) ~ num_doses_time,
    # Reset the origin to take account of tim
    data = df_survival %>%
      group_by(EAVE_LINKNO, num_doses_time) %>%
      mutate(
        offset = min(tstart),
        tstart = tstart - offset,
        tstop = tstop - offset
      ) %>%
      select(-offset) %>%
      ungroup() %>%
      droplevels()
  ) %>%
    tidy() %>%
    mutate(
      strata = gsub("num_doses_time=0", "0", strata),
      strata = gsub("num_doses_time=1", "1", strata),
      strata = gsub("num_doses_time=2", "2", strata),
      strata = gsub("num_doses_time=3", "3", strata),
      strata = gsub("num_doses_time=4", "4", strata),
      strata = gsub("num_doses_time=5", "5", strata)
    ) %>%
    ggplot(aes(x = time, y = 1 - estimate, colour = strata, fill = strata)) +
    geom_line() +
    geom_ribbon(aes(ymin = 1 - conf.high, ymax = 1 - conf.low), linetype = 2, alpha = 0.1) +
    xlab("Time (days)") +
    ylab("Cumulative incidence") +
    labs(color = "Number of doses") +
    labs(fill = "Number of doses") +
    theme_bw()

  ggsave(paste0(dir, "/cumulative_incidence.png"))

  # Cox model
  formula = as.formula(paste0("Surv(tstart, tstop, event) ~ ", paste(ind_vars, collapse = " + ")))
  model = coxph(formula, data = df_survival)

  # Results table
  write.csv(create_cox_results_table(df_survival, model, ind_vars), paste0(dir, "/results_table.csv"), row.names = FALSE)

  # Coefficients
  coef = data.frame("variable" = names(coef(model)), "coef" = coef(model), row.names = 1:length(coef(model)))
  write.csv(coef, paste0(dir, "/coef.csv"), row.names = FALSE)

  # Covariance matrix
  cov = vcov(model)
  write.csv(cov, paste0(dir, "/cov.csv"), row.names = FALSE)
}


# Creates dataframe that is ready for survival model to be fitted.
# event_col is a variable whose value is the name of the event column. This is binary and says whether person had an event.
# event_date_col is a variable whose value is the name of the event date column.
# df must both of these columns.
create_df_survival = function(df, event_col, calendar_days, study_start, study_end) {
  event_date_col = paste0(event_col, "_date")

  # Add in endpoints
  cohort_endpoints = select(endpoints, EAVE_LINKNO, !!sym(event_col), !!sym(event_date_col)) %>%
    # Only use events that happened in the study period
    filter(!!sym(event_date_col) >= study_start & !!sym(event_date_col) <= study_end) %>%
    # Get date of first event
    group_by(EAVE_LINKNO) %>%
    mutate(!!sym(event_date_col) := min(!!sym(event_date_col), na.rm = TRUE)) %>%
    ungroup() %>%
    distinct()

  df = left_join(df, cohort_endpoints) %>%
    mutate(surv_date = pmin(study_end, !!sym(event_date_col), na.rm = TRUE)) %>%
    mutate(surv_date = pmin(surv_date, NRS.Date.Death, na.rm = TRUE)) %>%
    # Can potentially censor them if they get a fifth dose
    # Or we could say that if they have a fifth dose, they remain fully vaccinated
    # mutate(surv_date = pmin(surv_date, date_vacc_5, na.rm = TRUE)) %>%
    mutate(event = case_when(
      surv_date == !!sym(event_date_col) ~ 1,
      TRUE ~ 0
    )) %>%
    mutate(surv_time = as.numeric(surv_date - (study_start))) %>%
    # Some people died before study start and will therefore have negative survival time. Remove them.
    filter(surv_time >= 0) %>%
    # If they had an event on day of vaccination, assume that the event happened shortly after vaccination
    # The exact amount of time you add on doesn't matter when fitting Cox model, only the ordering.
    mutate(surv_time = ifelse(surv_time == 0, 0.001, surv_time))

  df = tmerge(df, df, id = EAVE_LINKNO, event = event(surv_time, event))

  # dataframe of start times for different vaccination status
  df_vs = select(df, EAVE_LINKNO, date_vacc_1, date_vacc_2, date_vacc_3, date_vacc_4, date_vacc_5) %>%
    mutate(
      vacc_0 = -Inf,
      vacc_1 = as.numeric(date_vacc_1 - study_start),
      vacc_2 = as.numeric(date_vacc_2 - study_start),
      vacc_3 = as.numeric(date_vacc_3 - study_start),
      vacc_4 = as.numeric(date_vacc_4 - study_start),
      vacc_5 = as.numeric(date_vacc_5 - study_start)
    ) %>%
    select(-date_vacc_1, -date_vacc_2, -date_vacc_3, -date_vacc_4, -date_vacc_5)

  df_vs = pivot_longer(df_vs,
    cols = c("vacc_0", "vacc_1", "vacc_2", "vacc_3", "vacc_4", "vacc_5"),
    names_to = "vs", values_to = "time"
  ) %>%
    filter(!is.na(time))

  # Add calendar period as a time-dependent variable
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
    num_doses_time = strtoi(str_sub(vs_time, -1, -1))
  ) %>%
    mutate(num_doses_time = factor(num_doses_time, levels = c(0, 1, 2, 3, 4, 5))) %>%
    select(-vs_time)
}


# Creates a summary table of Cox model results, with number of events,
# person years and hazard ratios
create_cox_results_table = function(df_survival, model, ind_vars) {
  person_years = data.frame()

  for (var in ind_vars) {

    formula = as.formula(paste0("Surv(duration, event) ~ ", var))

    person_years = bind_rows(
      person_years,
      pyears(formula,
        data = df_survival,
        scale = 1,
        data.frame = TRUE
        )$data %>%
        mutate(Variable = var) %>%
        rename(Levels = !!sym(var)) %>%
        mutate(Levels = as.character(Levels),
               HR = case_when( Levels == as.character(levels( as.factor(pull(df_survival, !!sym(var))))[1])  ~ 'Ref',
               TRUE ~ NA_character_))
    )
  }

  person_years =
    # Person years in the table are measured in units of thousands of years
    mutate(person_years, pyears = pyears / (365.25 * 1000),
      term = paste0(Variable, Levels)
    ) %>%
    select(Variable, Levels, pyears, n, event, term, HR)

  cox_results = tidy(model) %>%
    filter(!grepl("calendar", term)) %>%
    mutate(
      ucl = estimate + 1.96 * std.error,
      lcl = estimate - 1.96 * std.error
    ) %>%
    select(term, estimate, lcl, ucl) %>%
    mutate_if(is.numeric, ~ formatC(round(exp(.), 2), format = "f", big.mark = ",", drop0trailing = TRUE)) %>%
    mutate(estimate = paste0(estimate, " (", lcl, ", ", ucl, ")")) %>%
    select(-lcl, -ucl) %>%
    rename(HR = estimate) %>%
    data.frame()

  results_table = left_join(person_years, cox_results, by = 'term') %>%
    mutate(
      HR = coalesce(HR.x, HR.y),
      #HR = replace_na(HR, "Ref"),
      Variable = names_map[Variable],
      Variable = case_when(
        !duplicated(Variable) ~ Variable,
        TRUE ~ ""
      )
    ) %>%
    select(-term, -HR.x, -HR.y) %>%
    rename(
      `Number of events` = event,
      `Hazard ratio` = HR,
      `Persons` = n,
      `Person years (thousands)` = pyears
    )
}




#### Undervaccination analysis

## Logistic regression
lr_dep_var = "fully_vaccinated"
lr_ind_vars = c("Sex", "age_gp_2", "ur6_2016_name", "simd2020_sc_quintile", "n_risk_gps")

lr_analysis(df_cohort, "5-15", lr_dep_var, lr_ind_vars)
lr_analysis(df_cohort, "16-74", lr_dep_var, lr_ind_vars)
# age_gp_2 doesn't have any subcategories for over 75, so remove it
lr_analysis(df_cohort, "75+", lr_dep_var, setdiff(lr_ind_vars, "age_gp_2"))


## Naive bayes model
nb_dep_var = "fully_vaccinated"
nb_ind_vars = c("Sex", "ur6_2016_name", "simd2020_sc_quintile", "n_risk_gps", "ageYear")

nb_analysis(df_cohort, "5-15", nb_dep_var, nb_ind_vars)
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







#### Severe outcomes analysis

# End dates
study_end_hosp = max(max(smr$ADMISSION_DATE, na.rm = TRUE),
  max(smr$DISCHARGE_DATE, na.rm = TRUE),
  na.rm = TRUE
)
study_end_hosp = min(study_end, study_end_hosp)

study_end_death = max(all_deaths$NRS.Date.Death, na.rm = TRUE)
study_end_death = min(study_end, study_end_death)

study_end_hosp_death = min(study_end_hosp, study_end_death, na.rm = TRUE)
study_end_hosp_death = min(study_end, study_end_hosp_death)

cox_ind_vars = c(
  "Sex",
  "age_gp_2",
  "ur6_2016_name",
  "simd2020_sc_quintile",
  "n_risk_gps",
  "last_positive_test_group",
  "shielding",
  "num_tests_6m_group",
  "covid_hosp_ever",
  "num_doses_time"
)

cox_analysis(df_cohort,
  age_group = "5-15", 
  dep_var = "covid_hosp_death", 
  ind_vars = cox_ind_vars,
  calendar_days = NA, 
  study_start = study_start, 
  study_end = study_end_hosp_death
)

cox_analysis(df_cohort,
  age_group = "16-74", 
  dep_var = "covid_hosp_death", 
  ind_vars = cox_ind_vars,
  calendar_days = NA, 
  study_start = study_start, 
  study_end = study_end_hosp_death
)

# There is only one age group in 75+, so remove age_gp_2 from cox_ind_vars
cox_analysis(df_cohort,
  age_group = "75+", 
  dep_var = "covid_hosp_death", 
  ind_vars = setdiff(cox_ind_vars, 'age_gp_2'),
  calendar_days = NA, 
  study_start = study_start, 
  study_end = study_end_hosp_death
)





#### Tests
# 
# df_cohort_backup = df_cohort
# 
# df_cohort = df_cohort %>%
#   filter(age_gp == "16-74")
# df_cohort = filter(df_cohort, EAVE_LINKNO %in% sample(df_cohort$EAVE_LINKNO, 10000)) %>%
#   droplevels()
# 
# df_survival = create_df_survival(df_cohort,
#   event_col = "covid_hosp_death", calendar_days = 0, study_start = study_start, study_end = study_end
# )
# 
# # Check everything is kosher - verify that bob is indeed one's uncle
# bob = select(
#   df_survival, EAVE_LINKNO, age_gp, date_vacc_1, date_vacc_2, date_vacc_3, date_vacc_4, date_vacc_5,
#   vacc_type_1, vacc_type_2, vacc_type_3, vacc_type_4, vacc_type_5,
#   num_doses_start, num_doses_recent, vacc_seq_start, vacc_seq_recent,
#   vs_start, vs_recent, fully_vaccinated, vs_mixed_start, event, surv_date,
#   tstart, tstop, duration, num_doses_time
# )
# 
# # Cox model
# formula = as.formula(paste0("Surv(tstart, tstop, event) ~ ", paste(cox_ind_vars, collapse = " + ")))
# model = coxph(formula, data = df_survival)
# 
# results_table = create_cox_results_table(df_cohort, model, cox_ind_vars)
# 
# cox_analysis(df_cohort,
#   age_group = "16-74", 
#   dep_var = "covid_hosp_death",
#   ind_vars = cox_ind_vars, 
#   calendar_days = 0, 
#   study_start, 
#   study_end
# )
# 
# df_cohort = df-cohort_backup
# 
# rm(df_cohort_backup)
