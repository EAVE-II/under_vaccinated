######################################################################
## Title: [Insert full title of paper]
## Code author: Steven Kerr steven.kerr@ed.ac.uk
## Description: Analysis of severe COVID-19 events
##              01_data_cleaning and 02_data_prep need to be run first.
##              names_map from  03_cohort_description is also needed.
######################################################################

library(survival)
library(lubridate)
library(broom)
library(survey)

output_dir = "./output/"

# Pre-requisites
#source("./code/01_data_cleaning.R")
#source("./code/02_data_prep.R")
# The only thing needed from 02_cohort_description is names_map
# source('./output/02_cohort_description.R')

# Add more entries to names_map to automatically generate display names
names_map["last_positive_test_group"] = "Last positive test"
names_map["num_tests_6m_group"] = "Number of tests in last 6 months"
names_map["shielding"] = "Ever shielding"
names_map["num_doses_time_2"] = "Number of doses"
names_map["age_17cat"] = "Age group"
names_map["covid_mcoa_hosp_ever"] = "Previous COVID-19 hospitalisation"
names_map["urban_rural_2cat"] = "Urban/rural classification"
names_map['competing_hosp'] = "Competing non-COVID-19 hospitalisation"

#### Functions

# Saves table of person years, cumulative incidence curve, fits cox model and
# saves results table, coefficients and covariance matrix
cox_analysis = function(df_cohort, age_group, dep_var, ind_vars, calendar_days, study_start, study_end, exclude_non_unit_weight) {

  # Hospitalisation or death
  dir = paste0(output_dir, "cox_", dep_var, "_", age_group)

  df_cohort = filter(df_cohort, age_3cat == age_group) %>%
    droplevels()

  df_survival = create_df_survival(df_cohort, dep_var, calendar_days = calendar_days, study_start = study_start, study_end = study_end)

  # Make factor
  if (age_group %in% c("5-17", "18-74")) {
    df_survival = mutate(
      df_survival,
      num_doses_time_2 = fct_relevel(num_doses_time_2, "0", "1", "2", "3+")
    )
  } else if (age_group == "75+") {
    df_survival = mutate(
      df_survival,
      num_doses_time_2 = fct_relevel(num_doses_time_2, "0", "1", "2", "3", "4+")
    )
  }

  formula = as.formula(paste0("Surv(tstart, tstop, event) ~ ", paste(ind_vars, collapse = " + ")))
  
  if (exclude_non_unit_weight){
    
    dir = paste0(dir, "_weight1")
    if (!dir.exists(dir)) {
      dir.create(dir)
    }
    
    
    df_survival = filter(df_survival, eave_weight == 1)
    
    # Robust standard errors are used by default if weights are not all equal to 1
    model = coxph(formula, data = df_survival)
    
    # This takes up too much space
    # tryCatch({
    #   saveRDS(model, paste0(dir, "/model.rds"))
    # }, 
    # error = function(e){cat("ERROR :",conditionMessage(e), "\n")}
    # )
    
  } else {
    
    if (!dir.exists(dir)) {
      dir.create(dir)
    }
    
    # Weighted Cox model
    # Throws errors - not clear why
    #
    # survey_design = svydesign(id = ~1,
    #                           weights = ~ eave_weight,
    #                           data = df_survival)
    #
    # model = svycoxph( formula,
    #                   design = survey_design,
    #                   data = df_survival)
    
    # Normalise weight to 1 so standard errors are not distorted
    df_survival = mutate(df_survival, eave_weight = eave_weight * nrow(df_survival) / sum(eave_weight))
    
    # Robust standard errors are used by default if weights are not all equal to 1
    model = coxph(formula, data = df_survival, weights = eave_weight)
    
    # This takes up too much space
    # tryCatch({
    #   saveRDS(model, paste0(dir, "/model.rds"))
    # }, 
    # error= function(e){cat("ERROR :",conditionMessage(e), "\n")}
    # )
  }
  
  # Results table
  write.csv(create_cox_results_table(df_survival, model, ind_vars), paste0(dir, "/results.csv"), row.names = FALSE)

  # Coefficients
  coef = data.frame("variable" = names(coef(model)), "coef" = coef(model), row.names = 1:length(coef(model)))
  write.csv(coef, paste0(dir, "/coef.csv"), row.names = FALSE)

  # Covariance matrix
  cov = vcov(model)
  write.csv(cov, paste0(dir, "/cov.csv"))
  
  # Cumulative incidence curve
  plot = survfit(Surv(tstart, tstop, event) ~ num_doses_time_2,
                 # Reset the origin to take account of tim
                 data = df_survival %>%
                   group_by(individual_id, num_doses_time_2) %>%
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
      strata = gsub("num_doses_time_2=0", "0", strata),
      strata = gsub("num_doses_time_2=1", "1", strata),
      strata = gsub("num_doses_time_2=2", "2", strata),
      strata = gsub("num_doses_time_2=3", "3", strata),
      strata = gsub("num_doses_time_2=4", "4", strata),
      strata = gsub("num_doses_time_2=5", "5", strata),
      strata = gsub("num_doses_time_2=3+", "3+", strata),
      strata = gsub("num_doses_time_2=3+", "4+", strata)
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
}


# Creates dataframe that is ready for survival model to be fitted.
# event_col is a variable whose value is the name of the event column. This is binary and says whether person had an event.
# event_date_col is a variable whose value is the name of the event date column.
create_df_survival = function(df, event_col, calendar_days, study_start, study_end) {
  
  event_date_col = paste0(event_col, "_date")

  df = df %>%
    mutate(surv_date = pmin(study_end, !!sym(event_date_col), na.rm = TRUE)) %>%
    mutate(surv_date = pmin(surv_date, death_date, na.rm = TRUE)) %>%
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

  df = tmerge(df, df, id = individual_id, event = event(surv_time, event))

  # dataframe of start times for different vaccination status
  df_vs = select(df, individual_id, date_vacc_1, date_vacc_2, date_vacc_3, date_vacc_4, date_vacc_5) %>%
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

  # Add in vaccination status as a time dependent variable
  df = tmerge(df, df_vs, id = individual_id, exposure = tdc(time, vs)) %>%
    rename(vs_time = exposure)
  
  # If event involves hospitalisation, put in competing hospitalisation events
  # as time-dependent variable
  if (str_detect(event_col, 'hosp')){
    
    # dataframe of start times for other hospitalisations
    df_competing_hosp = select(df, individual_id, !!sym(paste0('non_', event_date_col))) %>%
      mutate(
        `0` = -Inf,
        `1` = as.numeric(!!sym(paste0('non_', event_date_col)) - study_start),
      ) %>%
      select(-!!sym(paste0('non_', event_date_col)))
    
    df_competing_hosp = pivot_longer(df_competing_hosp,
                         cols = c("0", "1"),
                         names_to = "competing_hosp", values_to = "time"
    ) %>%
      filter(!is.na(time)) %>%
      mutate(competing_hosp = as.numeric(competing_hosp))
    
    # Add in vaccination status as a time dependent variable
    df = tmerge(df, df_competing_hosp, id = individual_id, exposure = tdc(time, competing_hosp)) %>%
      rename(competing_hosp = exposure)
  }
  

  # dataframe of start times for different calendar periods
  # df_calendar = select(df, individual_id)
  #
  # calendar_start = seq(study_start, max(df$surv_date), make_difftime(calendar_days * 24 * 60 * 60, 'days'))
  #
  # for (i in 1:length(calendar_start)){
  #   df_calendar[as.character(calendar_start[i])] = as.numeric( calendar_start[i] - study_start)
  # }
  #
  # df_calendar = pivot_longer(df_calendar, cols = as.character(calendar_start),
  #                            names_to = 'calendar', values_to = 'time')
  
  # Add in calendar period as a time dependent variable
  # df = tmerge(df, df_calendar, id = individual_id, exposure = tdc(time, calendar)) %>%
  #   rename(calendar = exposure)

  df = mutate(df,
    duration = tstop - tstart,
    num_doses_time = strtoi(str_sub(vs_time, -1, -1))
  ) %>%
    mutate(num_doses_time_2 = case_when(
      age_3cat %in% c("5-17", "18-74") & num_doses_time >= 3 ~ "3+",
      age_3cat == "75+" & num_doses_time >= 4 ~ "4+",
      TRUE ~ as.character(num_doses_time)
    )) %>%
    mutate(num_doses_time_2 = factor(num_doses_time_2)) %>%
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
        mutate(
          Levels = as.character(Levels),
          HR = case_when(
            Levels == as.character(levels(as.factor(pull(df_survival, !!sym(var))))[1]) ~ "Ref",
            TRUE ~ NA_character_
          )
        )
    )
  }

  person_years =
    # Person years in the table are measured in units of thousands of years
    mutate(person_years,
      pyears = pyears / (365.25 * 1000),
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

  results_table = left_join(person_years, cox_results, by = "term") %>%
    mutate(
      HR = coalesce(HR.x, HR.y),
      # HR = replace_na(HR, "Ref"),
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



#### Severe outcomes analysis

cox_ind_vars = c(
  "sex",
  "age_17cat",
  "urban_rural_2cat",
  "simd2020_sc_quintile",
  "n_risk_gps",
  "last_positive_test_group",
  "shielding",
  "num_tests_6m_group",
  "covid_mcoa_hosp_ever",
  "competing_hosp",
  "num_doses_time_2"
)

endpoint_names = c(
  "covid_death", 
  "covid_acoa_hosp",
  "covid_mcoa_hosp",
  "covid_mcoa_28_2_hosp",
  "covid_mcoa_14_2_hosp",
  "covid_mcoa_hosp_death")



for (dep_var in endpoint_names) {
  print(dep_var)
  
  for (age_group in c("5-17", "18-74", "75+")) {
    print(age_group)
    
    for (exclude_non_unit_weight in c(TRUE, FALSE)){
      print(exclude_non_unit_weight)
      
      if (dep_var == 'covid_death'){
        ind_vars = setdiff(cox_ind_vars, 'competing_hosp')
      } else {
        ind_vars = cox_ind_vars
      }
  
      cox_analysis(df_cohort,
        age_group = age_group,
        dep_var = dep_var,
        ind_vars = ind_vars,
        calendar_days = NA,
        study_start = study_start,
        study_end = study_end,
        exclude_non_unit_weight = exclude_non_unit_weight
    )
    }
  }
}


#### Tests
#
# df_cohort_test = df_cohort %>%
#   filter(age_3cat == "18-74")
# df_cohort_test = filter(df_cohort_test, individual_id %in% sample(df_cohort_test$individual_id, 100000)) %>%
#   droplevels()
# 
# df_survival = create_df_survival(df_cohort_test,
#   event_col = "covid_mcoa_hosp", calendar_days = 0, study_start = study_start, study_end = study_end
# )
# 
# # # Check everything is ok - verify that bob is indeed one's uncle
# bob = select(
#   df_survival, individual_id, age_3cat, date_vacc_1, date_vacc_2, date_vacc_3, date_vacc_4, date_vacc_5,
#   vacc_type_1, vacc_type_2, vacc_type_3, vacc_type_4, vacc_type_5,
#   num_doses_start, num_doses_recent, vacc_seq_start, vacc_seq_recent,
#   vs_start, vs_recent, fully_vaccinated, vs_mixed_start,
#   covid_death_date, covid_acoa_hosp_date, covid_mcoa_hosp_date,
#   non_covid_mcoa_hosp_date, covid_mcoa_28_2_hosp_date, covid_mcoa_14_2_hosp, covid_mcoa_hosp_death_date, 
#   non_covid_acoa_hosp_date, non_covid_mcoa_28_2_hosp_date,
#   non_covid_mcoa_14_2_hosp_date, competing_hosp, event, surv_date,
#   tstart, tstop, duration, num_doses_time
# )
# 
# # Weighted Cox model
# formula = as.formula(paste0("Surv(tstart, tstop, event) ~ ", paste(cox_ind_vars, collapse = " + ")))
# 
# # survey_design = svydesign(
# #   id = ~1,
# #   weights = ~eave_weight,
# #   data = df_survival
# # )
# #
# # model = svycoxph(formula,
# #   design = survey_design,
# #   data = df_survival
# # )
# 
# df_survival = mutate(df_survival, eave_weight = eave_weight * nrow(df_survival) / sum(eave_weight))
# 
# model = coxph(formula, data = df_survival, weights = eave_weight)
# 
# results_table = create_cox_results_table(df_survival, model, cox_ind_vars)
# 
# cox_analysis(df_cohort,
#   age_group = "18-74",
#   dep_var = "covid_mcoa_hosp_death",
#   ind_vars = cox_ind_vars,
#   calendar_days = 0,
#   study_start,
#   study_end,
#   exclude_non_unit_weight = TRUE
# )
