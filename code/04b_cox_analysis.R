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
library(survminer)
library(gridExtra)
library(car)

setwd("/conf/EAVE/GPanalysis/analyses/under_vaccinated")

output_dir = "./output/"

# Pre-requisites
#source("./code/01_data_cleaning.R")
#source("./code/02_data_prep.R")
# The only thing needed from 02_cohort_description is names_map
# source('./output/02_cohort_description.R')

# Read in
df_cohort = readRDS('./data/df_cohort.rds')

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

# Saves table of person years, cumulative incidence curve, fits cox model and
# saves results table, coefficients and covariance matrix
cox_analysis = function(df_cohort, age_group, dep_var, ind_vars, calendar_days, study_start, study_end, 
                        exclude_non_unit_weight, weight, name){
  
  # Hospitalisation or death
  dir = paste0(output_dir, "cox_", dep_var, "_", age_group, "_", name)
  
  # Remove inconsistent vaccination records
  df_cohort = df_cohort %>%
    filter(flag_incon_date == 0, incon_timing == 0)

  if (age_group != 'all'){
  df_cohort = filter(df_cohort, age_3cat == age_group) %>%
    droplevels()
  }

  df_survival = create_df_survival(
    df_cohort, 
    dep_var, 
    ind_vars, 
    calendar_days = calendar_days, 
    study_start = study_start, 
    study_end = study_end)

  # Relevel vaccination status factor
  if (age_group == "5-15" & 'num_doses_time_2' %in% ind_vars) {
    df_survival = df_survival %>% 
      mutate(
        num_doses_time_2 = fct_relevel(num_doses_time_2, "2+", "0", "1")
      )
  } else if (age_group == "16-74" & 'num_doses_time_2' %in% ind_vars) {
    df_survival = df_survival %>% 
      mutate(
        num_doses_time_2 = fct_relevel(num_doses_time_2, "3+", "0", "1", "2")
    )
  } else if (age_group == "75+" & 'num_doses_time_2' %in% ind_vars){
    df_survival = df_survival %>%
      mutate(
        num_doses_time_2 = fct_relevel(num_doses_time_2, "4+", "0", "1", "2", "3")
    )
  } else if (age_group == "all" & 'num_doses_time_3' %in% ind_vars){
    df_survival = df_survival %>%
      mutate(
        num_doses_time_3 = fct_relevel(num_doses_time_3, "4+", "0", "1", "2", "3")
    )
  }
  
  # Relevel age group for children and elderly
  if (age_group %in% c("5-15") & 'num_doses_time_2' %in% ind_vars) {
    df_survival = df_survival %>% 
      mutate(
        age_17cat = fct_relevel(age_17cat, "5-11")
      )
  } else if (age_group == "75+" & 'num_doses_time_2' %in% ind_vars){
    df_survival = df_survival %>% 
      mutate(
        age_17cat = fct_relevel(age_17cat, "75-79")
      )

  } 
  
  formula = as.formula(paste0("Surv(tstart, tstop, event) ~ ", paste(ind_vars, collapse = " + ") ))
  
  if ('under_vaccinated_time' %in% ind_vars){
    dose_time_var = "under_vaccinated_time"
    
    plot = labs(color = "Vaccination status", fill = "Vaccination status")
  } else if (age_group == 'all'){
    dose_time_var = "num_doses_time_3"
    
    plot = labs(color = "Number of doses", fill = "Number of doses") 
  } else {
    dose_time_var = "num_doses_time_2"
    
    plot = labs(color = "Number of doses", fill = "Number of doses") 
  }
  
  
  if (exclude_non_unit_weight){
    
    dir = paste0(dir, "_weight1")
    if (!dir.exists(dir)) {
      dir.create(dir)
    }
    
    df_survival = filter(df_survival, eave_weight == 1)
    
    model_adj = coxph(formula, data = df_survival)
    model_unadj = coxph(as.formula(paste0("Surv(tstart, tstop, event) ~ ", dose_time_var)), data = df_survival)
  } else if (weight){
    
    dir = paste0(dir, "_weighted")
    if (!dir.exists(dir)) {
      dir.create(dir)
    }
    
    # Weighted Cox model
    # Throws errors - not clear why
    #
    # survey_design = svydesign(
    #  id = ~1,
    #  weights = ~ eave_weight,
    #  data = df_survival)
    #
    # model = svycoxph( 
    #    formula,
    #    design = survey_design,
    #    data = df_survival)
    
    # Normalise weight to 1 so standard errors are not distorted
    df_survival = mutate(df_survival, eave_weight = eave_weight * nrow(df_survival) / sum(eave_weight))
    
    # Robust standard errors are used by default if weights are not all equal to 1
    model_adj = coxph(formula, data = df_survival, weights = eave_weight)
    model_unadj = coxph(as.formula(paste0("Surv(tstart, tstop, event) ~ ", dose_time_var)), data = df_survival, weights = eave_weight)
  } else {
    
    if (!dir.exists(dir)) {
      dir.create(dir)
    }
    
    model_adj = coxph(formula, data = df_survival)
    model_unadj = coxph(as.formula(paste0("Surv(tstart, tstop, event) ~ ", dose_time_var)), data = df_survival)
  }

  # Unadjusted
  # Results table
  write.csv(create_cox_results_table(df_survival, model_unadj, dose_time_var), paste0(dir, "/results_unadj.csv"), row.names = FALSE)
  
  # Adjusted
  # Results table
  write.csv(create_cox_results_table(df_survival, model_adj, ind_vars), paste0(dir, "/results_adj.csv"), row.names = FALSE)
  # Coefficients
  coef = data.frame("variable" = names(coef(model_adj)), "coef" = coef(model_adj), row.names = 1:length(coef(model_adj)))
  write.csv(coef, paste0(dir, "/coef.csv"), row.names = FALSE)
  # Covariance matrix
  cov = vcov(model_adj)
  write.csv(cov, paste0(dir, "/cov.csv"))
  
  # Schoenfeld residuals
  # In tryCatch because sometimes the smoothing algorithm doesn't converge
  tryCatch(
    {
      ph_test = cox.zph(model_adj)
      df = data.frame(
        print(ph_test )
      )
      write.csv(df, paste0(dir, "/ph_test.csv"))
      #pdf(paste0(dir, "/ggcoxzph.pdf"))
      
      #plot_ph = ggcoxzph(ph_test )
      #dev.off()
      #ggsave(paste0(dir, "/ggcoxzph.pdf"), arrangeGrob(grobs = plot_ph))
      
      resid = data.frame(
        trans_time = ph_test$x,
        time = ph_test$time,
        ph_test$y
      ) %>% 
        pivot_longer(
          cols = c(-trans_time, -time)
        ) 
      
      resid$name =names_map[resid$name]
      
      resid %>% 
        ggplot(aes(
          x = trans_time,
          y = value
        )) +
        facet_wrap(~name) +
        geom_smooth()

      ggsave(paste0(dir, "/schoenfeld_resid.png"))
    },
    error = function(e) {
      cat("ERROR :", conditionMessage(e), "\n")
    }
  )
 
  # Cumulative incidence curve
  if ('under_vaccinated_time' %in% ind_vars){
    CI_table = survfit(
      as.formula("Surv(duration, event) ~ under_vaccinated_time"),
      # Reset the origin to take account of time
      # That is, e.g. for dose 1 we start the clock at 0 on the day they received their dose
      # Then when they receive dose 2, we restrart the clock again at 0
      # Necessary to do it this way rather than just using the duration column
      # because they could have other time-dependent variables
      data = df_survival %>%
        group_by(individual_id, under_vaccinated_time) %>%
        mutate(
          offset = min(tstart),
          tstart = tstart - offset,
          tstop = tstop - offset
        ) %>%
        select(-offset) %>%
        ungroup() %>%
        droplevels()
    )
  } else {
    CI_table = survfit(
      as.formula("Surv(duration, event) ~ under_vaccinated_time"),
      # Reset the origin to take account of time
      # That is, e.g. for dose 1 we start the clock at 0 on the day they received their dose
      # Then when they receive dose 2, we restrart the clock again at 0
      # Necessary to do it this way rather than just using the duration column
      # because they could have other time-dependent variables
      data = df_survival %>%
        mutate(
          
          under_vaccinated_time = case_when(
            age_3cat == '5-15' & !!sym(dose_time_var) == '2+' ~ 'fully_vaccinated',
            age_3cat == '5-15' ~ 'under_vaccinated',
            age_3cat == '16-74' & !!sym(dose_time_var) == '3+' ~ 'fully_vaccinated',
            age_3cat == '16-74' ~ 'under_vaccinated',
            age_3cat == '75+' & !!sym(dose_time_var) == '4+' ~ 'fully_vaccinated',
            age_3cat == '75+' ~ 'under_vaccinated',
            
          )
        ) %>%
        group_by(individual_id, under_vaccinated_time) %>%
        mutate(
          offset = min(tstart),
          tstart = tstart - offset,
          tstop = tstop - offset
        ) %>%
        select(-offset) %>%
        ungroup() %>%
        droplevels()
    )
  }
  
  CI_table = CI_table %>%
    tidy() %>%
    mutate(
      strata = gsub("num_doses_time_2=0", "0", strata),
      strata = gsub("num_doses_time_2=1", "1", strata),
      strata = gsub("num_doses_time_2=2", "2", strata),
      strata = gsub("num_doses_time_2=3", "3", strata),
      strata = gsub("num_doses_time_2=4", "4", strata),
      strata = gsub("num_doses_time_2=5", "5", strata),
      strata = gsub("num_doses_time_2=3+", "3+", strata),
      strata = gsub("num_doses_time_2=3+", "4+", strata),
      
      strata = gsub("num_doses_time_3=0", "0", strata),
      strata = gsub("num_doses_time_3=1", "1", strata),
      strata = gsub("num_doses_time_3=2", "2", strata),
      strata = gsub("num_doses_time_3=3", "3", strata),
      strata = gsub("num_doses_time_3=4", "4", strata),
      strata = gsub("num_doses_time_3=5", "5", strata),
      strata = gsub("num_doses_time_3=3+", "3+", strata),
      strata = gsub("num_doses_time_3=3+", "4+", strata),
      
      strata = gsub("under_vaccinated_time=under_vaccinated", "Under-vaccinated", strata),
      strata = gsub("under_vaccinated_time=fully_vaccinated", "Fully-vaccinated", strata)
    ) 
  
  CI_table %>%
    ggplot(aes(x = time, y = 1 - estimate, colour = strata, fill = strata)) +
    geom_line() +
    geom_ribbon(aes(ymin = 1 - conf.high, ymax = 1 - conf.low), linetype = 2, alpha = 0.1) +
    xlab("Time (days)") +
    ylab("Cumulative incidence") +
    theme_bw() 
  
  ggsave(paste0(dir, "/cumulative_incidence.png"))
  
  CI_table = CI_table %>%
    mutate(
      week = floor(time/7)
    ) %>%
    group_by(week) %>%
    mutate(
      n.risk = sum(n.risk),
      n.event = sum(n.event)
    ) %>%
    select(
      week, strata, n.risk, n.event
    ) %>%
    unique()
  
  write.csv(CI_table, paste0(dir, "/CI_table.csv"))
  
  # OE analysis
  #write.csv(OE_analysis(model, df_survival, age_group), paste0(dir, "/OE_analysis.csv"))
  
  # Variance inflation factor, for testing collinearity
  tryCatch(
    {
      write.csv(vif(model_adj), paste0(dir, "/vif.csv"))
    },
    error = function(e) {
      cat("ERROR :", conditionMessage(e), "\n")
    }
  )
}


# Creates dataframe that is ready for survival model to be fitted.
# event_col is a variable whose value is the name of the event column. This is binary and says whether person had an event.
# event_date_col is a variable whose value is the name of the event date column.
create_df_survival = function(df, event_col, ind_vars, calendar_days, study_start, study_end) {
  
  event_date_col = paste0(event_col, "_date")

  df = df %>%
    mutate(surv_date = pmin(study_end, !!sym(event_date_col), death_date, date_transfer_out, na.rm = TRUE)) %>%
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

  # If event involves hospitalisation, put in competing hospitalisation events
  # as time-dependent variable
  if ("competing_hosp" %in% ind_vars){
    
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
    
    df = tmerge(df, df_competing_hosp, id = individual_id, exposure = tdc(time, competing_hosp)) %>%
      rename(competing_hosp = exposure)
  }
  
  # Add vaccination dose number as time-dependent variable
  if (any(c('num_doses_time_2', 'num_doses_time_3') %in% ind_vars)){
  
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
    
    df = tmerge(df, df_vs, id = individual_id, exposure = tdc(time, vs)) %>%
      rename(vs_time = exposure)
    
    df = df %>%
      mutate(num_doses_time = strtoi(str_sub(vs_time, -1, -1))) %>%
      mutate(
        num_doses_time_2 = case_when(
          age_3cat == "5-15" & num_doses_time >= 2 ~ "2+",
          age_3cat =="16-74" & num_doses_time >= 3 ~ "3+",
          age_3cat == "75+" & num_doses_time >= 4 ~ "4+",
          TRUE ~ as.character(num_doses_time)
        ),
        num_doses_time_3 = case_when(
          num_doses_time >= 4 ~ "4+",
          TRUE ~ as.character(num_doses_time)  
        )
      ) %>%
      mutate(num_doses_time_2 = factor(num_doses_time_2)) %>%
      select(-vs_time)
  }
  
  # Add fully/under vaccinated as time dependent variable
  if ('under_vaccinated_time' %in% ind_vars){
    
    # dataframe of start times for different vaccination status
    df_vs = select(df, individual_id, date_fully_vaccinated) %>%
      mutate(
        under_vaccinated = -Inf,
        fully_vaccinated = as.numeric(date_fully_vaccinated - study_start)
      ) %>%
      select(-date_fully_vaccinated)
    
    df_vs = pivot_longer(df_vs,
                         cols = c("under_vaccinated", "fully_vaccinated"),
                         names_to = "vs", values_to = "time"
    ) %>%
      filter(!is.na(time))
    
    df = tmerge(df, df_vs, id = individual_id, exposure = tdc(time, vs)) %>%
      rename(vs_time = exposure)
    
    df = df %>%
      rename(under_vaccinated_time = vs_time) %>%
      mutate(under_vaccinated_time = fct_relevel(under_vaccinated_time, "fully_vaccinated"))
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

  df = df %>%
    mutate(duration = tstop - tstart)
}


# Creates a summary table of Cox model results, with number of events,
# person years and hazard ratios
create_cox_results_table = function(df_survival, model, ind_vars) {
  person_years = data.frame()

  for (var in ind_vars) {
    formula = as.formula(paste0("Surv(duration, event) ~ ", var))

    new_df = pyears(formula,
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
    
    # For binary variables, only keep 1 row
    if (length(setdiff(new_df$Levels, c(0,1))) == 0){
      new_df = filter(new_df, Levels == 1) %>%
        mutate(Levels = '')
    }
    
    person_years = bind_rows(
      person_years,
      new_df)
  }

  person_years =
    # Person years in the table are measured in units of thousands of years
    mutate(person_years,
      pyears = pyears / (365.25 * 1000),
      term = paste0(Variable, Levels),
      rate = round(event/pyears, 2)
    ) %>%
    select(Variable, Levels, pyears, n, event, rate, term, HR)

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
      `Person years (thousands)` = pyears,
      `Event rate per thousand person years` = rate
    )
}


# Observed-expected calculations where expected is a fully vaccinated
# counterfactual
OE_analysis = function(model, df_survival, age_group){

  under_vacc = df_survival %>%
    filter(under_vaccinated == 1, !duplicated(individual_id))
  
  if (age_group %in% c('5-15', '16-74')){
    under_vacc = under_vacc %>%
      mutate(
        num_doses_time_2 = case_when(
          age_4cat == '5-11' ~ '1',
          age_4cat == '12-15' ~ '2',
          age_4cat == '16-74' ~ '3+'
        ) %>% factor(levels = c('0', '1', '2', '3+', '4+')),
        tstart = 0,
        tstop = 121
      )
  } else if (age_group %in% c('75+')){
    under_vacc = under_vacc %>%
      mutate(
        num_doses_time_2 = case_when(
          age_4cat == '5-11' ~ '1',
          age_4cat == '12-15' ~ '2',
          age_4cat == '16-74' ~ '3',
          age_4cat == '75+' ~ '4'
        ) %>% factor(levels = c('0', '1', '2', '3', '4+')),
        tstart = 0,
        tstop = 121
      )
    
  }
  
  observed = sum(under_vacc$event)

  fully_vacc_pred = predict(
    model,
    newdata = under_vacc, 
    type = 'expected', 
    se.fit = TRUE)
  
  expected_fully_vacc = sum(fully_vacc_pred[['fit']], na.rm = TRUE)
  expected_fully_vacc_sd = sqrt(sum(fully_vacc_pred[['se.fit']]^2, na.rm = TRUE))
  
  df = data.frame(
    reduction = observed - expected_fully_vacc,
    reduction_lower_95_ci = observed - (expected_fully_vacc + 1.96 * expected_fully_vacc_sd),
    reduction_upper_95_ci = observed - (expected_fully_vacc - 1.96 * expected_fully_vacc_sd)
  )
}


#### Severe outcomes analysis

cox_ind_vars_extended = c(
  "num_doses_time_2",
  "sex",
  "age_17cat",
  #"bmi_imputed_cat",
  "urban_rural_2cat",
  "simd2020_sc_quintile",
  "ethnicity_5cat",
  "n_risk_gps_6cat",
  "n_risk_gps_3cat",
  "health_board",
  "last_positive_test_group",
  "shielding",
  "other_household_shielding",
  "num_tests_6m_group",
  "num_pos_tests_6m_group",
  "covid_mcoa_hosp_ever",
  "competing_hosp"
  #"under_vaccinated_time"
)

cox_ind_vars_minimal = c(
  "num_doses_time_2",
  "sex",
  "age_17cat",
  "urban_rural_2cat",
  "simd2020_sc_quintile",
  "ethnicity_5cat",
  "n_risk_gps_6cat",
  "n_risk_gps_3cat"
  #"under_vaccinated_time"
)

var_list = list(
  "minimal" =  cox_ind_vars_minimal,
  "extended" = cox_ind_vars_extended)

endpoint_names = c(
  "covid_death", 
  
  #"covid_acoa_hosp",
  #"covid_mcoa_hosp",
  #"covid_mcoa_28_2_hosp",
  #"covid_mcoa_14_2_hosp",
  #"covid_mcoa_hosp_death",
  #"covid_acoa_hosp_death",
  
  #"covid_acoa_emerg_hosp",
  "covid_mcoa_emerg_hosp",
  # "covid_mcoa_28_2_emerg_hosp",
  # "covid_mcoa_14_2_emerg_hosp",
  "covid_mcoa_emerg_hosp_death"
  #"covid_acoa_emerg_hosp_death"
  )


for (dep_var in endpoint_names) {
  print(dep_var)
  
  for (age_group in c("5-15", "16-74", "75+")) {
    print(age_group)
    
    for (var_set in 1:length(var_list)){
      
      ind_vars = var_list[[var_set]]
      name = names(var_list)[var_set]
      
      print(name)
    
      if (dep_var == 'covid_death'){
        ind_vars = setdiff(ind_vars, 'competing_hosp')
      }
      
      if (age_group == 'all'){
        # num_doses_time_2 has levels 0, 1, 2, 3, 3+, 4+, depending on age group
        # num_doses_time_3 has levels 0, 1, 2, 3, 4+
        ind_vars = replace(ind_vars, ind_vars == "num_doses_time_2", "num_doses_time_3")
      } else if (age_group == '5-15'){
        ind_vars = setdiff(ind_vars, c('n_risk_gps_6cat', 'bmi_imputed_cat'))
      } else {
        ind_vars = setdiff(ind_vars, 'n_risk_gps_3cat')
      }
  
      cox_analysis(df_cohort,
        age_group = age_group,
        dep_var = dep_var,
        ind_vars = ind_vars,
        calendar_days = NA,
        study_start = study_start,
        study_end = study_end,
        exclude_non_unit_weight = TRUE,
        weight = FALSE,
        name = name
      )
    }
  }
}


#### Tests
# 
# df_cohort_test = df_cohort %>%
#   filter(age_3cat == "16-74", flag_incon_date == 0, incon_timing == 0)
# df_cohort_test = filter(df_cohort_test, individual_id %in% sample(df_cohort_test$individual_id, 100000)) %>%
#   droplevels()
# 
# df_survival = create_df_survival(df_cohort_test,
#   event_col = "covid_mcoa_hosp", ind_vars = cox_ind_vars_minimal, calendar_days = 0, study_start = study_start, study_end = study_end
# )
# 
# # Check everything is ok - verify that bob is indeed one's uncle
# bob = select(
#   df_survival, individual_id, age_3cat,
# 
#   # Vaccination variables
#   date_vacc_1, date_vacc_2, date_vacc_3, date_vacc_4, date_vacc_5,
#   vacc_type_1, vacc_type_2, vacc_type_3, vacc_type_4, vacc_type_5,
#   num_doses_start, num_doses_recent, vacc_seq_start, vacc_seq_recent,
#   vs_start, vs_recent, under_vaccinated, vs_mixed_start, date_transfer_out, hosp_date,
#   date_fully_vaccinated,
# 
#   # Event dates
#   death_date, covid_death_date, covid_acoa_hosp_date, covid_mcoa_hosp_date,
#   covid_mcoa_28_2_hosp_date, covid_mcoa_14_2_hosp_date, covid_mcoa_hosp_death_date,
#   covid_acoa_hosp_death_date,
# 
#   covid_acoa_emerg_hosp_date, covid_mcoa_emerg_hosp_date,
#   covid_mcoa_28_2_emerg_hosp_date, covid_mcoa_14_2_emerg_hosp_date, covid_mcoa_emerg_hosp_death_date,
#   covid_acoa_emerg_hosp_death_date,
# 
# 
#   non_covid_acoa_hosp_date, non_covid_mcoa_hosp_date,
#   non_covid_mcoa_28_2_hosp_date, non_covid_mcoa_14_2_hosp_date, non_covid_mcoa_hosp_death_date,
#   non_covid_acoa_hosp_death_date,
# 
#   non_covid_acoa_emerg_hosp_date, non_covid_mcoa_emerg_hosp_date,
#   non_covid_mcoa_28_2_emerg_hosp_date, non_covid_mcoa_14_2_emerg_hosp_date, non_covid_mcoa_emerg_hosp_death_date,
#   non_covid_acoa_emerg_hosp_death_date,
# 
# 
#   competing_hosp, event, surv_date,
#   tstart, tstop, duration,  num_doses_time_2
# )
# 
# formula = as.formula(paste0("Surv(tstart, tstop, event) ~ ", paste(cox_ind_vars_extended, collapse = " + ")))
# 
# # Weighted Cox model
# #
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
# #
# # df_survival = mutate(df_survival, eave_weight = eave_weight * nrow(df_survival) / sum(eave_weight))
# #
# # model = coxph(formula, data = df_survival, weights = eave_weight)
# 
# model = coxph(formula, data = df_survival, model = TRUE)
# 
# results_table = create_cox_results_table(df_survival, model, cox_ind_vars_full)
# 
# cox_analysis(df_cohort_test,
#   age_group = "16-74",
#   dep_var = "covid_mcoa_hosp",
#   ind_vars = cox_ind_vars_extended,
#   calendar_days = 0,
#   study_start,
#   study_end,
#   exclude_non_unit_weight = FALSE,
#   weight = FALSE,
#   name = "extended"
# )
