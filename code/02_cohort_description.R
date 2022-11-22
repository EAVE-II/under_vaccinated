######################################################################
## Title: [Insert full title of paper]
## Code author: Steven Kerr steven.kerr@ed.ac.uk
## Description: Create cohort description tables.
##              01_data_setup is sourced.
######################################################################

# Libraries
library(tidyverse)
library(lubridate)
library(finalfit)
library(dtplyr)
library(janitor)
library(survival)
library(broom)
library(scales)

setwd("/conf/EAVE/GPanalysis/analyses/under_vaccinated")

# source('./code/01_data_setup.R')

# Create output directory
output_dir ="./output/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# List of QCovid risk group variable names
rgs = setdiff(q_names, "Q_BMI")

# List of binary variables that will appear in summary tables
bin_vars = setdiff(rgs, c("Q_HOME_CAT", "Q_LEARN_CAT", "Q_DIAG_CKD_LEVEL", "Q_ETHNICITY"))

# Names of variables that will be displayed in table
var_names =c(
  "Total N (%)",
  "sex",
  "age",
  "age_group",
  "simd2020_sc_quintile",
  "urban_rural_class",
  "mean_household_age",
  "num_ppl_household",
  "n_risk_gps",
  "bmi_cat",
  rgs,
  "smoking_status",
  "blood_pressure"
)

# Display names for Qcovid groups in the table
q_display_names =c(
  "Ethnicity",
  "Housing category",
  "Learning disability/Down's syndrome",
  "Chronic Kidney Disease",
  "Atrial Fibrillation",
  "Asthma",
  "Blood cancer",
  "Congestive Cardiac Failure",
  "Cerebral Palsy",
  "Coronary heart disease",
  "Liver cirrhosis",
  "Congenital heart disease",
  "COPD",
  "Dementia",
  "Diabetes Type 1",
  "Diabetes Type 2",
  "Epilepsy",
  "Hip, wrist, spine, humerus fracture",
  "HIV/AIDS",
  "Severe combine immunodeficiency",
  "Neurological conditions",
  "Parkinson’s",
  "Pulmonary hypertension",
  "Cystic fibrosis, bronchiectasis or alveolitis",
  "Peripheral vascular disease",
  "SLE or rheumatoid arthritis",
  "Lung, oral cancer",
  "Severe mental illness",
  "Sickle cell disease or combined immune deficiency syndrome",
  "Stroke, transient ischaemic attack",
  "Venous thromboembolism"
)

# All display names for the table
display_names =c(
  "Total",
  "Sex",
  "Age",
  "Age group",
  "Socioeconomic status",
  "Urban rural classification",
  "Mean household age",
  "Number of people in household",
  "Number of risk groups",
  "BMI",
  q_display_names,
  "Smoking status",
  "Blood pressure"
)

# This will be used for changing variables names to display names
names_map = setNames(display_names, var_names)


# Fully vaccinated table
#
# summary_tbl_wt_chrt =summary_factorlist(
#   df_cohort %>%
#     mutate(fully_vaccinated = factor(fully_vaccinated,
#       labels = c("Under-vaccinated", "Fully-vaccinated")
#     )) %>%
#     # Put 1 as the first level, to make it easier to remove the level 0
#     # in the final table
#     mutate_at(bin_vars, ~ factor(., levels = c(1, 0))),
#   dependent = "fully_vaccinated",
#   explanatory = setdiff(var_names, "Total N (%)"),
#   add_col_totals = TRUE,
#   weights = "eave_weight",
#   na_include = TRUE
# )
# 
# 
# # Only display one level for binary variables
# summary_tbl_wt_chrt = summary_tbl_wt_chrt %>%
#   mutate(label_dup = label) %>%
#   mutate(label_dup = case_when(
#     label_dup == "" ~ NA_character_,
#     TRUE ~ label_dup
#   )) %>%
#   fill(label_dup, .direction = "down") %>%
#   filter(!(label_dup %in% bin_vars) | (label_dup %in% bin_vars & levels == "1")) %>%
#   select(-label_dup) %>%
#   mutate(levels = case_when(
#     levels == 1 & label %in% bin_vars ~ "",
#     TRUE ~ levels
#   ))
# 
# # Change to display names
# summary_tbl_wt_chrt$label =names_map[summary_tbl_wt_chrt$label]
# summary_tbl_wt_chrt$label =replace_na(summary_tbl_wt_chrt$label, "")
# 
# write.csv(summary_tbl_wt_chrt, paste0(output_dir, "/summary_table_weights_cohort.csv"), row.names = F)



# Summary table with more detailed vaccination status categories
summary_tbl_wt_chrt = summary_factorlist(df_cohort %>%
                                            # Put 1 as the first level, to make it easier to remove the level 0
                                            # in the final table
                                            mutate_at(bin_vars, ~factor(., levels = c(1,0)) ),
                                          dependent = "vs_mixed_start",
                                          explanatory = setdiff(var_names, 'Total N (%)'),
                                          add_col_totals = TRUE,
                                          weights = 'eave_weight',
                                          na_include = TRUE)


# Only display one level for binary variables
summary_tbl_wt_chrt$label[summary_tbl_wt_chrt$label == ''] = NA

summary_tbl_wt_chrt = summary_tbl_wt_chrt %>%
  mutate(label_dup = label) %>%
  mutate(label_dup = case_when(label_dup == '' ~ NA_character_,
                               TRUE ~ label_dup)) %>%
  fill(label_dup, .direction = 'down') %>%
  fill(label_dup, .direction = 'down') %>%
  filter(!(label_dup %in% bin_vars) |  (label_dup %in% bin_vars & levels == '1') ) %>%
  select(-label_dup) %>%
  mutate(levels = case_when(levels == 1 & label %in% bin_vars ~ '',
                            TRUE ~ levels))

# Change to display names
summary_tbl_wt_chrt$label = names_map[summary_tbl_wt_chrt$label]
summary_tbl_wt_chrt$label = replace_na(summary_tbl_wt_chrt$label, '')

write.csv(summary_tbl_wt_chrt ,paste0(output_dir, '/summary_table_weights_cohort.csv'), row.names = F)





#### Summary tables for each age group

# 5-11
summary_tbl_wt_chrt_5_11 =summary_factorlist(
  df_cohort %>%
    filter(age_group == "5-11") %>%
    # Put 1 as the first level, to make it easier to remove the level 0
    # in the final table
    mutate_at(bin_vars, ~ factor(., levels = c(1, 0))) %>%
    droplevels(),
  dependent = "num_doses_start",
  explanatory = setdiff(var_names, c("Total N (%)", "age_group")),
  add_col_totals = TRUE,
  weights = "eave_weight",
  na_include = TRUE
)


# Only display one level for binary variables
summary_tbl_wt_chrt_5_11$label[summary_tbl_wt_chrt_5_11$label == ""] =NA

summary_tbl_wt_chrt_5_11 =summary_tbl_wt_chrt_5_11 %>%
  mutate(label_dup = label) %>%
  mutate(label_dup = case_when(
    label_dup == "" ~ NA_character_,
    TRUE ~ label_dup
  )) %>%
  fill(label_dup, .direction = "down") %>%
  fill(label_dup, .direction = "down") %>%
  filter(!(label_dup %in% bin_vars) | (label_dup %in% bin_vars & levels == "1")) %>%
  select(-label_dup) %>%
  mutate(levels = case_when(
    levels == 1 & label %in% bin_vars ~ "",
    TRUE ~ levels
  ))

# Change to display names
summary_tbl_wt_chrt_5_11$label =names_map[summary_tbl_wt_chrt_5_11$label]
summary_tbl_wt_chrt_5_11$label =replace_na(summary_tbl_wt_chrt_5_11$label, "")

write.csv(summary_tbl_wt_chrt_5_11, paste0(output_dir, "/summary_table_weights_cohort_5_11.csv"), row.names = F)



# 12-15
summary_tbl_wt_chrt_12_15 =summary_factorlist(
  df_cohort %>%
    filter(age_group == "12-15") %>%
    # Put 1 as the first level, to make it easier to remove the level 0
    # in the final table
    mutate_at(bin_vars, ~ factor(., levels = c(1, 0))) %>%
    droplevels(),
  dependent = "num_doses_start",
  explanatory = setdiff(var_names, c("Total N (%)", "age_group")),
  add_col_totals = TRUE,
  weights = "eave_weight",
  na_include = TRUE
)


# Only display one level for binary variables
summary_tbl_wt_chrt_12_15$label[summary_tbl_wt_chrt_12_15$label == ""] =NA

summary_tbl_wt_chrt_12_15 =summary_tbl_wt_chrt_12_15 %>%
  mutate(label_dup = label) %>%
  mutate(label_dup = case_when(
    label_dup == "" ~ NA_character_,
    TRUE ~ label_dup
  )) %>%
  fill(label_dup, .direction = "down") %>%
  fill(label_dup, .direction = "down") %>%
  filter(!(label_dup %in% bin_vars) | (label_dup %in% bin_vars & levels == "1")) %>%
  select(-label_dup) %>%
  mutate(levels = case_when(
    levels == 1 & label %in% bin_vars ~ "",
    TRUE ~ levels
  ))

# Change to display names
summary_tbl_wt_chrt_12_15$label =names_map[summary_tbl_wt_chrt_12_15$label]
summary_tbl_wt_chrt_12_15$label =replace_na(summary_tbl_wt_chrt_12_15$label, "")

write.csv(summary_tbl_wt_chrt_12_15, paste0(output_dir, "/summary_table_weights_cohort_12_15.csv"), row.names = F)



# 16-74
summary_tbl_wt_chrt_16_74 =summary_factorlist(
  df_cohort %>%
    filter(age_group == "16-74") %>%
    # Put 1 as the first level, to make it easier to remove the level 0
    # in the final table
    mutate_at(bin_vars, ~ factor(., levels = c(1, 0))) %>%
    droplevels(),
  dependent = "num_doses_start",
  explanatory = setdiff(var_names, c("Total N (%)", "age_group")),
  add_col_totals = TRUE,
  weights = "eave_weight",
  na_include = TRUE
)


# Only display one level for binary variables
summary_tbl_wt_chrt_16_74$label[summary_tbl_wt_chrt_16_74$label == ""] =NA

summary_tbl_wt_chrt_16_74 = summary_tbl_wt_chrt_16_74 %>%
  mutate(label_dup = label) %>%
  mutate(label_dup = case_when(
    label_dup == "" ~ NA_character_,
    TRUE ~ label_dup
  )) %>%
  fill(label_dup, .direction = "down") %>%
  fill(label_dup, .direction = "down") %>%
  filter(!(label_dup %in% bin_vars) | (label_dup %in% bin_vars & levels == "1")) %>%
  select(-label_dup) %>%
  mutate(levels = case_when(
    levels == 1 & label %in% bin_vars ~ "",
    TRUE ~ levels
  ))

# Change to display names
summary_tbl_wt_chrt_16_74$label = names_map[summary_tbl_wt_chrt_16_74$label]
summary_tbl_wt_chrt_16_74$label = replace_na(summary_tbl_wt_chrt_16_74$label, "")

write.csv(summary_tbl_wt_chrt_16_74, paste0(output_dir, "/summary_table_weights_cohort_16_74.csv"), row.names = F)




# 75+
summary_tbl_wt_chrt_over_75 =summary_factorlist(
  df_cohort %>%
    filter(age_group == "75+") %>%
    # Put 1 as the first level, to make it easier to remove the level 0
    # in the final table
    mutate_at(bin_vars, ~ factor(., levels = c(1, 0))) %>%
    droplevels(),
  dependent = "num_doses_start",
  explanatory = setdiff(var_names, c("Total N (%)", "age_group")),
  add_col_totals = TRUE,
  weights = "eave_weight",
  na_include = TRUE
)


# Only display one level for binary variables
summary_tbl_wt_chrt_over_75 =summary_tbl_wt_chrt_over_75 %>%
  mutate(label_dup = label) %>%
  mutate(label_dup = case_when(
    label_dup == "" ~ NA_character_,
    TRUE ~ label_dup
  )) %>%
  fill(label_dup, .direction = "down") %>%
  filter(!(label_dup %in% bin_vars) | (label_dup %in% bin_vars & levels == "1")) %>%
  select(-label_dup) %>%
  mutate(levels = case_when(
    levels == 1 & label %in% bin_vars ~ "",
    TRUE ~ levels
  ))

# Change to display names
summary_tbl_wt_chrt_over_75$label = names_map[summary_tbl_wt_chrt_over_75$label]
summary_tbl_wt_chrt_over_75$label = replace_na(summary_tbl_wt_chrt_over_75$label, "")

write.csv(summary_tbl_wt_chrt_over_75, paste0(output_dir, "/summary_table_weights_cohort_over_75.csv"), row.names = F)



#### Vaccine uptake visualiations
# Stu's code

## Vaccinations by age group and week
d_study_weeks = seq(
  from = floor_date(study_start, "week"),
  to   = floor_date(study_end, "week"),
  by   = "1 week"
)

d_dummy_weeks = expand_grid(
  age_group  = sort(unique(df_cohort$age_group)),
  dose_num  = c("1", "2", "3", "4"),
  dose_week = d_study_weeks
)

d_vacc_week =
  df_cohort %>%
  select(
    age_group,
    `1` = date_vacc_1,
    `2` = date_vacc_2,
    `3` = date_vacc_3,
    `4` = date_vacc_4
  ) %>%
  pivot_longer(
    cols           = c('1', '2', '3', '4'),
    names_to       = "dose_num",
    values_to      = "dose_date",
    values_drop_na = TRUE
  ) %>%
  mutate(dose_week = floor_date(dose_date, "week")) %>%
  filter(dose_week <= study_end) %>%
  count(age_group, dose_num, dose_week) %>%
  full_join(d_dummy_weeks, by = c("age_group", "dose_num", "dose_week")) %>%
  mutate(n = replace_na(n, as.integer(0))) %>%
  arrange(age_group, dose_num, dose_week)


d_vacc_week %>%
  # Suppression
  # Alter this depending on what the rules for your TRE are
  mutate(n = if_else(1 <= n & n < 5, as.integer(3), n)) %>%
  ggplot(aes(x = dose_week, y = n, colour = dose_num)) +
  facet_wrap(~ age_group, ncol = 1, scales = "free_y") +
  geom_line() + 
  xlab('') + 
  ylab('') + 
  labs(color ='Dose') +
  scale_y_continuous(labels = scales::comma) +
  xlab('Date') + 
  ylab('Count')


ggsave(paste0(output_dir, 'vacc_count_by_week_dose_age_group.png'), height = 10)


## Time between vaccinations
d_vacc_diff =
  df_cohort %>%
  select(
    individual_id,
    age_group,
    dose1    = date_vacc_1,
    dose2    = date_vacc_2,
    dose3    = date_vacc_3,
    dose4    = date_vacc_4
  ) %>%
  # remove events after study window
  mutate(across(
    .cols = where(is.Date),
    .fns  = ~ if_else(.x <= study_end, .x, NA_Date_)
  )) %>%
  mutate(
    # weeks between doses
    `Dose 1 to 2` = interval(dose1, dose2) / dweeks(),
    `Dose 2 to 3` = interval(dose2, dose3) / dweeks(),
    `Dose 3 to 4` = interval(dose3, dose4) / dweeks(),
    # round off
    `Dose 1 to 2` = floor(`Dose 1 to 2`),
    `Dose 2 to 3` = floor(`Dose 2 to 3`),
    `Dose 3 to 4` = floor(`Dose 3 to 4`)
  ) %>%
  select(
    age_group,
    `Dose 1 to 2`,
    `Dose 2 to 3`,
    `Dose 3 to 4`
  ) %>%
  pivot_longer(
    cols           = -age_group,
    names_to       = "dose_to_dose",
    values_to      = "diff_weeks",
    values_drop_na = TRUE
  ) %>%
  count(age_group, dose_to_dose, diff_weeks)

d_vacc_diff %>%
  # Suppression
  # Alter this depending on what the rules for your TRE are
  mutate(n = if_else(1 <= n & n < 5, as.integer(3), n)) %>%
  ggplot(aes(
    x = diff_weeks,
    y = n,
  )) +
  facet_grid(age_group ~ dose_to_dose, scales = "free") +
  geom_col() +
  scale_y_continuous(labels = scales::comma) +
  xlab('Weeks') + 
  ylab('Count')

ggsave(paste0(output_dir, 'time_between_doses.png'), height = 10, width = 10)


### Uptake cumulative incidence
lkp_events = c(
  "Unvaccinated",
  "Dose 1",
  "Dose 2",
  "Dose 3",
  "Dose 4",
  "Death",
  "study_end"
)

# make multi-state date set for dose
d_mstate_vacc =
  df_cohort %>%
  mutate(
    Death = date_death,
    study_end = study_end
  ) %>%
  select(
    individual_id,
    age_group,
    `Dose 1`    = date_vacc_1,
    `Dose 2`    = date_vacc_2,
    `Dose 3`    = date_vacc_3,
    `Dose 4`    = date_vacc_4,
    Death,
    study_end
  ) %>%
  # any first day events just add half a day
  mutate(across(
    .cols = where(is.Date),
    .fns  = ~ if_else(.x == study_start, .x + ddays(0.5), as_datetime(.x))
  )) %>%
  # remove events after someone has either moved out, died, or reached the
  # end of the study window
  mutate(
    end_follow_up = pmin(Death, study_end, na.rm = TRUE)
  ) %>%
  mutate(across(
    .cols = where(is.POSIXct),
    .fns  = ~ if_else(.x <= end_follow_up, .x, NA_POSIXct_)
  )) %>%
  select(-end_follow_up) %>%
  pivot_longer(
    cols           = c(-individual_id, -age_group),
    names_to       = "event_name",
    values_to      = "event_date",
    values_drop_na = TRUE
  ) %>%
  # add row number for alf
  lazy_dt() %>%
  arrange(individual_id, event_date) %>%
  group_by(individual_id) %>%
  mutate(alf_seq = row_number()) %>%
  as_tibble() %>%
  # define survival columns
  mutate(
    tstart     = if_else(alf_seq == 1, as_datetime(study_start), lag(event_date)),
    tstart     = interval(as_datetime(study_start), tstart) / ddays(),
    state_from = if_else(alf_seq == 1, "Unvaccinated", lag(event_name)),
    state_from = factor(state_from, lkp_events) %>% fct_drop(),
    tstop      = event_date,
    tstop      = interval(as_datetime(study_start), tstop) / ddays(),
    state_to   = factor(event_name, lkp_events) %>% fct_drop()
  ) %>%
  select(
    -event_name,
    -event_date
  ) %>%
  # anyone who has event on last day, remove row for that transition from event to study end
  filter(!(tstart > 0 & tstart == tstop & state_to == "study_end")) %>%
  # any same day events that go from dose to death / move out just add half a day again
  mutate(
    tstop = if_else(
      condition = state_from %in% c("Dose 1", "Dose 2", "Dose 3", "Dose 4") &
        state_to  == "Death" &
        tstart == tstop,
      true = tstop + 0.5,
      false = tstop
    )
  ) %>%
  # finalise censored category
  mutate(
    state_to = state_to %>%
      fct_collapse("(censored)" = "study_end") %>%
      fct_relevel("(censored)")
  )

# explore
d_mstate_vacc %>% tabyl(state_from, state_to)

# fit mstate
mstate_vacc = survfit(
  formula = Surv(tstart, tstop, state_to) ~ age_group,
  data = d_mstate_vacc,
  id = individual_id
)


mstate_vacc %>%
  tidy() %>%
  mutate(
    strata = strata %>%
      str_replace("age_group=", "") %>%
      factor() %>%
      fct_relevel("5-11", "12-15", "16-74", "75+"),
    State = state %>%
      factor() %>%
      fct_recode("Unvaccinated" = "(s0)") %>%
      fct_relevel("Death")
  ) %>%
  # Convert time to date
  mutate(
  time = study_start + time
  ) %>%
  ggplot(aes(
    x = time,
    y = estimate,
    fill = State
  )) +
  facet_wrap(~strata, ncol = 1) +
  geom_area() +
  scale_fill_brewer( type = "qual") +
  ylab('Cumulative incidence') +
  xlab('Day')

ggsave(paste0(output_dir, 'uptake_cumulative_incidence.png'), height = 10)
