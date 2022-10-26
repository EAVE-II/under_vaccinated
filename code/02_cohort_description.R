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

setwd("/conf/EAVE/GPanalysis/analyses/under_vaccinated")

# source('./code/01_data_setup.R')

# Create output directory
output_dir ="./output/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# List of QCovid risk group variable names
rgs =setdiff(q_names, c("Q_BMI", "Q_ETHNICITY"))

# List of binary variables that will appear in summary tables
bin_vars =setdiff(rgs, c("Q_HOME_CAT", "Q_LEARN_CAT", "Q_DIAG_CKD_LEVEL"))

# Names of variables that will be displayed in table
var_names =c(
  "Total N (%)",
  "Sex",
  "ageYear",
  "age_gp",
  "simd2020_sc_quintile",
  "ur6_2016_name",
  "n_risk_gps",
  "ave_hh_age",
  "n_hh_gp",
  "bmi_cat",
  "EAVE_Smoke",
  "EAVE_BP",
  rgs
)

# Display names for Qcovid groups in the table
q_display_names =c(
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
  "Number of risk groups",
  "Average household age",
  "Number of people in household",
  "BMI",
  "Smoking status",
  "Blood pressure",
  q_display_names
)

# This will be used for changing variables names to display names
names_map =setNames(display_names, var_names)

summary_tbl_wt_chrt =summary_factorlist(
  df_cohort %>%
    mutate(fully_vaccinated = factor(fully_vaccinated,
      labels = c("Under-vaccinated", "Fully-vaccinated")
    )) %>%
    # Put 1 as the first level, to make it easier to remove the level 0
    # in the final table
    mutate_at(bin_vars, ~ factor(., levels = c(1, 0))),
  dependent = "fully_vaccinated",
  explanatory = setdiff(var_names, "Total N (%)"),
  add_col_totals = TRUE,
  weights = "eave_weight",
  na_include = TRUE
)


# Only display one level for binary variables
summary_tbl_wt_chrt = summary_tbl_wt_chrt %>%
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
summary_tbl_wt_chrt$label =names_map[summary_tbl_wt_chrt$label]
summary_tbl_wt_chrt$label =replace_na(summary_tbl_wt_chrt$label, "")

write.csv(summary_tbl_wt_chrt, paste0(output_dir, "/summary_table_weights_cohort.csv"), row.names = F)






# # Summary table with more detailed vaccination status categories
# summary_tbl_wt_chrt_2 = summary_factorlist(df_cohort %>%
#                                             # Put 1 as the first level, to make it easier to remove the level 0
#                                             # in the final table
#                                             mutate_at(bin_vars, ~factor(., levels = c(1,0)) ),
#                                           dependent = "vs_mixed_start",
#                                           explanatory = setdiff(var_names, 'Total N (%)'),
#                                           add_col_totals = TRUE,
#                                           weights = 'eave_weight',
#                                           na_include = TRUE)
#
#
# # Only display one level for binary variables
# summary_tbl_wt_chrt_2$label[summary_tbl_wt_chrt_2$label == ''] = NA
#
# summary_tbl_wt_chrt_2 = summary_tbl_wt_chrt_2 %>%
#   mutate(label_dup = label) %>%
#   mutate(label_dup = case_when(label_dup == '' ~ NA_character_,
#                                TRUE ~ label_dup)) %>%
#   fill(label_dup, .direction = 'down') %>%
#   fill(label_dup, .direction = 'down') %>%
#   filter(!(label_dup %in% bin_vars) |  (label_dup %in% bin_vars & levels == '1') ) %>%
#   select(-label_dup) %>%
#   mutate(levels = case_when(levels == 1 & label %in% bin_vars ~ '',
#                             TRUE ~ levels))
#
# # Change to display names
# summary_tbl_wt_chrt_2$label = names_map[summary_tbl_wt_chrt_2$label]
# summary_tbl_wt_chrt_2$label = replace_na(summary_tbl_wt_chrt_2$label, '')
#
# write.csv(summary_tbl_wt_chrt_2 ,paste0(output_dir, '/summary_table_weights_cohort_2.csv'), row.names = F)




# Summary tables by age

# 5-15
summary_tbl_wt_chrt_5_15 =summary_factorlist(
  df_cohort %>%
    filter(age_gp == "5-15") %>%
    # Put 1 as the first level, to make it easier to remove the level 0
    # in the final table
    mutate_at(bin_vars, ~ factor(., levels = c(1, 0))) %>%
    droplevels(),
  dependent = "num_doses_start",
  explanatory = setdiff(var_names, c("Total N (%)", "age_gp")),
  add_col_totals = TRUE,
  weights = "eave_weight",
  na_include = TRUE
)


# Only display one level for binary variables
summary_tbl_wt_chrt_5_15$label[summary_tbl_wt_chrt_5_15$label == ""] =NA

summary_tbl_wt_chrt_5_15 =summary_tbl_wt_chrt_5_15 %>%
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
summary_tbl_wt_chrt_5_15$label =names_map[summary_tbl_wt_chrt_5_15$label]
summary_tbl_wt_chrt_5_15$label =replace_na(summary_tbl_wt_chrt_5_15$label, "")

write.csv(summary_tbl_wt_chrt_5_15, paste0(output_dir, "/summary_table_weights_cohort_5_15.csv"), row.names = F)







# 16-74
summary_tbl_wt_chrt_16_74 =summary_factorlist(
  df_cohort %>%
    filter(age_gp == "16-74") %>%
    # Put 1 as the first level, to make it easier to remove the level 0
    # in the final table
    mutate_at(bin_vars, ~ factor(., levels = c(1, 0))) %>%
    droplevels(),
  dependent = "num_doses_start",
  explanatory = setdiff(var_names, c("Total N (%)", "age_gp")),
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
    filter(age_gp == "75+") %>%
    # Put 1 as the first level, to make it easier to remove the level 0
    # in the final table
    mutate_at(bin_vars, ~ factor(., levels = c(1, 0))) %>%
    droplevels(),
  dependent = "num_doses_start",
  explanatory = setdiff(var_names, c("Total N (%)", "age_gp")),
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
