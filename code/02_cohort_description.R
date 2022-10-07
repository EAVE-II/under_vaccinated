######################################################################
## Title: [Insert full title of paper]
## Code author: Steven Kerr steven.kerr@ed.ac.uk
## Description: Create cohort description tables.
##              01_data_setup is sourced.
######################################################################

#Libraries
library(tidyverse)
library(survival)
library(lubridate)
library(finalfit)

setwd('/conf/EAVE/GPanalysis/analyses/under_vaccinated')

#source('./code/01_data_setup.R')

# Create output directory
output_dir = "./output/"
if (!dir.exists(output_dir)) {dir.create(output_dir)}


rgs = setdiff(q_names, c('Q_BMI', 'Q_ETHNICITY'))
bin_vars = setdiff(rgs, c('Q_HOME_CAT', 'Q_LEARN_CAT', 'Q_DIAG_CKD_LEVEL'))

# Names of variables that will be displayed in table
var_names = c("Total N (%)",
              "ageYear", 
              "age_gp",  
              "ur6_2016_name", 
              "n_risk_gps",
              "ave_hh_age", 
              "n_hh_gp", 
              "bmi_cat", 
              'EAVE_Smoke',
              'EAVE_BP',
              rgs)

# Display names for table qcovid groups
q_display_names = c('Housing category',
                    'Learning disability/Down\'s syndrome',
                    'Chronic Kidney Disease', 
                    'Atrial Fibrillation', 
                    'Asthma', 
                    'Blood cancer', 
                    'Congestive Cardiac Failure', 
                    'Cerebral Palsy', 
                    'Coronary heart disease', 
                    'Liver cirrhosis', 
                    'Congenital heart disease', 
                    'COPD',
                    'Dementia', 
                    'Diabetes Type 1', 
                    'Diabetes Type 2', 
                    'Epilepsy', 
                    'Hip, wrist, spine, humerus fracture',
                    'HIV/AIDS',
                    'Severe combine immunodeficiency',
                    'Neurological conditions', 
                    'Parkinson’s', 
                    'Pulmonary hypertension', 
                    'Cystic fibrosis, bronchiectasis or alveolitis', 
                    'Peripheral vascular disease', 
                    'SLE or rheumatoid arthritis', 
                    'Lung, oral cancer', 
                    'Severe mental illness', 
                    'Sickle cell disease or combined immune deficiency syndrome', 
                    'Stroke, transient ischaemic attack', 
                    'Venous thromboembolism')

# Display names for table characteristics
display_names = c("Total",
                  "Age", 
                  "Age group",  
                  "Urban rural classification", 
                  "Number of risk groups",
                  "Average household age", 
                  "Number of people in household", 
                  "BMI", 
                  'Smoking status',
                  'Blood pressure',
                  q_display_names)

summary_tbl_wt_chrt <- summary_factorlist(df_cohort %>% 
                                            mutate(fully_vaccinated = factor(fully_vaccinated, 
                                              labels = c('Under-vaccinated', 'Fully-vaccinated'))) %>%
                                            # Put 1 as the first level, to make it easier to remove the level 0
                                            # in the final table
                                            mutate_at(bin_vars, ~factor(., levels = c(1,0)) ),
                                          dependent = "fully_vaccinated", 
                                          explanatory = setdiff(var_names, 'Total N (%)'),
                                          column = FALSE,
                                          add_col_totals = TRUE,
                                          weights = 'eave_weight')


# Only display one level for binary variables
summary_tbl_wt_chrt = filter(summary_tbl_wt_chrt, !(label =='' & levels == 0)) %>%
  mutate(levels = case_when(levels == 1 & label %in% bin_vars ~ '',
                            TRUE ~ levels))

# Change to display names
map = setNames(display_names, var_names)
summary_tbl_wt_chrt$label = map[summary_tbl_wt_chrt$label]

write.csv(summary_tbl_wt_chrt ,paste0(output_dir, '/summary_table_weights_cohort.csv'), row.names = F)




# Summary table with more detailed vaccination status categories
summary_tbl_wt_chrt_2 <- summary_factorlist(df_cohort %>%
                                            # Put 1 as the first level, to make it easier to remove the level 0
                                            # in the final table
                                            mutate_at(bin_vars, ~factor(., levels = c(1,0)) ),
                                          dependent = "vs_mixed", 
                                          explanatory = setdiff(var_names, 'Total N (%)'),
                                          column = FALSE,
                                          add_col_totals = TRUE,
                                          weights = 'eave_weight')


# Only display one level for binary variables
summary_tbl_wt_chrt_2 = filter(summary_tbl_wt_chrt_2, !(label =='' & levels == 0)) %>%
  mutate(levels = case_when(levels == 1 & label %in% bin_vars ~ '',
                            TRUE ~ levels))

# Change to display names
map = setNames(display_names, var_names)
summary_tbl_wt_chrt_2$label = map[summary_tbl_wt_chrt_2$label]

write.csv(summary_tbl_wt_chrt_2 ,paste0(output_dir, '/summary_table_weights_cohort_2.csv'), row.names = F)


