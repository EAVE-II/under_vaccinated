######################################################################
## Title: [Insert full title of paper]
## Code author: Steven Kerr steven.kerr@ed.ac.uk
## Description: Create study dataframe, adding derived variables
##              01_data_cleaning is sourced.
######################################################################

library(tidyverse)
library(lubridate)

Location <- "/conf/"

setwd("/conf/EAVE/GPanalysis/analyses/under_vaccinated")

#source("./code/01_data_cleaning.R")

# For reproducing any randomisation
set.seed(1234)

study_start = as.Date("2022-06-01")
study_end = as.Date("2022-09-30")

#### Create cohort dataframe with all required variables by merging togehter individual datasets

# Remove any EAVE_LINKNO which two old EAVE_LINKNOs are getting mapped onto
remove_EAVE_LINKNOs <- EAVE_LINKNO_refresh %>%
  filter(duplicated(EAVE_LINKNO)) %>%
  pull(EAVE_LINKNO)

df_cohort <- EAVE_LINKNO_refresh %>%
  #For testing
  #sample_n(10000) %>%
  select(EAVE_LINKNO_old, 
         EAVE_LINKNO,
         EAVE_LINKNO_change,
         chili_checked,
         unvalidatedCHI_flag,
         date_transfer_out = DATE_TRANSFER_OUT,
         date_transfer_in = DATE_TRANSFER_IN
  ) %>%
  filter(unvalidatedCHI_flag!= 1) %>%
  filter(!EAVE_LINKNO %in% remove_EAVE_LINKNOs) %>%
  filter(
    is.na(date_transfer_out) | date_transfer_out > study_start,
     is.na(date_transfer_in) | date_transfer_in <= study_start)
 
# Add deaths from any cause, and remove people who died before study_start
df_cohort <- df_cohort %>%
  left_join(all_deaths, by = "EAVE_LINKNO") %>%
  mutate(covid_death_date = ifelse(covid_death == 1, death_date, as.Date(NA)) %>%
           as.Date(origin = "1970-01-01")) %>%
  filter(is.na(death_date) | death_date >= study_start)

df_cohort <- df_cohort %>% 
  # Join by EAVE_LINKNO_old, since data pre-dates December 2020
  left_join(EAVE_cohort, by = c("EAVE_LINKNO_old" = "EAVE_LINKNO")) %>%
  filter(age >= 5) %>%
  mutate(
    age_3cat = cut(age,
      breaks = c(5, 16, 75, Inf),
      labels = c("5-15", "16-74", "75+"),
      right = FALSE
    ),
    age_4cat = cut(age,
                   breaks = c(5, 12, 16, 75, Inf),
                   labels = c("5-11", "12-15", "16-74", "75+"),
                   right = FALSE
    ),
    age_17cat = cut(age,
      breaks = c(5, 12, 16, 18, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, Inf),
      labels = c(
        "5-11", "12-15", "16-17", "18-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54",
        "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85+"
      ),
      right = FALSE
    ) %>%
      fct_relevel("18-24")
  ) %>%
  # Join by EAVE_LINKNO_old, since data pre-dates December 2020
  left_join(EAVE_Weights, by = c("EAVE_LINKNO_old" = "EAVE_LINKNO"))

# Add health board
df_cohort <- df_cohort %>%
  left_join(dz_hb_lookup, by = 'data_zone')

# Add imputed BMI
df_cohort <- df_cohort %>%
  left_join(imputation, by = "EAVE_LINKNO")

df_cohort$eave_weight[is.na(df_cohort$eave_weight)] <- mean(df_cohort$eave_weight, na.rm = T)

# Add risk groups
df_cohort <- qcovid_feb22 %>%
  select(-n_risk_gps) %>%
  right_join(df_cohort, by = "EAVE_LINKNO")

# Join by EAVE_LINKNO_old, since data pre-dates December 2021
df_cohort <- left_join(
  df_cohort,
  qcovid_old %>% 
    select(EAVE_LINKNO, Q_DIAG_DIABETES_1, Q_DIAG_DIABETES_2), 
  by = c("EAVE_LINKNO_old" = "EAVE_LINKNO"))

# If they aren't present in the new qcovid extract, take their data from the old qcovid extract  
df_cohort <- left_join(
  df_cohort,
  qcovid_old %>% 
    filter(!EAVE_LINKNO %in% qcovid_feb22$EAVE_LINKNO) %>%
    rename(EAVE_LINKNO_old = EAVE_LINKNO)
  )

df_cohort <- df_cohort %>%
  mutate(
    bmi_cat = cut(
      Q_BMI,
      breaks = c(10, 18.5, 25, 30, 35, 40, Inf),
      labels = c("10-18.5", "18.5-25", "25-30", "30-35", "35-40", "40+"),
      include.lowest = TRUE,
      right = FALSE
    ) %>%
      fct_relevel('25-30'),
    bmi_imputed_cat = cut(
      Q_BMI_imputed,
      breaks = c(10, 18.5, 25, 30, 35, 40, Inf),
      labels = c("10-18.5", "18.5-25", "25-30", "30-35", "35-40", "40+"),
      include.lowest = TRUE,
      right = FALSE
    )%>%
      fct_relevel('25-30')
  ) %>%
  mutate(
    Q_HOME_CAT = ifelse(is.na(Q_HOME_CAT), "Neither", Q_HOME_CAT) %>%
      factor() %>%
      fct_relevel("Neither"),
    Q_LEARN_CAT = ifelse(is.na(Q_LEARN_CAT), "Neither", Q_LEARN_CAT) %>%
      factor() %>%
      fct_relevel("Neither"),
    Q_DIAG_CKD_LEVEL = ifelse(is.na(Q_DIAG_CKD_LEVEL), "No CKD", Q_DIAG_CKD_LEVEL) %>%
      factor() %>%
      fct_relevel("No CKD")
  ) 
  
# Join by EAVE_LINKNO_old, since data pre-dates December 2021
df_cohort <- left_join(df_cohort, eave_rg, by = c("EAVE_LINKNO_old" = "EAVE_LINKNO"))

# Add ethnicity
df_cohort <- left_join(df_cohort, ethnicity, by = "EAVE_LINKNO") %>%
  mutate(
    ethnicity_5cat = replace_na(ethnicity_5cat, "Unknown"),
    ethnicity_18cat = replace_na(ethnicity_18cat, "Unknown")
    )

# Add Vaccinations, and vaccine related derived variables
df_cohort <- left_join(df_cohort, Vaccinations, by = "EAVE_LINKNO") %>%
  mutate(
    # If someone has e.g. a date_vacc_3, but no record of second dose, assume they got a second dose and set vacc_type_2 to Unknown
    vacc_type_1 = case_when(
      is.na(vacc_type_1) & (!is.na(date_vacc_2) | !is.na(date_vacc_3)| !is.na(date_vacc_4) | !is.na(date_vacc_5)) ~ 'Unk',
      TRUE ~ vacc_type_1),
    vacc_type_2 = case_when(
      is.na(vacc_type_2) & (!is.na(date_vacc_3) | !is.na(date_vacc_4) | !is.na(date_vacc_5)) ~ 'Unk',
      TRUE ~ vacc_type_2),
    vacc_type_3 = case_when(
      is.na(vacc_type_3) & (!is.na(date_vacc_4) | !is.na(date_vacc_5)) ~ 'Unk',
      TRUE ~ vacc_type_3),
    vacc_type_4 = case_when(
      is.na(vacc_type_4) &!is.na(date_vacc_5) ~ 'Unk',
      TRUE ~ vacc_type_4)
  ) %>%
  # Number of doses they have had according to most recent data, up to dose 5
  mutate(
    num_doses_recent = case_when(
      !is.na(date_vacc_5) ~ 5,
      !is.na(date_vacc_4) ~ 4,
      !is.na(date_vacc_3) ~ 3,
      !is.na(date_vacc_2) ~ 2,
      !is.na(date_vacc_1) ~ 1,
      TRUE ~ 0
    )
  ) %>%
  mutate(
    # Add winter booster dose as numbered vaccine dose
    # 19 day gap required between doses
    # Vaccine type
    vacc_type_5 = ifelse(is.na(date_vacc_5) & date_vacc_wb > date_vacc_4 + 19, vacc_type_wb, vacc_type_5),
    vacc_type_4 = ifelse(is.na(date_vacc_5) & is.na(date_vacc_4) & date_vacc_wb > date_vacc_3 + 19, vacc_type_wb, vacc_type_4),
    vacc_type_3 = ifelse(is.na(date_vacc_5) & is.na(date_vacc_4) & is.na(date_vacc_3) & date_vacc_wb > date_vacc_2 + 19,
                         vacc_type_wb, vacc_type_3),
    vacc_type_2 = ifelse(is.na(date_vacc_5) & is.na(date_vacc_4) & is.na(date_vacc_3) & is.na(date_vacc_2) & 
                           date_vacc_wb > date_vacc_1 + 19, vacc_type_wb, vacc_type_2),
    vacc_type_1 = ifelse(is.na(date_vacc_5) & is.na(date_vacc_4) & is.na(date_vacc_3) & is.na(date_vacc_2) & is.na(date_vacc_1) & 
                           !is.na(date_vacc_wb), vacc_type_wb, vacc_type_1),
    
    # Vaccine date
    date_vacc_5 = ifelse(is.na(date_vacc_5) & date_vacc_wb > date_vacc_4 + 19, date_vacc_wb, date_vacc_5) %>%
      as.Date(origin = "1970-01-01"),
    date_vacc_4 = ifelse(is.na(date_vacc_5) & is.na(date_vacc_4) & date_vacc_wb > date_vacc_3 + 19, date_vacc_wb, date_vacc_4) %>%
      as.Date(origin = "1970-01-01"),
    date_vacc_3 = ifelse(is.na(date_vacc_5) & is.na(date_vacc_4) & is.na(date_vacc_3) & date_vacc_wb > date_vacc_2 + 19,
      date_vacc_wb, date_vacc_3) %>%
      as.Date(origin = "1970-01-01"),
    date_vacc_2 = ifelse(is.na(date_vacc_5) & is.na(date_vacc_4) & is.na(date_vacc_3) & is.na(date_vacc_2) & 
      date_vacc_wb > date_vacc_1 + 19, date_vacc_wb, date_vacc_2) %>%
      as.Date(origin = "1970-01-01"),
    date_vacc_1 = ifelse(is.na(date_vacc_5) & is.na(date_vacc_4) & is.na(date_vacc_3) & is.na(date_vacc_2) & is.na(date_vacc_1) & 
      !is.na(date_vacc_wb), date_vacc_wb, date_vacc_1) %>%
      as.Date(origin = "1970-01-01")
  ) %>%
  mutate(
    # Number of doses at study_start, up to dose 5
    num_doses_start = case_when(
      date_vacc_5 <= study_start ~ 5,
      date_vacc_4 <= study_start ~ 4,
      date_vacc_3 <= study_start ~ 3,
      date_vacc_2 <= study_start ~ 2,
      date_vacc_1 <= study_start ~ 1,
      TRUE ~ 0
    ),
    # Sequence of vaccines according to most recent data, up to fifth dose
    vacc_seq_recent = paste0(
      ifelse(is.na(vacc_type_1), "", paste0(vacc_type_1, "_")),
      ifelse(is.na(vacc_type_2), "", paste0(vacc_type_2, "_")),
      ifelse(is.na(vacc_type_3), "", paste0(vacc_type_3, "_")),
      ifelse(is.na(vacc_type_4), "", paste0(vacc_type_4, "_")),
      ifelse(is.na(vacc_type_5), "", vacc_type_5)
    ),
    # Sequence of vaccines they had at study_start, up to fifth dose
    vacc_seq_start = paste0(
      ifelse(!is.na(date_vacc_1) & date_vacc_1 <= study_start, paste0(vacc_type_1, "_"), ""),
      ifelse(!is.na(date_vacc_2) & date_vacc_2 <= study_start, paste0(vacc_type_2, "_"), ""),
      ifelse(!is.na(date_vacc_3) & date_vacc_3 <= study_start, paste0(vacc_type_3, "_"), ""),
      ifelse(!is.na(date_vacc_4) & date_vacc_4 <= study_start, paste0(vacc_type_4, "_"), ""),
      ifelse(!is.na(vacc_type_5) & date_vacc_5 <= study_start, vacc_type_5, "")
    ),
    # Whether they had mixed vaccines according to most recent data, up to fifth dose
    mixed_vacc_recent = case_when(
      num_doses_recent == 0 ~ 0,
      num_doses_recent == 1 ~ 0,
      vacc_type_2 != vacc_type_1 ~ 1,
      vacc_type_3 != vacc_type_2 ~ 1,
      vacc_type_4 != vacc_type_3 ~ 1,
      vacc_type_5 != vacc_type_4 ~ 1,
      TRUE ~ 0
    ),
    # Whether they had mixed vaccines at study_start, up to fifth dose
    mixed_vacc_start = case_when(
      num_doses_start == 0 ~ 0,
      num_doses_start == 1 ~ 0,
      vacc_type_2 != vacc_type_1 & date_vacc_2 <= study_start ~ 1,
      vacc_type_3 != vacc_type_2 & date_vacc_3 <= study_start ~ 1,
      vacc_type_4 != vacc_type_3 & date_vacc_4 <= study_start ~ 1,
      vacc_type_5 != vacc_type_4 & date_vacc_5 <= study_start ~ 1,
      TRUE ~ 0
    )
  ) %>%
  mutate(
    vacc_seq_recent = trimws(vacc_seq_recent, whitespace = "_"),
    vacc_seq_start = trimws(vacc_seq_start, whitespace = "_")
  ) %>%
  mutate(
    vacc_seq_start = case_when(
      vacc_seq_start == "" ~ "uv",
      TRUE ~ vacc_seq_start
    ),
    vacc_seq_recent = case_when(
      vacc_seq_recent == "" ~ "uv",
      TRUE ~ vacc_seq_recent
    )
  )

df_cohort <- mutate(df_cohort,
  # Dose number and type of most recent vaccine dose, up to fifth dose
  vs_recent = case_when(
    !is.na(date_vacc_5) ~ paste0("v5_", vacc_type_5),
    !is.na(date_vacc_4) ~ paste0("v4_", vacc_type_4),
    !is.na(date_vacc_3) ~ paste0("v3_", vacc_type_3),
    !is.na(date_vacc_2) ~ paste0("v2_", vacc_type_2),
    !is.na(date_vacc_1) ~ paste0("v1_", vacc_type_1),
    TRUE ~ "uv"
  ),
  # Dose number and type of most recent vaccine dose as of study_start, up to fifth dose
  vs_start = case_when(
    date_vacc_5 <= study_start ~ paste0("v5_", vacc_type_5),
    date_vacc_4 <= study_start ~ paste0("v4_", vacc_type_4),
    date_vacc_3 <= study_start ~ paste0("v3_", vacc_type_3),
    date_vacc_2 <= study_start ~ paste0("v2_", vacc_type_2),
    date_vacc_1 <= study_start ~ paste0("v1_", vacc_type_1),
    TRUE ~ "uv"
  )
) %>%
  mutate(
    under_vaccinated = case_when(
      age_4cat == "5-11" & date_vacc_1 <= study_start ~ 0,
      age_4cat == "12-15" & date_vacc_2 <= study_start ~ 0,
      age_4cat == "16-74" & date_vacc_3 <= study_start ~ 0,
      age_4cat == "75+" & date_vacc_4 <= study_start ~ 0,
      TRUE ~ 1
    ),
    date_fully_vaccinated = case_when(
      age_4cat == "5-11" ~ date_vacc_1,
      age_4cat == "12-15" ~ date_vacc_2,
      age_4cat == "16-74" ~ date_vacc_3,
      age_4cat == "75+" ~ date_vacc_4,
      TRUE ~ as.Date(NA)
    )
  ) %>%
  mutate(under_vaccinated = factor(under_vaccinated))

df_cohort <- mutate(df_cohort,
  # Mixed vaccine status at study start
  vs_mixed_start = case_when(
    num_doses_start == 0 ~ "uv",
    num_doses_start == 1 ~ vacc_seq_start,
    num_doses_start == 2 & vacc_seq_start %in% c("AZ_AZ", "PB_PB", "Mo_Mo") ~ vacc_seq_start,
    num_doses_start == 2 & vacc_type_1 == "AZ" & vacc_type_2 %in% c("PB", "Mo") ~ "v2_AZ_mrna",
    num_doses_start == 2 & vacc_type_1 %in% c("PB", "Mo") & vacc_type_2 %in% c("PB", "Mo") ~ "v2_mixed_mrna",
    num_doses_start == 2 ~ "v2_other",
    num_doses_start == 3 & age_4cat %in% c("5-11", "12-15") ~ "3+",
    num_doses_start == 3 & vacc_seq_start %in%
      c("AZ_AZ_PB", "AZ_AZ_Mo", "AZ_AZ_AZ", "PB_PB_PB", "PB_PB_Mo", "Mo_Mo_PB", "Mo_Mo_Mo") ~ vacc_seq_start,
    num_doses_start == 3 & vacc_type_3 == "PB" ~ "other_mixed_2_dose_PB",
    num_doses_start == 3 & vacc_type_3 == "Mo" ~ "other_mixed_2_dose_Mo",
    num_doses_start == 3 ~ "v3_other",
    num_doses_start == 4 & age_4cat %in% c("5-11", "12-15", "16-74") ~ "4+",
    num_doses_start == 4 & vacc_type_1 == "AZ" & vacc_type_2 == "AZ" ~ "v4_AZ_AZ_any",
    num_doses_start == 4 & vacc_type_1 == "PB" & vacc_type_2 == "PB" ~ "v4_PB_PB_any",
    num_doses_start == 4 ~ "v4_other",
    num_doses_start > 4 ~ "5+"
  ),
  # Mixed vaccine status according to most recent data
  vs_mixed_recent = case_when(
    num_doses_recent == 0 ~ "uv",
    num_doses_recent == 1 ~ vacc_seq_recent,
    num_doses_recent == 2 & vacc_seq_recent %in% c("AZ_AZ", "PB_PB", "Mo_Mo") ~ vacc_seq_recent,
    num_doses_recent == 2 & vacc_type_1 == "AZ" & vacc_type_2 %in% c("PB", "Mo") ~ "v2_AZ_mrna",
    num_doses_recent == 2 & vacc_type_1 %in% c("PB", "Mo") & vacc_type_2 %in% c("PB", "Mo") ~ "v2_mixed_mrna",
    num_doses_recent == 2 ~ "v2_other",
    num_doses_recent == 3 & age_4cat %in% c("5-11", "12-15") ~ "3+",
    num_doses_recent == 3 & vacc_seq_recent %in%
      c("AZ_AZ_PB", "AZ_AZ_Mo", "AZ_AZ_AZ", "PB_PB_PB", "PB_PB_Mo", "Mo_Mo_PB", "Mo_Mo_Mo") ~ vacc_seq_recent,
    num_doses_recent == 3 & vacc_type_3 == "PB" ~ "other_mixed_2_dose_PB",
    num_doses_recent == 3 & vacc_type_3 == "Mo" ~ "other_mixed_2_dose_Mo",
    num_doses_recent == 3 ~ "v3_other",
    num_doses_recent == 4 & age_4cat %in% c("5-11", "12-15", "16-74") ~ "4+",
    num_doses_recent == 4 & vacc_type_1 == "AZ" & vacc_type_2 == "AZ" ~ "v4_AZ_AZ_any",
    num_doses_recent == 4 & vacc_type_1 == "PB" & vacc_type_2 == "PB" ~ "v4_PB_PB_any",
    num_doses_recent == 4 ~ "v4_other",
    num_doses_recent > 4 ~ "5+"
  )
)


q_names <- grep("Q", names(df_cohort), value = TRUE)

# Create a list of cols where NA means 0
cols = c(
  "covid_death", 
  "flag_incon_date",
  "incon_timing",
  "within_12wk_of_050922",
  "within_12wk_of_awb",
  "no_two_dose_b4_050922",
  setdiff(q_names, c("Q_BMI", "Q_ETHNICITY", "Q_HOME_CAT", "Q_LEARN_CAT", "Q_DIAG_CKD_LEVEL"))
  )

df_cohort <- mutate_at(
  df_cohort, 
  cols,
  ~ as.numeric(.))

df_cohort <- mutate_at(
  df_cohort, 
  cols, ~ 
    case_when(
      is.na(.) ~ 0,
      TRUE ~ .
    )
  )

# Add number of risk groups
# This has to be recalculated because we fused older and newer qcovid risk groups together
df_cohort <- df_cohort %>%
  mutate(Q_CKD_bin = ifelse(Q_DIAG_CKD_LEVEL == 'No CKD', 0, 1)
  ) 

df_cohort <- df_cohort %>%
  mutate(
    n_risk_gps = 
      select(df_cohort, 
        c('Q_CKD_bin', setdiff(q_names, c("Q_BMI", "Q_BMI_imputed", "Q_ETHNICITY", "Q_HOME_CAT", "Q_LEARN_CAT", "Q_DIAG_CKD_LEVEL")))
      ) %>%
      rowSums()
  ) %>%
  select(-Q_CKD_bin) %>%
  mutate(
    n_risk_gps_6cat =
        cut(n_risk_gps, 
        breaks = c(0, 1, 2, 3, 4, 5, Inf),
        labels = c('0', '1', '2', '3', '4', '5+'),
        right = FALSE
    ),
    n_risk_gps_3cat =
      cut(n_risk_gps, 
          breaks = c(0, 1, 2, Inf),
          labels = c('0', '1', '2+'),
          right = FALSE
      )
  )

# Add household characteristics
# Join by EAVE_LINKNO_old, since data pre-dates December 2021
df_cohort <- left_join(df_cohort, household, by = c("EAVE_LINKNO_old" = "EAVE_LINKNO"))

# Add shielding list
df_cohort <- left_join(df_cohort, ever_shielding, by = "EAVE_LINKNO") %>%
  mutate(shielding = replace_na(shielding, 0))

# Add whether or not another household member has ever been shielding
df_cohort <- df_cohort %>%
  group_by(household_id) %>%
  mutate(num_household_shielding = sum(shielding)) %>%
  ungroup() %>%
  mutate(other_household_shielding = ifelse(num_household_shielding - shielding >= 1, 1, 0))

# Add in time since last positive test and variant of last positive test
# Community testing stopped on 18th March, which could affect variant recording
df_cohort <- wgs %>%
  rename(
    last_positive_variant = lineage
  ) %>%
  mutate(
    last_positive_variant = case_when(
      last_positive_variant == "B.1.1.7" ~ "Alpha",
      last_positive_variant == "B.1.617.2" ~ "Delta",
      last_positive_variant %in% c("BA.1", "BA.2", "BA.3", "BA.4", "BA.5") ~ "Omicron",
      TRUE ~ "Other"
    )
  ) %>%
  arrange(EAVE_LINKNO, desc(specimen_date)) %>%
  filter(!duplicated(EAVE_LINKNO)) %>%
  right_join(df_cohort, by = "EAVE_LINKNO") %>%
  mutate(
    last_positive_variant = case_when(
      is.na(specimen_date) ~ "never_positive",
      is.na(last_positive_variant) ~ "not_sequenced",
      TRUE ~ last_positive_variant
    )
  ) %>%
  mutate(
    last_positive_variant = factor(last_positive_variant, c("never_positive", "not_sequenced", "Alpha", "Delta", "Omicron", "Other"))
  )


tests_before_study_start = filter(cdw, specimen_date < study_start)
pos_tests_before_study_start = filter(tests_before_study_start, test_result == "POSITIVE")

# Add date of last positive test
df_cohort <- pos_tests_before_study_start %>%
  rename(last_positive_test_date = specimen_date) %>%
  arrange(EAVE_LINKNO, desc(last_positive_test_date)) %>%
  filter(!duplicated(EAVE_LINKNO)) %>% 
  right_join(df_cohort) %>%
  mutate(
    last_positive_test_days = as.numeric(study_start - last_positive_test_date),
    last_positive_test_group = case_when(
      is.na(last_positive_test_days) ~ 'never_positive',
      last_positive_test_days >= 0 & last_positive_test_days < 92 ~ "0-13_weeks",
      last_positive_test_days >= 92 & last_positive_test_days < 183 ~ "14-26_weeks",
      last_positive_test_days >= 183 ~ "27+_weeks")
    ) %>%
  mutate(
    last_positive_test_group = factor(last_positive_test_group, c("never_positive", "0-13_weeks", "14-26_weeks", "27+_weeks"))
  ) %>%
  select(-last_positive_test_days)

# Add number of PCR tests in last 6 months before study_start
df_cohort <- tests_before_study_start %>%
  filter(specimen_date >= study_start %m-% months(6)) %>%
  count(EAVE_LINKNO) %>%
  rename(num_tests_6m = n) %>%
  dplyr::select(EAVE_LINKNO, num_tests_6m) %>%
  right_join(df_cohort, by = "EAVE_LINKNO") %>%
  mutate(num_tests_6m = replace_na(num_tests_6m, 0)) %>%
  mutate(num_tests_6m_group = cut(num_tests_6m,
    breaks = c(0, 1, 2, 3, 4, 5, Inf),
    labels = c("0", "1", "2", "3", "4", "5+"),
    right = FALSE
  ))

# Add number of positive PCR tests in last 6 months before study_start
df_cohort <- pos_tests_before_study_start %>%
  filter(specimen_date >= study_start %m-% months(6)) %>%
  count(EAVE_LINKNO) %>%
  rename(num_pos_tests_6m = n) %>%
  dplyr::select(EAVE_LINKNO, num_pos_tests_6m) %>%
  right_join(df_cohort, by = "EAVE_LINKNO") %>%
  mutate(num_pos_tests_6m = replace_na(num_pos_tests_6m, 0)) %>%
  mutate(num_pos_tests_6m_group = cut(num_pos_tests_6m,
    breaks = c(0, 1, 2, Inf),
    labels = c("0", "1", "2+"),
    right = FALSE
  ))

rm(tests_before_study_start, pos_tests_before_study_start)


## Add severe endpoints

study_hosps = filter(smr, hosp_date >= study_start, hosp_date <= study_end)

# acoa = any cause of admission
study_covid_acoa_hosps = study_hosps %>%
  filter(covid_main_diag_admit == 1 | covid_main_other_ep == 1) %>%
  left_join(all_hospitalisations, by = c("EAVE_LINKNO", "hosp_date"))

# mcoa = main cause of admission
study_covid_mcoa_hosps = study_covid_acoa_hosps %>%
  filter(covid_main_diag_admit == 1)

# Add covid hospitalistions with covid as any cause of admission
df_cohort = study_covid_acoa_hosps %>%
  arrange(EAVE_LINKNO, hosp_date) %>%
  filter(!duplicated(EAVE_LINKNO)) %>%
  dplyr::select(
    EAVE_LINKNO, 
    covid_acoa_hosp_date = hosp_date) %>%
  right_join(df_cohort, by = "EAVE_LINKNO") %>%
  mutate(
    covid_acoa_hosp = ifelse(!is.na(covid_acoa_hosp_date), 1, 0)
  )

# Add emergency covid hospitalistions with covid as any cause of admission
df_cohort = study_covid_acoa_hosps %>%
  filter(emergency == TRUE) %>%
  arrange(EAVE_LINKNO, hosp_date) %>%
  filter(!duplicated(EAVE_LINKNO)) %>%
  dplyr::select(
    EAVE_LINKNO, 
    covid_acoa_emerg_hosp_date = hosp_date) %>%
  right_join(df_cohort, by = "EAVE_LINKNO") %>%
  mutate(
    covid_acoa_emerg_hosp = ifelse(!is.na(covid_acoa_emerg_hosp_date), 1, 0)
  )

# Add covid hospitalistions with covid as main cause of admission
df_cohort = study_covid_mcoa_hosps %>%
  arrange(EAVE_LINKNO, hosp_date) %>%
  filter(!duplicated(EAVE_LINKNO)) %>%
  dplyr::select(
    EAVE_LINKNO, 
    covid_mcoa_hosp_date = hosp_date) %>%
  right_join(df_cohort, by = "EAVE_LINKNO") %>%
  mutate(
    covid_mcoa_hosp = ifelse(!is.na(covid_mcoa_hosp_date), 1, 0)
    )

# Add emergency covid hospitalistions with covid as main cause of admission
df_cohort = study_covid_mcoa_hosps %>%
  filter(emergency == TRUE) %>%
  arrange(EAVE_LINKNO, hosp_date) %>%
  filter(!duplicated(EAVE_LINKNO)) %>%
  dplyr::select(
    EAVE_LINKNO, 
    covid_mcoa_emerg_hosp_date = hosp_date) %>%
  right_join(df_cohort, by = "EAVE_LINKNO") %>%
  mutate(
    covid_mcoa_emerg_hosp = ifelse(!is.na(covid_mcoa_emerg_hosp_date), 1, 0)
  )

# Add covid hospitalistions with covid as main cause of admission and positive PCR test within 28 days prior to admission or 2 days after
df_cohort = study_covid_mcoa_hosps %>%
  left_join(cdw %>% filter(test_result == 'POSITIVE'), by = "EAVE_LINKNO") %>%
  mutate(
    specimen_to_hosp = hosp_date - specimen_date) %>%
  filter(specimen_to_hosp <= 28 & specimen_to_hosp >= -2) %>%
  arrange(EAVE_LINKNO, hosp_date) %>%
  filter(!duplicated(EAVE_LINKNO)) %>%
  dplyr::select(
    EAVE_LINKNO, 
    covid_mcoa_28_2_hosp_date = hosp_date) %>%
  right_join(df_cohort, by = "EAVE_LINKNO") %>%
  mutate(
    covid_mcoa_28_2_hosp = ifelse(!is.na(covid_mcoa_28_2_hosp_date), 1, 0))

# Add emergency covid hospitalistions with covid as main cause of admission and positive PCR test within 28 days prior to admission or 2 days after
df_cohort = study_covid_mcoa_hosps %>%
  filter(emergency == TRUE) %>%
  left_join(cdw %>% filter(test_result == 'POSITIVE'), by = "EAVE_LINKNO") %>%
  mutate(
    specimen_to_hosp = hosp_date - specimen_date) %>%
  filter(specimen_to_hosp <= 28, specimen_to_hosp >= -2) %>%
  arrange(EAVE_LINKNO, hosp_date) %>%
  filter(!duplicated(EAVE_LINKNO)) %>%
  dplyr::select(
    EAVE_LINKNO, 
    covid_mcoa_28_2_emerg_hosp_date = hosp_date) %>%
  right_join(df_cohort, by = "EAVE_LINKNO") %>%
  mutate(
    covid_mcoa_28_2_emerg_hosp = ifelse(!is.na(covid_mcoa_28_2_emerg_hosp_date), 1, 0))


# Add covid hospitalistions with covid as main cause of admission and positive PCR test within 14 days prior to admission or 2 days after
df_cohort = study_covid_mcoa_hosps %>%
  left_join(cdw %>% filter(test_result == 'POSITIVE'), by = "EAVE_LINKNO") %>%
  mutate(
    specimen_to_hosp = hosp_date - specimen_date) %>%
  filter(specimen_to_hosp <= 14 & specimen_to_hosp >= -2) %>%
  arrange(EAVE_LINKNO, hosp_date) %>%
  filter(!duplicated(EAVE_LINKNO)) %>%
  dplyr::select(
    EAVE_LINKNO, 
    covid_mcoa_14_2_hosp_date = hosp_date) %>%
  right_join(df_cohort, by = "EAVE_LINKNO") %>%
  mutate(
    covid_mcoa_14_2_hosp = ifelse(!is.na(covid_mcoa_14_2_hosp_date), 1, 0))


# Add emergency covid hospitalistions with covid as main cause of admission and positive PCR test within 14 days prior to admission or 2 days after
df_cohort = study_covid_mcoa_hosps %>%
  filter(emergency == TRUE) %>%
  left_join(cdw %>% filter(test_result == 'POSITIVE'), by = "EAVE_LINKNO") %>%
  mutate(
    specimen_to_hosp = hosp_date - specimen_date) %>%
  filter(specimen_to_hosp <= 14, specimen_to_hosp >= -2) %>%
  arrange(EAVE_LINKNO, hosp_date) %>%
  filter(!duplicated(EAVE_LINKNO)) %>%
  dplyr::select(
    EAVE_LINKNO, 
    covid_mcoa_14_2_emerg_hosp_date = hosp_date) %>%
  right_join(df_cohort, by = "EAVE_LINKNO") %>%
  mutate(
    covid_mcoa_14_2_emerg_hosp = ifelse(!is.na(covid_mcoa_14_2_emerg_hosp_date), 1, 0))

# Add date of first hospitalisation
df_cohort = study_hosps %>% 
  arrange(EAVE_LINKNO, hosp_date) %>%
  filter(!duplicated(EAVE_LINKNO)) %>%
  dplyr::select(
    EAVE_LINKNO, 
    hosp_date) %>%
  right_join(df_cohort, by = "EAVE_LINKNO")

# Add covid hospitalisation or death
df_cohort = df_cohort %>%
  mutate(
    covid_mcoa_hosp_death = ifelse(covid_mcoa_hosp == 1 | covid_death == 1, 1, 0),
    covid_mcoa_hosp_death_date = ifelse(covid_mcoa_hosp_death == 1, pmin(covid_mcoa_hosp_date, covid_death_date, na.rm = TRUE), as.Date(NA)) %>%
      as.Date(origin = "1970-01-01"),
    
    covid_acoa_hosp_death = ifelse(covid_acoa_hosp == 1 | covid_death == 1, 1, 0),
    covid_acoa_hosp_death_date = ifelse(covid_acoa_hosp_death == 1, pmin(covid_acoa_hosp_date, covid_death_date, na.rm = TRUE), as.Date(NA)) %>%
      as.Date(origin = "1970-01-01"),
    
    covid_acoa_emerg_hosp_death = ifelse(covid_acoa_emerg_hosp == 1 | covid_death == 1, 1, 0),
    covid_acoa_emerg_hosp_death_date = ifelse(covid_acoa_emerg_hosp_death == 1, pmin(covid_acoa_emerg_hosp_date, covid_death_date, na.rm = TRUE), as.Date(NA)) %>%
      as.Date(origin = "1970-01-01"),
    
    covid_mcoa_emerg_hosp_death = ifelse(covid_mcoa_emerg_hosp == 1 | covid_death == 1, 1, 0),
    covid_mcoa_emerg_hosp_death_date = ifelse(covid_mcoa_emerg_hosp_death == 1, pmin(covid_mcoa_emerg_hosp_date, covid_death_date, na.rm = TRUE), as.Date(NA)) %>%
      as.Date(origin = "1970-01-01")
  ) 

# Add competing non-COVID hospitalistions
df_cohort = df_cohort %>%
  mutate(
    non_covid_mcoa_hosp_date = ifelse(hosp_date < covid_mcoa_hosp_date | is.na(covid_mcoa_hosp_date), hosp_date, as.Date(NA)) %>%
      as.Date(origin = "1970-01-01"),
    non_covid_mcoa_emerg_hosp_date = ifelse(hosp_date < covid_mcoa_emerg_hosp_date | is.na(covid_mcoa_emerg_hosp_date), hosp_date, as.Date(NA)) %>%
      as.Date(origin = "1970-01-01"),
    
    non_covid_acoa_hosp_date = ifelse(hosp_date < covid_acoa_hosp_date | is.na(covid_acoa_hosp_date), hosp_date, as.Date(NA)) %>%
      as.Date(origin = "1970-01-01"),
    non_covid_acoa_emerg_hosp_date = ifelse(hosp_date < covid_acoa_emerg_hosp_date | is.na(covid_acoa_emerg_hosp_date), hosp_date, as.Date(NA)) %>%
      as.Date(origin = "1970-01-01"),
    
    non_covid_mcoa_28_2_hosp_date = ifelse(hosp_date < covid_mcoa_28_2_hosp_date | is.na(covid_mcoa_28_2_hosp_date), hosp_date, as.Date(NA)) %>%
      as.Date(origin = "1970-01-01"),
    non_covid_mcoa_28_2_emerg_hosp_date = ifelse(hosp_date < covid_mcoa_28_2_emerg_hosp_date | is.na(covid_mcoa_28_2_emerg_hosp_date), hosp_date, as.Date(NA)) %>%
      as.Date(origin = "1970-01-01"),
    
    non_covid_mcoa_14_2_hosp_date = ifelse(hosp_date < covid_mcoa_14_2_hosp_date | is.na(covid_mcoa_14_2_hosp_date), hosp_date, as.Date(NA)) %>%
      as.Date(origin = "1970-01-01"),
    non_covid_mcoa_14_2_emerg_hosp_date = ifelse(hosp_date < covid_mcoa_14_2_emerg_hosp_date | is.na(covid_mcoa_14_2_emerg_hosp_date), hosp_date, as.Date(NA)) %>%
      as.Date(origin = "1970-01-01"),
    
    non_covid_mcoa_hosp_death_date = ifelse(hosp_date < covid_mcoa_hosp_death_date | is.na(covid_mcoa_hosp_death_date), hosp_date, as.Date(NA)) %>%
      as.Date(origin = "1970-01-01"),
    non_covid_mcoa_emerg_hosp_death_date = ifelse(hosp_date < covid_mcoa_emerg_hosp_death_date | is.na(covid_mcoa_emerg_hosp_death_date), hosp_date, as.Date(NA)) %>%
      as.Date(origin = "1970-01-01"),
    
    non_covid_acoa_hosp_death_date = ifelse(hosp_date < covid_acoa_hosp_death_date | is.na(covid_acoa_hosp_death_date), hosp_date, as.Date(NA)) %>%
      as.Date(origin = "1970-01-01"),
    non_covid_acoa_emerg_hosp_death_date = ifelse(hosp_date < covid_acoa_emerg_hosp_death_date | is.na(covid_acoa_emerg_hosp_death_date), hosp_date, as.Date(NA)) %>%
      as.Date(origin = "1970-01-01")
  )
           
# Add previous covid hospitalisation
df_cohort = filter(smr, covid_main_diag_admit == 1, hosp_date < study_start, hosp_date >= as.Date("2020-03-01")) %>%
  mutate(covid_mcoa_hosp_ever = 1) %>%
  dplyr::select(EAVE_LINKNO, covid_mcoa_hosp_ever) %>%
  filter(!duplicated(EAVE_LINKNO)) %>%
  right_join(df_cohort, by = "EAVE_LINKNO") %>%
  mutate(covid_mcoa_hosp_ever = replace_na(covid_mcoa_hosp_ever, 0))
  
rm(study_hosps, study_covid_acoa_hosps, study_covid_mcoa_hosps) 


#### Re-weighting, taking into account contact with healthcare

# Load in other datasets that indicate contact with healthcare
bnf <- readRDS(paste0(Location, "EAVE/GPanalysis/data/BNF_paragraphs.rds"))
sicsag <- readRDS(paste0(Location, "EAVE/GPanalysis/data/SICSAG_episode_level_.rds"))
pis <- readRDS(paste0(Location, "EAVE/GPanalysis/data/PIS_2019_2021_EAVELink.rds"))
lft <- readRDS(paste0(Location, "EAVE/GPanalysis/data/lft_positives.rds"))

ids <- c(
  Vaccinations$EAVE_LINKNO,
  all_deaths$EAVE_LINKNO,
  all_hospitalisations$EAVE_LINKNO,
  bnf$EAVE_LINKNO,
  cdw$EAVE_LINKNO,
  sicsag$EAVE_LINKNO,
  smr$EAVE_LINKNO,
  lft$EAVE_LINKNO,
  pis$EAVE_LINKNO,
  wgs$EAVE_LINKNO
) %>%
  unique()

# This re-weighting no longer makes sense. The number of people who have
# had contact with healthcare now exceeds total weight in the cohort.
# This results in some people getting negative weight.
z_N <- round(sum(df_cohort$eave_weight))
z_k <- sum(df_cohort$EAVE_LINKNO %in% ids)
z_m <- round(sum(filter(df_cohort, (EAVE_LINKNO %in% ids))$eave_weight))

df_cohort <- df_cohort %>%
  mutate(eave_weight = if_else(EAVE_LINKNO %in% ids, 1, eave_weight * (z_N - z_k) / (z_N - z_m)))

# Remove data sets only used for identification of individuals who exist
rm(bnf, sicsag, pis, lft)

df_cohort <- data.frame(df_cohort)

df_cohort <- rename(df_cohort, individual_id = EAVE_LINKNO)

# Check NA
sapply(df_cohort, function(x) sum(is.na(x)))

# Save out
saveRDS(df_cohort, './data/df_cohort.rds')

#### df_cohort column descriptions:

# $individual_id
# [1] "character"
# it's what you think it is
# 
# $covid_mcoa_hosp_ever
# [1] "numeric"
# whether they have had a covid hospitalisation with covid as main cause of admission at any time before study_start
# 
# $death_date
# [1] "Date"
# 
# $covid_death
# [1] "numeric"
# whether or not they had a death with covid as any cause of death on death certificate
# 
# $hosp_date
# [1] "Date"
# date of their first, any-cause hospitalisation in the study period
# 
# $covid_mcoa_14_2_hosp_date
# [1] "Date"
# date of first hospitalisation with covid as main cause of admission and with a positive PCR test in the 14 days before admission or 2 days after
# 
# $covid_mcoa_28_2_hosp_date
# [1] "Date"
# date of first hospitalisation with covid as main cause of admission and with a positive PCR test in the 28 days before admission or 2 days after
# 
# $covid_mcoa_hosp_date
# [1] "Date"
# date of first hospitalisation that has covid as the main cause of admission
# 
# $covid_acoa_hosp_date
# [1] "Date"
# date of first hospitalisation with covid as any cause of admission
# 
# $num_pos_tests_6m
# [1] "factor"
# number of positive tests they have had in the last 6 months before study_start
# Levels are: 0, 1, 2+
# 
# $num_tests_6m
# [1] "integer"
# number of tests they have had in the last 6 months before study_start
# 
# $specimen_date
# [1] "Date"
# specimen date of last positive test
# 
# $last_positive_variant
# [1] "factor"
# 
# $sex
# [1] "factor"
# Levels are Female and Male
# 
# $age
# [1] "numeric"
# 
# $urban_rural_6cat
# [1] "factor"
# Levels are:
# Large Urban Areas
# Other Urban Areas
# Accessible Small Towns
# Remote Small Towns
# Accessible Rural
# Remote Rural
# 
# $simd2020_sc_quintile
# [1] "factor"
# quintile of the Scottish index of multiple deprivation
# Confusingly, 1 is the most deprived, 5 is the most affluent
# 
# $urban_rural_2cat
# [1] "character"
# Levels are:
# Urban, Rural
# 
# $age_3cat
# [1] "factor"
# Levels are:
# 5-15, 16-74, 75+
# 
# $age_4cat
# [1] "factor"
# Levels are:
# 5-11, 12-15, 16-74, 75+
# 
# $age_17cat
# [1] "factor"
# Levels are:
# 5-11, 12-15, 16-17, 18-24, 25-29, 30-34, 35-39, 40-44, 
# 45-49, 50-54, 55-59, 60-64, 65-69, 70-74, 75-79, 80-84, 85+
# 
# $eave_weight
# [1] "numeric"
# their weighting
# 
# $ethnicity_18cat
# [1] "factor"
# 
# $ethnicity_5cat
# [1] "factor"
# Levels are:
# White, Black, Asian, Mixed, Other


###### variables starting with Q are qcovid categories

# $Q_BMI
# [1] "numeric"
# 
# $Q_ETHNICITY
# [1] "character"
# 
# $Q_HOME_CAT
# [1] "factor"
# 
# $Q_LEARN_CAT
# [1] "factor"
# 
# $Q_DIAG_CKD_LEVEL
# [1] "factor"
# 
# $Q_DIAG_AF
# [1] "numeric"
# 
# $Q_DIAG_ASTHMA
# [1] "numeric"
# 
# $Q_DIAG_BLOOD_CANCER
# [1] "numeric"
# 
# $Q_DIAG_CCF
# [1] "numeric"
# 
# $Q_DIAG_CEREBRALPALSY
# [1] "numeric"
# 
# $Q_DIAG_CHD
# [1] "numeric"
# 
# $Q_DIAG_CIRRHOSIS
# [1] "numeric"
# 
# $Q_DIAG_CONGEN_HD
# [1] "numeric"
# 
# $Q_DIAG_COPD
# [1] "numeric"
# 
# $Q_DIAG_DEMENTIA
# [1] "numeric"
# 
# $Q_DIAG_EPILEPSY
# [1] "numeric"
# 
# $Q_DIAG_FRACTURE
# [1] "numeric"
# 
# $Q_DIAG_HIV_AIDS
# [1] "numeric"
# 
# $Q_DIAG_IMMU
# [1] "numeric"
# 
# $Q_DIAG_NEURO
# [1] "numeric"
# 
# $Q_DIAG_PARKINSONS
# [1] "numeric"
# 
# $Q_DIAG_PULM_HYPER
# [1] "numeric"
# 
# $Q_DIAG_PULM_RARE
# [1] "numeric"
# 
# $Q_DIAG_PVD
# [1] "numeric"
# 
# $Q_DIAG_RA_SLE
# [1] "numeric"
# 
# $Q_DIAG_RESP_CANCER
# [1] "numeric"
# 
# $Q_DIAG_SEV_MENT_ILL
# [1] "numeric"
# 
# $Q_DIAG_SICKLE_CELL
# [1] "numeric"
# 
# $Q_DIAG_STROKE
# [1] "numeric"
# 
# $Q_DIAG_VTE
# [1] "numeric"
#
# $n_risk_gps
# [1] "numeric"
# count of qcovid clinical risk groups. This does not include BMI, ethnicity,
# learning disability or housing category.
#
# $n_risk_gps_6cat
# [1] "factor"
# count of qcovid clinical risk groups. This does not include BMI, ethnicity,
# learning disability or housing category. Levels aer:
# 0, 1, 2, 3, 4, 5+
#
# $n_risk_gps_3cat
# [1] "factor"
# count of qcovid clinical risk groups. This does not include BMI, ethnicity,
# learning disability or housing category. Levels aer:
# 0, 1, 2+
# 
# $Q_DIAG_DIABETES_1
# [1] "numeric"
# 
# $Q_DIAG_DIABETES_2
# [1] "numeric"
# 
####

# $bmi_cat
# [1] "factor"
# body mass index category, and Levels are:
# 10-18.5, 18.5-25, 25-30, 30-35, 35-40, 40+
# 
# $smoking_status
# [1] "factor"
# Levels are:
# Ex-smoker, Non-smoker, Smoker
# 
# $blood_pressure
# [1] "factor"
# Levels are:
# Normal, Low, High, Very high
# 

###### date_vacc_x is the date of dose number x
###### vacc_type_x is the type of dose number x
###### wb stands for winter booster

# $vacc_type_1
# [1] "character"
# 
# $vacc_type_2
# [1] "character"
# 
# $date_vacc_1
# [1] "Date"
# 
# $date_vacc_2
# [1] "Date"
# 
# $vacc_type_3
# [1] "character"
# 
# $date_vacc_3
# [1] "Date"
# 
# $vacc_type_4
# [1] "character"
# 
# $date_vacc_4
# [1] "Date"
# 
# $vacc_type_5
# [1] "character"
# 
# $date_vacc_5
# [1] "Date"
# 
# $vacc_type_wb
# [1] "character"
# 
# $date_vacc_wb
# [1] "Date"

######

# $num_doses_recent
# [1] "numeric"
# number of vaccine doses they have had according to most recent data
# 
# $num_doses_start
# [1] "numeric"
# number of doses they had at study_start
# 
# $vacc_seq_recent
# [1] "character"
# sequence of vaccine types they have received according to most recent data
# 
# $vacc_seq_start
# [1] "character"
# sequence of vaccine types they had received at study_start
# 
# $mixed_vacc_recent
# [1] "numeric"
# whether they have had any mixed vaccine type according to most recent data
# 
# $mixed_vacc_start
# [1] "numeric"
# whether they have had any mixed vaccine type at study start
# 
# $vs_recent
# [1] "character"
# dose number and type of the last dose they had according to most recent data
# 
# $vs_start
# [1] "character"
# dose number and type of the last dose they had at study_start
# 
# $under_vaccinated
# [1] "factor"
# whether they have had the recommended
# vaccine schedule for their age group at study_start
# 
# $vs_mixed_start
# [1] "character"
# a more coarse_grained version of vacc_seq_start, with categories of vaccine sequences
# received by study_start. Details are in the draft tables document.
# 
# $vs_mixed_recent
# [1] "character"
# a more coarse_grained version of vacc_seq_start, with categories of vaccine sequences
# received according to most recent data. Details are in the draft tables document.
# 
# $mean_household_age
# [1] "numeric"
# 
# $num_ppl_household
# [1] "factor"
# number of people in household
# Levels are:
# 1, 2, 3-5, 6-10, 11-30, 31-100, 101+
# 
# $shielding
# [1] "numeric"
# whether they have ever been on the shielding list
# 
# $last_positive_test
# [1] "numeric"
# date of last positive test before study_start
# 
# $last_positive_test_group
# [1] "factor"
# Levels are:
# never_positive, 0-13_weeks, 14-26_weeks, 27+_weeks
# 
# last_positive_test_variant
# Levels are:
# never_positive, not_sequenced, Alpha, Delta, Omicron, Other
#
# $num_tests_6m_group
# [1] "factor"
# Levels are:
# 0, 1, 2, 3, 4-9, 10+
# 
# $covid_acoa_hosp
# [1] "numeric"
# whether they had a hospitalisation with covid as any cause of admission
# 
# $covid_mcoa_hosp
# [1] "numeric"
# whether they had a hospitalisation with covid as main cause of admission
# 
# $covid_mcoa_28_2_hosp
# [1] "numeric"
# whether they had a hospitalisation with covid as main cause of admission,
# and a positive PCR test in the 28 days prior to admission or 2 days after admission
# 
# $covid_mcoa_14_2_hosp
# [1] "numeric"
# whether they had a hospitalisation with covid as main cause of admission,
# and a positive PCR test in the 28 days prior to admission or 2 days after admission
# 
# $covid_death_date
# [1] "Date"
# 
# $covid_mcoa_hosp_death
# [1] "numeric"
# whether they had a hospitalisation with covid as main cause of admission
# or a covid death
# 
# $covid_mcoa_hosp_death_date
# [1] "Date"
# date of first hospitalisation with covid as main cause of admission
# or covid death
# 
# $non_covid_mcoa_hosp_date
# [1] "Date"
# date of first hospitalisation that is not a covid_mcoa_hosp
# 
# $non_covid_acoa_hosp_date
# [1] "Date"
# date of first hospitalisation that is not a covid_acoa_hosp
# 
# $non_covid_mcoa_28_2_hosp_date
# [1] "Date"
# date of first hospitalisation that is not a covid_28_2_mcoa_hosp
# 
# $non_covid_mcoa_14_2_hosp_date
# [1] "Date"
# date of first hospitalisation that is not a covid_14_2_mcoa_hosp
# 
# $non_covid_mcoa_hosp_death_date
# [1] "Date"
# date of first hospitalisation that is not a covid_mcoa_hosp_death

