######################################################################
## Title: [Insert full title of paper]
## Code author: Steven Kerr steven.kerr@ed.ac.uk
## Description: Load in all data and create a cohort dataframe
######################################################################

# Libraries
library(tidyverse)
library(lubridate)

Location = "/conf/"

setwd("/conf/EAVE/GPanalysis/analyses/under_vaccinated")

# study_start is the date where there under-vaccination status is determined
study_start = as.Date("2022-06-01")
# study_end is the start of follow-up for serious events
study_end = as.Date("2022-09-30")

# R does't have a built-in mode function
Mode = function(x) {
  x = x[!is.na(x)]

  if (length(x) == 0) {
    return(NA)
  } else {
    ux = unique(x)
    return(ux[which.max(tabulate(match(x, ux)))])
  }
}

#### 1 Load data

EAVE_cohort = readRDS(paste0(Location, "EAVE/GPanalysis/outputs/temp/Cohort_Demog_Endpoints_Dates2021-07-28.rds")) %>%
  filter(!duplicated(EAVE_LINKNO)) %>%
  # Original cohort is from approximately 2 years ago. Update their age.
  mutate(ageYear = ageYear + 2) %>%
  # remove all who have died before the beginning
  select(EAVE_LINKNO:ur6_2016_name)

EAVE_Weights = readRDS(paste0(Location, "EAVE/GPanalysis/outputs/temp/CR_Cohort_Weights.rds"))

# Get qcovid risk groups
qcovid_rg = readRDS(paste0(Location, "EAVE/GPanalysis/data/cleaned_data/QCOVID_feb22.rds")) %>%
  select(-Age, -Sex) %>%
  mutate(Q_BMI = as.numeric(Q_BMI)) %>%
  # Drop diabetes 1 and 2 cateogries because they are not useable
  # Ethnicity is not reliable
  select(-Q_DIAG_DIABETES_1, -Q_DIAG_DIABETES_2, -Q_ETHNICITY)

# Get old QCovid risk groups. We will use this for diabetes only, because diabetes grouping
# from more recent QCovid extract are not useable
qcovid_diabetes = readRDS("/conf/EAVE/GPanalysis/progs/CR/Vaccine/output/temp/Qcovid.rds") %>%
  select(EAVE_LINKNO, Q_DIAG_DIABETES_1, Q_DIAG_DIABETES_2)

qcovid_rg = full_join(qcovid_rg, qcovid_diabetes)

# Get EAVE_BP and EAVE_Smoke risk groups
eave_rg = readRDS(paste0(Location, "EAVE/GPanalysis/outputs/temp/CR_Cohort_RG_EAVE_BP_Smoke.rds")) %>%
  filter(!duplicated(EAVE_LINKNO)) %>%
  select(EAVE_LINKNO, EAVE_Smoking_Status_Worst, EAVE_BP) %>%
  rename(EAVE_Smoke = EAVE_Smoking_Status_Worst) %>%
  mutate(
    EAVE_Smoke = as.character(EAVE_Smoke),
    EAVE_BP = as.character(EAVE_BP)
  ) %>%
  mutate(
    EAVE_Smoke = case_when(
      EAVE_Smoke == "Unknown" ~ NA_character_,
      TRUE ~ EAVE_Smoke
    ),
    EAVE_BP = case_when(
      EAVE_BP == "No Investigation" ~ NA_character_,
      TRUE ~ EAVE_BP
    )
  ) %>%
  mutate(
    EAVE_Smoke = as.factor(EAVE_Smoke),
    EAVE_BP = as.factor(EAVE_BP)
  )

# Combine risk groups
rg = qcovid_rg %>%
  full_join(eave_rg, by = "EAVE_LINKNO") %>%
  filter(!duplicated(EAVE_LINKNO))

# Tests
# PCR
cdw = readRDS(paste0(Location, "EAVE/GPanalysis/data/CDW_full.rds")) %>%
  mutate(
    flag_covid_symptomatic = if_else(!is.na(flag_covid_symptomatic) & flag_covid_symptomatic == "true", 1L, 0L),
    date_ecoss_specimen = as.Date(date_ecoss_specimen)
  ) %>%
  filter(date_ecoss_specimen <= Sys.Date())

lft = readRDS(paste0(Location, "EAVE/GPanalysis/data/lft_positives.rds")) %>%
  arrange(EAVE_LINKNO, date_ecoss_specimen) %>%
  # Get one test per person per day
  filter(!duplicated(paste(EAVE_LINKNO, date_ecoss_specimen)))

# All deaths
all_deaths = readRDS(paste0(Location, "EAVE/GPanalysis/data/all_deaths.rds")) %>%
  rowwise() %>%
  mutate(covid_death = ifelse(rowSums(across(UNDERLYING_CAUSE_OF_DEATH:CAUSE_OF_DEATH_CODE_9, ~ .x %in% c("U071", "U072")), na.rm = T) > 0, 1, 0))

# All hospitalisations
all_hospitalisations = readRDS(paste0(Location, "EAVE/GPanalysis/data/automated_any_hospitalisation_post_01022020.rds"))

# SMR
smr = readRDS(paste0(Location, "EAVE/GPanalysis/data/SMR01_allstays.rds")) %>%
  mutate(
    ADMISSION_DATE = as.Date(ADMISSION_DATE),
    DISCHARGE_DATE = as.Date(DISCHARGE_DATE)
  )

# Covid hospitalisation is defined as hospitalisation with positive PCR test in 28 days prior, or 2 days after admission
# with covid as main cause of admission
covid_hospitalisations = smr %>%
  filter(covid_main_diag_admit == 1 | covid_main_other_ep == 1) %>%
  left_join(cdw) %>%
  mutate(specimen_to_hosp = ADMISSION_DATE - date_ecoss_specimen) %>%
  filter(specimen_to_hosp <= 28 & specimen_to_hosp >= -2) %>%
  select(EAVE_LINKNO, ADMISSION_DATE)

# Import vaccination data
# source('/conf/EAVE/GPanalysis/progs/Data_Cleaning/00_Read_DV_Vaccinations_Dose4.R')
Vaccinations = readRDS(paste0(Location, "/EAVE/GPanalysis/data/cleaned_data/C19vaccine_dvprod_cleaned.rds")) %>%
  select(
    EAVE_LINKNO, first_dose, second_dose, third_dose, fourth_dose, fifth_dose,
    vacc_dose_number, vacc_product_name
  ) %>%
  mutate(
    vacc_product_name = gsub("Covid-19 Vaccine AstraZeneca", "AZ", vacc_product_name),
    vacc_product_name = gsub("Covid-19 mRNA Vaccine Pfizer", "PB", vacc_product_name),
    vacc_product_name = gsub("Covid-19 mRNA Vaccine Moderna", "Mo", vacc_product_name),
    vacc_product_name = gsub("Covid-19 mRNA Vaccine Spikevax Bivalent Moderna", "Mo", vacc_product_name),
    vacc_product_name = gsub("Covid-19 mRNA Vaccine Comirnaty Bivalent Pfizer", "PB", vacc_product_name),
    vacc_product_name = gsub("Covid-19 Vaccine Novavax ", "No", vacc_product_name)
  ) %>%
  # 4th and 5th dose vaccinations with the moderna bivalent vaccine are not currently picked up
  mutate(vacc_dose_number = case_when(
    !is.na(fifth_dose) & is.na(vacc_dose_number) ~ 5,
    !is.na(fourth_dose) & is.na(vacc_dose_number) ~ 4,
    TRUE ~ vacc_dose_number
  )) %>%
  pivot_wider(
    id_cols = c(
      "EAVE_LINKNO", "first_dose", "second_dose", "third_dose", "fourth_dose",
      "fifth_dose"
    ),
    names_prefix = "vacc_type_", names_from = vacc_dose_number,
    values_from = vacc_product_name, values_fill = NA, values_fn = Mode
  ) %>%
  rename(
    date_vacc_1 = first_dose,
    date_vacc_2 = second_dose,
    date_vacc_3 = third_dose,
    date_vacc_4 = fourth_dose,
    date_vacc_5 = fifth_dose
  ) %>%
  data.frame()

# Clean rows with duplicated EAVE_LINKNO
duplicate_vacc_ids = Vaccinations$EAVE_LINKNO[duplicated(Vaccinations$EAVE_LINKNO)]

duplicate_vaccs = filter(Vaccinations, EAVE_LINKNO %in% duplicate_vacc_ids) %>%
  group_by(EAVE_LINKNO) %>%
  mutate_at(setdiff(names(Vaccinations), "EAVE_LINKNO"), Mode) %>%
  distinct()

Vaccinations = bind_rows(filter(Vaccinations, !(EAVE_LINKNO %in% duplicate_vacc_ids)), duplicate_vaccs)

# Shielding is no longer available in vaccination data. Load up an historical version
# to see if they have ever been shielding
# shielding = readRDS(paste0(Location,"EAVE/GPanalysis/data/Shielding_list.rds")) %>%
#   mutate(shielding = 1)
ever_shielding = readRDS(paste0(Location, "/EAVE/GPanalysis/data/cleaned_data/C19vaccine_dvprod_cleaned_incl_cohorts_20220926.rds")) %>%
  select(EAVE_LINKNO, shielding) %>%
  arrange(EAVE_LINKNO, desc(shielding)) %>%
  filter(!duplicated(EAVE_LINKNO))

# Household information (from Sept 2020)
Cohort_Household = readRDS(paste0(Location, "EAVE/GPanalysis/outputs/temp/Cohort_Household.rds")) %>%
  mutate(n_hh_gp = cut(n_hh,
    breaks = c(0, 1, 2, 5, 10, 30, 100, max(n_hh)),
    labels = c("1", "2", "3-5", "6-10", "11-30", "31-100", "101+")
  )) %>%
  mutate(ave_hh_age = if_else(is.na(ave_hh_age), mean(ave_hh_age, na.rm = T), ave_hh_age))

# endpoints
# This only includes covid hospitalisations/deaths when they are the listed as the cause
# Should I also take deaths/hospitalisations within 28 days of positive test?
# endpoints = readRDS(paste0(Location,'/EAVE/GPanalysis/outputs/temp/severe_endpoints2022-06-23.rds')) %>%
#   mutate(covid_hosp = case_when( covid_mcoa_hosp == 1 & emergency == 1 ~ 1,
#                                  TRUE ~ 0),
#          covid_death = case_when( covid_cod == 1 ~ 1,
#                                   TRUE ~ 0),
#          covid_hosp_death = case_when(covid_hosp == 1 | covid_death == 1 ~ 1,
#                                       TRUE ~ 0) ) %>%
#   mutate(covid_hosp_date = case_when(covid_hosp == 1 ~ hosp_admit_date,
#                                      TRUE ~ as.Date(NA)),
#          covid_death_date = case_when(covid_death == 1 ~ NRS.Date.Death,
#                                      TRUE ~ as.Date(NA)),
#          covid_hosp_death_date = case_when(covid_hosp_death == 1 ~ pmin(covid_hosp_date, covid_death_date, na.rm = TRUE),
#                                       TRUE ~ as.Date(NA)) ) %>%
#     select(EAVE_LINKNO, hosp_admit_date, NRS.Date.Death, covid_hosp, covid_hosp_date,
#            covid_death, covid_death_date,
#            covid_hosp_death, covid_hosp_death_date)

endpoints = covid_hospitalisations %>%
  rename(covid_hosp_date = ADMISSION_DATE) %>%
  mutate(covid_hosp = 1) %>%
  full_join(all_deaths %>%
    filter(covid_death == 1) %>%
    rename(covid_death_date = NRS.Date.Death) %>%
    select(EAVE_LINKNO, covid_death, covid_death_date)) %>%
  mutate(covid_hosp_death = case_when(
    covid_hosp == 1 | covid_death == 1 ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(covid_hosp_death_date = case_when(
    covid_hosp_death == 1 ~ pmin(covid_hosp_date, covid_death_date, na.rm = TRUE),
    TRUE ~ as.Date(NA)
  )) %>%
  distinct()



#### 2 Create cohort dataframe

df_cohort = EAVE_cohort %>%
  # It has been a while since the EAVE cohort extract, add 2 to age.
  mutate(ageYear = ageYear + 2) %>%
  filter(ageYear >= 5) %>%
  mutate(
    age_gp = cut(ageYear,
      breaks = c(5, 16, 75, Inf),
      labels = c("5-15", "16-74", "75+"),
      right = FALSE
    ),
    age_gp_2 = cut(ageYear,
      breaks = c(5, 12, 16, 18, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, Inf),
      labels = c(
        "5-11", "12-15", "16-17", "18-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54",
        "55-59", "60-64", "65-69", "70-74", "75+"
      ),
      right = FALSE
    )
  ) %>%
  left_join(EAVE_Weights, by = "EAVE_LINKNO")

df_cohort$eave_weight[is.na(df_cohort$eave_weight)] = mean(df_cohort$eave_weight, na.rm = T)

# Add risk groups
df_cohort = left_join(df_cohort, rg, by = "EAVE_LINKNO")

# Im taking 10 to be the lowest feasible human BMI, and 100 the largest
# There are lots of odd values for BMI, which probably happens because issues
# with units that calculation is done in.
df_cohort$Q_BMI = ifelse(df_cohort$Q_BMI < 10 | df_cohort$Q_BMI > 100, NA_real_, df_cohort$Q_BMI)

df_cohort = mutate(df_cohort, 
  bmi_cat = cut(Q_BMI,
    breaks = c(10, 18.5, 25, 30, 35, 40, Inf),
    labels = c("10-18.5", "18.5-25", "25-30", "30-35", "35-40", "40+"),
    include.lowest = TRUE,
    right = FALSE
))


# Add vaccinations
df_cohort = left_join(df_cohort, Vaccinations, by = "EAVE_LINKNO") %>%
  # Number of doses they have had according to most recent data
  mutate(
    num_doses_recent = case_when(
      !is.na(date_vacc_5) ~ 5,
      !is.na(date_vacc_4) ~ 4,
      !is.na(date_vacc_3) ~ 3,
      !is.na(date_vacc_2) ~ 2,
      !is.na(date_vacc_1) ~ 1,
      TRUE ~ 0
    ),
    # Number of doses at study_start
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
    mixed_vacc_start = case_when(
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
    num_doses_recent = factor(num_doses_recent),
    num_doses_start = factor(num_doses_start),
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

# Count all combinations of vaccines that anyone has received
vacc_seq_start_count = count(df_cohort, vacc_seq_start)
vacc_seq_start_count_pruned = filter(vacc_seq_start_count, n >= 1000)

# More vaccine related derived variables
df_cohort = mutate(df_cohort,
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
    fully_vaccinated = case_when(
      age_gp == "5-15" & date_vacc_2 <= study_start ~ 1,
      age_gp == "16-74" & date_vacc_3 <= study_start ~ 1,
      age_gp == "75+" & date_vacc_4 <= study_start ~ 1,
      TRUE ~ 0
    ),
    # Mixed vaccine status
    vs_mixed_start = case_when(
      mixed_vacc_start == 1 & !(vacc_seq_start %in% vacc_seq_start_count_pruned$vacc_seq) ~ paste0("v", num_doses_start, "_mixed"),
      TRUE ~ vacc_seq_start
    )
  ) %>%
  mutate(fully_vaccinated = factor(fully_vaccinated))

# Add household characteristics
df_cohort = Cohort_Household %>%
  select(EAVE_LINKNO, n_hh_gp, ave_hh_age) %>%
  right_join(df_cohort)

# Add deaths from any cause
df_cohort = all_deaths %>%
  select(EAVE_LINKNO, NRS.Date.Death, covid_death) %>%
  right_join(df_cohort) %>%
  mutate(covid_death = replace_na(covid_death, 0))

# Add in covid hospitalisations
df_cohort = endpoints %>%
  group_by(EAVE_LINKNO) %>%
  summarise(covid_hosps = sum(covid_hosp)) %>%
  mutate(covid_hosps = replace_na(covid_hosps, 0)) %>%
  ungroup() %>%
  right_join(df_cohort) %>%
  mutate(covid_hosp_ever = case_when(
    covid_hosps >= 1 ~ 1,
    TRUE ~ 0
  ))

# Add most recent positive test
df_cohort = filter(cdw, test_result == "POSITIVE") %>%
  select(EAVE_LINKNO, date_ecoss_specimen) %>%
  arrange(EAVE_LINKNO, desc(date_ecoss_specimen)) %>%
  filter(!duplicated(EAVE_LINKNO)) %>%
  mutate(last_positive_test = as.numeric(study_start - date_ecoss_specimen)) %>%
  select(EAVE_LINKNO, last_positive_test) %>%
  right_join(df_cohort) %>%
  mutate(last_positive_test_group = cut(last_positive_test,
    breaks = c(0, 92, 183, Inf),
    labels = c("0-13_weeks", "14-26_weeks", "27+_weeks"),
    right = FALSE
  )) %>%
  # Make NA a level for last positive test group
  mutate(last_positive_test_group = addNA(last_positive_test_group)) %>%
  mutate(last_positive_test_group = fct_relevel(last_positive_test_group, NA))

# Add number of PCR tests in last 6 months
df_cohort = filter(cdw, date_ecoss_specimen >= study_start %m-% months(6)) %>%
  count(EAVE_LINKNO) %>%
  rename(num_tests = n) %>%
  select(EAVE_LINKNO, num_tests) %>%
  right_join(df_cohort) %>%
  mutate(num_tests = replace_na(num_tests, 0)) %>%
  mutate(num_tests_6m_group = cut(num_tests,
    breaks = c(0, 1, 2, 3, 4, 10, Inf),
    labels = c("0", "1", "2", "3", "4-9", "10+"),
    right = FALSE
  ))

# Add shielding list
df_cohort = ever_shielding %>%
  right_join(df_cohort) %>%
  mutate(shielding = replace_na(shielding, 0))

# Re-weight, taking into account contact with healthcare
bnf = readRDS(paste0(Location, "EAVE/GPanalysis/data/BNF_paragraphs.rds"))
sicsag = readRDS(paste0(Location, "EAVE/GPanalysis/data/SICSAG_episode_level_.rds"))
pis = readRDS(paste0(Location, "EAVE/GPanalysis/data/PIS_2019_2021_EAVELink.rds"))

z_ids = c(
  Vaccinations$EAVE_LINKNO,
  all_deaths$EAVE_LINKNO,
  all_hospitalisations$EAVE_LINKNO,
  bnf$EAVE_LINKNO,
  cdw$EAVE_LINKNO,
  sicsag$EAVE_LINKNO,
  smr$EAVE_LINKNO,
  lft$EAVE_LINKNO,
  pis$EAVE_LINKNO
) %>%
  unique()

z_N = round(sum(df_cohort$eave_weight))
z_k = sum(df_cohort$EAVE_LINKNO %in% z_ids)
z_m = round(sum(filter(df_cohort, (EAVE_LINKNO %in% z_ids))$eave_weight))

df_cohort = df_cohort %>%
  mutate(eave_weight = if_else(EAVE_LINKNO %in% z_ids, 1, eave_weight * (z_N - z_k) / (z_N - z_m)))

# Remove data sets only used for identification of patients who exist
rm(bnf, sicsag, pis)

# Create a list of cols where NA means 0
q_names = grep("Q", names(df_cohort), value = TRUE)

df_cohort = mutate_at(df_cohort, q_names, ~ as.numeric(.))

cols = setdiff(q_names, c("Q_BMI", "Q_ETHNICITY"))
cols = c("covid_hosps", "covid_death", cols)

df_cohort = mutate_at(df_cohort, cols, ~ case_when(
  is.na(.) ~ 0,
  TRUE ~ .
)) %>%
  mutate(
    Q_HOME_CAT = case_when(
      Q_HOME_CAT == 0 ~ "Neither",
      Q_HOME_CAT == 1 ~ "Care home",
      Q_HOME_CAT == 2 ~ "Homeless"
    ),
    Q_LEARN_CAT = case_when(
      Q_LEARN_CAT == 0 ~ "Neither",
      Q_LEARN_CAT == 1 ~ "Learning disability",
      Q_LEARN_CAT == 2 ~ "Down's syndrome"
    ),
    Q_DIAG_CKD_LEVEL = case_when(
      Q_DIAG_CKD_LEVEL == 0 ~ "No CKD",
      Q_DIAG_CKD_LEVEL == 3 ~ "CKD 3",
      Q_DIAG_CKD_LEVEL == 4 ~ "CKD 4",
      Q_DIAG_CKD_LEVEL == 5 ~ "CKD 5"
    )
  ) %>%
  mutate_at(setdiff(q_names, "Q_BMI"), ~ as.factor(.))

df_cohort = mutate(df_cohort,
  n_risk_gps = as.character(n_risk_gps)
) %>%
  mutate(
    n_risk_gps = case_when(
      is.na(n_risk_gps) ~ "0",
      TRUE ~ n_risk_gps
    )
  ) %>%
  mutate(n_risk_gps = as.factor(n_risk_gps))


df_cohort = mutate_at(df_cohort, c("EAVE_Smoke", "EAVE_BP"), ~ as.character(.)) %>%
  mutate(
    EAVE_Smoke = case_when(
      EAVE_Smoke == "Ex Smoker" ~ "Ex-smoker",
      EAVE_Smoke == "Non Smoker" ~ "Non-smoker",
      TRUE ~ EAVE_Smoke
    ),
    EAVE_BP = case_when(
      EAVE_BP == "Very High" ~ "Very high",
      TRUE ~ EAVE_BP
    ),
    ur6_2016_name = case_when(
      ur6_2016_name == "1 Large Urban Areas" ~ "Large Urban Areas",
      ur6_2016_name == "2 Other Urban Areas" ~ "Other Urban Areas",
      ur6_2016_name == "3 Accessible Small Towns" ~ "Accessible Small Towns",
      ur6_2016_name == "4 Remote Small Towns" ~ "Remote Small Towns",
      ur6_2016_name == "5 Accessible Rural" ~ "Accessible Rural",
      ur6_2016_name == "6 Remote Rural" ~ "Remote Rural"
    )
  ) %>%
  mutate_at(c("EAVE_Smoke", "EAVE_BP"), ~ as.factor(.)) %>%
  mutate(
    EAVE_Smoke = fct_relevel(EAVE_Smoke, "Non-smoker"),
    EAVE_BP = fct_relevel(EAVE_BP, "Normal", "Low", "High", "Very high"),
    Q_HOME_CAT = fct_relevel(Q_HOME_CAT, "Neither"),
    Q_LEARN_CAT = fct_relevel(Q_LEARN_CAT, "Neither"),
    Q_DIAG_CKD_LEVEL = fct_relevel(Q_DIAG_CKD_LEVEL, "No CKD"),
    ur6_2016_name = fct_relevel(
      ur6_2016_name, "Large Urban Areas",
      "Other Urban Areas",
      "Accessible Small Towns",
      "Remote Small Towns",
      "Accessible Rural",
      "Remote Rural"
    ),
    simd2020_sc_quintile = as.factor(as.character(simd2020_sc_quintile))
  ) %>%
  data.frame()

# Check NA
sapply(df_cohort, function(x) sum(is.na(x)))

df_cohort = as.data.frame(df_cohort)

# df_cohort has the folloiwng columns:
#
# [1] "EAVE_LINKNO"              "shielding"                "num_tests"                "last_positive_test"      
# [5] "covid_hosps"              "NRS.Date.Death"           "covid_death"              "n_hh_gp"                 
# [9] "ave_hh_age"               "Sex"                      "ageYear"                  "simd2020_sc_quintile"    
# [13] "DataZone"                 "ur6_2016_name"            "age_gp"                   "age_gp_2"                
# [17] "eave_weight"              "Q_BMI"                    "Q_HOME_CAT"               "Q_LEARN_CAT"             
# [21] "Q_DIAG_CKD_LEVEL"         "Q_DIAG_AF"                "Q_DIAG_ASTHMA"            "Q_DIAG_BLOOD_CANCER"     
# [25] "Q_DIAG_CCF"               "Q_DIAG_CEREBRALPALSY"     "Q_DIAG_CHD"               "Q_DIAG_CIRRHOSIS"        
# [29] "Q_DIAG_CONGEN_HD"         "Q_DIAG_COPD"              "Q_DIAG_DEMENTIA"          "Q_DIAG_EPILEPSY"         
# [33] "Q_DIAG_FRACTURE"          "Q_DIAG_HIV_AIDS"          "Q_DIAG_IMMU"              "Q_DIAG_NEURO"            
# [37] "Q_DIAG_PARKINSONS"        "Q_DIAG_PULM_HYPER"        "Q_DIAG_PULM_RARE"         "Q_DIAG_PVD"              
# [41] "Q_DIAG_RA_SLE"            "Q_DIAG_RESP_CANCER"       "Q_DIAG_SEV_MENT_ILL"      "Q_DIAG_SICKLE_CELL"      
# [45] "Q_DIAG_STROKE"            "Q_DIAG_VTE"               "n_risk_gps"               "Q_DIAG_DIABETES_1"       
# [49] "Q_DIAG_DIABETES_2"        "EAVE_Smoke"               "EAVE_BP"                  "bmi_cat"                 
# [53] "date_vacc_1"              "date_vacc_2"              "date_vacc_3"              "date_vacc_4"             
# [57] "date_vacc_5"              "vacc_type_1"              "vacc_type_2"              "vacc_type_3"             
# [61] "vacc_type_4"              "vacc_type_5"              "vacc_type_NA"             "vacc_type_6"             
# [65] "vacc_type_7"              "num_doses_recent"         "num_doses_start"          "vacc_seq_recent"         
# [69] "vacc_seq_start"           "mixed_vacc_start"         "vs_recent"                "vs_start"                
# [73] "fully_vaccinated"         "vs_mixed_start"           "covid_hosp_ever"          "last_positive_test_group"
# [77] "num_tests_6m_group"
