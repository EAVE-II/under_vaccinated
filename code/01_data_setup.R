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

# study_start is the date where their under-vaccination status is determined
study_start = as.Date("2022-06-01")
# study_end is the end of follow-up for serious events
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
  select(EAVE_LINKNO:ur6_2016_name) %>%
  rename(sex = Sex,
         age = ageYear,
         urban_rural_class = ur6_2016_name)

EAVE_Weights = readRDS(paste0(Location, "EAVE/GPanalysis/outputs/temp/CR_Cohort_Weights.rds"))

# Get qcovid risk groups
qcovid_rg = readRDS(paste0(Location, "EAVE/GPanalysis/data/cleaned_data/QCOVID_feb22.rds")) %>%
  select(-Age, -Sex) %>%
  mutate(
    Q_BMI = as.numeric(Q_BMI)
  ) %>%
  # Drop diabetes 1 and 2 cateogries because they are not useable
  select(-Q_DIAG_DIABETES_1, -Q_DIAG_DIABETES_2)

# Get old QCovid risk groups. We will use this for diabetes only, because diabetes grouping
# from more recent QCovid extract are not useable
qcovid_diabetes = readRDS("/conf/EAVE/GPanalysis/progs/CR/Vaccine/output/temp/Qcovid.rds") %>%
  select(EAVE_LINKNO, Q_DIAG_DIABETES_1, Q_DIAG_DIABETES_2)

qcovid_rg = full_join(qcovid_rg, qcovid_diabetes)

# Get blood pressure and smoking risk groups
eave_rg = readRDS(paste0(Location, "EAVE/GPanalysis/outputs/temp/CR_Cohort_RG_EAVE_BP_Smoke.rds")) %>%
  filter(!duplicated(EAVE_LINKNO)) %>%
  select(EAVE_LINKNO, EAVE_Smoking_Status_Worst, EAVE_BP) %>%
  rename(smoking_status = EAVE_Smoking_Status_Worst,
         blood_pressure = EAVE_BP) %>%
  mutate(
    smoking_status = as.character(smoking_status),
    blood_pressure = as.character(blood_pressure)
  ) %>%
  mutate(
    smoking_status = case_when(
      smoking_status == "Unknown" ~ NA_character_,
      TRUE ~ smoking_status
    ),
    blood_pressure = case_when(
      blood_pressure == "No Investigation" ~ NA_character_,
      TRUE ~ blood_pressure
    )
  ) %>%
  mutate(
    smoking_status = as.factor(smoking_status),
    blood_pressure = as.factor(blood_pressure)
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
  rename(date_death = NRS.Date.Death) %>%
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
  select(EAVE_LINKNO, first_dose, second_dose, third_dose, fourth_dose, fifth_dose,
    vacc_dose_number, vacc_product_name
  ) %>%
  mutate(
    vacc_product_name = gsub("Covid-19 Vaccine AstraZeneca", "AZ", vacc_product_name),
    vacc_product_name = gsub("Covid-19 mRNA Vaccine Pfizer", "PB", vacc_product_name),
    vacc_product_name = gsub("Covid-19 mRNA Vaccine Moderna", "Mo", vacc_product_name),
    vacc_product_name = gsub("Covid-19 mRNA Vaccine Spikevax Bivalent Moderna", "Mo", vacc_product_name),
    vacc_product_name = gsub("Covid-19 mRNA Vaccine Comirnaty Bivalent Pfizer", "PB", vacc_product_name),
    vacc_product_name = gsub("Covid-19 Vaccine Novavax", "No", vacc_product_name)
  ) %>%
  # 4th and 5th dose vaccinations with the moderna bivalent vaccine are not currently picked up
  mutate(vacc_dose_number = case_when(
    !is.na(fifth_dose) & is.na(vacc_dose_number) ~ 5,
    !is.na(fourth_dose) & is.na(vacc_dose_number) ~ 4,
    TRUE ~ vacc_dose_number
  )) %>%
  pivot_wider(
    id_cols = c(
      "EAVE_LINKNO", "first_dose", "second_dose", "third_dose", "fourth_dose", "fifth_dose"
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
  rename(mean_household_age = ave_hh_age) %>%
  mutate(num_ppl_household = cut(n_hh,
    breaks = c(0, 1, 2, 5, 10, 30, 100, max(n_hh)),
    labels = c("1", "2", "3-5", "6-10", "11-30", "31-100", "101+")
  )) %>%
  mutate(mean_household_age = if_else(is.na(mean_household_age), mean(mean_household_age, na.rm = T), mean_household_age))

# Whole genome sequencing
wgs = readRDS(paste0(Location, "EAVE/GPanalysis/data/WGS_latest.rds")) 

# Endpoints
endpoints = covid_hospitalisations %>%
  rename(covid_hosp_date = ADMISSION_DATE) %>%
  mutate(covid_hosp = 1) %>%
  full_join(all_deaths %>%
    filter(covid_death == 1) %>%
    rename(covid_death_date = date_death) %>%
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
  # Original cohort is from approximately 2 years ago. Update their age.
  mutate(age = age + 2) %>%
  filter(age >= 5) %>%
  mutate(
    age_group = cut(age,
      breaks = c(5, 12, 16, 75, Inf),
      labels = c("5-11", "12-15", "16-74", "75+"),
      right = FALSE
    ),
    age_group_2 = cut(age,
                    breaks = c(5, 18, 75, Inf),
                    labels = c("5-17", "18-74", "75+"),
                    right = FALSE
    ),
    age_group_3 = cut(age,
      breaks = c(5, 12, 16, 18, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, Inf),
      labels = c(
        "5-11", "12-15", "16-17", "18-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54",
        "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85+"
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

# Count all combinations of vaccines that anyone has received
vacc_seq_start_count = count(df_cohort, vacc_seq_start)

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
      age_group == "5-11" & date_vacc_1 <= study_start ~ 1,
      age_group == "12-15" & date_vacc_2 <= study_start ~ 1,
      age_group == "16-74" & date_vacc_3 <= study_start ~ 1,
      age_group == "75+" & date_vacc_4 <= study_start ~ 1,
      TRUE ~ 0
    )
  ) %>%
  mutate(fully_vaccinated = factor(fully_vaccinated))

df_cohort = mutate(df_cohort,
    # Mixed vaccine status at study start
    vs_mixed_start = case_when(
      num_doses_start == 0 ~ 'uv',
      num_doses_start == 1 ~ vacc_seq_start,
      num_doses_start == 2 & vacc_seq_start %in% c('AZ_AZ', 'PB_PB', 'Mo_Mo') ~ vacc_seq_start,
      num_doses_start == 2 & vacc_type_1 == 'AZ' & vacc_type_2 %in% c('PB', 'Mo') ~ 'v2_AZ_mrna',
      num_doses_start == 2 & vacc_type_1 %in% c('PB', 'Mo') & vacc_type_2 %in% c('PB', 'Mo') ~ 'v2_mixed_mrna',
      num_doses_start == 2 ~ 'v2_other',
      num_doses_start == 3 & age_group %in% c('5-11', '12-15') ~ '3+',
      num_doses_start == 3 & vacc_seq_start %in% 
        c('AZ_AZ_PB', 'AZ_AZ_Mo', 'AZ_AZ_AZ', 'PB_PB_PB', 'PB_PB_Mo', 'Mo_Mo_PB', 'Mo_Mo_Mo') ~ vacc_seq_start,
      num_doses_start == 3 & vacc_type_3 == 'PB' ~ 'other_mixed_2_dose_PB',
      num_doses_start == 3 & vacc_type_3 == 'Mo' ~ 'other_mixed_2_dose_Mo',
      num_doses_start == 3 ~ 'v3_other',
      num_doses_start == 4 & age_group %in% c('5-11', '12-15', '16-74') ~ '4+',
      num_doses_start == 4 & vacc_type_1 == 'AZ' & vacc_type_2 == 'AZ' ~ 'v4_AZ_AZ_any',
      num_doses_start == 4 & vacc_type_1 == 'PB' & vacc_type_2 == 'PB' ~ 'v4_PB_PB_any',
      num_doses_start == 4 ~ 'v4_other',
      num_doses_start > 4 ~ '5+'
    ),
    # Mixed vaccine status according to most recent data
    vs_mixed_recent = case_when(
      num_doses_recent == 0 ~ 'uv',
      num_doses_recent == 1 ~ vacc_seq_recent,
      num_doses_recent == 2 & vacc_seq_recent %in% c('AZ_AZ', 'PB_PB', 'Mo_Mo') ~ vacc_seq_recent,
      num_doses_recent == 2 & vacc_type_1 == 'AZ' & vacc_type_2 %in% c('PB', 'Mo') ~ 'v2_AZ_mrna',
      num_doses_recent == 2 & vacc_type_1 %in% c('PB', 'Mo') & vacc_type_2 %in% c('PB', 'Mo') ~ 'v2_mixed_mrna',
      num_doses_recent == 2 ~ 'v2_other',
      num_doses_recent == 3 & age_group %in% c('5-11', '12-15') ~ '3+',
      num_doses_recent == 3 & vacc_seq_recent %in% 
        c('AZ_AZ_PB', 'AZ_AZ_Mo', 'AZ_AZ_AZ', 'PB_PB_PB', 'PB_PB_Mo', 'Mo_Mo_PB', 'Mo_Mo_Mo') ~ vacc_seq_recent,
      num_doses_recent == 3 & vacc_type_3 == 'PB' ~ 'other_mixed_2_dose_PB',
      num_doses_recent == 3 & vacc_type_3 == 'Mo' ~ 'other_mixed_2_dose_Mo',
      num_doses_recent == 3 ~ 'v3_other',
      num_doses_recent == 4 & age_group %in% c('5-11', '12-15', '16-74') ~ '4+',
      num_doses_recent == 4 & vacc_type_1 == 'AZ' & vacc_type_2 == 'AZ' ~ 'v4_AZ_AZ_any',
      num_doses_recent == 4 & vacc_type_1 == 'PB' & vacc_type_2 == 'PB' ~ 'v4_PB_PB_any',
      num_doses_recent == 4 ~ 'v4_other',
      num_doses_recent > 4 ~ '5+'
    )
  ) 


# Add household characteristics
df_cohort = Cohort_Household %>%
  select(EAVE_LINKNO, num_ppl_household, mean_household_age) %>%
  right_join(df_cohort)

# Add deaths from any cause
df_cohort = all_deaths %>%
  select(EAVE_LINKNO, date_death, covid_death) %>%
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


# This step is acquiring duplicates!
# Add in time since last positive test and variant of last positive test
df_cohort = filter(cdw, date_ecoss_specimen < study_start & test_result == "POSITIVE") %>%
  select(EAVE_LINKNO, date_ecoss_specimen) %>%
  left_join(wgs %>% 
    select(
      EAVE_LINKNO,
      date_ecoss_specimen = Collection_Date,
      last_positive_variant = lineage) %>%
    mutate(
      last_positive_variant = case_when(
        last_positive_variant == "B.1.1.7" ~ "Alpha",
        last_positive_variant == "B.1.617.2" ~ "Delta",
        last_positive_variant %in% c('BA.1', 'BA.2', 'BA.3', 'BA.4', 'BA.5') ~ "Omicron",
        TRUE ~ "Other"))
    ) %>%
  arrange(EAVE_LINKNO, desc(date_ecoss_specimen)) %>%
  filter(!duplicated(EAVE_LINKNO)) %>%
  right_join(df_cohort) %>%
  mutate(
    last_positive_test = as.numeric(study_start - date_ecoss_specimen),
    last_positive_variant = case_when( 
      is.na(date_ecoss_specimen) ~ 'never_positive',
      is.na(last_positive_variant) ~ "not_sequenced",
      TRUE ~ last_positive_variant)
) %>%
  mutate(last_positive_test_group = as.character(cut(last_positive_test,
    breaks = c(0, 92, 183, Inf),
    labels = c("0-13_weeks", "14-26_weeks", "27+_weeks"),
    right = FALSE
  ))) %>%
  mutate(last_positive_test_group = ifelse(is.na(last_positive_test_group), 'never_positive', last_positive_test_group)) %>%
  mutate(last_positive_test_group = factor(last_positive_test_group, c("never_positive", "0-13_weeks", "14-26_weeks", "27+_weeks")),
         last_positive_variant = factor(last_positive_variant, c('never_positive', 'not_sequenced', 'Alpha', 'Delta', 'Omicron')))

# Add number of PCR tests in last 6 months before study_start
df_cohort = filter(cdw, date_ecoss_specimen >= study_start %m-% months(6) & date_ecoss_specimen < study_start) %>%
  count(EAVE_LINKNO) %>%
  rename(num_tests_6m = n) %>%
  select(EAVE_LINKNO, num_tests_6m) %>%
  right_join(df_cohort) %>%
  mutate(num_tests_6m = replace_na(num_tests_6m, 0)) %>%
  mutate(num_tests_6m_group = cut(num_tests_6m,
    breaks = c(0, 1, 2, 3, 4, 10, Inf),
    labels = c("0", "1", "2", "3", "4-9", "10+"),
    right = FALSE
  ))


# Add number of positive PCR tests in last 6 months before study_start
df_cohort = filter(cdw, date_ecoss_specimen >= study_start %m-% months(6) & date_ecoss_specimen < study_start & test_result == 'POSITIVE') %>%
  count(EAVE_LINKNO) %>%
  rename(num_pos_tests_6m = n) %>%
  select(EAVE_LINKNO, num_pos_tests_6m) %>%
  right_join(df_cohort) %>%
  mutate(num_pos_tests_6m = replace_na(num_pos_tests_6m, 0)) %>%
  mutate(num_pos_tests_6m = cut(num_pos_tests_6m,
                                  breaks = c(0, 1, 2,  Inf),
                                  labels = c("0", "1", "2+"),
                                  right = FALSE
  ))


# Add shielding list
df_cohort = ever_shielding %>%
  right_join(df_cohort) %>%
  mutate(shielding = replace_na(shielding, 0))


# Re-weight, taking into account contact with healthcare
# Load in other datasets that indicate contact with healthcare
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
  pis$EAVE_LINKNO,
  wgs$EAVE_LINKNO
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

cols = setdiff(q_names, c("Q_BMI", "Q_ETHNICITY"))
cols = c("covid_hosps", "covid_death", cols)

df_cohort = mutate_at(df_cohort, cols, ~ as.numeric(.))

df_cohort = mutate_at(df_cohort, cols, ~ case_when(
  is.na(.) ~ 0,
  TRUE ~ .
)) %>%
  mutate(
    # Ethnicity coding taken from QCovid algorithm
    Q_ETHNICITY = case_when(
      Q_ETHNICITY %in% c(1, 2, 3) ~ 'White',
      Q_ETHNICITY %in% c(4, 5, 6, 7) ~ 'Mixed',
      Q_ETHNICITY %in% c(8, 9, 10, 11, 15) ~ 'Asian',
      Q_ETHNICITY %in% c(12, 13, 14) ~ 'Black',
      Q_ETHNICITY == 16 ~ 'Other',
      TRUE ~ NA_character_),
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
  mutate(n_risk_gps = as.factor(n_risk_gps),
         num_doses_start = factor(num_doses_start),
         num_doses_recent = factor(num_doses_recent))

df_cohort = mutate_at(df_cohort, c("smoking_status", "blood_pressure"), ~ as.character(.)) %>%
  mutate(
    smoking_status = case_when(
      smoking_status == "Ex Smoker" ~ "Ex-smoker",
      smoking_status == "Non Smoker" ~ "Non-smoker",
      TRUE ~ smoking_status
    ),
    blood_pressure = case_when(
      blood_pressure == "Very High" ~ "Very high",
      TRUE ~ blood_pressure
    ),
    urban_rural_class = case_when(
      urban_rural_class == "1 Large Urban Areas" ~ "Large Urban Areas",
      urban_rural_class == "2 Other Urban Areas" ~ "Other Urban Areas",
      urban_rural_class == "3 Accessible Small Towns" ~ "Accessible Small Towns",
      urban_rural_class == "4 Remote Small Towns" ~ "Remote Small Towns",
      urban_rural_class == "5 Accessible Rural" ~ "Accessible Rural",
      urban_rural_class == "6 Remote Rural" ~ "Remote Rural"
    )
  ) %>%
  mutate_at(c("smoking_status", "blood_pressure"), ~ as.factor(.)) %>%
  mutate(
    smoking_status = fct_relevel(smoking_status, "Non-smoker"),
    blood_pressure = fct_relevel(blood_pressure, "Normal", "Low", "High", "Very high"),
    Q_HOME_CAT = fct_relevel(Q_HOME_CAT, "Neither"),
    Q_LEARN_CAT = fct_relevel(Q_LEARN_CAT, "Neither"),
    Q_DIAG_CKD_LEVEL = fct_relevel(Q_DIAG_CKD_LEVEL, "No CKD"),
    urban_rural_class = fct_relevel(
      urban_rural_class, "Large Urban Areas",
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

#### df_cohort has the following columns:

# [1] "EAVE_LINKNO"              "shielding"                "num_pos_tests_6m"         "num_tests_6m"             "date_ecoss_specimen"     
# [6] "last_positive_variant"    "covid_hosps"              "date_death"               "covid_death"              "num_ppl_household"       
# [11] "mean_household_age"       "sex"                      "age"                      "simd2020_sc_quintile"     "DataZone"                
# [16] "urban_rural_class"        "age_group"                "age_group_2"              "age_group_3"              "eave_weight"             
# [21] "Q_BMI"                    "Q_ETHNICITY"              "Q_HOME_CAT"               "Q_LEARN_CAT"              "Q_DIAG_CKD_LEVEL"        
# [26] "Q_DIAG_AF"                "Q_DIAG_ASTHMA"            "Q_DIAG_BLOOD_CANCER"      "Q_DIAG_CCF"               "Q_DIAG_CEREBRALPALSY"    
# [31] "Q_DIAG_CHD"               "Q_DIAG_CIRRHOSIS"         "Q_DIAG_CONGEN_HD"         "Q_DIAG_COPD"              "Q_DIAG_DEMENTIA"         
# [36] "Q_DIAG_EPILEPSY"          "Q_DIAG_FRACTURE"          "Q_DIAG_HIV_AIDS"          "Q_DIAG_IMMU"              "Q_DIAG_NEURO"            
# [41] "Q_DIAG_PARKINSONS"        "Q_DIAG_PULM_HYPER"        "Q_DIAG_PULM_RARE"         "Q_DIAG_PVD"               "Q_DIAG_RA_SLE"           
# [46] "Q_DIAG_RESP_CANCER"       "Q_DIAG_SEV_MENT_ILL"      "Q_DIAG_SICKLE_CELL"       "Q_DIAG_STROKE"            "Q_DIAG_VTE"              
# [51] "n_risk_gps"               "Q_DIAG_DIABETES_1"        "Q_DIAG_DIABETES_2"        "smoking_status"           "blood_pressure"          
# [56] "bmi_cat"                  "date_vacc_1"              "date_vacc_2"              "date_vacc_3"              "date_vacc_4"             
# [61] "date_vacc_5"              "vacc_type_1"              "vacc_type_2"              "vacc_type_3"              "vacc_type_4"             
# [66] "vacc_type_5"              "vacc_type_NA"             "vacc_type_6"              "vacc_type_7"              "num_doses_recent"        
# [71] "num_doses_start"          "vacc_seq_recent"          "vacc_seq_start"           "mixed_vacc_recent"        "mixed_vacc_start"        
# [76] "vs_recent"                "vs_start"                 "fully_vaccinated"         "vs_mixed_start"           "vs_mixed_recent"         
# [81] "covid_hosp_ever"          "last_positive_test"       "last_positive_test_group" "num_tests_6m_group" 


#### Variable descriptions

# EAVE_LINKNO is an individual identifier

# shielding is whether they have ever been on the shielding list

# num_pos_tests_6m is the number of positive tests they have had in the last 6 months before study_start
# Its levels are: 0, 1, 2+

# num_tests_6m is integer - number of tests they have had in the last 6 months before study_start

# num_tests_6m_group has the following levels:
# 0, 1, 2, 3, 4-9, 10+

# last_positive_test_group has the following levels:
# never_positive, 0-13_weeks, 14-26_weeks, 27+_weeks

# covid_hosps is integer - number of covid hospitalisations they have ever had

# simd2020_sc_quintile is their quintile of the Scottish index of multiple deprivation
# 1 is the most deprived, 5 is the least deprived

# urban_rural_class has the following levels:
# Large Urban Areas
# Other Urban Areas
# Accessible Small Towns
# Remote Small Towns
# Accessible Rural
# Remote Rural

# bmi_cat is their body mass index cateogry, and has the following levels:
# 10-18.5, 18.5-25, 25-30, 30-35, 35-40, 40+

# age_group has the following levels:
# 5-11, 12-15, 16-74, 75+

# age_group_2 has the following levels:
# 5-17, 18-74, 75+

# age_group_3 has the following levels:
# 5-11, 12-15, 16-17, 18-24, 25-29, 30-34, 35-39, 40-44, 
# 45-49, 50-54, 55-59, 60-64, 65-69, 70-74, 75-79, 80-84, 85+

# num_doses_recent is the number of vaccine doses they have had according to most recent data
# num_doses_start is the number of doses they had at study_start

# vacc_seq_recent is the sequence of vaccine types they have received according to most recent data
# vacc_seq_start is the sequence of vaccine types they had received at study_start

# mixed_vacc_start is a binary variable indicating whether they have had any mixed vaccine type at study start

# vs_recent is the dose number and type of the last dose they had according to most recent data
# vs_start is the dose number and type of the last dose they had at study_start

# fully_vaccinated is a binary variable indicating whether they have had the recommended
# vaccine schedule for their age group at study_start

# vs_mixed_start is a more coarse_grained version of vacc_seq_start, with categories of vaccine sequences
# detailed in the draft tables document
# vs_mixed_recent is similar but for vacc_seq_recent

# covid_hosp_ever is whether they have ever had a covid hospitalisation

