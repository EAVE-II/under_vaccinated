######################################################################
## Title: [Insert full title of paper]
## Code author: Steven Kerr steven.kerr@ed.ac.uk
## Description: Load in all data and clean it
######################################################################

library(tidyverse)
library(lubridate)

Location = "/conf/"

setwd("/conf/EAVE/GPanalysis/analyses/under_vaccinated")

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

## Demographics
EAVE_cohort <- readRDS(paste0(Location, "EAVE/GPanalysis/outputs/temp/Cohort_Demog_Endpoints_Dates2021-07-28.rds")) %>%
  select(
    EAVE_LINKNO,
    sex = Sex,
    age = ageYear,
    urban_rural_6cat = ur6_2016_name,
    simd2020_sc_quintile) %>%
  filter(!duplicated(EAVE_LINKNO)) %>%
  mutate(
    sex = recode_factor(sex, F = "Female", M = "Male"),
    urban_rural_6cat = case_when(
      urban_rural_6cat == "1 Large Urban Areas" ~ "Large Urban Areas",
      urban_rural_6cat == "2 Other Urban Areas" ~ "Other Urban Areas",
      urban_rural_6cat == "3 Accessible Small Towns" ~ "Accessible Small Towns",
      urban_rural_6cat == "4 Remote Small Towns" ~ "Remote Small Towns",
      urban_rural_6cat == "5 Accessible Rural" ~ "Accessible Rural",
      urban_rural_6cat == "6 Remote Rural" ~ "Remote Rural"
    ),
    urban_rural_6cat = fct_relevel(
      urban_rural_6cat, 
      "Large Urban Areas",
      "Other Urban Areas",
      "Accessible Small Towns",
      "Remote Small Towns",
      "Accessible Rural",
      "Remote Rural"
    ),
    urban_rural_2cat = case_when(
      urban_rural_6cat %in% c("Large Urban Areas", "Other Urban Areas") ~ "Urban",
      urban_rural_6cat %in% c("Accessible Small Towns", "Remote Small Towns", "Accessible Rural", "Remote Rural") ~ "Rural"
    ),
    simd2020_sc_quintile = as.factor(as.character(simd2020_sc_quintile))
  ) %>%
  # Original extract was in 2020 - add 2 to age
  mutate(age = age + 2)

EAVE_Weights <- readRDS(paste0(Location, "EAVE/GPanalysis/outputs/temp/CR_Cohort_Weights.rds"))

# Household information
household <- readRDS(paste0(Location, "EAVE/GPanalysis/outputs/temp/Cohort_Household.rds")) %>%
  rename(mean_household_age = ave_hh_age) %>%
  mutate(
    num_ppl_household = cut(
      n_hh,
      breaks = c(0, 1, 2, 5, 10, 30, 100, max(n_hh)),
      labels = c("1", "2", "3-5", "6-10", "11-30", "31-100", "101+")
    )
  ) %>%
  mutate(
    mean_household_age = if_else(is.na(mean_household_age), mean(mean_household_age, na.rm = T), mean_household_age)
  ) %>%
  select(
    EAVE_LINKNO,
    mean_household_age,
    num_ppl_household)


## Clinical characteristics

# QCovid risk groups
qcovid_rg <- readRDS(paste0(Location, "EAVE/GPanalysis/data/cleaned_data/QCOVID_feb22.rds")) %>%
  select(-Age, -Sex) %>%
  mutate(Q_BMI = as.numeric(Q_BMI)) %>%
  # Im taking 10 to be the lowest feasible human BMI, and 100 the largest
  # There are lots of odd values for BMI, which probably happens because issues
  # with units that calculation is done in.
  mutate(
    Q_BMI = ifelse(Q_BMI < 10 | Q_BMI > 100, NA_real_, Q_BMI),
    # Ethnicity coding taken from QCovid algorithm
    Q_ETHNICITY = case_when(
      Q_ETHNICITY %in% c(1, 2, 3) ~ "White",
      Q_ETHNICITY %in% c(4, 5, 6, 7) ~ "Mixed",
      Q_ETHNICITY %in% c(8, 9, 10, 11, 15) ~ "Asian",
      Q_ETHNICITY %in% c(12, 13, 14) ~ "Black",
      Q_ETHNICITY == 16 ~ "Other",
      TRUE ~ NA_character_
    ),
    Q_HOME_CAT = case_when(
      Q_HOME_CAT == 1 ~ "Care home",
      Q_HOME_CAT == 2 ~ "Homeless"
    ),
    Q_LEARN_CAT = case_when(
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
  # Drop diabetes 1 and 2 cateogries because they are not useable
  select(-Q_DIAG_DIABETES_1, -Q_DIAG_DIABETES_2)

# Get old QCovid risk groups. We will use this for diabetes only, because diabetes grouping
# from more recent QCovid extract are not useable
qcovid_diabetes <- readRDS("/conf/EAVE/GPanalysis/progs/CR/Vaccine/output/temp/Qcovid.rds") %>%
  select(EAVE_LINKNO, Q_DIAG_DIABETES_1, Q_DIAG_DIABETES_2)

qcovid_rg <- full_join(qcovid_rg, qcovid_diabetes)

# Blood pressure and smoking status
eave_rg <- readRDS(paste0(Location, "EAVE/GPanalysis/outputs/temp/CR_Cohort_RG_EAVE_BP_Smoke.rds")) %>%
  filter(!duplicated(EAVE_LINKNO)) %>%
  select(
    EAVE_LINKNO,
    smoking_status = EAVE_Smoking_Status_Worst,
    blood_pressure = EAVE_BP
  ) %>%
  mutate(
    smoking_status = recode_factor(
      smoking_status,
      "Non-smoker" = "Non-smoker",
      "Ex Smoker" = "Ex-smoker",
      "Non Smoker" = "Non-smoker",
      "Unknown" = NA_character_
    ),
    blood_pressure = recode_factor(blood_pressure, "Very High" = "Very high", "No Investigation" = NA_character_),
    blood_pressure = fct_relevel(blood_pressure, "Normal", "Low", "High", "Very high")
  )


## Tests
# PCR
cdw <- readRDS(paste0(Location, "EAVE/GPanalysis/data/CDW_full.rds")) %>%
  select(
    EAVE_LINKNO,
    specimen_date = date_ecoss_specimen,
    test_result) %>%
  mutate(
    specimen_date = as.Date(specimen_date)
  ) %>%
  filter(specimen_date <= Sys.Date())


# All deaths
all_deaths <- readRDS(paste0(Location, "EAVE/GPanalysis/data/all_deaths.rds")) %>%
  rowwise() %>%
  mutate(covid_death = ifelse(rowSums(across(UNDERLYING_CAUSE_OF_DEATH:CAUSE_OF_DEATH_CODE_9, ~ .x %in% c("U071", "U072")), na.rm = T) > 0, 1, 0)) %>%
  select(
    EAVE_LINKNO,
    death_date = NRS.Date.Death,
    covid_death) 


## Hospitalisations
all_hospitalisations <- readRDS(paste0(Location, "EAVE/GPanalysis/data/automated_any_hospitalisation_post_01022020.rds")) %>%
  select(
    EAVE_LINKNO,
    hosp_date = admission_date)


# SMR - Scottish morbidity record
smr <- readRDS(paste0(Location, "EAVE/GPanalysis/data/SMR01_allstays.rds")) %>%
  select(
    EAVE_LINKNO,
    hosp_date = ADMISSION_DATE,
    covid_main_diag_admit,
    covid_main_other_ep) %>%
  mutate(
    hosp_date = as.Date(hosp_date)
  )

## Vaccinations
Vaccinations <- readRDS(paste0(Location, "/EAVE/GPanalysis/data/cleaned_data/C19vaccine_dvprod_cleaned.rds"))

Vaccinations = Vaccinations %>%
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


## Shielding
ever_shielding <- readRDS(paste0(Location, "/EAVE/GPanalysis/data/cleaned_data/C19vaccine_dvprod_cleaned_incl_cohorts_20220926.rds")) %>%
  select(EAVE_LINKNO, shielding) %>%
  arrange(EAVE_LINKNO, desc(shielding)) %>%
  filter(!duplicated(EAVE_LINKNO))


## Whole genome sequencing
wgs <- readRDS(paste0(Location, "EAVE/GPanalysis/data/WGS_latest.rds")) %>%
  select(
     EAVE_LINKNO,
     specimen_date = Collection_Date,
     lineage)
