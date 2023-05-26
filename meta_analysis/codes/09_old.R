#old code cox read in

#add in the ni estimates
z <- ni_ests %>%  
  mutate(variance=std.error^2) %>% 
  mutate(country="ni",
         index=paste0(xvar,xlbl),
         title=z_title, model=z_model) %>% 
  dplyr::select(country,index,xvar,xlbl,estimate,variance, title, model) %>% 
  dplyr::rename(coef=estimate) %>% 
  mutate(xvar=case_when(xvar=="mdm_quintile" ~ "wimd2019",
                        xvar=="n_risk_groups" ~ "qc_count_cat",
                        TRUE ~ xvar)  ) %>% 
  mutate(xlbl=case_when(xlbl=="M" ~ "male",
                        xlbl=="5 - least deprived" ~ "5_least",
                        xlbl=="85+" ~ "85plus",
                        xlbl=="Rural" ~ "rural",
                        xlbl=="5+" ~ "5plus",
                        TRUE~xlbl)) %>% 
  mutate(xlbl=gsub("-","_",xlbl))
if (z_title == "Hospital Death Age 5-15") {z <- z %>% filter(!(xvar=="age_5y_cat" & xlbl=="16_24")) %>% 
  filter(!(xvar=="qc_count_cat" & xlbl %in% c("3","4","5plus"))) %>%
  mutate(xlbl=if_else(xvar=="qc_count_cat" & xlbl=="2", "2plus", xlbl)) %>% 
  mutate(xvar=if_else(xvar=="qc_count_cat","qc_count_cat3", xvar))
}
df_coefs <- bind_rows(df_coefs, z)

##########################################

#old cox table read in when 

if ("england_1" %in% names(tables)) {
  z1 <- tables[["england"]]
  z2 <- tables[["england_1"]]
  z <- z1 %>% full_join(z2, by=c("xvar","xlbl"), suffix=c("_1","_2"))
  z <- z %>% mutate(across(total_n_1:covid_mcoa_em_hosp_death_p_2, ~replace_na(.,0))) %>% 
    mutate(total_n = total_n_1 + total_n_2,
           covid_mcoa_em_hosp_death_n = covid_mcoa_em_hosp_death_n_1+covid_mcoa_em_hosp_death_n_2) %>% 
    dplyr::select(xvar, xlbl, total_n, covid_mcoa_em_hosp_death_n)
  tables[["england"]] <- z
  tables[["england_1"]] <- NULL
}

######################################################
#02 table plots cox

z0 <- z_df %>% filter(label=="Number of risk groups") %>% 
  mutate(levels=if_else(levels %in% c("0","1"), levels, "2+")) %>% 
  group_by(country, label,levels) %>% 
  dplyr::summarise(across(pyears:events, ~sum(., na.rm=T))) %>% ungroup() 
z_df <- z_df %>% filter(label !="Number of risk groups") %>% 
  bind_rows(z0) %>% 
  arrange(country, label, levels)


#impute the suppressed events
z <- df_table %>% filter(is.na(events))
z_notw <- df_table %>% filter(label %in% unique(z$label) & levels %in% unique(z$levels) & !is.na(events))
z_notw <- z_notw %>% group_by(label, levels) %>% dplyr::summarise(across(pyears:events, ~sum(.))) %>% 
  ungroup() %>% 
  mutate(rate = events/pyears)
z  <- z %>%  left_join( dplyr::select(z_notw, label, levels, rate )) %>% mutate(events = pyears*rate) %>% 
  mutate(events = case_when(events >= 9 ~ 9,
                            events <=1 ~ 1,
                            TRUE ~ round(events))) %>% 
  dplyr::select(-rate)

df_table_imp <- df_table %>% filter(!is.na(events)) %>% 
  bind_rows(z) %>% 
  arrange(country,label,levels)

df_table_imp <- df_table %>% filter(!(label %in% c("BMI", "Number of people in household", "Non Covid Hospitalisation")))

df_table_imp <- df_table_imp %>% mutate(rate_1000 = events/pyears*1000) %>% 
  group_by(country,label) %>% mutate(total_p = total_n/sum(total_n)*100,
                                     py_p = pyears/sum(pyears)*100) %>% 
  ungroup()


###########################
# 01 read data uv lr

#add in the ni estimates
z_age_gp <- gsub("Undervaccination Age ", "", z_title)
if (z_age_gp=="5-11") z_age_gp <- paste0("0",z_age_gp)
z_model_ni <- case_when(z_model=="min_adj" ~ "madj",
                        TRUE ~ "eadj")
z <- ni_ests %>% filter(age_group==z_age_gp & est_type==z_model_ni) %>% 
  mutate(variance=std.error^2) %>% 
  mutate(country="ni",
         index=paste0(xvar,xlbl),
         title=z_title, model=z_model) %>% 
  dplyr::select(country,index,xvar,xlbl,estimate,variance, title, model) %>% 
  dplyr::rename(coef=estimate) %>% 
  mutate(xvar=case_when(xvar=="mdm_quintile" ~ "wimd2019",
                        xvar=="n_risk_groups" ~ "qc_count_cat",
                        TRUE ~ xvar)  ) %>% 
  mutate(xlbl=case_when(xlbl=="M" ~ "male",
                        xlbl=="5 - least deprived" ~ "5_least",
                        xlbl=="85+" ~ "85plus",
                        xlbl=="80-84" ~ "80_84",
                        xlbl=="Rural" ~ "rural",
                        xlbl=="5+" ~ "5plus",
                        TRUE~xlbl)) %>% 
  mutate(xlbl=gsub("-","_",xlbl))
df_coefs <- bind_rows(df_coefs, z)

##########################################

#changes to labels - should be temporary
if (z_title == "Undervaccination Age 5-11" | z_title == "Undervaccination Age 12-15") {
  z <- df_coefs %>% filter(!(country=="wales" & xvar=="qc_count_cat" & xlbl %in% c("04","03","05plus"))) %>% 
    mutate(xvar=if_else(xvar=="qc_count_cat","qc_count_cat3", xvar)) %>% 
    mutate(xlbl=if_else(country=="wales" & xvar=="qc_count_cat3" & xlbl == "02", "2plus", xlbl)) %>% 
    mutate(xlbl=if_else(country=="wales" & xvar=="qc_count_cat3" & xlbl == "01", "1", xlbl))
  df_coefs <- z
  z <- df_coefs %>% filter(!(country=="ni" & xvar=="qc_count_cat3" & xlbl %in% c("4","3","5plus"))) %>% 
    mutate(xlbl=if_else(country=="ni" & xvar=="qc_count_cat3" & xlbl == "2", "2plus", xlbl))
  df_coefs <- z
}

############################################
# 01 read data demog

if(i =="scotland") {
  z_variables <- c("Total","Sex","Ethnicity","SIMD quintile","Urban rural classification",
                   "BMI","Number of people in household", "Number of risk groups")
  z <- z_df
  while(any(z$label=="")) z <- z %>% mutate(label = if_else(label=="", lag(label), label))
  z <- z %>% filter(label %in% z_variables)
  z <- z %>% pivot_longer(cols=starts_with("X"))
  z <- z %>% mutate(location = regexpr(" ", value)) %>% 
    mutate(value = as.numeric(substring(value,1,location-1)))
  if (z_title == "Undervaccination Age 75+") z <- z %>% mutate(vacc = if_else(name %in% c("X4","X5"), "full_v","under_v"))
  if (z_title == "Undervaccination Age 16-74") z <- z %>% mutate(vacc = if_else(name %in% c("X3","X4","X5"), "full_v","under_v"))
  if (z_title == "Undervaccination Age 12-15") z <- z %>% mutate(vacc = if_else(name %in% c("X2","X3","X4","X5"), "full_v","under_v"))
  if (z_title == "Undervaccination Age 5-11") z <- z %>% mutate(vacc = if_else(name %in% c("X1","X2","X3","X4","X5"), "full_v","under_v"))
  z <- z %>% group_by(country, label, levels, vacc) %>% 
    dplyr::summarise(N=sum(value)) %>% ungroup()
  z <- z %>% pivot_wider(id_cols = country:levels, names_from=vacc, values_from=N)
  z <- z %>% mutate(total_n = full_v + under_v) %>% 
    dplyr::rename(undervac_n = under_v) %>% 
    mutate(p_undervac=undervac_n/total_n*100) %>% 
    group_by(label) %>% mutate(total_p = total_n/sum(total_n)*100) %>% ungroup()
  z_df <- z
  #for household and urban rural aggregate
  z <- z_df %>% filter(label %in% c("Number of people in household","Urban rural classification")) %>% 
    mutate(levels=if_else(levels %in% c("101+","11-30","31-100"), "11plus",levels),
           levels=if_else(levels %in% c("Accessible Rural","Remote Rural"), "Rural",levels),
           levels=if_else(levels %in% c("Accessible Small Towns","Large Urban Areas",
                                        "Other Urban Areas", "Remote Small Towns"), "Urban",levels)) %>% 
    group_by(country, label, levels) %>% 
    dplyr::summarise(across(full_v:total_p, ~ sum(.))) %>% ungroup() %>% 
    mutate(p_undervac = undervac_n/total_n*100)
  
  z_df <- z_df %>% 
    filter(!(label %in% c("Number of people in household","Urban rural classification"))) %>%
    bind_rows(z)
}
############
#for household and urban rural aggregate
z <- z_df %>% filter(label %in% c("Number of people in household","Urban rural classification")) %>% 
  mutate(levels=if_else(levels %in% c("101+","11-30","31-100"), "11plus",levels),
         levels=if_else(levels %in% c("Accessible Rural","Remote Rural"), "Rural",levels),
         levels=if_else(levels %in% c("Accessible Small Towns","Large Urban Areas",
                                      "Other Urban Areas", "Remote Small Towns"), "Urban",levels)) %>% 
  group_by(country, label, levels) %>% 
  dplyr::summarise(across(full_v:total_p, ~ sum(.))) %>% ungroup() %>% 
  mutate(p_undervac = undervac_n/total_n*100)

z_df <- z_df %>% 
  filter(!(label %in% c("Number of people in household","Urban rural classification"))) %>%
  bind_rows(z)

#######################################################################

#02 table plots

if (g_age %in% c("5_11","12_15")  )  {z_df <- df_table %>% filter(label !="BMI") %>% filter(label!= "Age Group")  
z0 <- z_df %>% filter(label=="Number of risk groups") %>% 
  mutate(levels=if_else(levels %in% c("0","1"), levels, "2+")) %>% 
  group_by(country, label,levels) %>% 
  dplyr::summarise(across(full_v:total_n, ~sum(., na.rm=T))) %>% ungroup() %>% 
  mutate(p_undervac = undervac_n/total_n*100) %>% 
  group_by(country) %>% mutate(total_p = total_n/sum(total_n)*100) %>% ungroup()
z_df <- z_df %>% filter(label !="Number of risk groups") %>% 
  bind_rows(z0) %>% 
  arrange(country, label, levels)
z_df <- z_df %>% filter(!(label=="Number of people in household" & levels=="1"))
} 
if (g_age %in% c("16_74")  )  {z_df <- df_table %>%  
  filter(!(country == "ni" & label == "Age Group" & levels %in% c("05_11","12_15","75_79","80_84","85+")))
z_bmi_miss_p <- z_df %>% filter(country=="scotland" & label=="BMI" & levels=="Missing") %>% pull(total_p)
z_df <- z_df %>% mutate(total_p = if_else(country=="scotland" & label=="BMI" & levels!="Missing", total_p/(1-z_bmi_miss_p/100), total_p))
z_bmi_miss_p <- z_df %>% filter(country=="england" & label=="BMI" & levels=="Missing") %>% pull(total_p)
z_df <- z_df %>% mutate(total_p = if_else(country=="england" & label=="BMI" & levels!="Missing", total_p/(1-z_bmi_miss_p/100), total_p))
}

if (g_age %in% c("75plus")  )  {z_df <- df_table %>%  
  filter(!(country == "ni" & label == "Age Group" & levels %in% c("05_11","12_15","16_17", "18_24", "25_29","30_34","35_39",
                                                                  "40_44","45_49","50_54","55_59","60_64","65_69","70_74")))
z_bmi_miss_p <- z_df %>% filter(country=="scotland" & label=="BMI" & levels=="Missing") %>% pull(total_p)
z_df <- z_df %>% mutate(total_p = if_else(country=="scotland" & label=="BMI" & levels!="Missing", total_p/(1-z_bmi_miss_p/100), total_p))
z_bmi_miss_p <- z_df %>% filter(country=="england" & label=="BMI" & levels=="Missing") %>% pull(total_p)
z_df <- z_df %>% mutate(total_p = if_else(country=="england" & label=="BMI" & levels!="Missing", total_p/(1-z_bmi_miss_p/100), total_p))
}


  
  
)