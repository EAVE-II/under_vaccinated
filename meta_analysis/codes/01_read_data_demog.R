#read in basic demographic tables and amalgamate
library(tidyverse)
library(readr)
library(readxl)
library(metafor)

file_locations <- read.csv("data/Filenames.csv")
vars_labs_lkup <- read_xlsx("data/demographics_label_lookup.xlsx")
#ni_data <- read.csv("data/ni/t_undervac_overall_np.csv")

z_title <- "Undervaccination Age 5-11"  #5-11, 12-15, 16-74, 75+
z_model <- "demographics"
g_age <- gsub("Undervaccination Age ","",z_title)
g_age <- gsub("-","_",g_age)
g_age <- gsub("\\+","plus",g_age)
g_model <- substring(z_model,1,5)
z_files <- file_locations %>% filter(title==z_title & model==z_model & filename != "")
#if (g_age %in% c("5_11", "12_15")) ni_data <- filter(ni_data, xvar != "age_5y_cat")

tables<-list()

#read in the data sets and put in two lists
for (i in 1:nrow(z_files)) {
  #i <- 1
  z_type <- z_files[i,"type"]
  z_country <- z_files[i,"country"]
  z_fname <- z_files[i,"filename"]
  z_f <- if (z_country=="scotland") read.csv(z_fname) else read.csv(paste0("data/",z_country,"/",z_fname))
  print(c(z_country,  z_type, dim(z_f)))
  print(colnames(z_f))
  tables[[z_country]] <- z_f
}
#tables[["ni"]] <- ni_data


#combine tables and relabel
if (exists("df_table")) rm(df_table)
for (i in names(tables)){
  #i <- names(tables)[[3]]
  z_df <- as.data.frame(tables[[i]]) %>% mutate(country=i)
  if(i=="wales") {z_df <- z_df %>% 
    mutate(full_v = total_n-undervac_n) %>% 
    mutate(p_undervac= undervac_n/total_n*100) %>% 
    dplyr::select(country, xvar, xlbl, full_v, undervac_n, total_n, p_undervac, total_p) %>% 
    dplyr::rename(label=xvar, levels=xlbl) %>%
    filter(!(label %in% c("pcr_last_positive_cat", "health_board"))) 
    }
  if(i =="scotland") {
    z_variables <- c("Total","Age group", "Sex","Ethnicity","SIMD quintile","Urban/rural classification",
    "BMI","Number of people in household", "Number of risk groups")
    z <- z_df
    while(any(z$label=="")) z <- z %>% mutate(label = if_else(label=="", lag(label), label))
    z <- z %>% filter(label %in% z_variables)
    z_position_1 <- regexpr("\\(", z$Fully.vaccinated)-2
    z_position_2 <- regexpr("\\(", z$Sub.optimally.vaccinated)-2
    z <- z %>% mutate(Fully.vaccinated = as.numeric(substring(Fully.vaccinated, 1,z_position_1)),
                      Sub.optimally.vaccinated = as.numeric(substring(Sub.optimally.vaccinated, 1,z_position_2)) )
    z <- z %>% mutate(total_n = Fully.vaccinated + Sub.optimally.vaccinated) %>% 
        dplyr::rename(undervac_n = Sub.optimally.vaccinated) %>% 
        mutate(p_undervac=undervac_n/total_n*100) %>% 
        group_by(label) %>% mutate(total_p = total_n/sum(total_n)*100) %>% ungroup()
    z_df <- z %>% dplyr::select(-Fully.vaccinated)
   }
   if(i =="england") {z_df <- z_df %>% 
     mutate(full_v = total_n-undervac_n) %>% 
     mutate(p_undervac= undervac_n/total_n*100) %>% 
     dplyr::select(country, xvar, xlbl, full_v, undervac_n, total_n, p_undervac, total_p) %>% 
     dplyr::rename(label=xvar, levels=xlbl) %>% 
     filter(!(label %in% c("num_doses_start", "under_vaccinated","age_4cat", ""))) %>% 
     mutate(levels = if_else(label=="bmi_cat" & is.na(levels), "missing", levels))
     if (g_age %in% c("5_11","12_15")) z_df <- filter(z_df, label != "qc_count_cat")
     if (g_age %in% c("16_74","75plus")) z_df <-  filter(z_df, label != "qc_count_3cat")

   }
  if(i =="ni") {z_df <- z_df %>% 
    mutate(full_v = total_n-undervac_n) %>% 
    mutate(p_undervac= undervac_n/total_n*100) %>% 
    dplyr::select(country, xvar, xlbl, full_v, undervac_n, total_n, p_undervac, total_p) %>% 
    dplyr::rename(label=xvar, levels=xlbl) %>% 
    filter(!(label %in% c("pcr_last_positive_cat" ,"lgd","hh_n_cat")))
  }
  
  if (exists("df_table")) df_table <- bind_rows(df_table,z_df) else df_table <- z_df
}


#changes to labels and levels
#labels
z_labs <- vars_labs_lkup %>% filter(!duplicated(paste(country,orig_label))) %>%
  dplyr::select(country, orig_label,final_label)
z_df <- df_table %>% left_join(z_labs, by=c("country","label"="orig_label"))
#levels
z_labs <- vars_labs_lkup %>% 
  filter(!(country=="scotland" & orig_label =="n_risk_gps_6cat")) %>% 
  dplyr::select(-orig_label)
z_df <- z_df %>% left_join(z_labs, by=c("country","final_label", "levels"="orig_level"))
z_df <- z_df %>% dplyr::select(-label, -levels) %>% 
  dplyr::rename(label=final_label, levels=final_level) %>% 
  dplyr::relocate(label:levels, .after=country)
z_df <- z_df %>%  arrange(country, label, levels)
df_table <- z_df

if (g_age %in% c("5_11","12_15")) df_table <- df_table %>% mutate(levels = if_else(label=="Number of risk groups" & levels=="2","2+",levels))

remove(list=ls(pa="^z"))
rm(i)  
