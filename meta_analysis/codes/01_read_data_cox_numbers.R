#read in basic demographic tables and amalgamate
library(tidyverse)
library(readr)
library(readxl)
library(metafor)

file_locations <- read.csv("data/Filenames.csv")
vars_labs_lkup <- read_xlsx("data/demographics_label_lookup.xlsx")

z_title <- "Hospital Death Age 5-15" #5-15, 16-74, 75+
z_model <- "demographics"
g_age <- gsub("Hospital Death Age ","",z_title)
g_age <- gsub("-","_",g_age)
g_age <- gsub("\\+","plus",g_age)
g_model <- substring(z_model,1,5)
z_files <- file_locations %>% filter(title==z_title & model==z_model & filename != "")

tables<-list()

#read in the data sets and put in two lists
for (i in 1:nrow(z_files)) {
  #i <- 1
  z_type <- z_files[i,"type"]
  z_country <- z_files[i,"country"]
  z_fname <- z_files[i,"filename"]
  z_f <- if (z_country=="scotland") read.csv(z_fname) else { 
     if (z_country=="england_1") read.csv(paste0("data/",gsub("_1","",z_country),"/",z_fname)) else read.csv(paste0("data/",z_country,"/",z_fname))
    }
  print(c(z_country,  z_type, dim(z_f)))
  print(colnames(z_f))
  tables[[z_country]] <- z_f
}

#combine tables and relabel
if (exists("df_table")) rm(df_table)
for (i in names(tables)){
  #i <- names(tables)[[3]]
  z_df <- as.data.frame(tables[[i]]) %>% mutate(country=i)
  if(i=="wales") {z_df <- z_df %>% 
    dplyr::select(country, xvar, xlbl, total_n, pyears, events) %>% 
    dplyr::rename(Variable=xvar, Levels=xlbl) %>%
    filter(!(Variable %in% c("pcr_last_positive_cat", "health_board", "hosp_other_cat", "hh_n_cat", "bmi_cat"))) 
    }
  if(i =="scotland") {
    z_variables <- c("Total","Sex","Age group", "Ethnicity","SIMD quintile","Urban/rural classification",
    "Number of risk groups", "Vaccination status")
    z <- z_df
    while(any(z$Variable=="")) z <- z %>% mutate(Variable = if_else(Variable=="", lag(Variable), Variable))
    z <- z %>% filter(Variable %in% z_variables)
    z <- z %>% dplyr::select(-Hazard.ratio, -Event.rate.per.thousand.person.years)
    z <- z %>% dplyr::rename(total_n=Persons, events=Number.of.events, pyears=Person.years..thousands.)
    #get the total
    z_s <- z %>% filter(Variable=="Sex") %>% group_by(Variable, country) %>% 
      dplyr::summarise(across(pyears:events, ~sum(.))) %>% ungroup() %>% 
      mutate(Levels="Total", Variable="Total")
    z <- z %>% bind_rows(z_s) %>% mutate(pyears=pyears*1000)  
    z_df <- z %>% 
      mutate(Variable = if_else(Variable=="Urban/rural classification","Urban rural classification", Variable)
              )
   }
   if(i =="england") {z_df <- z_df %>%
     dplyr::rename(total_n = total_n_rounded5, events=events_rounded5, pyears=pyears_rounded1) %>% 
     dplyr::select(country, xvar, xlbl, total_n,  pyears, events) %>% 
     dplyr::rename(Variable=xvar, Levels=xlbl) %>%
      filter(!(Variable %in% c("region")))
   }
  if(i=="ni") {z_df <- z_df %>% 
    dplyr::select(country, xvar, xlbl, total_n, pyears, events) %>% 
    dplyr::rename(Variable=xvar, Levels=xlbl) %>%
    filter(!(Variable %in% c("pcr_last_positive_cat", "lgd", "hosp_other_cat", "hh_n_cat"))) 
  }
  
  if (exists("df_table")) df_table <- bind_rows(df_table,z_df) else df_table <- z_df
}


#changes to labels and levels
#labels
z_labs <- vars_labs_lkup %>% filter(!duplicated(paste(country,orig_label))) %>% 
  dplyr::select(country, orig_label,final_label)
z_df <- df_table %>% left_join(z_labs, by=c("country","Variable"="orig_label"))
#levels
z_labs <- vars_labs_lkup %>% filter(!(country=="scotland" & orig_label=="n_risk_gps_6cat")) %>%   #two risk group classifications for scotland
  filter(!duplicated(paste(country,final_label, orig_level))) %>%
  dplyr::select(-orig_label) 
z_df <- z_df %>% left_join(z_labs, by=c("country","final_label", "Levels"="orig_level"))
z_df <- z_df %>% dplyr::select(-Variable, -Levels) %>% 
  dplyr::rename(label=final_label, levels=final_level) %>% 
  dplyr::relocate(label:levels, .after=country)
z_df <- z_df %>%  arrange(country, label, levels)
df_table <- z_df

remove(list=ls(pa="^z"))
rm(i)  
