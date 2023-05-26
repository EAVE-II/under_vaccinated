#reading in the coefficients and variances for the under vaccinated logistic regression
library(tidyverse)
library(readr)
library(metafor)

file_locations <- read.csv("data/Filenames.csv")
vars_labs_lkup <- read.csv("data/uv_vars_labels_lookup.csv")
#ni_ests <- read.csv("data/ni/d_undervac_coef.csv")

z_title <- "Undervaccination Age 5-11"  #5-11, 12-15, 16-74, 75+
z_model <- "min_adj"

z_files <- file_locations %>% filter(title==z_title & model==z_model)

coefs<-list()
vcovs<-list()

#read in the data sets and put in two lists
for (i in 1:nrow(z_files)) {
  #i <- 1
  z_type <- z_files[i,"type"]
  z_country <- z_files[i,"country"]
  z_fname <- z_files[i,"filename"]
  z_f <- if (z_country=="scotland") read.csv(z_fname) else read.csv(paste0("data/",z_country,"/",z_fname))
  print(c(z_country,  z_type, dim(z_f)))
  print(colnames(z_f))
  if (z_type == "coef") coefs[[z_country]] <- z_f
  if (z_type == "vcov") vcovs[[z_country]] <- z_f
}

#coefs + se's
if (exists("df_coefs")) rm(df_coefs)
for (i in names(coefs)){
  #i <- names(coefs)[[4]]
  z_df <- as.data.frame(coefs[[i]]) %>% mutate(country=i)
  if(i=="wales") z_df <- z_df %>% 
      mutate(index =paste0(xvar,xlbl)) %>% 
      left_join(vars_labs_lkup, by=c("xvar","xlbl")) %>% 
      dplyr::select(country,index,xvar,xlbl,estimate) %>%
      dplyr::rename(coef=estimate) 
  if(i =="scotland") z_df <- z_df %>% 
      left_join(vars_labs_lkup, by=c("variable"="xvar_sc")) %>%
      dplyr::rename(index=variable) %>% 
      dplyr::select(country,index,xvar,xlbl,coef) %>% 
      filter(!grepl("health_board", index))
  if(i =="england") z_df <- z_df %>% 
      left_join(vars_labs_lkup, by=c("variable"="xvar_en")) %>%
      dplyr::rename(index=variable) %>% 
      dplyr::select(country,index,xvar,xlbl,coef)
  if(i=="ni") z_df <- z_df %>% 
      mutate(index =paste0(xvar,xlbl)) %>% dplyr::select(-xvar, -xlbl) %>% 
      left_join(vars_labs_lkup, by=c("index"="xvar_ni")) %>% 
      dplyr::select(country,index,xvar,xlbl,estimate) %>%
      dplyr::rename(coef=estimate)
  if (exists("df_coefs")) df_coefs <- bind_rows(df_coefs,z_df) else df_coefs <- z_df
}
df_coefs <- df_coefs %>%
  mutate(index =if_else(country=="wales" & grepl("Intercept",xvar), xvar, index)) %>% 
  mutate(index =if_else(country=="ni" & grepl("Intercept",xvar), xvar, index))

if (exists("df_ses")) rm(df_ses)
for (i in names(vcovs)){
  #i <- names(vcovs)[[2]]
  z_df <- as.data.frame(vcovs[[i]])
    if(i %in% c("england", "scotland")) {z_names <- z_df$X
       z_df <- z_df %>% dplyr::select(-X)
       names(z_df) <- z_names}
  z_se <- z_df %>% as.matrix() %>% diag()
  z_se <- data.frame(index=names(z_df), variance=z_se) %>% 
    mutate(country=i)
  if (i=="scotland") z_se <- z_se %>% filter(!grepl("health_board", index))
  if (exists("df_ses")) df_ses <- bind_rows(df_ses,z_se) else df_ses <- z_se
}
df_ses <- df_ses %>%
  mutate(index =if_else(country=="wales" & grepl("Intercept",index), "(Intercept)", index)) %>% 
  mutate(index =if_else(country=="ni" & grepl("Intercept",index), "(Intercept)", index)) %>% 
  mutate(index =if_else(country=="ni" & grepl("\\.$",index), gsub("\\.$","+",index), index)) %>% 
  mutate(index =if_else(country=="ni" & grepl("\\.\\.\\.",index), gsub("\\.\\.\\."," - ",index), index)) %>% 
  mutate(index =if_else(country=="ni" & grepl("age_5y_cat",index), gsub("\\.","-",index), index)) %>% 
  mutate(index =if_else(country=="ni" & grepl("mdm_quintile",index), gsub("\\."," ",index), index))

df_coefs <- df_coefs %>% left_join(df_ses, by=c("country","index"))
df_coefs <- df_coefs %>% mutate(title=z_title, model=z_model)


#changes to labels - should be temporary
if (z_title == "Undervaccination Age 5-11" | z_title == "Undervaccination Age 12-15") {
  z <- df_coefs %>%  
    mutate(xvar=if_else(xvar=="qc_count_cat","qc_count_cat3", xvar)) %>% 
    mutate(xlbl=if_else(xvar=="qc_count_cat3" & xlbl == "01", "1", xlbl)) %>% 
    mutate(xlbl=if_else(xvar=="qc_count_cat3" & xlbl == "02", "2plus", xlbl)) %>% 
    mutate(xlbl=if_else(xvar=="qc_count_cat3" & xlbl == "2", "2plus", xlbl)) 
  df_coefs <- z
  }


z <- df_coefs %>% 
  mutate(xlbl = if_else(country=="wales" & xvar=="qc_count_cat" & xlbl %in% paste0(0,1:4), gsub("0","",xlbl), xlbl)) %>% 
  mutate(xlbl = if_else(country=="wales" & xvar=="qc_count_cat" & xlbl == "05", "5plus", xlbl)) 
#if (z_title == "Undervaccination Age 16-74" & z_model == "min_adj") { #wales has one age group less
#  z1 <- z %>% filter(country=="wales" & xvar=="age_5y_cat" & xlbl == "25_29") %>%
#    mutate(xlbl="18_24")
#  z <- rbind.data.frame(z,z1)
#}
z <- z %>%  arrange(country, xvar, xlbl)
df_coefs <- z

remove(list=ls(pa="^z"))
rm(i,df_ses)  
