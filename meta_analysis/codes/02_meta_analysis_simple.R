#meta analysis
#uses df_coefs from 01_read_data.R

z_vars <- df_coefs %>% pull(xvar) %>% unique()  # get unique variables names
z_vars <- z_vars[z_vars!="(Intercept)"]  #drop the intercept for the logisitc regression
print(z_vars)  #check the order as z_refs needs to be in the same order
z_title <- unique(df_coefs$title)
if (grepl("Hospital Death",z_title)) z_title_short <- "hosp_death"
if (grepl("Undervaccination",z_title)) z_title_short <- "under_vacc"
z_model <- unique(df_coefs$model)
if (grepl("Undervaccination", z_title)) z_y_label <- "Odds Ratio of Under Vaccination"
if (grepl("Hospital", z_title) | grepl("Death", z_title)) z_y_label <- "Hazard Ratio"

if (z_title=="Undervaccination Age 75+") z_refs <- c("75_79","White","0","female","urban","1_most")
if (z_title=="Undervaccination Age 16-74") z_refs <- c("18_24","White","0","female","urban","1_most")
if (z_title=="Undervaccination Age 12-15") z_refs <- c("White","0","female","urban","1_most")
if (z_title=="Undervaccination Age 5-11") z_refs <- c("White","0","female","urban","1_most")

if (z_title=="Hospital Death Age 75+") z_refs <- c("75_79","White","0","female","urban","full","1_most")
if (z_title=="Hospital Death Age 16-74") z_refs <- c("18_24","White","0","female","urban","full","1_most")
if (z_title=="Hospital Death Age 5-15") z_refs <- c("05_11","White","0","female","urban","full","1_most")

names(z_refs) <- z_vars
z_age_gp <- gsub("Undervaccination Age ", "", z_title)
z_age_gp <- gsub("Hospital Death Age ", "", z_title)
z_age_gp <- gsub("\\+","plus",z_age_gp)
z_age_gp <- gsub("-","_",z_age_gp)

if (exists("z_out")) rm(z_out)

for (i in 1:length(z_vars)) {
#i <-2
z_xvar <- z_vars[i]
z_ref <- z_refs[z_xvar]
z_df <- df_coefs %>% filter(xvar==z_xvar) %>%  
  mutate(xlbl=factor(xlbl, levels=c(z_ref,unique(xlbl)))) %>% 
  mutate(se=sqrt(variance)) %>% 
  mutate(lower=coef - 1.96*se, upper=coef + 1.96*se)
z_df_ma <- if (z_xvar %in% c("qc_count_cat","qc_count_3cat","qc_count_cat3")) filter(z_df, country != "ni") else z_df
z_mv <- metafor::rma.mv(yi=coef,V=variance, mods= ~ -1 + xlbl,  data=z_df_ma)
z_res <- data.frame(country="All", xvar = z_xvar, xlbl=dimnames(z_mv$beta)[[1]],
                    coef = z_mv$beta, se = z_mv$se, lower = z_mv$ci.lb, upper = z_mv$ci.ub) %>%
  mutate(xlbl = gsub("xlbl","", xlbl))
rownames(z_res) <- NULL
if (nrow(z_res)==1) z_res$xlbl <- unique(z_df$xlbl)

z<- z_df %>% dplyr::select(country, xvar, xlbl, coef, lower, upper, se)
z_res <- bind_rows(z,z_res)

z <- z_res %>% filter(!duplicated(country)) %>% mutate(xlbl=z_ref) %>% 
  mutate(across(coef:upper, ~ 0)) %>% 
  mutate(se=0)

z_res <- bind_rows(z_res,z) %>% 
  arrange(country, xlbl) %>% 
  mutate(across(coef:upper, ~exp(.)))

z_res <- z_res %>% mutate(age_gp=z_age_gp) %>% 
  dplyr::relocate(age_gp)

if (exists("z_out")) z_out <- bind_rows(z_out,z_res) else z_out <- z_res
}

#link in labels for plotting and tables
plot_labels_lookup <- read_xlsx("data/plot_label_lookup.xlsx")
z <- z_out %>% left_join(plot_labels_lookup, by=c("xvar","xlbl")) %>% 
  dplyr::select(-xvar,-xlbl) %>% 
  dplyr::rename(xvar=plot_xvar, xlbl=plot_xlbl) %>% 
  mutate(country=factor(country, levels = c("All", "england","ni","scotland","wales"), 
                        labels= c("All", "England","NI","Scotland","Wales") ) )
z <- z %>% filter(!(xvar=="Ethnicity" & xlbl=="Unknown"))
z_out <- z

write_csv(z_out, paste0("output/ma_",z_title_short,"_",z_model,"_",z_age_gp,".csv"))


z_out <- z_out %>% 
  mutate(coef =  if_else(!is.finite(coef) | coef >= 100, NA_real_,coef),
         lower = if_else(!is.finite(lower) | lower >= 100, NA_real_,lower),
         upper = if_else(!is.finite(upper) | upper >= 100, NA_real_,upper)) %>% 
  mutate(coef =  if_else(!is.finite(coef) | coef <= 0.01, NA_real_,coef),
         lower = if_else(!is.finite(lower) | lower <= 0.01, NA_real_,lower),
         upper = if_else(!is.finite(upper) | upper <= 0.01, NA_real_,upper))

z_out %>% filter(country == "All") %>% mutate(country="UK - pooled") %>% 
  ggplot(aes(
  x = xlbl,
  y = coef,
  ymin = lower,
  ymax = upper)) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  scale_y_log10() +
  coord_flip() +
#  facet_wrap(~country, ncol=4) +
  facet_grid(rows=vars(xvar), cols=vars(country), space = "free_y", scales = "free_y", switch = "y") +
  labs(title = "", x="",  y=z_y_label, colour="") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  theme_grey(base_size = 12) +
  theme(
    axis.title.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.background.y = element_blank(), strip.background.x = element_blank(),
    strip.text.y.left = element_text(angle= 0, face = "bold"), strip.text.x.top = element_text(angle= 0, face = "bold")
  )
ggsave(paste0("output/",z_title_short,"_",z_model,"_",z_age_gp,".png"), height=20, width=25, units="cm")

z_out %>% filter(country != "All") %>% #mutate(xlbl=paste(xvar,xlbl)) %>% 
  ggplot(aes(
    x = xlbl,
    y = coef,
    ymin = lower,
    ymax = upper)) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  scale_y_log10() +
  coord_flip() +
  #facet_wrap(~country, ncol=4) +
  facet_grid(rows=vars(xvar), cols=vars(country), space = "free_y", scales = "free_y", switch = "y") +
  labs(title = "", x="",  y=z_y_label, colour="") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  theme_grey(base_size = 12) +
  theme(
    axis.title.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.background.y = element_blank(), strip.background.x = element_blank(),
    strip.text.y.left = element_text(angle= 0, face = "bold"), strip.text.x.top = element_text(angle= 0, face = "bold")
  )
ggsave(paste0("output/",z_title_short,"_",z_model,"_countries_",z_age_gp,".png"), height=20, width=25, units="cm")
