
if (g_age %in% c("5_11","12_15")  )  {z_df <- df_table %>% filter(label !="BMI") %>% filter(label!= "Age Group")  
  z_df <- z_df %>% filter(!(label=="Number of people in household" & levels=="1"))
} 
if (g_age %in% c("16_74", "75plus")  )  {z_df <- df_table 
  z_bmi_miss_p <- z_df %>% filter(country=="scotland" & label=="BMI" & levels=="Missing") %>% pull(total_p)
  z_df <- z_df %>% mutate(total_p = if_else(country=="scotland" & label=="BMI" & levels!="Missing", total_p/(1-z_bmi_miss_p/100), total_p))
  z_bmi_miss_p <- z_df %>% filter(country=="england" & label=="BMI" & levels=="Missing") %>% pull(total_p)
  z_df <- z_df %>% mutate(total_p = if_else(country=="england" & label=="BMI" & levels!="Missing", total_p/(1-z_bmi_miss_p/100), total_p))
}

#take the unknown Ethnicity person years and spread throughout the other categories
z <- z_df %>% filter(label== "Ethnicity" & levels=="Unknown") %>% dplyr::select(country, label, total_p)
z_df <- z_df %>% left_join(z, by=c("country", "label"), suffix=c("","_unk")) %>% 
  mutate(total_p = case_when(label=="Ethnicity" & levels!="Unknown" ~ 100*total_p/(100-total_p_unk),
                          TRUE ~ total_p)) %>% 
  filter(!(label=="Ethnicity" & levels=="Unknown")) %>% 
  dplyr::select(-total_p_unk)


z_df <- z_df %>% filter(levels!="Missing")

z_df %>% filter(label != "Total") %>% ggplot(aes(x=levels, y=total_p , fill=country)) + geom_col(position="dodge") +
  facet_wrap(~label, scales="free_x") +
  theme(axis.text.x = element_text(angle = 45, vjust=0.5)) + 
  labs(x="", y="Percentage")
ggsave(paste0("output/",g_model,"_",g_age,".png"), height=20, width=25, units="cm")


z_df %>%  ggplot(aes(x=levels, y=p_undervac , fill=country)) + geom_col(position="dodge") +
  facet_wrap(~label, scales="free_x") +
  theme(axis.text.x = element_text(angle = 45, vjust=0.5)) + 
  labs(x="", y="Percentage Undervaccinated")
ggsave(paste0("output/perc_under_vacc_",g_model,"_",g_age,".png"), height=20, width=25, units="cm")


#overall figures

z <- z_df %>% filter(!(country=="ni" & label=="Number of risk groups")) %>% 
  group_by(label,levels) %>% dplyr::summarise(across(all_of(c("total_n", "undervac_n")) , ~sum(.))) %>% 
  mutate(p_undervac = undervac_n/total_n*100) %>% 
  group_by(label) %>% mutate(total_p = total_n/sum(total_n)*100) %>% 
  ungroup()
z_ci <- Hmisc::binconf(z$undervac_n, n=z$total_n, method="wilson")
z <- z %>% mutate(ucl = z_ci[,3]*100,
                        lcl = z_ci[,2]*100)

z %>% ggplot(aes(x=levels, y=p_undervac)) + geom_point(aes(size=total_p)) + geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.2) + 
  facet_wrap(~label, scales="free") +
  theme(axis.text.x = element_text(angle = 45, vjust=0.5))  +
  labs(x="", y="Percentage", title="Under Vaccination", size="Percentage")
ggsave(paste0("output/perc_under_vacc_point_",g_model,"_",g_age,".png"), height=20, width=25, units="cm")

write_csv(z, paste0("output/perc_under_vacc_numbers_",g_model,"_",g_age,".csv"))
