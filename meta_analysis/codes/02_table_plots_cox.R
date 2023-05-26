if (g_age %in% c("5_15")  )  {z_df <- df_table %>% filter(label !="BMI")  
z_df <- z_df %>% filter(!(label=="Number of people in household" & levels=="1"))
z_df <- z_df %>% mutate(levels = if_else(label=="Number of risk groups" & levels=="2", "2+", levels))
df_table <- z_df
}

#impute the suppressed events
#z <- df_table %>% filter(is.na(events))
#z_notw <- df_table %>% filter(label %in% unique(z$label) & levels %in% unique(z$levels) & !is.na(events))
z_notw <- df_table %>% filter(!is.na(events)) %>% group_by(label, levels) %>% 
  dplyr::summarise(across(pyears:events, ~sum(.))) %>% 
  ungroup() %>% 
  mutate(rate = events/pyears)
z  <- df_table %>%  left_join( dplyr::select(z_notw, label, levels, rate )) %>% 
  mutate(events_imp = as.integer(round(pyears*rate,0))) %>% 
  mutate(events = if_else(is.na(events), events_imp, events)) %>% 
  dplyr::select(-rate, -events_imp)

df_table_imp <- z %>% 
  arrange(country,label,levels)

df_table_imp <- df_table_imp %>% mutate(rate_1000 = events/pyears*1000) %>% 
  group_by(country,label) %>% mutate(total_p = total_n/sum(total_n)*100,
                                     py_p = pyears/sum(pyears)*100) %>% 
  ungroup()

#take the unknown Ethnicity person years and spread throughout the other categories
z <- df_table_imp %>% filter(label== "Ethnicity" & levels=="Unknown") %>% dplyr::select(country, label, py_p)
df_table_imp <- df_table_imp %>% left_join(z, by=c("country", "label"), suffix=c("","_unk")) %>% 
  mutate(py_p = case_when(label=="Ethnicity" & levels!="Unknown" ~ 100*py_p/(100-py_p_unk),
                          TRUE ~ py_p)) %>% 
  filter(!(label=="Ethnicity" & levels=="Unknown")) %>% 
  dplyr::select(-py_p_unk)

df_table_imp %>% filter(label != "Total") %>% ggplot(aes(x=levels, y=py_p , fill=country)) + geom_col(position="dodge") +
  facet_wrap(~label, scales="free_x") +
  theme(axis.text.x = element_text(angle = 45, vjust=0.5)) +
  labs(x="", y="Person Years - Percentage")
ggsave(paste0("output/py_perc_",g_model,"_",g_age,".png"), height=20, width=25, units="cm")


df_table_imp %>%  ggplot(aes(x=levels, y=rate_1000 , fill=country)) + geom_col(position="dodge") +
  facet_wrap(~label, scales="free") +
  theme(axis.text.x = element_text(angle = 45, vjust=0.5))  +
  labs(x="", y="Event rate per 1000 person years", title="Covid Hospitalisation or Death")
ggsave(paste0("output/hosp_death_",g_model,"_",g_age,".png"), height=20, width=25, units="cm")


z_df <- df_table_imp %>% filter(!(country=="ni" & label=="Number of risk groups")) %>% 
  group_by(label,levels) %>% dplyr::summarise(across(all_of(c("pyears", "events")) , ~sum(.))) %>% 
  mutate(rate_1000 = events/pyears*1000) %>% 
  group_by(label) %>% mutate(py_p = pyears/sum(pyears)*100) %>% 
  ungroup()
z <- survival::cipoisson(z_df$events, time=z_df$pyears/1000, method="anscombe")
z_df <- z_df %>% mutate(ucl = z[,2],
         lcl = z[,1])

z_df %>% ggplot(aes(x=levels, y=rate_1000)) + geom_point(aes(size=py_p)) + geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.2) + 
  facet_wrap(~label, scales="free") +
  theme(axis.text.x = element_text(angle = 45, vjust=0.5))  +
  labs(x="", y="Event rate per 1000 person years", title="Covid Hospitalisation or Death", size="Percentage")
ggsave(paste0("output/hosp_death_point_",g_model,"_",g_age,".png"), height=20, width=25, units="cm")

write_csv(z_df, paste0("output/hosp_death_numbers_",g_model,"_",g_age,".csv"))
