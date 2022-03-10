#### optional script to repair names in metadata
#
# NB in the dea model, group names must be single string without spaces 


df_samples$group %>% 
  as.factor() %>% 
  fct_recode(
    group_1 = "Group 1",
    group_2 = "Group 2" 
  ) -> df_samples$group

contrast_name <- "group_2-group_1"
