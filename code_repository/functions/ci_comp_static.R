ci_comp_static <- function(df){
  
  # process cue. if detect _i, suggesting individual cue. 
  if(stringr::str_detect(df$cue, "-i")){cue_1 = paste0(gsub("*-i", "", df$cue), "1")}
  if(stringr::str_detect(df$cue, "-i", negate = T)){cue_1 = df$cue}
  
  if(stringr::str_detect(df$cue_2, "-i")){cue_2 = paste0(gsub("*-i", "", df$cue_2), "2")}
  if(stringr::str_detect(df$cue_2, "-i", negate = T)){cue_2 = df$cue_2}
  
  if(df$cue == "I-i+Ig-i"){cue_1 = "I1+Ig1"}
  if(df$cue_2 == "I-i+Ig-i"){cue_2 = "I2+Ig2"}
  
  # process log10
  if(stringr::str_detect(df$log, "log")){log_1 = "log10"}
  if(stringr::str_detect(df$log, "none")){log_1 = "none"}
  
  if(stringr::str_detect(df$log_2, "log")){log_2 = "log10"}
  if(stringr::str_detect(df$log_2, "none")){log_2 = "none"}
  
  # get fitness difference between strain and strain 2 
  # simulated fort 30 days for transmission resolution (infectino may not resolve but transmission accumulation has halted)
  # past test suggests very little difference between 20 vs 30 days.
  dyn <- chabaudi_ci_clean(parameters_cr_1 = c(df$var1, df$var2, df$var3, df$var4),
                           parameters_cr_2 = c(df$var1_2, df$var2_2, df$var3_2, df$var4_2),
                           immunity = "tsukushi",
                           parameters = parameters_tsukushi,
                           cue_1 = cue_1,
                           cue_2 = cue_2,
                           cue_range_1 = seq(df$low, df$high, by = df$by),
                           cue_range_2 = seq(df$low_2, df$high_2, by = df$by_2),
                           log_cue_1 = log_1,
                           log_cue_2 = log_2,
                           solver = "vode",
                           time_range = seq(0, 30, 0.001),
                           dyn = T)
  
  dyn2 <- cbind(id_1 = rep(df$id, nrow(dyn)), id_2 = rep(df$id_2, nrow(dyn)),
                dyn)
  filename <- paste0(df$id, "_", df$id_2, "_static.parquet")
  write_parquet(dyn2, here(paste0("data/ci_static/", filename)))
  
  # get fitness
  fitness_1 <- max(dyn2 %>% filter(variable == "tau_cum1") %>% select(value), na.rm = T)
  fitness_2 <- max(dyn2 %>% filter(variable == "tau_cum2") %>% select(value), na.rm = T)
  
  # return list
  return(list(df$id, df$id_2, fitness_1, fitness_2))
}