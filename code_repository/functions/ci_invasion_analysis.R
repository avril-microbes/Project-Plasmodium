ci_invasion_analysis <- function(df){
  
  # get cue
  if(stringr::str_detect(df$cue, "-i")){mut_cue = paste0(gsub("*-i", "", df$cue), "1")} ## mutant cue
  if(stringr::str_detect(df$cue, "-i", negate = T)){mut_cue = df$cue} ## mutant cue
  if(df$cue == "I-i+Ig-i"){mut_cue = "I1+Ig1"}
  
  if(stringr::str_detect(df$log, "log")){mut_log = "log10"} ## mutant log
  if(stringr::str_detect(df$log, "none")){mut_log = "none"} ## mutant log
  mut_cue_range <- seq(df$low, df$high, by = df$by) ## mutant cue range
  
  # get resident info (being invaded)
  res_par <- c(df$var1_2, df$var2_2, df$var3_2, df$var4_2) ## resident parameter
  if(stringr::str_detect(df$cue_2, "-i")){res_cue = paste0(gsub("*-i", "", df$cue_2), "2")} ## resident cue
  if(stringr::str_detect(df$cue_2, "-i", negate = T)){res_cue = df$cue_2} ## resident cue
  if(df$cue_2 == "I-i+Ig-i"){res_cue = "I2+Ig2"}
  if(stringr::str_detect(df$log_2, "log")){res_log = "log10"} ## resident log
  if(stringr::str_detect(df$log_2, "none")){res_log = "none"} ## resident log
  res_cue_range <- seq(df$low_2, df$high_2, by = df$by_2) ## resident cue range
  
  # get ids
  mut_id <- df$id
  res_id <- df$id_2
  
  # one way optimization
  ## start cluster
  cl <- makeCluster(detectCores()); setDefaultCluster(cl = cl)
  
  # start optimization LFBGS
  res_1 <- optimParallel(
    par = rep(0.5,4), # start at 0.5x4
    fn = chabaudi_ci_clean, 
    control = list(trace = 6, fnscale = -1),
    parameters_cr_2 = res_par,
    immunity = "tsukushi",
    parameters = parameters_tsukushi,
    time_range = seq(0, 20, by = 1e-3),
    cue_range_1 = mut_cue_range,
    cue_range_2 = res_cue_range,
    cue_1 = mut_cue,
    cue_2 = res_cue,
    log_cue_1 = mut_log,
    log_cue_2 = res_log,
    solver = "vode")
  
  res_2 <- optimParallel(
    par = c(df$var1, df$var2, df$var3, df$var4), # start at previous optimal parameter
    fn = chabaudi_ci_clean, 
    control = list(trace = 6, fnscale = -1),
    parameters_cr_2 = res_par,
    immunity = "tsukushi",
    parameters = parameters_tsukushi,
    time_range = seq(0, 20, by = 1e-3),
    cue_range_1 = mut_cue_range,
    cue_range_2 = res_cue_range,
    cue_1 = mut_cue,
    cue_2 = res_cue,
    log_cue_1 = mut_log,
    log_cue_2 = res_log,
    solver = "vode")
  
  # close cluster
  stopCluster(cl)
  
  # assign final result to which starting point that yields highest fitness
  if(res_1$value>res_2$value){final_res <- res_1}
  if(res_2$value>=res_1$value){final_res <- res_2}
  
  # get model output
  fitness <- final_res$value ## fitness difference between mutant and residence
  mut_opt_par <- final_res$par ## optimized parameters of mutant
  
  # produce output
  output <- cbind.data.frame(mut_id = mut_id, res_id = res_id, 
                             fitness = fitness, 
                             mut_var1 = df$var1, mut_var2 = df$var2, mut_var3 = df$var3, mut_var4 = df$var4,
                             res_var1 = df$var1_2, res_var2 = df$var2_2, res_var3 = df$var3_2, res_var4 = df$var4_2,
                             mut_var1_opt = mut_opt_par[1], mut_var2_opt = mut_opt_par[2], mut_var3_opt = mut_opt_par[3], mut_var4_opt = mut_opt_par[4])
  
  write.csv(output, here(paste0("data/ci_invasion/", mut_id, "_", res_id, ".csv")))
  return(output)
}