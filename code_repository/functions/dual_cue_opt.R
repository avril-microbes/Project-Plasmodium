dual_cue_opt <- function(df){
  # process cues
  cue <- df$cue
  cue_b <- df$cue_b
  
  # process log
  log <- ifelse(str_detect("log", df$log), "log10", "none")
  log_b <- ifelse(str_detect("log", df$log_b), "log10", "none")
  
  # process cue_range. cannot use by to ensure that both cue ranges are of the same length
  cue_range <- seq(df$low, df$high, length.out = 500)
  cue_range_b <- seq(df$low_b, df$high_b, length.out = 500)
  
  # optimization
  cl <- makeCluster(detectCores()); setDefaultCluster(cl = cl)
  
  # start optimization LFBGS with inital value of 0.5x5
  res <- optimParallel(
    par = rep(0.5, 9),
    fn = chabaudi_si_clean, 
    control = list(trace = 6, fnscale = -1),
    immunity = "tsukushi",
    parameters = parameters_tsukushi,
    time_range = seq(0, 20, 0.01),
    cue = cue,
    cue_b = cue_b,
    cue_range = cue_range,
    cue_range_b = cue_range_b,
    log_cue = log,
    log_cue_b = log_b,
    solver = "vode",
    gam = "te")
  
  
  
  # close cluster
  stopCluster(cl)
  
  # get model output
  fitness <- res$value ## fitness difference between mutant and residence
  par <- res$par ## optimized parameters of mutant
  
  # produce output
  output <- cbind.data.frame(id = df$Var1, id_b = df$Var2, 
                             label = df$label, label_b = df$label_b,
                             fitness = fitness,
                             par1 = par[1], par2 = par[2], par3 = par[3], par4 = par[4], par5 = par[5], par6 = par[6], par7 = par[7], par8 = par[8], par9 = par[9])
  
  write.csv(output, here(paste0("code_repository/data/dual_cue_opt/", df$Var1, "_", df$Var2, ".csv")))
  return(output)
  
}