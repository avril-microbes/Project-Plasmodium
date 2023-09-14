# function to perform differential evolution optimziation (global optimization) of the dual-cue model

dual_cue_opt_global <- function(df, itermax = 500, steptol = 50){
  # process cues
  cue <- df$cue
  cue_b <- df$cue_b
  
  # process log
  log <- ifelse(stringr::str_detect("log", df$log), "log10", "none")
  log_b <- ifelse(stringr::str_detect("log", df$log_b), "log10", "none")
  
  # process cue_range. cannot use by to ensure that both cue ranges are of the same length
  cue_range <- seq(df$low, df$high, length.out = 500)
  cue_range_b <- seq(df$low_b, df$high_b, length.out = 500)
  

  # DE optimization
  global <- DEoptim::DEoptim(
    fn = chabaudi_si_clean, 
    control = list(trace = 1, parallelType = "parallel", itermax = itermax,
                    NP = 90, F = 0.8, CR = 0.9,
                   packages = c("crone", "splines2", "mgcv", "deSolve", "tidyr", "stringr", "dplyr")), # NP, F, CR values set according to rec by storn et al 2006
    lower = c(-10, -500, -1000, -1000, -250, -500, -1000, -500, -250), # lower and upper bound values are derived empirically based on spline space
    upper = c(10, 500, 1000, 1000, 250, 500, 1000, 500, 250),
    immunity = "tsukushi",
    parameters = parameters_tsukushi,
    time_range = seq(0, 20, 0.05), # this is the largest time step we can do without convergence issues
    cue = cue,
    cue_b = cue_b,
    cue_range = cue_range,
    cue_range_b = cue_range_b,
    log_cue = log,
    log_cue_b = log_b,
    solver = "vode",
    gam = "te",
    neg = T)
  
  # local optimization
  cl <- makeCluster(8); setDefaultCluster(cl = cl)
  res <- optimParallel(
    par = c(global$optim$bestmem), # starting value is the optimal from global 
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
  
  write.csv(output, here(paste0("code_repository/data/dual_cue_opt_global/", df$Var1, "_", df$Var2, ".csv")))
  return(output)
  
}