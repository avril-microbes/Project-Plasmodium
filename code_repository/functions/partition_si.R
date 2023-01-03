partition_si <- function(df){
  # source all functions
  source(here("functions/chabaudi_si_clean.R"))
  source(here("functions/chabaudi_si_clean_R.R"))
  source(here("functions/chabaudi_si_clean_N.R"))
  source(here("functions/chabaudi_si_clean_W.R"))
  
  # get parameters
  parameters_tsukushi <- c(R1 = 8.89*(10^6), # slightly higher
                           lambda = 3.7*(10^5),
                           mu = 0.025, 
                           p = 8*(10^-6), # doubled form original
                           alpha = 1, 
                           alphag = 2, 
                           beta = 5.721, 
                           mum = 48, 
                           mug = 4, 
                           I0 = 43.85965, 
                           Ig0 = 0, 
                           a = 150, 
                           b = 100, 
                           sp = 1,
                           psin = 16.69234,
                           psiw = 0.8431785,
                           phin = 0.03520591, 
                           phiw = 550.842,
                           iota = 2.18*(10^6),
                           rho = 0.2627156)
  
  # get all characteristics
  cue <- df$cue
  log <- ifelse(df$log == "log", "log10", "none")
  cue_range <- seq(df$low, df$high, by = df$by)
  id <- df$id
  
  # run optimization with static RBC
  cl <- makeCluster(detectCores()); setDefaultCluster(cl = cl)
  static_R <- optimParallel(
    par = rep(0.5,4), # start at 0.5x4
    fn = chabaudi_si_clean_R, 
    control = list(trace = 6, fnscale = -1),
    immunity = "tsukushi",
    parameters = parameters_tsukushi,
    time_range = seq(0, 20, by = 1e-3),
    cue_range =  cue_range,
    cue = cue,
    log_cue = log,
    solver = "vode")
  
  
  # run optimization with no indiscriminant immunity
  no_N <- optimParallel(
    par = rep(0.5,4), # start at 0.5x4
    fn = chabaudi_si_clean_N, 
    control = list(trace = 6, fnscale = -1),
    immunity = "tsukushi",
    parameters = parameters_tsukushi,
    time_range = seq(0, 20, by = 1e-3),
    cue_range =  cue_range,
    cue = cue,
    log_cue = log,
    solver = "vode")
  
  # run simulation with no targeted immunity (W)
  no_W <- optimParallel(
    par = rep(0.5,4), # start at 0.5x4
    fn = chabaudi_si_clean_W, 
    control = list(trace = 6, fnscale = -1),
    immunity = "tsukushi",
    parameters = parameters_tsukushi,
    time_range = seq(0, 20, by = 1e-3),
    cue_range =  cue_range,
    cue = cue,
    log_cue = log,
    solver = "vode")
  
  stopCluster(cl)
  
  # join together all data
  res <- cbind(id = id,
               cue = cue,
               log = log,
               fitness_R = static_R$value,
               fitness_N = no_N$value,
               fitness_W = no_W$value,
               var_R1 = static_R$par[[1]],
               var_R2 = static_R$par[[2]], 
               var_R3 = static_R$par[[3]],
               var_R4 = static_R$par[[4]],
               var_N1 = no_N$par[[1]],
               var_N2 = no_N$par[[2]], 
               var_N3 = no_N$par[[3]],
               var_N4 = no_N$par[[4]],
               var_W1 = no_W$par[[1]],
               var_W2 = no_W$par[[2]], 
               var_W3 = no_W$par[[3]],
               var_W4 = no_W$par[[4]])
  
  write.csv(res, paste0(here("data/partition/"), id, "_partition.csv"))
}