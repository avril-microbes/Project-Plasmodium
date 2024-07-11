# This function executes single inefction simulation with only 1 parameter varying. Note for dual cue models!!!
mc_single_dual <- function(par, cue, cue_b, cue_range, cue_range_b, log, log_b, rand_df){
  
  # source single infection function
  source(here("code_repository/functions/chabaudi_si_clean.R"))
  
  # get parameter values form rand_df
  rho_rand <- rand_df$rho # proportion of RBC recovered
  beta_rand <- rand_df$burst # burst size
  psin_rand <- rand_df$iota_N1 # activation strengths for indiscriminate immunity
  psiw_rand <- rand_df$iota_N2 # activation strengths for targeted immunity
  phin_rand <- rand_df$phi_N1 # half-life for indiscriminate immunity
  phiw_rand <- rand_df$phi_N2 # half-life for targeted immunity
  
  # id of the run (cue)
  id <- gsub("log10", "log", paste0(cue, "_", log, "-", cue_b, "_", log_b))
  
  # iter of the run
  iter <- rand_df$iter
  
  # assign parameter lists where for each one, only one parameter varies
  parameters_rho <- c(R1 = 8.89*10^6, 
                      lambda = 3.7*10^5,
                      mu = 0.025, 
                      p = 8*10^-6, 
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
                      rho = rho_rand)
  
  parameters_beta <- c(R1 = 8.89*10^6, 
                       lambda = 3.7*10^5,
                       mu = 0.025, 
                       p = 8*10^-6, 
                       alpha = 1, 
                       alphag = 2, 
                       beta = beta_rand, 
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
  
  parameters_psin <- c(R1 = 8.89*10^6, 
                       lambda = 3.7*10^5,
                       mu = 0.025, 
                       p = 8*10^-6, 
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
                       psin = psin_rand,
                       psiw = 0.8431785,
                       phin = 0.03520591, 
                       phiw = 550.842,
                       iota = 2.18*(10^6),
                       rho = 0.2627156)
  
  parameters_psiw <- c(R1 = 8.89*10^6, 
                       lambda = 3.7*10^5,
                       mu = 0.025, 
                       p = 8*10^-6, 
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
                       psiw = psiw_rand,
                       phin = 0.03520591, 
                       phiw = 550.842,
                       iota = 2.18*(10^6),
                       rho = 0.2627156)
  
  parameters_phin <- c(R1 = 8.89*10^6, 
                       lambda = 3.7*10^5,
                       mu = 0.025, 
                       p = 8*10^-6, 
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
                       phin = phin_rand, 
                       phiw = 550.842,
                       iota = 2.18*(10^6),
                       rho = 0.2627156)
  
  parameters_phiw <- c(R1 = 8.89*10^6, 
                       lambda = 3.7*10^5,
                       mu = 0.025, 
                       p = 8*10^-6, 
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
                       phiw = phiw_rand,
                       iota = 2.18*(10^6),
                       rho = 0.2627156)
  
  # time
  time_range <- seq(0, 20, by = 1e-2)
  
  # get fitness when only one parameter is randomized
  dyn_rho <- chabaudi_si_clean(
    parameters_cr = par, 
    parameters = parameters_rho, 
    immunity = "tsukushi",
    time_range = time_range, 
    cue = cue, 
    cue_b = cue_b,
    cue_range = cue_range, 
    cue_range_b = cue_range_b,
    log_cue = log,
    log_cue_b = log_b,
    solver = "vode",
    gam = "te",
    dyn = T)
  
  dyn_beta <- chabaudi_si_clean(
    parameters_cr = par, 
    parameters = parameters_beta, 
    immunity = "tsukushi",
    time_range = time_range, 
    cue = cue, 
    cue_b = cue_b,
    cue_range = cue_range, 
    cue_range_b = cue_range_b,
    log_cue = log,
    log_cue_b = log_b,
    solver = "vode",
    gam = "te",
    dyn = T)
  
  dyn_psin <- chabaudi_si_clean(
    parameters_cr = par, 
    parameters = parameters_psin, 
    immunity = "tsukushi",
    time_range = time_range, 
    cue = cue, 
    cue_b = cue_b,
    cue_range = cue_range, 
    cue_range_b = cue_range_b,
    log_cue = log,
    log_cue_b = log_b,
    solver = "vode",
    gam = "te",
    dyn = T)
  
  dyn_psiw <- chabaudi_si_clean(
    parameters_cr = par, 
    parameters = parameters_psiw, 
    immunity = "tsukushi",
    time_range = time_range, 
    cue = cue, 
    cue_b = cue_b,
    cue_range = cue_range, 
    cue_range_b = cue_range_b,
    log_cue = log,
    log_cue_b = log_b,
    solver = "vode",
    gam = "te",
    dyn = T)
  
  dyn_phin <- chabaudi_si_clean(
    parameters_cr = par, 
    parameters = parameters_phin, 
    immunity = "tsukushi",
    time_range = time_range, 
    cue = cue, 
    cue_b = cue_b,
    cue_range = cue_range, 
    cue_range_b = cue_range_b,
    log_cue = log,
    log_cue_b = log_b,
    solver = "vode",
    gam = "te",
    dyn = T)
  
  dyn_phiw <- chabaudi_si_clean(
    parameters_cr = par, 
    parameters = parameters_phiw, 
    immunity = "tsukushi",
    time_range = time_range, 
    cue = cue, 
    cue_b = cue_b,
    cue_range = cue_range, 
    cue_range_b = cue_range_b,
    log_cue = log,
    log_cue_b = log_b,
    solver = "vode",
    gam = "te",dyn = T)
  
  # get fitness
  fitness_rho <- dyn_rho %>% dplyr::filter(variable == "tau_cum") %>% 
    dplyr::summarise(max_fitness_rho = max(value,  na.rm = T))
  
  fitness_beta <- dyn_beta %>% dplyr::filter(variable == "tau_cum") %>% 
    dplyr::summarise(max_fitness_beta = max(value,  na.rm = T))
  
  fitness_psin <- dyn_psin %>% dplyr::filter(variable == "tau_cum") %>% 
    dplyr::summarise(max_fitness_psin = max(value,  na.rm = T))
  
  fitness_psiw <- dyn_psiw %>% dplyr::filter(variable == "tau_cum") %>% 
    dplyr::summarise(max_fitness_psiw = max(value,  na.rm = T))
  
  fitness_phin <- dyn_phin %>% dplyr::filter(variable == "tau_cum") %>% 
    dplyr::summarise(max_fitness_phin = max(value,  na.rm = T))
  
  fitness_phiw <- dyn_phiw %>% dplyr::filter(variable == "tau_cum") %>% 
    dplyr::summarise(max_fitness_phiw = max(value,  na.rm = T))
  
  # arrange results
  res <- cbind.data.frame(
    rand_df,
    fitness_rho = fitness_rho,
    fitness_beta = fitness_beta,
    fitness_psin = fitness_psin,
    fitness_psiw = fitness_psiw,
    fitness_phin = fitness_phin,
    fitness_phiw = fitness_phiw
  ) %>% mutate(id = id)
  
  # add associated id and iter to dynamics
  ## get list of dynamics
  dyn.ls <- list(dyn_rho, dyn_beta, dyn_psin, dyn_psiw, dyn_phin, dyn_phiw)
  ## attach id and iter
  dyn_p.ls <- mapply(function(x, y){
    df_p <- x %>% 
      mutate(id = id, iter = iter, parameter = y)
    ## write df
    arrow::write_parquet(x = df_p, sink = here(paste0("code_repository/data/mc_single_dyn/", id, "_", iter, "_", y, ".parquet")))
  },
  dyn.ls, c("rho", "beta", "psin", "psiw", "phin", "phiw"))
  
  # write fitness 
  write.csv(res, here(paste0("code_repository/data/mc_single_fitness2/", id, "_", iter, "_single.csv")))
}