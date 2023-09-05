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
  
  # id of the run
  id <- rand_df$id
  
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
  fitness_rho <- chabaudi_si_clean(
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
    gam = "te")
  
  fitness_beta <- chabaudi_si_clean(
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
    gam = "te")
  
  fitness_psin <- chabaudi_si_clean(
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
    gam = "te")
  
  fitness_psiw <- chabaudi_si_clean(
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
    gam = "te")
  
  fitness_phin <- chabaudi_si_clean(
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
    gam = "te")
  
  fitness_phiw <- chabaudi_si_clean(
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
    gam = "te")
  
  # arrange results
  res <- cbind.data.frame(
    cue = cue,
    cue_b = cue_b,
    log = log,
    log_b = log_b,
    rand_df,
    fitness_rho = fitness_rho,
    fitness_beta = fitness_beta,
    fitness_psin = fitness_psin,
    fitness_psiw = fitness_psiw,
    fitness_phin = fitness_phin,
    fitness_phiw = fitness_phiw
  )
  
  # write results
  write.csv(res, here(paste0("code_repository/data/mc_single_fitness/", cue, "_", log, "-", cue_b, "_", log_b, "_", id, "_single.csv")))
}