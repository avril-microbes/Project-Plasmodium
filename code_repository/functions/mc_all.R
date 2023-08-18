# code that runs the semi-stochastic model with all six parameters varying at the same time
# The function takes a list of coefficients that determines the parasite's strategy (par), cue choice (cue),
# cue range, whether the cue is logged or not, and rand_df (dataframe containing parameter values drawn from the
# posterior distribution). 

mc_all <- function(par, cue, cue_range, log, rand_df){
  
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
  
  # other parameters
  parameters_rand <- c(R1 = 8.89*10^6, 
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
                           psin = psin_rand,
                           psiw = psiw_rand,
                           phin = phin_rand, 
                           phiw = phiw_rand,
                           iota = 2.18*(10^6),
                           rho = rho_rand)
  
  # run for 20 days with the larger time step of 0.01. note that this does not really change the overall fitness
  # whenever comparisons are made between the randomized run and non-randomized run, it must be compared to the
  # reference run that is only simualted with 0.01 step size (id = 1)
  time_range <- seq(0, 20, by = 1e-2)
  
  # run dynamics
  dyn <- chabaudi_si_clean(
    parameters_cr = par, 
    parameters = parameters_rand, 
    immunity = "tsukushi",
    time_range = time_range, 
    cue = cue, 
    cue_range = cue_range, 
    log_cue = log,
    solver = "vode",
    dyn = TRUE)
  
  # get fitness
  fitness <- dyn %>% dplyr::filter(variable == "tau_cum") %>% 
    dplyr::summarise(max_fitness = max(value,  na.rm = T))
  
  # attach fitness to associated meta data and write and return it
  fitness.df <- cbind(fitness  = fitness,
                      cue = cue,
                      log = log,
                      rand_df)
  
  # keep only conversion rate and transmission potential and write it into MC_all_dyn
  dyn_f <- dyn %>% 
    dplyr::filter(variable %in% c("cr", "tau")) %>% 
    dplyr::group_by(variable) %>% 
    dplyr::arrange(time) %>% 
   # filter(row_number() %% 10 == 1) %>% # keep only the 0.01th 
    dplyr::mutate(cue = cue, 
           log = log) # add cue and log

  # write files
  write.csv(fitness.df, here(paste0("code_repository/data/mc_all_fitness/", cue, "_", log, "_", id, ".csv")))
  arrow::write_parquet(x = dyn_f, sink = here(paste0("code_repository/data/mc_all_dyn/", cue, "_", log, "_", id, ".parquet")))
}