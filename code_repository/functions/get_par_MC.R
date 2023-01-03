get_par_MC <- function(par, cue, cue_range, log){
  # run MC with only rho
  res <- mc_five2(par = par,
                  cue = cue,
                  cue_range = cue_range,
                  log = log,
                  rand_var = "rho",
                  rand_mean = 0.2627156,
                  rand_sd = 0.2579136)
  
  fitness_name <- paste0("mc_par_", cue, "_", log, "_rho_fitness.csv")
  write.csv(res[[1]], here("data/MC_partitioned", fitness_name))
  
  dyn_name <- paste0("mc_par_", cue, "_", log, "_rho_dyn.parquet")
  write_parquet(res[[3]], here("data/MC_partitioned", dyn_name))
  
  # beta
  res <- mc_five2(par = par,
                  cue = cue,
                  cue_range = cue_range,
                  log = log,
                  rand_var = "beta",
                  rand_mean = 5.721,
                  rand_sd = 0.1722868)
  
  fitness_name <- paste0("mc_par_", cue, "_", log, "_beta_fitness.csv")
  write.csv(res[[1]], here("data/MC_partitioned", fitness_name))
  
  dyn_name <- paste0("mc_par_", cue, "_", log, "_beta_dyn.parquet")
  write_parquet(res[[3]], here("data/MC_partitioned", dyn_name))
  
  # psin
  res <- mc_five2(par = par,
                  cue = cue,
                  cue_range = cue_range,
                  log = log,
                  rand_var = "psin",
                  rand_mean = 16.69234,
                  rand_sd = 0.5778196)
  
  fitness_name <- paste0("mc_par_", cue, "_", log, "_psin_fitness.csv")
  write.csv(res[[1]], here("data/MC_partitioned", fitness_name))
  
  dyn_name <- paste0("mc_par_", cue, "_", log, "_psin_dyn.parquet")
  write_parquet(res[[3]], here("data/MC_partitioned", dyn_name))
  
  # psiw
  res <- mc_five2(par = par,
                  cue = cue,
                  cue_range = cue_range,
                  log = log,
                  rand_var = "psiw",
                  rand_mean = 0.8431785,
                  rand_sd = 0.2355804)
  
  fitness_name <- paste0("mc_par_", cue, "_", log, "_psiw_fitness.csv")
  write.csv(res[[1]], here("data/MC_partitioned", fitness_name))
  
  dyn_name <- paste0("mc_par_", cue, "_", log, "_psiw_dyn.parquet")
  write_parquet(res[[3]], here("data/MC_partitioned", dyn_name))
  
  # phin
  res <- mc_five2(par = par,
                  cue = cue,
                  cue_range = cue_range,
                  log = log,
                  rand_var = "phin",
                  rand_mean = 0.03520591,
                  rand_sd = 0.02609495)
  
  fitness_name <- paste0("mc_par_", cue, "_", log, "_phin_fitness.csv")
  write.csv(res[[1]], here("data/MC_partitioned", fitness_name))
  
  dyn_name <- paste0("mc_par_", cue, "_", log, "_phin_dyn.parquet")
  write_parquet(res[[3]], here("data/MC_partitioned", dyn_name))
  
  # phiw
  res <- mc_five2(par = par,
                  cue = cue,
                  cue_range = cue_range,
                  log = log,
                  rand_var = "phiw",
                  rand_mean = 550.842,
                  rand_sd = 0.8286213)
  
  fitness_name <- paste0("mc_par_", cue, "_", log, "_phiw_fitness.csv")
  write.csv(res[[1]], here("data/MC_partitioned", fitness_name))
  
  dyn_name <- paste0("mc_par_", cue, "_", log, "_phiw_dyn.parquet")
  write_parquet(res[[3]], here("data/MC_partitioned", dyn_name))
}