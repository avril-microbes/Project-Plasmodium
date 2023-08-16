# Base code for running 1000 iterations of semi-stochastic modelling (all six parameters vary). You need to split the 5000 runs
# into 5 iterations or else R will crash!

# par = coefficients for strategy
# cue = cue used
# cue range
# is log true or false?
# set seed for reproducibility
# rand_df = dataframe containing 

mc_run <- function(par, cue, cue_range, log, seed){
  # run 5 iterations each with 1000
  cl <- makeCluster(detectCores()); setDefaultCluster(cl = cl)
  registerDoParallel(cl)
  registerDoRNG(seed)
  mc_res <- foreach(i= 1:1000, .packages = c("doParallel", "doRNG", "deSolve", "splines2", "stringr", "dplyr", "tidyr", "crone", "here")) %dorng% {
    
    # source function
    source(here("functions/chabaudi_si_clean.R"))
    
    # get parameter values that are stochastic
    rho_rand <- 0.2627156*exp(rnorm(n = 1, mean = 0, sd = 0.2579136)) # proportion of RBC recovered
    beta_rand <- 5.721*exp(rnorm(n = 1, mean = 0, sd = 0.1722868)) # burst size
    psin_rand <- 16.69234*exp(rnorm(n = 1, mean = 0, sd = 0.5778196))
    psiw_rand <- 0.8431785*exp(rnorm(n = 1, mean = 0, sd = 0.2355804))
    phin_rand <- 0.03520591*exp(rnorm(n = 1, mean = 0, sd = 0.02609495))
    phiw_rand <- 550.842*exp(rnorm(n = 1, mean = 0, sd = 0.8286213))
    
    parameters_tsukushi <- c(R1 = 8.89*10^6, # slightly higher
                             lambda = 3.7*10^5,
                             mu = 0.025,
                             p = 8*10^-6, # doubled form original
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
    
    time_range <- seq(0, 30, by = 1e-3)
    
    # run dynamics
    chabaudi_si_clean(
      parameters_cr = par, 
      parameters = parameters_tsukushi, 
      immunity = "tsukushi",
      time_range = time_range, 
      cue = cue, 
      cue_range = cue_range, 
      log_cue = log,
      solver = "vode",
      dyn = TRUE)
    
  }
  
  stopCluster(cl)
  
  # process data
  fitness.ls <- mclapply(mc_res, function(x){
    fitness <- x %>% 
      dplyr::filter(variable == "tau_cum") %>% 
      dplyr::summarise(max_fitness = max(value))
    
    return(fitness)
  })
  
  ## get fitness.df
  fitness.df <- do.call(rbind, fitness.ls)
  
  # get cue vs dynamics
  ## produce wide df for plotting cue vs cr graph
  df_wide.ls <- mclapply(mc_res, function(x){
    ## convert to wide
    wide <- tidyr::pivot_wider(x, names_from = variable, values_from = value, id_cols = c(time))
    ## get every 10th row
    wide_filtered <- wide %>% dplyr::filter(row_number() %% 10 == 1)
    return(wide_filtered)
  })
  
  # get wide list
  df_wide.df <- dplyr::bind_rows(df_wide.ls, .id = "id")
  
  # return output
  return(list(fitness.df, df_wide.df))
}