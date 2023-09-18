dual_cue_opt_beta <- function(df, default = F){
  
  # process inputs
  ifelse(default == F, par <- c(df$par1, df$par2,df$par3,df$par4,df$par5,df$par6,df$par7,df$par8,df$par9),
         par <- rep(0.5, 9))
  id <- df$id
  id_b <- df$id_b
  beta <- df$beta
  cue <- df$cue
  cue_b <- df$cue_b
  
  # process log
  ifelse(df$log == "log", log <- "log10", log <- "none")
  ifelse(df$log_b == "log", log_b <- "log10", log_b <- "none")
  
  # process cue_range. cannot use by to ensure that both cue ranges are of the same length
  cue_range <- seq(df$low, df$high, length.out = 500)
  cue_range_b <- seq(df$low_b, df$high_b, length.out = 500)
  
  # parameter sets
  parameters_tsukushi_beta <- c(R1 = 8.89*10^6, 
  lambda = 3.7*10^5,
  mu = 0.025, 
  p = 8*10^-6, 
  alpha = 1, 
  alphag = 2, 
  beta = beta, ## instead of 5.721, the beta value is substituted
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
  
  # optimization
  cl <- makeCluster(detectCores()); setDefaultCluster(cl = cl)
  
  # start optimization LFBGS with inital value of 0.5x5
  res <- optimParallel(
    par = par,
    fn = chabaudi_si_clean, 
    control = list(trace = 6, fnscale = -1),
    immunity = "tsukushi",
    parameters = parameters_tsukushi_beta,
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
  output <- cbind.data.frame(id = paste0(id, "+", id_b), 
                             cue = paste(paste(cue, "+", cue_b)),
                             fitness = fitness,
                             par1 = par[1], par2 = par[2], par3 = par[3], par4 = par[4], par5 = par[5], par6 = par[6], par7 = par[7], par8 = par[8], par9 = par[9],
                             default = ifelse(default == F, F, T))
  
  # write output
  beta_label <- gsub("\\.", "-", beta) ### code to substitute the decimal point of beta with (-) for file naming purposes
  
  if(default == F){write.csv(output, here(paste0("code_repository/data/mc_burst_opt/", id, "_", id_b, "_", beta_label, ".csv")))
  } else{
    write.csv(output, here(paste0("code_repository/data/mc_burst_opt/",  id, "_", id_b, "_", beta_label, "_default.csv")))
  }
  
}