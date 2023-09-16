# this function performs local optimization with parameter sets with different burst size stipulation
# burst size parameter of 5.721. 

## note if default = F, the optimization will start with an initiap parameter set = optimized parameter when beta = 5.721

opt_local_beta <- function(df, default = F){
  # process input
  beta <- df$beta ## beta values to replace 
  ifelse(default == F, par <- c(df$var1, df$var2, df$var3, df$var4), par <- rep(0.5,4)) ## the original optimized parameter set (with beta = 5.721) is used as the starting point
  time_range <- seq(0, 20, by = 1e-3) ## default time range values
  cue <- df$cue
  ifelse(df$log == "log", log <- "log10", log <- "none")
  cue_range <- seq(df$low, df$high, by = df$by)
  id <- df$id
  
  
  #  obtained altered parmeter set with beta value replaced
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
  
  
  # start cluster
  cl <- makeCluster(detectCores()); setDefaultCluster(cl = cl)
  
  # start optimization LFBGS
  res <- optimParallel(
    par = par,
    fn = chabaudi_si_clean, 
    control = list(trace = 6, fnscale = -1),
    immunity = "tsukushi",
    parameters = parameters_tsukushi_beta,
    time_range = time_range,
    cue = cue,
    log_cue = log,
    cue_range = cue_range,
    solver = "vode")
  
  # close cluster
  stopCluster(cl)
  
  # get model output
  fitness <- res$value
  par <- res$par 
  
  # produce output
  output <- cbind.data.frame(id = id,
                             cue = cue, 
                             log = log,
                             beta = beta,
                             fitness = fitness,
                             par1 = par[1], par2 = par[2], par3 = par[3], par4 = par[4],
                             default = default)
  
  # write output
  beta_label <- gsub("\\.", "-", beta) ### code to substitute the decimal point of beta with (-) for file naming purposes
  
  gsub("5.721", "\\.", "_")
  ifelse(default == F, 
         write.csv(output, here(paste0("code_repository/data/mc_burst_opt/", id, "_", beta_label, ".csv"))),
         write.csv(output, here(paste0("code_repository/data/mc_burst_opt/", id, "_", beta_label, "_default.csv"))))
}