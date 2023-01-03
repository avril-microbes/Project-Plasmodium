validate_si <- function(df){
  # cluster initialization
  cl <- makeCluster(detectCores()); setDefaultCluster(cl = cl)
  registerDoParallel(cl)
  registerDoRNG(137)
  
  res <- foreach(i= 1:1000, .verbose = T, .packages = c("doParallel", "doRNG", "deSolve", "splines2", "stringr", "dplyr", "tidyr", "crone", "here")) %dorng% {
    # generate random par, bounds for parameters based on previous determined bounds
    par1 <- runif(1, -1, 1)
    par2 <- runif(1, -100, 100)
    par3 <- runif(1, -1000, 10000)
    par4 <- runif(1, -5000, 5000)
    
    # get random parameters
    par_rand <- c(par1, par2, par3, par4)
    
    # get log
    if(df$log == "log"){log = "log10"}
    if(df$log == "none"){log = "none"}
    
    # get cue_range
    cue_range <- seq(df$low, df$high, by = df$by)
    
    # get fitness. 20 days simulation to cut down on computation time
    fitness <- chabaudi_si_clean(parameters_cr = par_rand,
                                 parameters= parameters_tsukushi,
                                 immunity = "tsukushi",
                                 time_range = seq(0,30,0.001),
                                 cue = df$cue,
                                 log_cue = log,
                                 cue_range = cue_range,
                                 dyn = F)
    
    print(fitness)
  }
  
  ## get fitness.df
  fitness.df <- do.call("rbind", res)
  ## rbind
  fitness.df <- as.data.frame(cbind(fitness.df, 
                                    cue = rep(df$cue, nrow(fitness.df)),
                                    log = rep(df$log, nrow(fitness.df)),
                                    id = rep(df$id, nrow(fitness.df))))
  
  filename <- paste0("si_", df$id, "_cf.csv")
  
  ## write.csv
  write.csv(fitness.df, here(paste0("data/si_validation/", filename)))
}