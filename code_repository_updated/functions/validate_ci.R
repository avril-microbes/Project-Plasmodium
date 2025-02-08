validate_ci <- function(df, n = 1000){
  
  # get label
  label <- df$label
  
  # cluster initialization
  cl <- makeCluster(detectCores()); setDefaultCluster(cl = cl)
  registerDoParallel(cl)
  registerDoRNG(137)
  
  res <- foreach(i= 1:n) %dorng% {
    # source function
    source(here::here("functions/chabaudi_ci_clean.R"))
    
    # parameters 
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
    
    # read in df info
    if(stringr::str_detect(df$cue, "-i") && df$cue != "I-i+Ig-i"){cue_1 = paste0(gsub("*-i", "", df$cue), "1")}
    if(stringr::str_detect(df$cue, "-i", negate = T)){cue_1 = df$cue}
    if(stringr::str_detect(df$cue, "-i") &&  df$cue != "I-i+Ig-i"){cue_2 = paste0(gsub("*-i", "", df$cue), "2")}
    if(stringr::str_detect(df$cue, "-i", negate = T)){cue_2 = df$cue}
    if(df$id == "I-i+Ig-i_none" | df$id == "I-i+Ig-i_log"){
      cue_1 = "I1+Ig1"
      cue_2 = "I2+Ig2"
    }
    if(stringr::str_detect(df$log, "log")){log_1 = "log10"}
    if(stringr::str_detect(df$log, "none")){log_1 = "none"}
    if(stringr::str_detect(df$log, "log")){log_2 = "log10"}
    if(stringr::str_detect(df$log, "none")){log_2 = "none"}
    par <- c(df$var1, df$var2, df$var3, df$var4)
    cue_range <- seq(df$low, df$high, by = df$by)
    
    # generate random par, bounds for parameters based on previous determined bounds
    par1 <- runif(1, -1, 1)
    par2 <- runif(1, -100, 100)
    par3 <- runif(1, -1000, 10000)
    par4 <- runif(1, -5000, 5000)
    
    par_rand <- c(par1, par2, par3, par4)
    
    # get fitness. 20 days simulation to cut down on computation time
    fitness <- chabaudi_ci_clean(parameters_cr_1 = par,
                                 parameters_cr_2 = par_rand,
                                 immunity = "tsukushi",
                                 parameters = parameters_tsukushi,
                                 cue_1 = cue_1,
                                 cue_2 = cue_2,
                                 cue_range_1 = cue_range,
                                 cue_range_2 = cue_range,
                                 log_cue_1 = log_1,
                                 log_cue_2 = log_2,
                                 solver = "vode",
                                 time_range = seq(0, 20, by = 1e-3),
                                 dyn = F)
    
  }
  stopCluster(cl)
  
  ## get fitness.df
  fitness.df <- do.call("rbind", res)
  
  ## rbind
  fitness.df <- as.data.frame(cbind(fitness.df, label = rep(label, nrow(fitness.df))))
  
  filename <- paste0(df$id, "_validation.csv")
  write.csv(fitness.df, paste0(here("data/ci_validation/"), filename))
  return(fitness.df)
}