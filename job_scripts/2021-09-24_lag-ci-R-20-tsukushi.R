library(deSolve)
library(splines)
library(optimParallel)
library(magrittr)
library(here)


source(here("functions/par_to_df.R"))
source(here("functions/chabaudi_ci_opt_lag.R"))
source(here("functions/co_infection_opt.R"))

parameters_tsukushi <- c(R1 = 8.89*10^6, # slightly higher
                         lambda = 3.7*10^5,
                         mu = 0.025, 
                         p = 8*10^-6, # doubled form original
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

time_range <- seq(0, 20, 1e-3)
R_range_log10 <- seq(0, log10(8.89*10^6), by = log10(8.89*10^6)/10000)

cl <- makeCluster(detectCores()); setDefaultCluster(cl = cl)
lag_ci_I_20.opt_tsukushi_log10_R <- co_infection_opt(parameters_cr = rep(0.5, 4),
                                                     limit = 0.01, 
                                                     model = chabaudi_ci_opt_lag,
                                                     immunity = "tsukushi",
                                                     parameters = parameters_tsukushi,
                                                     time_range = time_range,
                                                     df = 3,
                                                     cue_1 = "R",
                                                     cue_2 = "R",
                                                     cue_range_1 = R_range_log10,
                                                     cue_range_2 = R_range_log10,
                                                     solver = "vode",
                                                     log_cue_1 = "log10",
                                                     log_cue_2 = "log10")
stopCluster(cl)

print(lag_ci_I_20.opt_tsukushi_log10_R)

