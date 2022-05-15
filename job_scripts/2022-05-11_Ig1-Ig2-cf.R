# compute canada code for running validation 

# load packages
library(doParallel)
library(here)
library(doRNG)

source(here("chabaudi_ci_clean.R"))

# parameters 
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

# cluster initialization
ncores = Sys.getenv("SLURM_CPUS_PER_TASK")
registerDoParallel(cores=ncores)# Shows the number of Parallel Workers to be used
print(ncores) # this how many cores are available, and how many you have requested.
getDoParWorkers()# you can compare with the number of actual workers
registerDoRNG(137)

res <- foreach(i= 1:1000, .verbose = T, .packages = c("doParallel", "doRNG", "deSolve", "splines2", "stringr", "dplyr", "tidyr", "crone", "here")) %dorng% {
  # generate random par, bounds for parameters based on previous determined bounds
  par1 <- runif(1, -1, 1)
  par2 <- runif(1, -100, 100)
  par3 <- runif(1, -1000, 10000)
  par4 <- runif(1, -5000, 5000)
  
  par_rand <- c(par1, par2, par3, par4)
  
  # get fitness. 20 days simulation to cut down on computation time
  fitness <- chabaudi_ci_clean(parameters_cr_1 = c(1.000000,  -32.831645,  378.108931,-4999.999969),
                               parameters_cr_2 = par_rand,
                               immunity = "tsukushi",
                               parameters = parameters_tsukushi,
                               cue_1 = "Ig1+Ig2",
                               cue_2 = "Ig1+Ig2",
                               cue_range_1 = seq(0, 6*(10^6), by = (6*(10^6))/5000),
                               cue_range_2 = seq(0, 6*(10^6), by = (6*(10^6))/5000),
                               log_cue_1 = "none",
                               log_cue_2 = "none",
                               solver = "vode",
                               time_range = seq(0, 30, by = 1e-3),
                               dyn = F)
  
  print(fitness)
}

## get fitness.df
fitness.df <- do.call("rbind", res)

## rbind
fitness.df <- as.data.frame(cbind(fitness.df, label = rep("Total sexual iRBC", nrow(fitness.df))))
write.csv(fitness.df, here("validation/Ig1_Ig2_cf.csv"))

