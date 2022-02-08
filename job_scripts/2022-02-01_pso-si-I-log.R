# script to submit to compute canada parallelized ppso runs

# load libraries
library(Rmpi)
library(ppso)
library(here)

#--------------------#
# configure rmpi
#---------------------#
sprintf("TEST mpi.universe.size() =  %i", mpi.universe.size())
ns <- mpi.universe.size() - 1
sprintf("TEST attempt to spawn %i slaves", ns)
mpi.spawn.Rslaves(nslaves=ns)
mpi.remote.exec(paste("I am",mpi.comm.rank(),"of",mpi.comm.size()))
mpi.remote.exec(paste(mpi.comm.get.parent()))

#-----------------------------#
# load in parameters
#-----------------------------#
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

time_range <- seq(0, 20, by = 1e-3)

I_range <- seq(0, log10(10^8), by = log10(10^8)/5000)

parameter_bounds <- cbind(c(-5, -100, -500, -1000),
                          c(5, 100, 500, 1000))
#-------------#
# run optim_ppso parallelized
#-------------#
# defining function to send to mpi slave
opt_ppso <- function(){optim_ppso_robust(
  objective_function = function(X)
    chabaudi_si_clean(
      parameters_cr = X, 
      parameters = parameters_tsukushi, 
      time_range = time_range, 
      cue = "I", 
      cue_range = I_range, 
      log_cue = "log10",
      immunity = "tsukushi",
      solver = "vode",
      neg = TRUE),
  number_of_parameters = 4,
  number_of_particles = 70, # on the higher end to https://doi.org/10.1016/j.swevo.2020.100718
  max_number_of_iterations = 1000,
  parameter_bounds = parameter_bounds,
  tryCall = TRUE,
  logfile = "2022-02-01_pso-si-I-log.log"
)}

# run parallelized ppso
set.seed(137)
ppso_res <- mpi.remote.exec(opt_ppso)

# end run
mpi.close.Rslaves()
mpi.quit()
