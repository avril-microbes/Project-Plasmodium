# R script to submit parallelized genetic algorithm optimzied runs.

# load libraries
library(doParallel)
library(GA)
library(here)
library(doRNG)

source(here::here("functions/ga_verbose.R"), local = T)
source(here::here("functions/chabaudi_si_clean.R"), local = T)
# load function
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

R_range <- seq(log10(10^6), log10(10^8), (log10(10^8)-log10(10^6))/5000)

#-----------------------------#
# Begin parallelized code
#----------------------------#


# Create an array from the NODESLIST environnement variable
nodeslist = unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))

# Create the cluster with the nodes name. One process per count of node name.
# nodeslist = node1 node1 node2 node2, means we are starting 2 processes on node1, likewise on node2.
cl = makeCluster(nodeslist, type = "PSOCK") 
registerDoParallel(cl)
clusterExport(cl,c("ga_verbose", "time_range", "parameters_tsukushi", "R_range", "chabaudi_si_clean")) 
clusterCall(cl, library, package = "mclust", character.only = TRUE)
ga_res <- ga_verbose(type = "real-valued", 
                     function(x)
                       chabaudi_si_clean(
                         parameters_cr = c(x[1], x[2], x[3], x[4]), 
                         parameters = parameters_tsukushi, 
                         time_range = time_range, 
                         cue = "R", 
                         cue_range = R_range, 
                         log_cue = "log10",
                         immunity = "tsukushi",
                         solver = "vode"),
                     lower = c(-1, -100, -1000, -5000), # range determined that would alter shape of spline
                     upper = c(1, 100, 1000, 5000),  
                     popSize = 250, 
                     maxiter = 1000, # change to 10 for testing purpose 
                     pmutation = 0.3,
                     keepBest = TRUE,
                     run = 50,
                     parallel = cl,
                     seed = 137,
                     monitor = TRUE,
                     id = "2022-03-05_si-R-log-emp")
stopCluster(cl)

print(list(ga_res@bestSol, ga_res@fitnessValue))


