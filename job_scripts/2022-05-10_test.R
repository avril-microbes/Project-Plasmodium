# R script to submit parallelized genetic algorithm optimzied runs.

# load libraries
library(doParallel)
library(optimParallel)
library(here)
library(doRNG)

source(here::here("chabaudi_ci_clean.R"), local = T)
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

I_range <- seq(0, log10(10^7), log10(10^7)/5000)

#-----------------------------#
# Begin parallelized code
#----------------------------#


# Create an array from the NODESLIST environnement variable
nodeslist = unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))

# Create the cluster with the nodes name. One process per count of node name.
# nodeslist = node1 node1 node2 node2, means we are starting 2 processes on node1, likewise on node2.
cl = makeCluster(nodeslist, type = "PSOCK") 
registerDoParallel(cl)
clusterExport(cl,c("time_range", "parameters_tsukushi", "I_range", "chabaudi_ci_clean")) 
clusterCall(cl, library, package = "mclust", character.only = TRUE)
res <- optimParallel::optimParallel(
  par = rep(0.5,4),
  fn = chabaudi_ci_clean, 
  control = list(trace = 6, fnscale = -1),
  parameters_cr_2 = rep(1,4),
  immunity = "tsukushi",
  parameters = parameters_tsukushi,
  time_range = seq(0, 20, by = 1e-3),
  cue_range_1 = seq(0, (6*10^6), by = (6*10^6)/5000),
  cue_range_2 = seq(0, (6*10^6), by = (6*10^6)/5000),
  cue_1 = "I1",
  cue_2 = "Ig2",
  log_cue_1 = "log10",
  log_cue_2 = "none",
  solver = "vode")
stopCluster(cl)

print(res$par)


