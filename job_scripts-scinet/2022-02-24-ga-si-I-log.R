# R script to submit parallelized genetic algorithm optimzied runs.

# load libraries
library(doParallel)
library(parallel)
library(GA)
library(here)
library(doRNG)

# load functions
source(here::here("functions/ga_verbose.R"), local = T)
source(here::here("functions/chabaudi_si_clean.R"), local = T)

######################################################
#
#  worker code
#

# first define the function which will be run on all the cluster nodes.  This is just a test function.  
# Put your real worker code here.
testfunc <- function(a) {
  
  # this part is just to waste time
  b <- 0
  for (i in 1:100000000) {
    b <- b + 1
  }
  
  s <- Sys.info()['nodename']
  return(paste0(s, " ", a[1], " ", a[2]))
  
}


######################################################
#
#  head node code
#

# Create a bunch of index pairs to feed to the worker function.  These could be parameters,
# or whatever your code needs to vary across jobs.  Note that the worker function only 
# takes a single argument; each entry in the list must contain all the information 
# that the function needs to run.  In this example, each entry contains a list which
# contains two pieces of information, a pair of indices.
indexlist <- list()
index <- 1
for (i in 1:10) {
  for (j in 1:10) {
    indexlist[index] <- list(c(i,j))
    index <- index +1
  }
}


# Now set up the cluster.

# First load the parallel library.
library(parallel)

# Next find all the nodes which the scheduler has given to us.
# These are given by the SLURM_JOB_NODELIST environment variable.
nodelist <- Sys.getenv("SLURM_JOB_NODELIST")
# Get your SCRATCH directory.
my.scratch <- Sys.getenv("SCRATCH")
node_ids <- unlist(strsplit(nodelist,split="[^a-z0-9-]"))[-1]

if (length(node_ids)>0) {
  expanded_ids <- lapply(node_ids, function (id) {
    ranges <- as.numeric(
      unlist(strsplit(id, split="[-]"))
    )
    if (length(ranges)>1) seq(ranges[1], ranges[2], by=1) else ranges
  })
  
  nodelist <- sprintf("nia%04d", unlist(expanded_ids))
}

## We now have the nodelist, but we need to repeat it for each core in the nodes.
nodelist <- rep(nodelist, 40)

# Now launch the cluster, using the list of nodes and our Rscript
# wrapper.  This assumes that your MyRscript.sh is at the base of your SCRATCH directory.
cl <- makePSOCKcluster(names = nodelist, rscript = paste(my.scratch, "/job_scripts/2022-02-24-ga-si-I-log.sh", sep = ""))

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

#----------------------#
# run simulation
#-----------------------#
registerDoParallel(cl)
clusterExport(cl, c("ga_verbose", "time_range", "parameters_tsukushi", "I_range", "chabaudi_si_clean")) 
clusterApplyLB(cl, indexlist, library, package = "mclust", character.only = TRUE)
ga_res <- ga_verbose(type = "real-valued", 
                     function(x)
                       chabaudi_si_clean(
                         parameters_cr = c(x[1], x[2], x[3], x[4]), 
                         parameters = parameters_tsukushi, 
                         time_range = time_range, 
                         cue = "I", 
                         cue_range = I_range, 
                         log_cue = "log10",
                         immunity = "tsukushi",
                         solver = "vode"),
                     lower = c(-5, -100, -500, -1000), # range determined that would alter shape of spline
                     upper = c(5,100,500, 1000),  
                     popSize = 250, 
                     maxiter = 1000, # change to 10 for testing purpose 
                     pmutation = 0.3,
                     keepBest = TRUE,
                     run =50,
                     parallel = cl,
                     seed = 137,
                     monitor = TRUE,
                     id = "2022-02-22_ga-si-I-log")
stopCluster(cl)

print(list(ga_res@bestSol, ga_res@fitnessValue))


