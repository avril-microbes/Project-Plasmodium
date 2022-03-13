# load library
library(doParallel) 
library(here)
library(doRNG)
library(arrow)

# source functions
source(here::here("chabaudi_si_sto.R"), local = T)
source(here::here("monte_carlo.R"), local = T)

# Create an array from the NODESLIST environnement variable
nodeslist = unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))

# Create the cluster with the nodes name. One process per count of node name.
# nodeslist = node1 node1 node2 node2, means we are starting 2 processes on node1, likewise on node2.
cl = makeCluster(nodeslist, type = "PSOCK") 
registerDoParallel(cl)

# you can compare with the number of actual workers
registerDoRNG(137) # set seed for reproducibility

# load parameters
#### CHANGE THIS BETWEEN RUNS!!!!
cue_range <- seq(0, log10(10^7), by = log10(10^7)/5000)
cue <- "I"
log <- "log10"

#### stays the same between simulations
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

# run monte carlo simulation. 10000 to be exact
res <- foreach(i= 1:10, .packages = c("doParallel", "doRNG", "deSolve", "splines2", "stringr", "dplyr", "tidyr", "crone")) %dorng% {
  monte_carlo(
    parameters_cr = c(6.437168, 0.774665, -19.588253 , 5.618792), 
    parameters = parameters_tsukushi, 
    time_range = time_range, 
    cue = cue, 
    cue_range = cue_range, 
    log_cue = log,
    rho_sd = 0.2579136, 
    beta_sd = 0.1722868, 
    psin_sd = 0.5778196,
    psiw_sd = 0.2355804, 
    phin_sd = 0.02609495, 
    phiw_sd = 0.8286213)}
stopCluster(cl)

res.df <- dplyr::bind_rows(res, .id = "id")

# write the result file in the MC folder
write_parquet(res.df, here("MC/test.parquet"))

