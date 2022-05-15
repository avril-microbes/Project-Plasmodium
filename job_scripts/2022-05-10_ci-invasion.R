# R script to run invasion analysis in cedar
# load libraries
library(doParallel)
library(here)
library(doRNG)
library(surveillance)
library(optimParallel)

source(here::here("chabaudi_ci_clean.R"), local = T)

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


# load df
ci_cue_pair.df3 <- read.csv(here("ci_cue_par.csv"))
ci_cue_pair.ls3 <- split(ci_cue_pair.df3, seq(nrow(ci_cue_pair.df3)))


#-------Function-----------#

# load function
ci_invasion_analysis <- function(df, default = F){
  
  # get mutant info (invader)
  if(default == F){mut_par <- c(df$var1, df$var2, df$var3, df$var4)} ## mutant parameter. if default is false, start optimization with optimized parameter set for mutant
  if(default == T){mut_par <- rep(0.5, 4)} ## if default is true, start optimization with 0.5x4
  if(stringr::str_detect(df$cue, "-i")){mut_cue = paste0(gsub("*-i", "", df$cue), "1")} ## mutant cue
  if(stringr::str_detect(df$cue, "-i", negate = T)){mut_cue = df$cue} ## mutant cue
  if(stringr::str_detect(df$log, "log")){mut_log = "log10"} ## mutant log
  if(stringr::str_detect(df$log, "none")){mut_log = "none"} ## mutant log
  mut_cue_range <- seq(df$low, df$high, by = df$by) ## mutant cue range
  
  # get resident info (being invaded)
  res_par <- c(df$var1_2, df$var2_2, df$var3_2, df$var4_2) ## resident parameter
  if(stringr::str_detect(df$cue_2, "-i")){res_cue = paste0(gsub("*-i", "", df$cue_2), "2")} ## resident cue
  if(stringr::str_detect(df$cue_2, "-i", negate = T)){res_cue = df$cue_2} ## resident cue
  if(stringr::str_detect(df$log_2, "log")){res_log = "log10"} ## resident log
  if(stringr::str_detect(df$log_2, "none")){res_log = "none"} ## resident log
  res_cue_range <- seq(df$low_2, df$high_2, by = df$by_2) ## resident cue range
  
  # get ids
  mut_id <- df$id
  res_id <- df$id_2
  
  # one way optimization
  # start optimization LFBGS
  res <- optimParallel::optimParallel(
    par = mut_par,
    fn = chabaudi_ci_clean, 
    control = list(trace = 6, fnscale = -1),
    parameters_cr_2 = res_par,
    immunity = "tsukushi",
    parameters = parameters_tsukushi,
    time_range = seq(0, 20, by = 1e-3),
    cue_range_1 = mut_cue_range,
    cue_range_2 = res_cue_range,
    cue_1 = mut_cue,
    cue_2 = res_cue,
    log_cue_1 = mut_log,
    log_cue_2 = res_log,
    solver = "vode",
    parallel = list(cluster = cl))
  
  # get model output
  fitness <- res$value ## fitness difference between mutant and residence
  mut_opt_par <- res$par ## optimized parameters of mutant
  
  # produce output
  output <- cbind.data.frame(mut_id = mut_id, res_id = res_id, 
                             fitness = fitness, 
                             mut_var1 = df$var1, mut_var2 = df$var2, mut_var3 = df$var3, mut_var4 = df$var4,
                             res_var1 = df$var1_2, res_var2 = df$var2_2, res_var3 = df$var3_2, res_var4 = df$var4_2,
                             mut_var1_opt = mut_opt_par[1], mut_var2_opt = mut_opt_par[2], mut_var3_opt = mut_opt_par[3], mut_var4_opt = mut_opt_par[4])
  write.csv(output, here(paste0("ci_invasion/", mut_id, "_", res_id, ".csv")))
  return(output)
}

#-----------------------------#
# Begin parallelized code
#----------------------------#
nodeslist = unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))
cl = makeCluster(nodeslist, type = "PSOCK") 
registerDoParallel(cl)
clusterExport(cl, c("ci_cue_pair.ls3", "time_range", "parameters_tsukushi", "chabaudi_ci_clean", "ci_invasion_analysis")) 
clusterCall(cl, library, package = "mclust", character.only = TRUE)

surveillance::plapply(ci_cue_pair.ls3, ci_invasion_analysis)

stopCluster(cl)

