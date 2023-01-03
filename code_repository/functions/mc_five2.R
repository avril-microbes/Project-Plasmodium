mc_five2 <- function(par, cue, cue_range, log, rand_var, rand_mean, rand_sd){
  mc.1 <- partition_var(par = par,
                        cue = cue,
                        cue_range = cue_range,
                        log = log,
                        rand_var = rand_var,
                        rand_mean = rand_mean,
                        rand_sd = rand_sd, 
                        seed = 137)
  mc.2 <- partition_var(par = par,
                        cue = cue,
                        cue_range = cue_range,
                        log = log,
                        rand_var = rand_var,
                        rand_mean = rand_mean,
                        rand_sd = rand_sd, 
                        seed = 138)
  mc.3 <- partition_var(par = par,
                        cue = cue,
                        cue_range = cue_range,
                        log = log,
                        rand_var = rand_var,
                        rand_mean = rand_mean,
                        rand_sd = rand_sd, 
                        seed = 139)
  mc.4 <- partition_var(par = par,
                        cue = cue,
                        cue_range = cue_range,
                        log = log,
                        rand_var = rand_var,
                        rand_mean = rand_mean,
                        rand_sd = rand_sd, 
                        seed = 140)
  mc.5 <- partition_var(par = par,
                        cue = cue,
                        cue_range = cue_range,
                        log = log,
                        rand_var = rand_var,
                        rand_mean = rand_mean,
                        rand_sd = rand_sd, 
                        seed = 141)
  
  
  # get fitness
  fitness <- rbind(mc.1[[1]],mc.2[[1]],mc.3[[1]],mc.4[[1]],mc.5[[1]])
  
  # produce output for convergence of fitness
  ## get cum sum of fitness
  fitness.sum <- cumsum(fitness)
  ## get average fitness across iteration
  fitness.avg <- fitness.sum/(1:length(fitness.sum))
  
  # get wide data
  dyn <- rbind(
    cbind(mc.1[[2]], run = 1),
    cbind(mc.2[[2]], run = 2),
    cbind(mc.3[[2]], run = 3),
    cbind(mc.4[[2]], run = 4),
    cbind(mc.5[[2]], run = 5))
  
  # get unique id
  dyn <- dyn %>% 
    dplyr::mutate(run_id = paste0(run,"_", id, "_", rand_var))
  
  return(list(fitness, fitness.avg, dyn))
}