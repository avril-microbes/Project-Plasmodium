mc_five <- function(par, cue, cue_range, log){
  mc.1 <- mc_run(par = par,
                 cue = cue,
                 cue_range = cue_range,
                 log = log,
                 seed = 137)
  mc.2 <- mc_run(par = par,
                 cue = cue,
                 cue_range = cue_range,
                 log = log,
                 seed = 138)
  mc.3 <- mc_run(par = par,
                 cue = cue,
                 cue_range = cue_range,
                 log = log,
                 seed = 139)
  mc.4 <- mc_run(par = par,
                 cue = cue,
                 cue_range = cue_range,
                 log = log,
                 seed = 140)
  mc.5 <- mc_run(par = par,
                 cue = cue,
                 cue_range = cue_range,
                 log = log,
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
    dplyr::mutate(run_id = paste0(run,"_", id))
  
  return(list(fitness, fitness.avg, dyn))
}