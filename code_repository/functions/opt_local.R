opt_local <- function(
  par, # initial parameter sets
  cue, #
  log,
  cue_range,
  dyn = F
) {
  # start cluster
  cl <- makeCluster(detectCores()); setDefaultCluster(cl = cl)
  
  # start optimization LFBGS
  res <- optimParallel(
    par = par,
    fn = chabaudi_si_clean, 
    control = list(trace = 6, fnscale = -1),
    immunity = "tsukushi",
    parameters = parameters_tsukushi,
    time_range = time_range,
    cue = cue,
    log_cue = log,
    cue_range = cue_range,
    solver = "vode")
  
  # close cluster
  stopCluster(cl)
  
  return (res)
}