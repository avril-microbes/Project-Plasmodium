# function to perfor global optimization of the single infection model. Hence, the global
# equilevant of opt_local.R script. This uses the DEoptim package which uses the differential
# evolution algorithm. I have success with using the algorithm over other global optimization algo
# such as particle swarming and genetics algorthim. DE (differential evolution) works in similar ways
# than genetics algorithm. DE is generally suitable for  non-smooth solution spaces (whereby most
# derivative-based optimization algo work best across smooth surfaces) or solution spaces with many
# local optimums. For our purposes, we are cognscent of the probability of local optimums so 
# this might be a good approach.

opt_global <- function(cue, #
                       log,
                       cue_range,
                       dyn = F, # upper bound of solution space
                       ...){
  
  ## for time based optimization, need to match cue range to time range so overrides existing input cue_range
  if(cue == "t"){cue_range_global <- seq(0, 20, by = 1e-2)} 
  
  # start cluster
  cl <- makeCluster(8); setDefaultCluster(cl = cl)
  
  # run optimization. note that this is ran with a larger time step to increase speed
  
  global <- DEoptim::DEoptim(
    fn = chabaudi_si_clean,
    control = list(trace = 1, parallelType = "parallel", itermax = 500, steptol = 50),
    lower = c(-5, -500, -500, -500), # lower and upper bound values are derived empirically based on spline space
    upper = c(5, 500, 500, 500),
    immunity = "tsukushi",
    parameters = parameters_tsukushi,
    time_range = seq(0, 20, by = 1e-2), # increase speed
    cue = cue,
    log_cue = log,
    cue_range = cue_range_global,
    solver = "vode",
    neg = T # DE Is minimization function so must return negative fitness
  )
  
  # second round of optimzation with 0.001 time step
  res <- optimParallel(
    par = c(global$optim$bestmem), # starting value is the optimal from global 
    fn = chabaudi_si_clean, 
    control = list(trace = 6, fnscale = -1),
    immunity = "tsukushi",
    parameters = parameters_tsukushi,
    time_range = seq(0, 20, by = 1e-3), # small time step
    cue = cue,
    log_cue = log,
    cue_range = cue_range,
    solver = "vode")
  
  # close cluster
  stopCluster(cl)
    
  # save output
  return(res)
}