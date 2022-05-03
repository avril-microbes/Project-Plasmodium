# Function to perform MC simulation using the chabaudi_si_sto.R model
# 


monte_carlo <- function(parameters_cr, # parameter sets for strategy (use optimal)
                        parameters, # parameter set for model
                        time_range, # time range
                        cue, # cue to use
                        cue_range, # cue range
                        log_cue, # whether to log cue or not
                        rho_sd, # sd for rho
                        beta_sd, # sd for beta (burst size)
                        psin_sd, # sd for psin (activation strength for general RBC removal)
                        psiw_sd, # sd for psiw (activation strength for targeted iRBC removal)
                        phin_sd, # sd for general RBC removal half-life
                        phiw_sd, # sd for targeted RBC removal half-lfie 
                        factor = 1 # factor to divide sd with
){
  
  # forcing parameters
  force(parameters_cr)
  force(parameters)
  force(log_cue)
  force(cue)
  force(cue_range)
  force(rho_sd)
  force(beta_sd)
  force(psin_sd)
  force(psiw_sd)
  force(phin_sd)
  force(phiw_sd)
  
  ## define the random list of parameter values
  rho_rand <- parameters["rho"]*exp(rnorm(n = length(time_range) + 1, mean = 0, sd = rho_sd/factor)) # proportion of RBC recovered
  beta_rand <- parameters["beta"]*exp(rnorm(n = length(time_range) + 1, mean = 0, sd = beta_sd/factor)) # burst size
  psin_rand <- parameters["psin"]*exp(rnorm(n = length(time_range) + 1, mean = 0, sd = psin_sd/factor))
  psiw_rand <- parameters["psiw"]*exp(rnorm(n = length(time_range) + 1, mean = 0, sd = psiw_sd/factor))
  phin_rand <- parameters["phin"]*exp(rnorm(n = length(time_range) + 1, mean = 0, sd = phin_sd/factor))
  phiw_rand <- parameters["phiw"]*exp(rnorm(n = length(time_range) + 1, mean = 0, sd = phiw_sd/factor))
  
  ## run the numerical simulation
  sim <- chabaudi_si_sto(parameters_cr = parameters_cr,
                             immunity = "tsukushi",
                             parameters = parameters,
                             time_range = time_range,
                             cue = cue,
                             cue_range = cue_range, 
                             solver = "vode",
                             log_cue = log_cue, 
                             rho_rand = rho_rand,
                             beta_rand = beta_rand,
                             psin_rand = psin_rand,
                             psiw_rand = psiw_rand,
                             phin_rand = phin_rand,
                             phiw_rand = phiw_rand,
                             dyn = TRUE)
  
  ## return simulation results
  return(sim)
}