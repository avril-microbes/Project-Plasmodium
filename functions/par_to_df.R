# function to convert the output form an "opt" function to a dataframe of cue range vs conversion rate
par_to_df <- function(mod_opt, cue_range){
  # get parameter values
  par <- mod_opt$par
  
  # get df from length of parameter values
  df <- length(par)-1

  # define basic spline model
  basis_mod <- splines::bs(cue_range, df)
  
  # iterate through parameter values and basic spline model to create conversion rate values
  model_part <- list()
  for(i in 1:length(par)){
    model_part[[i]] <- par[i]*basis_mod[,i-1]
  }
  ## sum the conversion rate values up
  model_part2 <- rowSums(do.call(cbind, model_part))
  
  ## need to add par[1] given it is not included. Double exp to limit conversion rate between 0 and 1
  model <- exp(-exp(par[1]+model_part2))
  
  # Create dataframe from model
  s_ni_t.df <- data.frame(cue_range, model)
}


