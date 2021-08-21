# function to convert parameter from optimal parameter values to a dataframe of cue range vs conversion rate
par_to_df <- function(mod_opt, cue_range){
  # get parameter values
  par <- mod_opt$par
  
  # get df from length of parameter values
  df <- length(par)-1

  # define basic spline model
  basis_mod <- splines::bs(cue_range, df)
  
  # define function that creates spline function
  ## function creates individual parts of spline function
  model_part <- list()
  for(i in 1:length(par)){
    model_part[[i]] <- par[i]*basis_mod[,i-1]
  }
  ## sum them up
  model_part2 <- rowSums(do.call(cbind, model_part))
  
  ## need to add par[1] given it is not included. Double exp to limit between 0 and 1
  model <- exp(-exp(par[1]+model_part2))
  
  # Create dataframe from model
  s_ni_t.df <- data.frame(cue_range, model)
}


