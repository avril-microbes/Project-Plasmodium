# function to convert the output form an "opt" function to a data frame of 
# cue range vs conversion rate

# to plot the resulting data frame, using the following script:
## library(ggplot)
## ggplot() +
  # geom_line(data = df, aes(x = cue_range, y = cr))


par_to_df <- function(par, cue_range, max, transformation = "exp"){
  
  # ensure correct model arguments
  if(transformation != "norm" && transformation != "exp" && transformation != "logit"){
    stop("Transformation must be either 'norm' or 'exp' or 'logit'")
  }
  
  # heaviside transformation function
  heaviside_trans <- function(cue_range, max){
    res <- crone::heaviside(cue_range)*(cue_range)+(crone::heaviside(cue_range-max)*(max-cue_range))
    return(res)
  }
  
  # get  heaviside transformed cue_range
  cue_range_h <- heaviside_trans(cue_range, max = max)
  
  # get df from length of parameter values
  df <- length(par)-1

  # define basic spline model
  basis_mod <- splines2::bSpline(x = cue_range_h, df = df)
  
  # iterate through parameter values and basic spline model to create conversion rate values
  model_part <- list()
  for(i in 1:length(par)){
    model_part[[i]] <- par[i]*basis_mod[,i-1]
  }
  ## sum the conversion rate values up
  model_part2 <- rowSums(do.call(cbind, model_part))
  
  ## need to add par[1] given it is not included. Double exp to limit conversion rate between 0 and 1
  #cr <- exp(-exp(par[1]+model_part2))
  
  # if transformation is specified to be normalization, perform normalization
  if(transformation == "norm"){
  cr_pre <- par[1]+model_part2
  cr <- (cr_pre-min(cr_pre, na.rm = TRUE))/(max(cr_pre, na.rm = TRUE)-min(cr_pre, na.rm = TRUE))
  } 
  
  if(transformation == "exp"){
    cr <- exp(-exp(par[1]+model_part2))
  } else{
    cr <- 1/(1+exp(-(par[1]+model_part2)))
  }
  
  # Create dataframe from model
  cr.df <- data.frame(cue_range, cr)
}


