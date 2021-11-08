# function to convert dual cue to conversion rate df for plotting

dual_to_df <- function(par, cue_range, cue_range_b, df){
  
  # make dummy y values
  dummy_y.vals <- rep(0, length(cue_range)) 
  dummy_cr.data <- as.data.frame(cbind(cue_range, dummy_y.vals))
  
  # make model
  model <- lm(dummy_y.vals ~ splines2::bSpline(x = cue_range, degree = df, df = df)*
       splines2::bSpline(x = cue_range_b, degree = df, df = df))
  
  # assign coefficients
  model$coefficients <- par
  
  # create expand grid for cue ranges
  cr_grid <- expand.grid(cue_range, cue_range_b)
  names(cr_grid) <- c("cue_range", "cue_range_b")
  
  # predict
  cr_fit <- exp(-exp(predict(model, newdata = cr_grid)))
  cr_df <- data.frame(cr_grid, cr_fit)
  return(cr_df)
}