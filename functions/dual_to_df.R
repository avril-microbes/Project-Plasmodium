# function to convert dual cue to conversion rate df for plotting

dual_to_df <- function(par, cue_range, cue_range_b){
  
  # make dummy y values
  cr_grid <- expand.grid(cue_range, cue_range_b)
  ## rename
  names(cr_grid) <- c("cue_range", "cue_range_b")
  ## create dummy y
  dummy_y <- runif(length(cue_range_b), 0, 1)
  ## put together df
  dummy_df <- data.frame(cue_range, cue_range_b, dummy_y)
  
  dummy_cr.mod <- mgcv::gam(dummy_y ~ ti(cue_range, cue_range_b, 
                                         k = c(3,3)), 
                            data = dummy_df, 
                            method = "REML")
  
  # assign parameters
  dummy_cr.mod$coefficients <- parameters_cr

  cr_fit <- exp(-exp(predict(dummy_cr.mod, newdata = cr_grid)))
  cr_df <- data.frame(cr_grid, cr_fit)
  return(cr_df)
}
