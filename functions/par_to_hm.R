# function to convert coefficients into a heatmap

par_to_hm <- function(par, cue_range, cue_range_b, plot = FALSE, dyn= F){
  
  # get cr grid
  cr_grid <- expand.grid(cue_range, cue_range_b)
  
  ## rename
  names(cr_grid) <- c("cue_range", "cue_range_b")
  
  ## create dummy y
  dummy_y <- runif(length(cue_range_b), 0, 1)
  
  ## put together df
  dummy_df <- data.frame(cue_range, cue_range_b, dummy_y)
  
  ## gam model
  dummy_cr.mod <- mgcv::gam(dummy_y ~ ti(cue_range, cue_range_b, 
                                         k = c(3,3)), 
                            data = dummy_df)
  
  ## assign parameters
  dummy_cr.mod$coefficients <- par
  
  # predict
  cr_res <- exp(-exp(mgcv::predict.gam(dummy_cr.mod, 
                                       newdata = cr_grid)))
  
  cr_res2 <- cbind(cr_grid, cr = cr_res)
  
  # plot
  if(isTRUE(plot)){
    plot <- ggplot2::ggplot() +
      geom_raster(data = cr_res2, aes(x = cue_range, y = cue_range_b, fill = cr)) +
      scale_fill_viridis_c() +
      theme_bw()
    print(plot)
  }
  
  if(!isTRUE(plot)){return(cr_res2)}
  
}