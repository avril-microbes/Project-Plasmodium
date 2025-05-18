get_si_rn <- function(df){
  ## read in the parameter sets
  par <- c(df$var1, df$var2, df$var3, df$var4)
  ## get cue range
  cue_range <- seq(df$low, df$high, by = df$by)
  ## nested function to convert parameter set into basis spline function
  rn <- par_to_df(par = par, cue_range = cue_range)
  
  # if parasite is sensing logged cue, we have to exponentiate it back so they are on the same scale! 
  rn2 <- rn
  if(stringr::str_detect(df$log, "log")){rn2$cue_range <- 10^(rn2$cue_range)}
  
  # append label to reaction norm, dyn, and rug
  rn2 <- data.frame(rn2, id = df$id)
  
  return(rn2)
}