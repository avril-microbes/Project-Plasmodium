# function to produce a parameter set for dual cues that mimics that of the best single cue
dual_to_si_glb <- function(df, par_init = rep(0.5,9)){
  ## if first cue is better than second cue
  if(df$fitness_si > df$fitness_si_b){
    # get cue_range_si (single cue model)
    cue_range_si <- seq(df$low_si, df$high_si, length.out = 250)
    
    # get par_si (single cue model parameters)
    par_si <- c(df$var1,df$var2,df$var3,df$var4)
    
    # get cue_range_to_keep (first one)
    cue_range_to_keep <- "a"
  }
  
  ## if second cue is better than second cue
  if(df$fitness_si < df$fitness_si_b){
    # get cue_range_si (single cue model)
    cue_range_si <- seq(df$low_si_b, df$high_si_b, length.out = 250)
    
    # get par_si (single cue model parameters)
    par_si <- c(df$var1_b,df$var2_b,df$var3_b,df$var4_b)
    
    # get cue_range_to_keep (first one)
    cue_range_to_keep <- "b"
  }
  
  # get dual cue ranges
  cue_range <- seq(df$low, df$high, length.out = 250)
  cue_range_b <- seq(df$low_b, df$high_b, length.out = 250)
  
  ## nested objective function to generate a difference metric between the dual and single cue model
  compare_rn <- function(par, cue_range_si, cue_range, cue_range_b, par_si, cue_range_to_keep){
    ## import library and functions
    library(magrittr)
    source(here::here("code_repository/functions/par_to_df.R"))
    source(here::here("code_repository/functions/par_to_hm_te.R"))
    
    # generate single cue reaction norm
    si_rn <- par_to_df(par = par_si, cue_range = cue_range_si, max = max(cue_range_si))
    
    # generate double cue reaction norm
    dual_rn <- par_to_hm_te(par = par,
                            cue_range = cue_range,
                            cue_range_b = cue_range_b)
    
    # left join. If the best performing single cue is a (first cue range), join by first cue range
    if(cue_range_to_keep =="a"){
      si_dual_rn <- dual_rn %>% dplyr::left_join(dplyr::select(si_rn, cue_range, cr_si = cr), by = "cue_range") %>% 
        dplyr::select(cue_range_opt = cue_range, cue_range_sec = cue_range_b, cr, cr_si)
    }
    if(cue_range_to_keep =="b"){
      si_dual_rn <- dual_rn %>% dplyr::left_join(dplyr::select(si_rn, cue_range, cr_si = cr), by = c("cue_range_b" = "cue_range")) %>% 
        dplyr::select(cue_range_opt = cue_range_b, cue_range_sec = cue_range, cr, cr_si)
    }
    
    # calculate metric
    si_dual_rn.diff <-  si_dual_rn %>% 
      dplyr::group_by(cue_range_opt, cr_si) %>% 
      dplyr::summarise(mean_cr = mean(cr),
                       sd_cr = sd(cr)) %>%  # calculate cr variation and mean
      dplyr::mutate(cr_diff = abs(mean_cr - cr_si))
    
    # the metric is simply a sum of the difference in cr and standard deviation
    metric <- sum(si_dual_rn.diff$cr_diff) + sum(si_dual_rn.diff$sd_cr)
    return(metric)
  }
  
  # start optimize to minimize this metric
  cl <- makeCluster(detectCores()); setDefaultCluster(cl = cl)
  
  # start global optimization DEoptim with inital value of 0.5x9
  global <- DEoptim::DEoptim(
    control = list(trace = 1, parallelType = "parallel",
                   NP = 90,itermax = 500, steptol = 50,
                   packages = c("dplyr")), # NP, F, CR values set according to rec by storn et al 2006
    lower = c(-10, -500, -1000, -1000, -250, -500, -1000, -500, -250), # lower and upper bound values are derived empirically based on spline space
    upper = c(10, 500, 1000, 1000, 250, 500, 1000, 500, 250),
    fn = compare_rn, 
    cue_range_si = cue_range_si, 
    cue_range = cue_range, 
    cue_range_b = cue_range_b, 
    par_si = par_si, 
    cue_range_to_keep = cue_range_to_keep)
  
  # secondary local opt
  cl <- makeCluster(detectCores()); setDefaultCluster(cl = cl)
  
  # start optimization LFBGS with inital value of 0.5x5
  res <- optimParallel(
    par = c(global$optim$bestmem),
    fn = compare_rn, 
    cue_range_si = cue_range_si, 
    cue_range = cue_range, 
    cue_range_b = cue_range_b, 
    par_si = par_si, 
    cue_range_to_keep = cue_range_to_keep,
    control = list(trace = 6))
  
  # close cluster
  stopCluster(cl)
  
  # get df
  res_df <- cbind.data.frame(
    df,
    value = res$value,
    par1_init = res$par[[1]],
    par2_init = res$par[[2]],
    par3_init = res$par[[3]],
    par4_init = res$par[[4]],
    par5_init = res$par[[5]],
    par6_init = res$par[[6]],
    par7_init = res$par[[7]],
    par8_init = res$par[[8]],
    par9_init = res$par[[9]]
  )
  
  return(res_df)
}