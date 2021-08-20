# The following function will produce the optimal conversion rate strategy of Plasmodium chabaudi in a single infection scenario. Users must specify:
## 1. Initial parameter (for conversion rate) search space
## 2.whether immunity is present or not
## 3. parameter setting
## 4. time period that infection run
## 5. degrees of freedom for conversion rate spline function
## 6. the "cue" that conversion rate is dependent on
## 7. The range of cue

chabaudi_si_opt <- function(parameters_cr, immunity, parameter, time, df, cue, cue_range){
  #-------------------------#
  # Define initial condition
  #------------------------#
  state <- c(R = 8.5*10^6,
             I = 43.85965,
             Ig = 0,
             M = 0,
             G = 0)
  
  #-------------------------#
  # Function to describe population 
  # structure of initial inoculum
  #------------------------#
  pulseBeta <- function(I0, sp, t){ 
    res = rep(NA, length(t))
    for (num in 1:length(t)){
      res[num] = I0*(dbeta(t[num], sp, sp))
    }
    return(res)}
  
  #-------------------------#
  # Define conversion rate function
  #------------------------#
  
  
}
