# Deprecated version of chabaudi_si_opt_fast function. Still runs but lacks customizability and is very slow.


# The following function will produce the optimal conversion rate strategy of Plasmodium chabaudi in a single infection scenario. Users must specify:
## 1. Initial parameter (for conversion rate) search space
## 2. whether immunity is present or not ("ni" for no immunity and "si" for saturated immunity)
## 3. parameter setting
## 4. time period that infection run
## 5. degrees of freedom for conversion rate spline function
## 6. the "cue" that conversion rate is dependent on
## 7. The range of cue

## To run the function, please see "plasticity_model_journals.Rmd" for reference". Function best run 
## when paired with optimParallel to allow for parallel computing.

chabaudi_si_opt <- function(parameters_cr, immunity, parameters, time_range, df, cue, cue_range){
  #-------------------------#
  # Ensure values we inputted 
  # are available in environment
  #------------------------#
  force(parameters_cr)
  force(immunity)
  force(parameters)
  force(time_range)
  force(df)
  force(cue)
  force(cue_range)
  
  #-------------------------#
  # Ensure inputs are correct
  #------------------------#
  ## Ensure length of initial parameter search space matches with df
  if (length(parameters_cr) != df+1 ) {
    stop("Conversion rate parameters must match degrees of freedom")
  }
  ## Ensure immunity input is correct
  if (immunity != "ni" && immunity != "i") {
    stop("Immunity must be either 'ni' for no immunity or 'i' for saturating immunity")
  }
  
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
  ## Define dummy data of conversion rate
  dummy_y.vals <- rep(0, length(cue_range)) 
  dummy_cr.data <- as.data.frame(cbind(cue_range, dummy_y.vals))
  
  ## Define initial conversion rate as basic-Spline
  dummy_cr.mod <- lm(dummy_y.vals ~ splines::bs(cue_range, df), data = dummy_cr.data)
  dummy_cr.mod$data <- dummy_cr.data
    
  ## Assign coefficient to be optimized to the dummy conversion rate function
  dummy_cr.mod$coefficients <- parameters_cr
  
  ## Function to get conversion rates at various cue values (this speeds up fitting process)
  spline.fun <- function(x){
    cr.int <- exp(-exp(predict(dummy_cr.mod, newdata =
                                 data.frame(cue_range = x))))
    return(cr.int[[1]])}
  
  ## Put the conversion rate into a list
  cr.int2 <- rep(NA, length(cue_range)) # define empty list
  cr.int2 <- sapply(cue_range, spline.fun) # loop over function to get fitted value
  
  ## Produce a table with time value and conversion rate
  cr <- splinefun(cbind(cue_range, cr.int2))

  
  
  
  #-------------------------#
  # Define single-infection model
  #------------------------#
  single_infection.fun <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {

      ## Define the lag terms. lag[1] = R, lag[2] = I, lag[3] = Ig, lag[4] = M, lag[5] = G
      if(t>alpha){lag1 = deSolve::lagvalue(t-alpha)} # lag state for asexual development
      if(t>alphag){lag2 = deSolve::lagvalue(t-alphag)} # lag state for gametocyte development
      
      ## get lag term index given cue
      lag.i <- match(cue, names(state))
      
      ## convert cue to variable in state
      cue_state <- get(cue)
      
      ## Define K, carrying capacity of RBC
      K <- lambda*R1/(lambda-mu*R1)
      
      ## Define survival functions
      ### Survival of infected asexual RBC
      if(t>alpha && immunity == "ni"){
        S <- exp(-mu*alpha)
      } else if(t>alpha && immunity =="i"){
        integrand <- function(x) {mu+a/(b+I)}
        integrate_val <- integrate(Vectorize(integrand), lower = t-alpha, upper = t)
        S <- exp(-1*integrate_val$value)
      } else if(t<=alpha && immunity == "ni"){
        S <- exp(-mu*t)
      } else{
        integrand <- function(x) {mu+a/(b+I)}
        integrate_val <- integrate(Vectorize(integrand), lower = 0, upper = t)
        S <- exp(-1*integrate_val$value)
      }

      ### Survival of gametocytes. We assume that infected
      ### RBC with gametocyte is not removed by immune response
      if(t<=alphag){
        Sg <- exp(-mu*t)
      } else{
      Sg <- exp(-mu*alphag)}
      
      ## Define the models without lag terms. 
      dR <- lambda*(1-R/K)-mu*R-p*R*M # change in susceptible RBC
      if(immunity == "ni"){
        dI_nolag <- (1-cr(cue_state))*p*R*M-mu*I # change in infected RBC density
      } else{
        dI_nolag <- (1-cr(cue_state))*p*R*M-mu*I-(a*I)/(b+I) # change in infected RBC density with immunity
      }
      dIg_nolag <- cr(cue_state)*p*R*M-mu*Ig 
      dM_nolag <- -mum*M-p*R*M
      dG_nolag <- -mug*G
      
      ## Track states in initial cohort of infection
      if(t<=alpha){
        dI <- dI_nolag-pulseBeta(I0, sp, t)*S 
        dM <- dM_nolag+beta*pulseBeta(I0, sp, t)*S 
      }
      
      if(t<=alphag){
        dIg <- dIg_nolag-pulseBeta(Ig0, sp, t-1)*Sg 
        dG <- dG_nolag+pulseBeta(Ig0, sp, t-1)*Sg 
      }
      
      ## Track states after delay 
      if(t>alpha){
        dI <- dI_nolag-(1-cr(lag1[lag.i]))*p*lag1[1]*lag1[4]*S 
        dM <- dM_nolag+beta*(1-cr(lag1[lag.i]))*p*lag1[1]*lag1[4]*S 
      }
      
      if(t>alphag){
        dIg <- dIg_nolag-cr(lag2[lag.i])*p*lag2[1]*lag2[4]*Sg 
        dG <- dG_nolag+cr(lag2[lag.i])*p*lag2[1]*lag2[4]*Sg
      }
      
      ## Return the states
      return(list(c(dR, dI, dIg, dM, dG)))
    })
  }
  
  #-------------------------#
  # Run single-infection model
  #------------------------#
  single_infection.df <- as.data.frame(deSolve::dede(y = state,
                                                     times = time_range,
                                                     func = single_infection.fun,
                                                     p = parameters,
                                                     control=list(mxhist = 1e6)))
  
  #-------------------------#
  # Calculate fitness
  #------------------------#
  ## Get Gametocyte density time series data
  gam <- single_infection.df$G
  gam[gam<0] <- 0 # Assign negative gametocyte density to 0
  
  ## Get timeseries interval. Simplify first time after t=0
  int <- single_infection.df[2,1]
  
  ## Define the fitness parameter values
  aval <- -12.69
  bval <- 3.6
  dens <- log10(gam)
  
  ## Calculate the transmission potential at each time t
  tau.ls <- (exp(aval+bval*dens))/(1+exp(aval+bval*dens))
  
  ## Get approximation of cumulative transmission potential
  tau.sum <- sum(tau.ls*int)
  
  # return cumulative transmission potential. Turn negative to maximize
  return(tau.sum*-1) 
}


