# Single strain infection model with Plasmodium chabaudi
# By: Avril Wang. Code adapted from Greischar et al., 2016 Predicting optimal transmission investment in malaria parasites

# The following function will produce the optimal conversion rate strategy of Plasmodium chabaudi in a single infection scenario. Users must specify:
## 1. Initial parameter (for conversion rate) search space
## 2. whether immunity is present or not ("ni" for no immunity and "si" for saturated immunity)
## 3. parameter setting
## 4. time period that infection run
## 5. degrees of freedom for conversion rate spline function
## 6. the "cue" that conversion rate is dependent on
## 7. The range of cue
## 8. Choice of solver used by dede. lsoda is chosen as default.

# Function best run when paired with optimParallel to allow for parallel computing.
# To run the function, use the following code:
## library(optimParallel)
## source("path to this file")
## parameters <- c()
## time_range <- seq(0, max time, by = 1e-3)
## cue_range <- seq(0, max cue, by = ...)
## cl <- makeCluster(detectCores()); setDefaultCluster(cl = cl)
## mod.opt <- optimParallel(par = c(...), 
                                  # fn = chabaudi_si_opt_cpp, 
                                  # control = list(trace = 6),
                                  # immunity = "i",
                                  # parameters = parameters,
                                  # time_range = time_range,
                                  # df = 3,
                                  # cue = "t",
                                  # cue_range = time_range)
## stopCluster(cl) 

chabaudi_si_opt_cpp <- function(parameters_cr, immunity, parameters, time_range, df, cue, cue_range, solver = "lsoda"){
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
  # Define initial condition
  #------------------------#
  state <- c(R = 8.5*10^6,
             I = 43.85965,
             Ig = 0,
             M = 0,
             G = 0)
  
  #-------------------------#
  # Ensure inputs are correct
  #------------------------#
  ## Ensure length of initial parameter search space matches with df
  if (length(parameters_cr) != df+1) {
    stop("Conversion rate parameters must match degrees of freedom")
  }
  ## Ensure immunity input is correct
  if (immunity != "ni" && immunity != "i") {
    stop("Immunity must be either 'ni' for no immunity or 'i' for saturating immunity")
  }
  ## Ensure cue is correct
  if (!(cue %in% names(state)) && !(cue %in% paste0("d", names(state))) && cue != "t") {
    stop("Cue must be one of the states, derivative of states, or time")
  }
  ## Ensure that time_range is used as cue_range when t is used
  if(cue == "t" && !isTRUE(all.equal(cue_range, time_range))){
    stop("Time is chosen as cue. Cue_range must equal to time_range")
  }
  
  #-------------------------#
  # Function to describe population 
  # structure of initial inoculum
  #------------------------#
  pulseBeta_fun <- function(I0, sp, t){ 
    res = rep(NA, length(t))
    res = I0*(dbeta(t, sp, sp))
  }
  
  #-------------------------#
  # Define conversion rate function. 
  # Simplified to increase performance. 
  #------------------------#
  ## Define dummy data of conversion rate
  dummy_y.vals <- rep(0, length(cue_range)) 
  dummy_cr.data <- as.data.frame(cbind(cue_range, dummy_y.vals))
  
  ## fit basic cubic spline with no internal knots. Here, increasing df increases knots.
  ## degree defines shape of basis spline. This shouldn't have any impact on ultimate result
  ## What matter is the number of knots, or df. Increasing df should give us more complex reaction norm
  dummy_cr.mod <- lm(dummy_y.vals ~ splines2::bSpline(cue_range, degree = 3, df = df))
  dummy_cr.mod$data <- dummy_cr.data
  
  ## Assign coefficient to be optimized to the dummy conversion rate function
  dummy_cr.mod$coefficients <- parameters_cr
  
  ## use spline function to predict cr 
  cr_fit <- exp(-exp(predict(dummy_cr.mod, newdata = data.frame(cue_range))))
  
  ## Get spline function where cr ~ cue
  cr <- splinefun(cbind(cue_range, cr_fit))
  
  #-------------------------#
  # Define trapezoidal estimation 
  # for survival function approximation. Use this if speeding up model
  #------------------------#
  trapezoid <- function(f, lower, upper) {
    if (is.function(f) == FALSE) {
      stop('f must be a function with one parameter (variable)')
    }
    h <- upper - lower
    fxdx <- (h / 2) * (f(lower) + f(upper))
    return(fxdx)}
  
  #-------------------------#
  # Define single-infection model
  #------------------------#
  single_infection.fun <- function(t, state, parameters) {
    
    ## Rename parameters for cleaner code. With.list not used to speed up computation
    R1 <- parameters["R1"]
    lambda <- parameters["lambda"]
    mu <- parameters["mu"]
    p <- parameters["p"]
    alpha <- parameters["alpha"]
    alphag <- parameters["alphag"]
    beta <- parameters["beta"]
    mum <- parameters["mum"]
    mug <- parameters["mug"]
    I0 <- parameters["I0"]
    Ig0 <- parameters["Ig0"]
    a <- parameters["a"]
    b <- parameters["b"]
    sp <- parameters["sp"]
    
    # rename states for cleaner code
    R <- state["R"]
    I <- state["I"]
    Ig <- state["Ig"]
    M <- state["M"]
    G <- state["G"]
    
     ## Defining Pulse beta function based on current time
     pulseBeta <- pulseBeta_fun(I0, sp, t)
      
     ## Define the lag terms. lag[1] = R, lag[2] = I, lag[3] = Ig, lag[4] = M, lag[5] = G
     if(t>alpha){lag1 = deSolve::lagvalue(t-alpha)} # lag state for asexual development
     if(t>alphag){lag2 = deSolve::lagvalue(t-alphag)} # lag state for gametocyte development
      
     ## get lag term index given cue
        ### Only get lag index when it is a state-based cue
        if(cue != "t") {
          lag.i <- match(cue, names(state))}
        
        ### define lagged cue. Lag1 = alpha times ago, lag2 = alphag times ago
        if(t>alpha && cue == "t"){
          cue_lag1 <- t-alpha} 
        
        if(t>alpha && cue != "t"){
          cue_lag1 <- lag1[lag.i]} 
        
        if(t>alphag && cue == "t") {
          cue_lag2 <- t-alphag} 
        
        if(t>alphag && cue != "t") {
          cue_lag2 <- lag2[lag.i]}
      
    ## convert cue to variable in state or just time
    if(cue == "t"){
      cue_state <- t}
    else{
      cue_state <- state[cue]}
      
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
        dI <- dI_nolag-pulseBeta*S 
        dM <- dM_nolag+beta*pulseBeta*S 
      }
      
      if(t<=alphag){
        dIg <- dIg_nolag ## should have no cells form initial cohort
        dG <- 0
      }
      
      ## Track states after delay 
      if(t>alpha){
        dI <- dI_nolag-(1-cr(cue_lag1))*p*lag1[1]*lag1[4]*S 
        dM <- dM_nolag+beta*(1-cr(cue_lag1))*p*lag1[1]*lag1[4]*S 
      }
      
      if(t>alphag){
        dIg <- dIg_nolag-cr(cue_lag2)*p*lag2[1]*lag2[4]*Sg 
        dG <- dG_nolag+cr(cue_lag2)*p*lag2[1]*lag2[4]*Sg
      }
      
      ## Return the states
      return(list(c(dR, dI, dIg, dM, dG)))
  }
  
  #-------------------------#
  # Run single-infection model
  #------------------------#
  single_infection.df <- as.data.frame(deSolve::dede(y = state,
                                                     times = time_range,
                                                     func = single_infection.fun,
                                                     p = parameters,
                                                     method = solver,
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


