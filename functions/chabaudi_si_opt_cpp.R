# The following function will produce the optimal conversion rate strategy of Plasmodium chabaudi in a single infection scenario. Users must specify:
## 1. Initial parameter (for conversion rate) search space
## 2. whether immunity is present or not ("ni" for no immunity and "si" for saturated immunity)
## 3. parameter setting
## 4. time period that infection run
## 5. degrees of freedom for conversion rate spline function
## 6. the "cue" that conversion rate is dependent on
## 7. The range of cue
## 8. Choice of solver used by dede. lsoda is chosen as default.

## To run the function, please see "plasticity_model_journals.Rmd" for reference". Function best run 
## when paired with optimParallel to allow for parallel computing.

# the following function is also optimized with Rcpp to increase computation speed

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
    
    ## Defining Pulse beta function based on current time
    pulseBeta <- pulseBeta_fun(parameters["I0"], parameters["sp"], t)
      
      ## Define the lag terms. lag[1] = R, lag[2] = I, lag[3] = Ig, lag[4] = M, lag[5] = G
      if(t>parameters["alpha"]){lag1 = deSolve::lagvalue(t-parameters["alpha"])} # lag state for asexual development
      if(t>parameters["alphag"]){lag2 = deSolve::lagvalue(t-parameters["alphag"])} # lag state for gametocyte development
      
      ## get lag term index given cue
    ### Only get lag index when it is a state-based cue
    if(cue != "t") {
      lag.i <- match(cue, names(state))}
    
    ### define lagged cue. Lag1 = alpha times ago, lag2 = alphag times ago
    if(t>parameters["alpha"] && cue == "t"){
      cue_lag1 <- t-parameters["alpha"]} 
    
    if(t>parameters["alpha"] && cue != "t"){
      cue_lag1 <- lag1[lag.i]} 
    
    if(t>parameters["alphag"] && cue == "t") {
      cue_lag2 <- t-parameters["alphag"]} 
    
    if(t>parameters["alphag"] && cue != "t") {
      cue_lag2 <- lag2[lag.i]}
      
      ## convert cue to variable in state or just time
    if(cue == "t"){
      cue_state <- t}
    else{
      cue_state <- state[cue]}
      
      ## Define K, carrying capacity of RBC
      K <- parameters["lambda"]*parameters["R1"]/(parameters["lambda"]-parameters["mu"]*parameters["R1"])
      
      ## Define survival functions
      ### Survival of infected asexual RBC
      if(t>parameters["alpha"] && immunity == "ni"){
        S <- exp(-parameters["mu"]*parameters["alpha"])
      } else if(t>parameters["alpha"] && immunity =="i"){
        integrand <- function(x) {parameters["mu"]+parameters["a"]/(parameters["b"]+state["I"])}
        integrate_val <- integrate(Vectorize(integrand), lower = t-parameters["alpha"], upper = t)
        S <- exp(-1*integrate_val$value)
      } else if(t<=parameters["alpha"] && immunity == "ni"){
        S <- exp(-parameters["mu"]*t)
      } else{
        integrand <- function(x) {parameters["mu"]+parameters["a"]/(parameters["b"]+state["I"])}
        integrate_val <- integrate(Vectorize(integrand), lower = 0, upper = t)
        S <- exp(-1*integrate_val$value)
      }
      
      ### Survival of gametocytes. We assume that infected
      ### RBC with gametocyte is not removed by immune response
      if(t<=parameters["alphag"]){
        Sg <- exp(-parameters["mu"]*t)
      } else{
        Sg <- exp(-parameters["mu"]*parameters["alphag"])}
      
      ## Define the models without lag terms. 
      dR <- parameters["lambda"]*(1-state["R"]/K)-parameters["mu"]*state["R"]-parameters["p"]*state["R"]*state["M"] # change in susceptible RBC
     
       if(immunity == "ni"){
        dI_nolag <- (1-cr(cue_state))*parameters["p"]*state["R"]*state["M"]-parameters["mu"]*state["I"] # change in infected RBC density
      } else{
        dI_nolag <- (1-cr(cue_state))*parameters["p"]*state["R"]*state["M"]-parameters["mu"]*state["I"]-(parameters["a"]*state["I"])/(parameters["b"]+state["I"]) # change in infected RBC density with immunity
      }
      
      dIg_nolag <- cr(cue_state)*parameters["p"]*state["R"]*state["M"]-parameters["mu"]*state["Ig"] 
      dM_nolag <- -parameters["mum"]*state["M"]-parameters["p"]*state["R"]*state["M"]
      dG_nolag <- -parameters["mug"]*state["G"]
      
      ## Track states in initial cohort of infection
      if(t<=parameters["alpha"]){
        dI <- dI_nolag-pulseBeta*S 
        dM <- dM_nolag+parameters["beta"]*pulseBeta*S 
      }
      
      if(t<=parameters["alphag"]){
        dIg <- dIg_nolag ## should have no cells form initial cohort
        dG <- 0
      }
      
      ## Track states after delay 
      if(t>parameters["alpha"]){
        dI <- dI_nolag-(1-cr(cue_lag1))*parameters["p"]*lag1[1]*lag1[4]*S 
        dM <- dM_nolag+parameters["beta"]*(1-cr(cue_lag1))*parameters["p"]*lag1[1]*lag1[4]*S 
      }
      
      if(t>parameters["alphag"]){
        dIg <- dIg_nolag-cr(cue_lag2)*parameters["p"]*lag2[1]*lag2[4]*Sg 
        dG <- dG_nolag+cr(cue_lag2)*parameters["p"]*lag2[1]*lag2[4]*Sg
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


