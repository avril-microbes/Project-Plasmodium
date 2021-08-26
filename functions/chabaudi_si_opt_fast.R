# Single strain infection model with Plasmodium chabaudi
# By: Avril Wang. Code adapted from Greischar et al., 2016 Predicting optimal transmission investment in malaria parasites

#-------------------------#
# Getting optimal conversion 
# rate strategy
#------------------------#
# The following function will produce the optimal conversion rate strategy of Plasmodium chabaudi in a single infection scenario. Users must specify:
## 1. Initial parameter (for conversion rate) search space
## 2. whether immunity is present or not ("ni" for no immunity and "si" for saturated immunity)
## 3. parameter setting
## 4. time period that infection run
## 5. degrees of freedom for conversion rate spline function
## 6. the "cue" that conversion rate is dependent on
## 7. The range of cue
## 8. Choice of solver used by dede. lsoda is chosen as default.
## 9. Whether infection dynamics is simulated. If set as TRUE, rather than performing optimization, the function will
## perform dede on the model using the parameter provided

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
                                  # immunity = "...", ## "i" or "ni"
                                  # parameters = parameters, # list of time points. Use seq(lower time, upper time, by = time interval)
                                  # time_range = time_range,
                                  # df = ..., # give a number
                                  # cue = "...", # see state for choice of cue
                                  # cue_range = time_range) # list of cue points. Use seq(lower cue, upper cue, by = cue interval)
## stopCluster(cl) 

#-------------------------#
# Simulating infection dynamics
# Given a set of conversion rate
# parameters
#------------------------#
# Function can also be used to simulate infection dynamics (track states with time) if dyn is set to TRUE
# To perform infection dynamics simulation, using the following script
## library(ggplot)
## source("path to this file")
## parameters <- c()
## time_range <- seq(0, max time, by = 1e-3)
## cue_range <- seq(0, max cue, by = ...)
## mod.dyn <- chabaudi_si_opt_fast(parameters_cr = ...,
                                # immunity = "...",
                                # parameters = ...,
                                # time_range = ...,
                                # df = ...,
                                # cue = "...",
                                # cue_range = ...,
                                # solver = "...",
                                # dyn = TRUE)
## ggplot(mod.dyn, aes(x = time, y = value)) + # plot infection dynamics 
  # geom_line() +
  # facet_wrap(~variable, scales = "free") +
  # scale_y_continuous(labels = scales::scientific) +
  # theme_bw()                  

chabaudi_si_opt_fast <- function(parameters_cr, immunity, parameters, time_range, df, cue, cue_range, solver = "lsoda", integration = "integrate", dyn = FALSE) {
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
  force(solver)
  force(integration)
  
  #-------------------------#
  # Define initial condition
  #------------------------#
  if(immunity == "ni" || immunity == "i"){
  state <- c(R = 8.5*10^6,
             I = 43.85965,
             Ig = 0,
             M = 0,
             G = 0)
  }
  
  ## Input additional state parameters if Kochin's 
  ## innate immunity model is chosen
  if(immunity == "kochin"){
    state <- c(R = 8.5*10^6,
               I = 43.85965,
               Ig = 0,
               M = 0,
               G = 0,
               E = 0)} # targeted RBC removal (proportion activated)
  
  ## Input additional state parameters for Tsukushi's model
  if(immunity == "tsukushi"){
    state <- c(R = 8.5*10^6,
               I = 43.85965,
               Ig = 0,
               M = 0,
               G = 0,
               N = 0, # general RBC removal
               W = 0) # targeted RBC removal
  }
  
  #-------------------------#
  # Ensure inputs are correct
  #------------------------#
  ## Ensure length of initial parameter search space matches with df
  if (length(parameters_cr) != df+1) {
    stop("Conversion rate parameters must match degrees of freedom")
  }
  ## Ensure immunity input is correct
  if (immunity != "ni" && immunity != "i" && immunity != "kochin" && immunity != "tsukushi") {
    stop("Immunity must be either 'ni', 'i', 'kochin,' or 'tsukushi'")
  }
  ## Ensure cue is correct
  if (!(cue %in% names(state)) && !(cue %in% paste0("d", names(state))) && cue != "t") {
    stop("Cue must be one of the states, derivative of states, or time")
  }
  ## Ensure that time_range is used as cue_range when t is used
  if(cue == "t" && !isTRUE(all.equal(cue_range, time_range))){
    stop("Time is chosen as cue. Cue_range must equal to time_range")
  }
  ## Ensure integration is entered correctly
  if(integration != "integrate" && integration != "trapezoid" && integration != "simpson"){
    stop("Please enter the correct integration method. Must be 'integrate', trapezoid', or 'simpson'")
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
  dummy_cr.mod <- lm(dummy_y.vals ~ splines2::bSpline(x = cue_range, degree = 3, df = df))
  dummy_cr.mod$data <- dummy_cr.data
  
  ## Assign coefficient to be optimized to the dummy conversion rate function
  dummy_cr.mod$coefficients <- parameters_cr
  
  ## use spline function to predict cr 
  cr_fit <- exp(-exp(predict(dummy_cr.mod, newdata = data.frame(cue_range))))
  
  ## Get spline function where cr ~ cue
  cr <- splinefun(cbind(cue_range, cr_fit))
  
  #-------------------------#
  # Define integration method. Default to integrate. 
  #------------------------#
  if(integration == "integrate"){
    integrate_fun <- stats::integrate
    } else if (integration == "trapezoid"){
    integrate_fun <- function(f, lower, upper) {
      if (is.function(f) == FALSE) {
        stop('f must be a function with one parameter (variable)')}
      h <- upper - lower
      fxdx <- (h / 2) * (f(lower) + f(upper))
      return(fxdx)}
    } else {
    integrate_fun <- function(f, lower, upper) {
      if (is.function(f) == FALSE) {
        stop('f must be a function with one parameter (variable)')}
      h <- (upper - lower) / 2
      x0 <- lower
      x1 <- lower + h
      x2 <- upper
      s <- (h / 3) * (f(x0) + 4 * f(x1) + f(x2))
      return(s)
    }
  }
  
  #-------------------------#
  # Define single-infection model
  #------------------------#
  chabaudi_si_model <- function(t, state, parameters) {
    
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
    
    ## Additional parameters in Kochin
    if (immunity == "kochin"){ 
      sigma <- parameters["sigma"] # Probability of activating immune cell upon contact
      mue <- parameters["mue"] # inactivation of immune cells
      gamma <- parameters["gamma"] # Maximum removal rate of iRBC
    }
    
    ## Additional parameters of Tsukushi model
    if (immunity == "tsukushi") {
      psin <- parameters["psin"] # activation strength for general RBC removal
      psiw <- parameters["psiw"] # activation strength for targeted RBC removal
      phin <- parameters["phin"] # half life for general RBC removal
      phiw <- parameters["phiw"] # halflife for targeted RBC removal
      iota <- parameters["iota"] # cue strength (infected iRBC)
    }
    
    # rename states for cleaner code
    R <- state["R"]
    I <- state["I"]
    Ig <- state["Ig"]
    M <- state["M"]
    G <- state["G"]
    if (immunity == "kochin") {E <- state["E"]}
    if (immunity == "tsukushi"){
      N <- state["N"]
      W <- state["W"]
    }
    
    ## Defining Pulse beta function based on current time
    pulseBeta <- pulseBeta_fun(I0, sp, t)
    
    ## Define the lag terms. lag[1] = R, lag[2] = I, lag[3] = Ig, lag[4] = M, lag[5] = G
    if(t>alpha){lag1 = deSolve::lagvalue(t-alpha)} # lag state for asexual development
    if(t>alphag){lag2 = deSolve::lagvalue(t-alphag)} # lag state for gametocyte development
    
    ## get lag term index given cue. Cannot use else if given that during simulation, multiple iterations of cue_lag is used
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
    
    if(cue != "t"){ cue_state <- state[cue]}
    
    ## Define K, carrying capacity of RBC
    K <- lambda*R1/(lambda-mu*R1)
    
    ## Define survival functions
    ### Survival of infected asexual RBC
    if(t>alpha && immunity == "ni"){
      S <- exp(-mu*alpha)} 
    
    if(t>alpha && immunity =="i"){
      integrand <- function(x) {mu+a/(b+I)}
      integrate_val <- integrate_fun(Vectorize(integrand), lower = t-alpha, upper = t)
      if(integration == "integrate"){integrate_val <- integrate_val$value}
      S <- exp(-1*integrate_val)
    }  
    
    if(t>alpha && immunity == "kochin"){
      integrand <- function(x) {mu+gamma*E}
      integrate_val <- integrate_fun(Vectorize(integrand), lower = t-alpha, upper = t)
      if(integration == "integrate"){integrate_val <- integrate_val$value}
      S <- exp(-1*integrate_val)} 
    
    if(t>alpha && immunity == "tsukushi"){
      integrand <- function(x) {mu-log(1-N)-log(1-W)}
      integrate_val <- integrate_fun(Vectorize(integrand), lower = t-alpha, upper = t)
      if(integration == "integrate"){integrate_val <- integrate_val$value}
      S <- exp(-1*integrate_val)} 
    
    ################################
    
    if(t<=alpha && immunity == "ni"){
      S <- exp(-mu*t)} 
    
    if(t<=alpha && immunity == "i"){
      integrand <- function(x) {mu+a/(b+I)}
      integrate_val <- integrate_fun(Vectorize(integrand), lower = 0, upper = t)
      if(integration == "integrate"){integrate_val <- integrate_val$value}
      S <- exp(-1*integrate_val)
    }
    
    if(t<=alpha && immunity == "kochin"){
      integrand <- function(x) {mu+gamma*E}
      integrate_val <- integrate_fun(Vectorize(integrand), lower = 0, upper = t)
      if(integration == "integrate"){integrate_val <- integrate_val$value}
      S <- exp(-1*integrate_val)}
    
    if(t<=alpha && immunity == "tsukushi"){
      integrand <- function(x) {mu-log(1-N)-log(1-W)}
      integrate_val <- integrate_fun(Vectorize(integrand), lower = 0, upper = t)
      if(integration == "integrate"){integrate_val <- integrate_val$value}
      S <- exp(-1*integrate_val)}
    
    ### Survival of gametocytes. We assume that infected
    ### RBC with gametocyte is not removed by immune response
    if(t<=alphag){
      Sg <- exp(-mu*t)
    } 
    
    if(t>alphag){
      Sg <- exp(-mu*alphag)}
    
    ## Define the models without lag terms. 
    if(immunity != "tsukushi"){
      dR <- lambda*(1-R/K)-mu*R-p*R*M # change in susceptible RBC
    } 
    
    if(immunity == "tsukushi"){ #Tsukushi exclusive ODEs
      dR <- lambda*(1-R/K)-(mu-log(1-N))*R-p*R*M
      dI_nolag <- (1-cr(cue_state))*p*R*M-mu*I-(-log(1-N)-log(1-W))*I
      dN <- psin*iota*(1-N)-N/phin
      dW <- psiw*iota*(1-W)-W/phiw
    }
    
    if(immunity =="kochin"){
      dE <- sigma*I*(1-E)-mue*E # change in innate immune strength
      dI_nolag <- (1-cr(cue_state))*p*R*M-mu*I-gamma*E*I
    }
    
    if(immunity == "ni"){
      dI_nolag <- (1-cr(cue_state))*p*R*M-mu*I # change in infected RBC density
    }
    
    if(immunity == "i") {
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
    if (immunity == "ni" || immunity == "i") {return(list(c(dR, dI, dIg, dM, dG)))}
    
    if (immunity == "kochin") {return(list(c(dR, dI, dIg, dM, dG, dE)))}
    
    if (immunity == "tsukushi") {return(list(c(dR, dI, dIg, dM, dG, dN, dW)))}
  }
  
  #-------------------------#
  # Run single-infection model
  #------------------------#
  chabaudi_si.df <- as.data.frame(deSolve::dede(y = state,
                                                     times = time_range,
                                                     func = chabaudi_si_model,
                                                     p = parameters,
                                                     method = solver,
                                                     control=list(mxhist = 1e6)))
  
  #-------------------------#
  # Calculate fitness
  #------------------------#
  ## Get Gametocyte density time series data
  gam <- chabaudi_si.df$G
  gam[gam<0] <- 0 # Assign negative gametocyte density to 0
  
  ## Get timeseries interval. Simplify first time after t=0
  int <- chabaudi_si.df[2,1]
  
  ## Define the fitness parameter values
  aval <- -12.69
  bval <- 3.6
  dens <- log10(gam)
  
  ## Calculate the transmission potential at each time t
  tau.ls <- (exp(aval+bval*dens))/(1+exp(aval+bval*dens))
  
  ## Get approximation of cumulative transmission potential
  tau.sum <- sum(tau.ls*int)
  
  # return cumulative transmission potential. Turn negative to maximize
  if(dyn == FALSE){return(tau.sum*-1)}
  
  #-------------------------#
  # Simulating infection dynamics if Dyn == TRUE
  #------------------------# 
  if(dyn == TRUE) {
    ### calculate cumulative transmission potential gain
    tau_cum.ls <- cumsum(tau.ls*int)
    
    ### cbind results
    chabaudi_si.df$tau <- tau.ls
    chabaudi_si.df$tau_cum <- tau_cum.ls
    
    ### calculate CR based on cue
    if(cue != "t"){cr.ls <- cr(chabaudi_si.df[, cue])}
    if(cue == "t"){cr.ls <- cr(time_range)}
    chabaudi_si.df$cr <- cr.ls

    ### processing df for plotting
    chabaudi_si.df2 <- chabaudi_si.df %>% tidyr::gather(key = "variable", value = "value", -time)
    
    return(chabaudi_si.df2)
  }
}


