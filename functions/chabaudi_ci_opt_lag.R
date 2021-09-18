# co-infection strain infection model with Plasmodium chabaudi
# By: Avril Wang. Code adapted from Greischar et al., 2016 Predicting optimal transmission investment in malaria parasites
# the following script is altered to reflect next cycle conversion, where sexual commitment occurs in the previous
# transmission cycle. Upon receiving a cue, infected RBC does not decide, on site, whether to produce merozoite
# or gametocyte. Instead, all iRBC will produce either sexually committed merozoite (Mg) or asexual merozoite after
# 1 day of development time. Subsequent RBC that is infected by Mg or M then becomes sexually committed iRBC (Ig) or
# asexual iRBC (I). I assumed that all injected iRBC are asexual, hence, the earliest date at which gametocytes
# are produced is day 3.

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

chabaudi_ci_opt_lag <- function(parameters_cr_1,
                                parameters_cr_2,
                                immunity,
                                parameters,
                                time_range,
                                df,
                                cue,
                                cue_range,
                                solver = "lsoda",
                                #integration = "integrate",
                                transformation = "exp",
                                adaptive = FALSE,
                                dyn = FALSE,
                                log_cue = FALSE,
                                gamete_immune = TRUE) {
  #-------------------------#
  # Ensure values we inputted 
  # are available in environment
  #------------------------#
  force(parameters_cr_1)
  force(parameters_cr_2)
  force(immunity)
  force(parameters)
  force(time_range)
  force(df)
  force(cue)
  force(cue_range)
  force(solver)
  #force(integration)
  force(adaptive)
  
  #-------------------------#
  # Define initial condition
  #------------------------#
  if(immunity == "ni" || immunity == "i"){
    state <- c(R = parameters[["R1"]],
               M1 = 0,
               M2 = 0,
               Mg1 = 0,
               Mg2 = 0,
               ID = 0,
               I1 = parameters[["I0"]]/2,
               I2 = parameters[["I0"]]/2,
               Ig1 = 0,
               Ig2 = 0,
               G1 = 0, 
               G2 = 0,
               A = 0) # survival function. Just for tracking
  } else if(immunity == "kochin"){
    state <- c(R = parameters[["R1"]],
               M1 = 0,
               M2 = 0,
               Mg1 = 0,
               Mg2 = 0,
               ID = 0,
               I1 = parameters[["I0"]]/2,
               I2 = parameters[["I0"]]/2,
               Ig1 = 0,
               Ig2 = 0,
               G1 = 0,
               G2 = 0,
               E = 0,
               A = 0)
  } else{
    state <- c(R = parameters[["R1"]],
               M1 = 0,
               M2 = 0,
               Mg1 = 0,
               Mg2 = 0,
               ID = 0,
               I1 = parameters[["I0"]]/2,
               I2 = parameters[["I0"]]/2,
               Ig1 = 0,
               Ig2 = 0,
               G1 = 0,
               G2 = 0,
               N = 0, # general RBC removal
               W = 0,
               A = 0) # targeted RBC removal
  }
  
  #-------------------------#
  # Ensure inputs are correct
  #------------------------#
  ## Ensure length of initial parameter search space matches with df. 
  if (length(parameters_cr_1) != df+1 | length(parameters_cr_2) != df+1) {
    stop("Conversion rate parameters must match degrees of freedom")
  }
  ## Ensure immunity input is correct
  if (immunity != "ni" && immunity != "i" && immunity != "kochin" && immunity != "tsukushi") {
    stop("Immunity must be either 'ni', 'i', 'kochin,' or 'tsukushi'")
  }
  ## Ensure cue is correct
  ## Ensure cue is correct
  if (!(unlist(stringr::str_split(cue, "\\+|\\-|\\*|\\/")) %in% names(state)) && 
      !(cue %in% paste0("d", names(state))) && cue != "t") {
    stop("Cue must be one of the states, derivative of states, or time")
  }
  ## Ensure that time_range is used as cue_range when t is used
  if(cue == "t" && !isTRUE(all.equal(cue_range, time_range))){
    stop("Time is chosen as cue. Cue_range must equal to time_range")
  }
  ## Ensure integration is entered correctly. Deprecated
  #if(integration != "integrate" && integration != "trapezoid" && integration != "simpson"){
  #  stop("Please enter the correct integration method. Must be 'integrate', trapezoid', or 'simpson'")
  #}
  ## Ensure spline transformation is entered correct
  if(transformation != "norm" && transformation != "exp" && transformation != "logit"){
    stop("Transformation must be either 'norm' or 'exp' or 'logit'")
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
  dummy_cr.mod_1 <- lm(dummy_y.vals ~ splines2::bSpline(x = cue_range, degree = 3, df = df))
  dummy_cr.mod_1$data <- dummy_cr.data
  
  dummy_cr.mod_2 <- lm(dummy_y.vals ~ splines2::bSpline(x = cue_range, degree = 3, df = df))
  dummy_cr.mod_2$data <- dummy_cr.data
  
  ## Assign coefficient to be optimized to the dummy conversion rate function
  dummy_cr.mod_1$coefficients <- parameters_cr_1
  dummy_cr.mod_2$coefficients <- parameters_cr_2
  
  ## use spline function to predict cr 
  if(transformation == "norm"){
    cr_fit_1 <- predict(dummy_cr.mod_1, newdata = data.frame(cue_range))
    cr_fit_t_1 <- (cr_fit_1-min(cr_fit_1))/(max(cr_fit_1-min(cr_fit_1)))
    
    cr_fit_2 <- predict(dummy_cr.mod_2, newdata = data.frame(cue_range))
    cr_fit_t_2 <- (cr_fit_2-min(cr_fit_2))/(max(cr_fit_2-min(cr_fit_2)))
    
  } else if(transformation == "exp"){
    cr_fit_t_1 <- exp(-exp(predict(dummy_cr.mod_1, newdata = data.frame(cue_range))))
    cr_fit_t_2 <- exp(-exp(predict(dummy_cr.mod_2, newdata = data.frame(cue_range))))
  } else{
    cr_fit_1 <- predict(dummy_cr.mod_1, newdata = data.frame(cue_range))
    cr_fit_t_1 <- 1 / (1 + exp(-cr_fit_1))
    
    cr_fit_2 <- predict(dummy_cr.mod_2, newdata = data.frame(cue_range))
    cr_fit_t_2 <- 1 / (1 + exp(-cr_fit_2))
  }
  
  ## Get spline function where cr ~ cue
  cr_1 <- splinefun(cbind(cue_range, cr_fit_t_1))
  cr_2 <- splinefun(cbind(cue_range, cr_fit_t_2))
  
  #-------------------------#
  # Define integration method. Deprecated
  #------------------------#
  #if(integration == "integrate"){
  #  integrate_fun <- stats::integrate
  #} else if (integration == "trapezoid"){
  #  integrate_fun <- function(f, lower, upper) {
  #    if (is.function(f) == FALSE) {
  #      stop('f must be a function with one parameter (variable)')}
  #    h <- upper - lower
  #    fxdx <- (h / 2) * (f(lower) + f(upper))
  #    return(fxdx)}
  #} else {
  #  integrate_fun <- function(f, lower, upper) {
  #    if (is.function(f) == FALSE) {
  #      stop('f must be a function with one parameter (variable)')}
  #    h <- (upper - lower) / 2
  #    x0 <- lower
  #    x1 <- lower + h
  #    x2 <- upper
  #    s <- (h / 3) * (f(x0) + 4 * f(x1) + f(x2))
  #    return(s)
  #  }
  #}
  
  #-------------------------#
  # Define single-infection model
  #------------------------#
  chabaudi_ci_model_lag <- function(t, state, parameters) {
    
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
      rho <- parameters["rho"] #  proportion of the deviation from the homeostatic equilibrium restored by the host per day 
    }
    
    ## Additional parameters if adaptive immunity is incorporated
    if(adaptive == TRUE){
      psia <- parameters["psia"] # activation strength for adaptive immunity
      epsilon <- parameters["epsilon"] # time lag between immune response and parasetemia
      phia <- parameters["phia"] # decay rate for adaptive immunity
      theta <- parameters["theta"] # developement time for adaptive immunity
    }
    
    # rename states for cleaner code
    R <- state["R"]
    I1 <- state["I1"]
    I2 <- state["I2"]
    Ig1 <- state["Ig1"]
    Ig2 <- state["Ig2"]
    ID <- state["ID"]
    M1 <- state["M1"]
    M2 <- state["M2"]
    G1 <- state["G1"]
    G2 <- state["G2"]
    A <- state["A"]
    Mg1 <- state["Mg1"]
    Mg2 <- state["Mg2"]
    if (immunity == "kochin") {E <- state["E"]}
    if (immunity == "tsukushi"){
      N <- state["N"]
      W <- state["W"]
    }
    
    ## Defining Pulse beta function based on current time
    pulseBeta <- pulseBeta_fun(I0/2, sp, t)
    
    ## Define the lag terms. lag[1] = R, lag[2] = I, lag[3] = Ig, lag[4] = M, lag[5] = G
    if(t>alpha){lag1 = deSolve::lagvalue(t-alpha)} # lag state for asexual development
    if(t>alphag){lag2 = deSolve::lagvalue(t-alphag)} # lag state for gametocyte development
    ### extra lag term for adaptive immunity
    if(adaptive){
      if(t>epsilon){lag3 = deSolve::lagvalue(t-epsilon)}
    }
    
    ## get lag term index given cue. Cannot use else if given that during simulation, multiple iterations of cue_lag is used
    ### Only get lag index when it is a state-based cue. Multiple indexes are returned
    ### if multiple cues are given, multiple indexes are returned
    if(cue != "t") {
      lag.i <- match(unlist(stringr::str_split(cue, "\\+|\\-|\\*|\\/")), names(state))
    }
    
    ### define lagged cue. Lag1 = alpha times ago, lag2 = alphag times ago
    ### For simple cues (if it does not contain special characters)
    if(stringr::str_detect(cue, "\\+|\\-|\\*|\\/", negate = TRUE)){
      if(t>alpha && cue == "t"){
        cue_lag1 <- t-alpha} 
      
      if(t>alpha && cue != "t"){
        cue_lag1 <- lag1[lag.i]}
      
      if(t>alphag && cue == "t") {
        cue_lag2 <- t-alphag} 
      
      if(t>alphag && cue != "t") {
        cue_lag2 <- lag2[lag.i]}
      
      ### convert cue to time if time-based conversion rate strategy is used
      if(cue == "t"){
        cue_state <- t}
      
      ### get cue_state if it is state-based
      if(cue != "t"){
        cue_state <- state[cue]}
    } else{### manually create lag values if cues contain special characters
      if(stringr::str_detect(cue, "\\+")){ # if it contains plus
        if(t>alpha && cue != "t"){cue_lag1 <- lag1[lag.i[1]]+lag1[lag.i[2]]}
        if(t>alphag && cue != "t") {cue_lag2 <- lag2[lag.i[1]]+lag2[lag.i[2]]}
      }
      if(stringr::str_detect(cue, "\\-")){ # if it contains -
        if(t>alpha && cue != "t"){cue_lag1 <- lag1[lag.i[1]]-lag1[lag.i[2]]}
        if(t>alphag && cue != "t") {cue_lag2 <- lag2[lag.i[1]]-lag2[lag.i[2]]}
      }
      if(stringr::str_detect(cue, "\\*")){ # if it contains multiplication
        if(t>alpha && cue != "t"){cue_lag1 <- lag1[lag.i[1]]*lag1[lag.i[2]]}
        if(t>alphag && cue != "t") {cue_lag2 <- lag2[lag.i[1]]*lag2[lag.i[2]]}
      }
      if(stringr::str_detect(cue, "\\/")){ # if it contains division
        if(t>alpha && cue != "t"){cue_lag1 <- lag1[lag.i[1]]/lag1[lag.i[2]]}
        if(t>alphag && cue != "t") {cue_lag2 <- lag2[lag.i[1]]/lag2[lag.i[2]]}
      }
      ### get present states
      if(cue != "t"){cue_state <- eval(parse(text = cue))}
    }
    
    ## Define K, carrying capacity of RBC
    K <- lambda*R1/(lambda-mu*R1)
    
    ## Define adaptive immunity
    if(adaptive == FALSE){
      dA <- 0}
    
    if(adaptive == TRUE){
      if(t<theta){dA <- 0}
      if(t>=theta){
        dA <- psia*((lag3[2]+lag3[3])/iota)*(0.85-A)-(A/phia)}
      #dA <- psia*((lag3[2]+lag3[3])/iota)*(0.85-A)} # assume no decay
    }
    
    ## Define survival functions
    ### Survival of infected asexual RBC
    if(t>alpha && immunity == "ni"){
      S <- exp(-mu*alpha)} 
    
    if(t>alpha && immunity != "ni"){
      S <- exp(-ID + lag1[6])
    }
    
    #if(t>alpha && immunity =="i"){
    #integrand <- function(x) {mu+a/(b+I)}
    #integrate_val <- integrate_fun(Vectorize(integrand), lower = t-alpha, upper = t)
    #if(integration == "integrate"){integrate_val <- integrate_val$value}
    #S <- exp(-1*integrate_val)
    # S <- exp(-ID + lag1[4])
    #}  
    
    #if(t>alpha && immunity == "kochin"){
    #  integrand <- function(x) {mu+gamma*E}
    #  integrate_val <- integrate_fun(Vectorize(integrand), lower = t-alpha, upper = t)
    #  if(integration == "integrate"){integrate_val <- integrate_val$value}
    #  S <- exp(-1*integrate_val)} 
    
    # if(t>alpha && immunity == "tsukushi"){
    #  integrand <- function(x) {mu-log(1-N)-log(1-W)-log(1-A)}
    # integrate_val <- integrate_fun(Vectorize(integrand), lower = t-alpha, upper = t)
    #if(integration == "integrate"){integrate_val <- integrate_val$value}
    #S <- exp(-1*integrate_val)} 
    
    ################################
    
    if(t<=alpha && immunity == "ni"){
      S <- exp(-mu*t)} 
    
    if(t<=alpha && immunity != "ni"){
      S <- exp(-ID)} 
    
    #if(t<=alpha && immunity == "i"){
    #integrand <- function(x) {mu+a/(b+I)}
    #integrate_val <- integrate_fun(Vectorize(integrand), lower = 0, upper = t)
    #if(integration == "integrate"){integrate_val <- integrate_val$value}
    #S <- exp(-1*integrate_val)
    # S <- exp(-ID)
    #}
    
    # if(t<=alpha && immunity == "kochin"){
    # integrand <- function(x) {mu+gamma*E}
    #  integrate_val <- integrate_fun(Vectorize(integrand), lower = 0, upper = t)
    #  if(integration == "integrate"){integrate_val <- integrate_val$value}
    #  S <- exp(-1*integrate_val)}
    
    # if(t<=alpha && immunity == "tsukushi"){
    # integrand <- function(x) {mu-log(1-N)-log(1-W)-log(1-A)} # assume targetted removal is half as effective
    #  integrate_val <- integrate_fun(Vectorize(integrand), lower = 0, upper = t)
    # if(integration == "integrate"){integrate_val <- integrate_val$value}
    # S <- exp(-1*integrate_val)}
    
    ### Survival of gametocytes. We assume that infected
    ### RBC with gametocyte is removed by immune response for Tsukushi's model.
    if(immunity != "tsukushi"){
      if(t<=alpha){
        Sg <- 0 # not relevent
      } 
      
      if(t>alpha && t<=alpha+alphag){
        Sg <- exp(-mu*t+mu*alpha)} # not relevent
      
      if(t>alpha+alphag){
        Sg <- exp(-mu*alphag)
      }
    }
    
    if(immunity == "tsukushi"){
      if(t<=alpha){
        Sg <- 0 # not relevent
      }
      
      if(t>alpha && t<=alpha+alphag){
        Sg <- 0 # not relevent
      }
      
      if(t>alpha+alphag){
        Sg <- exp(-ID + lag2[6])
      }
      
      #if(t>alpha && t<=alpha+alphag){
      # #integrand <- function(x) {mu-log(1-N)}
      #  integrand <- function(x) {mu-log(1-N)-log(1-W)-log(1-A)} 
      # integrate_val <- integrate_fun(Vectorize(integrand), lower = alpha, upper = t)
      #  if(integration == "integrate"){integrate_val <- integrate_val$value}
      # Sg <- exp(-1*integrate_val)
      #} 
      
      # if(t>alpha+alphag){
      #  #integrand <- function(x) {mu-log(1-N)}
      # integrand <- function(x) {mu-log(1-N)-log(1-W)-log(1-A)}
      #  integrate_val <- integrate_fun(Vectorize(integrand), lower = t-alphag, upper = t)
      #  if(integration == "integrate"){integrate_val <- integrate_val$value}
      #  Sg <- exp(-1*integrate_val)}
    }
    
    ## Define the models without lag terms. 
    if(immunity != "tsukushi"){
      dR <- lambda*(1-(R/K))-(mu*R)-(p*R*M1)-(p*R*M2)-(p*R*Mg1)-(p*R*Mg2) # change in susceptible RBC
    } 
    
    if(immunity == "tsukushi"){ #Tsukushi exclusive ODEs
      #dR <- lambda*(1-R/K)-mu*R-p*R*M-(mu-log(1-N))*R
      dR <- R1*mu+rho*(R1-R)-(mu-log(1-N))*R-(p*R*M1)-(p*R*M2)-(p*R*Mg1)-(p*R*Mg2)
      dI1_nolag <- p*R*M1-mu*I1-(-log(1-N)-log(1-W)-log(1-A))*I1
      dI2_nolag <- p*R*M2-mu*I2-(-log(1-N)-log(1-W)-log(1-A))*I2
      #dIg_nolag <- cr(cue_state)*p*R*M-mu*Ig-(-log(1-N))*Ig #assume no targeted clearance
      if(gamete_immune){ # assuming if gametocyte triggers and is removed by targeted immunity
        dIg1_nolag <- p*R*Mg1-mu*Ig1-(-log(1-N)-log(1-W)-log(1-A))*Ig1
        dIg2_nolag <- p*R*Mg2-mu*Ig2-(-log(1-N)-log(1-W)-log(1-A))*Ig2
        dN <- psin*((I1+I2+Ig1+Ig2)/iota)*(1-N)-(N/phin)
        dW <- psiw*((I1+I2+Ig1+Ig2)/iota)*(1-W)-(W/phiw)}
      else { # assuming if gametocyte does not trigger and is not removed by targeted immunity
        dIg1_nolag <- p*R*Mg1-mu*Ig1-(-log(1-N)-log(1-A))*Ig1
        dIg2_nolag <- p*R*Mg2-mu*Ig2-(-log(1-N)-log(1-A))*Ig2
        dN <- psin*((I1+I2)/iota)*(1-N)-(N/phin)
        dW <- psiw*((I1+I2)/iota)*(1-W)-(W/phiw)}
      
      #dN <- psin*(I/iota)*(1-N)-(N/phin) # assume Ig does not elicit strong immune response. Not included in cue
      
      #dW <- psiw*(I/iota)*(1-W)-(W/phiw)
      
      dM1_nolag <- (-mum*M1)-(p*R*M1)
      dM2_nolag <- (-mum*M2)-(p*R*M2)
      dMg1_nolag <- (-mum*Mg1)-(p*R*Mg1)
      dMg2_nolag <- (-mum*Mg2)-(p*R*Mg2)
      dG1_nolag <- -mug*G1
      dG2_nolag <- -mug*G2
      dID <- mu-log(1-N)-log(1-W)-log(1-A)
    }
    
    if(immunity =="kochin"){
      if(gamete_immunity){
        dE <- sigma*(I1+I2+Ig1+Ig2)*(1-E)-mue*E
        dIg1_nolag <- p*R*Mg1-mu*Ig1-gamma*E*Ig1
        dIg2_nolag <- p*R*Mg2-mu*Ig2-gamma*E*Ig2} 
      else{
        dE <- sigma*(I1+I2)*(1-E)-mue*E
        dIg1_nolag <- p*R*Mg1-mu*Ig1
        dIg2_nolag <- p*R*Mg2-mu*Ig2
      }# change in innate immune strength
      dI1_nolag <- p*R*M1-mu*I1-gamma*E*I1
      dI2_nolag <- p*R*M2-mu*I2-gamma*E*I2
      dM1_nolag <- -mum*M1-p*R*M1
      dM2_nolag <- -mum*M2-p*R*M2
      dMg1_nolag <- -mum*Mg1-p*R*Mg1
      dMg2_nolag <- -mum*Mg2-p*R*Mg2
      dG1_nolag <- -mug*G1
      dG2_nolag <- -mug*G2
      dID <- mu+gamma*E
    }
    
    if(immunity == "ni"){
      dI1_nolag <- p*R*M1-mu*I1 # change in infected RBC density
      dI2_nolag <- p*R*M2-mu*I2
      dIg1_nolag <- p*R*Mg1-mu*Ig1
      dIg2_nolag <- p*R*Mg2-mu*Ig2
      dM1_nolag <- -mum*M1-p*R*M1
      dM2_nolag <- -mum*M2-p*R*M2
      dMg1_nolag <- -mum*Mg1-p*R*Mg1
      dMg2_nolag <- -mum*Mg2-p*R*Mg2
      dG1_nolag <- -mug*G1
      dG2_nolag <- -mug*G2
    }
    
    if(immunity == "i") {
      if(gamete_immune){
        dI1_nolag <- p*R*M1-mu*I1-(a*I1)/(b+I1+I2+Ig1+Ig2) # change in infected RBC density with immunity
        dI2_nolag <- p*R*M2-mu*I2-(a*I2)/(b+I1+I2+Ig1+Ig2)
        dIg1_nolag <- p*R*Mg1-mu*Ig1-(a*Ig1)/(b+I1+I2+Ig1+Ig2)
        dIg2_nolag <- p*R*Mg2-mu*Ig2-(a*Ig2)/(b+I1+I2+Ig1+Ig2)
        dID <- mu+a/(b+I1+I2+Ig1+Ig2)
      } else{
        dI1_nolag <- p*R*M1-mu*I1-(a*I1)/(b+I1+I2)
        dI2_nolag <- p*R*M2-mu*I2-(a*I2)/(b+I1+I2)
        dIg1_nolag <- p*R*Mg1-mu*Ig1
        dIg2_nolag <- p*R*Mg2-mu*Ig2
        dID <- mu+a/(b+I1+I2)
      }
      
      dM1_nolag <- -mum*M1-p*R*M1
      dM2_nolag <- -mum*M2-p*R*M2
      dMg1_nolag <- -mum*Mg1-p*R*Mg1
      dMg2_nolag <- -mum*Mg2-p*R*Mg2
      dG1_nolag <- -mug*G1
      dG2_nolag <- -mug*G2

    } 
    
    ## Track states in initial cohort of infection
    if(t<=alpha){
      dI1 <- dI1_nolag-pulseBeta*S 
      dI2 <- dI2_nolag-pulseBeta*S 
      
      dM1 <- dM1_nolag+beta*pulseBeta*S # all of them are asexual merozoite
      dM2 <- dM2_nolag+beta*pulseBeta*S
      dMg1 <- 0 # should have no Mg before day 1
      dMg2 <- 0
      dIg1 <- 0 #first wave starts on day alpha
      dIg2 <- 0
      dG1 <- 0 # first wave starts on day alpha+alphag
      dG2 <- 0
    }
    
    if(t<=alpha+alphag && t>alpha){
      dG1 <- 0
      dG2 <- 0
      dIg1 <- dIg1_nolag
      dIg2 <- dIg2_nolag
    }
    
    ## Track states after delay 
    if(t>alpha){
      dI1 <- dI1_nolag-p*lag1[1]*lag1[2]*S 
      dI2 <- dI2_nolag-p*lag1[1]*lag1[3]*S 
      if(log_cue){
        dM1 <- dM1_nolag+beta*(1-cr_1(log(cue_lag1)))*p*lag1[1]*lag1[2]*S
        dM2 <- dM2_nolag+beta*(1-cr_2(log(cue_lag1)))*p*lag1[1]*lag1[3]*S
        dMg1 <- dMg1_nolag+beta*cr_1(log(cue_lag1))*p*lag1[1]*lag1[2]*S
        dMg2 <- dMg2_nolag+beta*cr_2(log(cue_lag1))*p*lag1[1]*lag1[3]*S}
      else{
        dM1 <- dM1_nolag+beta*(1-cr_1(cue_lag1))*p*lag1[1]*lag1[2]*S 
        dM2 <- dM2_nolag+beta*(1-cr_2(cue_lag1))*p*lag1[1]*lag1[3]*S
        dMg1 <- dMg1_nolag+beta*cr_1(cue_lag1)*p*lag1[1]*lag1[2]*S
        dMg2 <- dMg2_nolag+beta*cr_2(cue_lag1)*p*lag1[1]*lag1[3]*S}
    }
    
    if(t>alpha+alphag){
      dG1 <- dG1_nolag+p*lag2[1]*lag2[4]*Sg
      dG2 <- dG2_nolag+p*lag2[1]*lag2[5]*Sg
      dIg1 <- dIg1_nolag-p*lag2[1]*lag2[4]*Sg
      dIg2 <- dIg2_nolag-p*lag2[1]*lag2[5]*Sg
    }
    
    ## Return the states. Must be in the same order as states!
    if (immunity == "ni" || immunity == "i") {return(list(c(dR, dM1, dM2, dMg1, dMg2, dID, dI1, dI2, dIg1, dIg2, dG1, dG2, dA)))}
    
    if (immunity == "kochin") {return(list(c(dR, dM1, dM2, dMg1, dMg2, dID, dI1, dI2, dIg1, dIg2, dG1, dG2, dE, dA)))}
    
    if (immunity == "tsukushi") {return(list(c(dR, dM1, dM2, dMg1, dMg2, dID, dI1, dI2, dIg1, dIg2, dG1, dG2, dN, dW, dA)))}
  }
  
  #-------------------------#
  # Run single-infection model
  #------------------------#
  chabaudi_ci.df <- as.data.frame(deSolve::dede(y = state,
                                                times = time_range,
                                                func = chabaudi_ci_model_lag,
                                                p = parameters,
                                                method = solver,
                                                control=list(mxhist = 1e6)))
  
  #-------------------------#
  # Calculate fitness
  #------------------------#
  ## Get Gametocyte density time series data
  gam_1 <- chabaudi_ci.df$G1
  gam_2 <- chabaudi_ci.df$G2
  gam_1[gam_1<0] <- 0 # Assign negative gametocyte density to 0
  gam_2[gam_2<0] <- 0
  
  ## Get timeseries interval. Simplify first time after t=0
  int <- 1e-3
  
  ## Define the fitness parameter values
  aval <- -12.69
  bval <- 3.6
  dens <- log10(gam_1+gam_2)
  ## Calculate the transmission potential at each time t
  tau1.ls <- (exp(aval+bval*dens))/(1+exp(aval+bval*dens))*(gam_1/(gam_1+gam_2))
  tau2.ls <- (exp(aval+bval*dens))/(1+exp(aval+bval*dens))*(gam_2/(gam_1+gam_2))
                                                            
  ## Get approximation of cumulative transmission potential
  tau1.sum <- sum(tau1.ls*int, na.rm = TRUE)
  tau2.sum <- sum(tau2.ls*int, na.rm = TRUE)
  
  # return cumulative transmission potential. Turn negative to maximize
  if(dyn == FALSE){return(tau1.sum)} # let tau1 be the optimizing strain
  
  #-------------------------#
  # Simulating infection dynamics if Dyn == TRUE
  #------------------------# 
  if(dyn == TRUE) {
    ### calculate cumulative transmission potential gain. NAs or became 0
    tau_cum1.ls <- cumsum(ifelse(is.na(tau1.ls*int), 0, tau1.ls*int)) + tau1.ls*0
    tau_cum2.ls <- cumsum(ifelse(is.na(tau2.ls*int), 0, tau2.ls*int)) + tau2.ls*0

    ### cbind results
    chabaudi_ci.df$tau1 <- tau1.ls
    chabaudi_ci.df$tau2 <- tau2.ls
    
    chabaudi_ci.df$tau_cum1 <- tau_cum1.ls
    chabaudi_ci.df$tau_cum2 <- tau_cum2.ls
    
    ### calculate CR based on cue
    if(cue != "t"){
      cue_for_cr.df <- chabaudi_ci.df %>% dplyr::mutate(cue_state = eval(parse(text = cue)))
      cue_for_cr <- cue_for_cr.df$cue_state
      
      if(log_cue){
        cr1.ls <- cr_1(log(cue_for_cr))
        cr2.ls <- cr_2(log(cue_for_cr))}

      else{
        cr1.ls <- cr_1(cue_for_cr)
        cr2.ls <- cr_2(cue_for_cr)}
    }
    
    if(cue == "t"){
      cr1.ls <- cr_1(time_range)
      cr2.ls <- cr_2(time_range)
      }
    
    chabaudi_ci.df$cr1 <- cr1.ls
    chabaudi_ci.df$cr2 <- cr2.ls
    
    ### processing df for plotting
    #### If no adaptive immunity, filter out adaptive immunity
    if(!adaptive){chabaudi_ci.df2 <- chabaudi_ci.df %>% dplyr::select(-A)}
    chabaudi_ci.df3 <- chabaudi_ci.df2 %>% tidyr::gather(key = "variable", value = "value", -time)
    #### processing to separate stran from variable
    chabaudi_ci.df4 <- chabaudi_ci.df3 %>% 
      dplyr::mutate(variable_alt = gsub("[[:digit:]]","", variable),
                    strain = gsub("[^[:digit:]]", "", variable))
    
    return(chabaudi_ci.df4)
  }
}


