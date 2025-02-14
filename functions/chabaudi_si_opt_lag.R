# Single strain infection model with Plasmodium chabaudi
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

chabaudi_si_opt_lag <- function(parameters_cr,
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
                                 log_cue = "none",
                                 delay = 0,
                                 drug = 0,
                                 admin = 0,
                                 ratio = 1,
                                 lag_deriv = FALSE,
                                 lag_smooth = 0) {
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
  #force(integration)
  force(adaptive)
  force(dyn)
  force(log_cue)
  force(delay)
  force(drug)
  force(admin)
  force(ratio)
  force(lag_deriv)
  force(lag_smooth)
  
  #-------------------------#
  # Define initial condition
  #------------------------#
  if(immunity == "ni" || immunity == "i"){
    state <- c(R = parameters[["R1"]],
               M = 0,
               Mg = 0,
               ID = 0,
               S = 0,
               I = 0,
               Ig = 0,
               G = 0, 
               A = 0) # survival function. Just for tracking
  } else if(immunity == "kochin"){
    state <- c(R = parameters[["R1"]],
               M = 0,
               Mg = 0,
               ID = 0,
               S = 0,
               I = 0,
               Ig = 0,
               G = 0,
               E = 0,
               A = 0)
  } else{
    state <- c(R = parameters[["R1"]],
               M = 0,
               Mg = 0,
               ID = 0,
               S = 0,
               I = 0,
               Ig = 0,
               G = 0,
               N = 0, # general RBC removal
               W = 0,
               A = 0) # targeted RBC removal
  }
  
  #-------------------------#
  # Ensure inputs are correct
  #------------------------#
  ## Ensure length of initial parameter search space matches with df. 
  if (length(parameters_cr) != df+1) {
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
  ## Ensure that cue transformation is entered correctly
  if(log_cue != "none" && log_cue != "log" && log_cue != "log10"){
    stop("log_cue must be either 'norne' or 'log' or 'log10'")
  }
  ## Ensure administration time is above 0 if drug is administered
  if(drug > 0 && admin <= 0){
    stop("Drug administration date must be above 0!")
  }
  ## for now, lag smooth only implemented for derivative cue
  if(lag_smooth > 0 && lag_deriv == FALSE){
    stop("lag smoothing only available for derivative based cue")
  }
  ## stop if lag smooth is negative
  if(lag_smooth < 0){
    stop("lag smooth must be 0 or positive!")
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
  if(transformation == "norm"){
    cr_fit <- predict(dummy_cr.mod, newdata = data.frame(cue_range))
    cr_fit_2 <- (cr_fit-min(cr_fit))/(max(cr_fit-min(cr_fit)))
  } else if(transformation == "exp"){
    cr_fit_2 <- exp(-exp(predict(dummy_cr.mod, newdata = data.frame(cue_range))))
  } else{
    cr_fit <- predict(dummy_cr.mod, newdata = data.frame(cue_range))
    cr_fit_2 <- 1 / (1 + exp(-cr_fit))
  }
  
  ## Get spline function where cr ~ cue
  cr <- splinefun(cbind(cue_range, cr_fit_2))
  
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
  chabaudi_si_model_lag <- function(t, state, parameters) {
    
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
    if(drug > 0){mud <- parameters["mud"]} # if drug action is included, add drug length
    if(drug == 0){mud <- 0} # assign no drug induced death if no drugs
    
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
    I <- state["I"]
    Ig <- state["Ig"]
    ID <- state["ID"]
    S <- state["S"]
    M <- state["M"]
    G <- state["G"]
    A <- state["A"]
    Mg <- state["Mg"]
    if (immunity == "kochin") {E <- state["E"]}
    if (immunity == "tsukushi"){
      N <- state["N"]
      W <- state["W"]
    }
    # threshold or not
    
    ## Defining Pulse beta function based on current time
    pulseBeta <- pulseBeta_fun(I0*ratio, sp, t-delay)
    
    ## Define the lag terms. lag[1] = R, lag[2] = I, lag[3] = Ig, lag[4] = M, lag[5] = G
    if(t>alpha+delay){
      lag1 = deSolve::lagvalue(t-alpha)
      dlag1 = deSolve::lagderiv(t-alpha)
      dlagsmooth = deSolve::lagvalue(t-lag_smooth) # for derivative based cue smoothing
      
    } # lag state for asexual development
    if(t>alphag+delay){
      lag2 = deSolve::lagvalue(t-alphag)
      dlag2 = deSolve::lagderiv(t-alphag)
    } # lag state for gametocyte development
    ### extra lag term for adaptive immunity
    if(adaptive == TRUE){
      if(t>epsilon+delay){lag3 = deSolve::lagvalue(t-epsilon)}
    }
    
    ## get lag term index given cue. Cannot use else if given that during simulation, multiple iterations of cue_lag is used
    ### Only get lag index when it is a state-based cue. Multiple indexes are returned
    ### if multiple cues are given, multiple indexes are returned
    if(cue != "t") {
      lag.i <- match(unlist(stringr::str_split(cue, "\\+|\\-|\\*|\\/")), names(state))
    }
    
    
    
    ### convert cue to time if time-based conversion rate strategy is used
    if(cue == "t"){
      cue_state <- t}
    
    ### get cue_state if it is state-based
    if(cue != "t"){
      cue_state <- state[cue]}
    
    ### define lagged cue. Lag1 = alpha times ago, lag2 = alphag times ago
    ### For simple cues (if it does not contain special characters)
    if(stringr::str_detect(cue, "\\+|\\-|\\*|\\/", negate = TRUE)){
      if(t>alpha+delay && cue == "t"){
        cue_lag1 <- t-alpha} 
      
      if(t>alpha+delay && cue != "t"){
        if(lag_deriv == FALSE){cue_lag1 <- lag1[lag.i]}
        if(lag_deriv == TRUE){
          if(lag_smooth == 0){
            cue_lag1 <- dlag1[lag.i]} # cue is perceived instatenous deriv alpha days ago
          if(lag_smooth > 0){ # cue is the average deriv over smoothing period. Recall average of derivative is integral/period. 
            cue_lag1 <- (cue_state-dlagsmooth[lag.i])/lag_smooth
          }
        }
      }
      
      if(t>alphag+delay && cue == "t") {
        cue_lag2 <- t-alphag
      } 
      
      if(t>alphag+delay && cue != "t") {
        if(lag_deriv == FALSE){cue_lag2 <- lag2[lag.i]}
        if(lag_deriv == TRUE){cue_lag2 <- dlag2[lag.i]}
      }
      
    } else{### manually create lag values if cues contain special characters
      if(stringr::str_detect(cue, "\\+")){ # if it contains plus
        if(lag_deriv == FALSE){
          if(t>alpha+delay && cue != "t"){cue_lag1 <- lag1[lag.i[1]]+lag1[lag.i[2]]}
          if(t>alphag+delay && cue != "t") {cue_lag2 <- lag2[lag.i[1]]+lag2[lag.i[2]]}}
        if(lag_deriv == TRUE){
          if(t>alpha+delay && cue != "t"){cue_lag1 <- dlag1[lag.i[1]]+dlag1[lag.i[2]]}
          if(t>alphag+delay && cue != "t") {cue_lag2 <- dlag2[lag.i[1]]+dlag2[lag.i[2]]}}
      }
      if(stringr::str_detect(cue, "\\-")){ # if it contains -
        if(lag_deriv == FALSE){
          if(t>alpha+delay && cue != "t"){cue_lag1 <- lag1[lag.i[1]]-lag1[lag.i[2]]}
          if(t>alphag+delay && cue != "t") {cue_lag2 <- lag2[lag.i[1]]-lag2[lag.i[2]]}}
        if(lag_deriv == TRUE){
          if(t>alpha+delay && cue != "t"){cue_lag1 <- dlag1[lag.i[1]]-dlag1[lag.i[2]]}
          if(t>alphag+delay && cue != "t") {cue_lag2 <- dlag2[lag.i[1]]-dlag2[lag.i[2]]}}
      }
      if(stringr::str_detect(cue, "\\*")){ # if it contains multiplication
        if(lag_deriv == FALSE){
          if(t>alpha+delay && cue != "t"){cue_lag1 <- lag1[lag.i[1]]*lag1[lag.i[2]]}
          if(t>alphag+delay && cue != "t") {cue_lag2 <- lag2[lag.i[1]]*lag2[lag.i[2]]}}
        if(lag_deriv == TRUE){
          if(t>alpha+delay && cue != "t"){cue_lag1 <- dlag1[lag.i[1]]*dlag1[lag.i[2]]}
          if(t>alphag+delay && cue != "t") {cue_lag2 <- dlag2[lag.i[1]]*dlag2[lag.i[2]]}}
      }
      if(stringr::str_detect(cue, "\\/")){ # if it contains division
        if(lag_deriv == FALSE){
          if(t>alpha+delay && cue != "t"){cue_lag1 <- lag1[lag.i[1]]/lag1[lag.i[2]]}
          if(t>alphag+delay && cue != "t") {cue_lag2 <- lag2[lag.i[1]]/lag2[lag.i[2]]}}
        if(lag_deriv == TRUE){
          if(t>alpha+delay && cue != "t"){cue_lag1 <- dlag1[lag.i[1]]/dlag1[lag.i[2]]}
          if(t>alphag+delay && cue != "t") {cue_lag2 <- dlag2[lag.i[1]]/dlag2[lag.i[2]]}}
      }
      ### get present states
      if(cue != "t"){cue_state <- eval(parse(text = cue))}
    }
    
    ## Define K, carrying capacity of RBC
    K <- lambda*R1/(lambda-mu*R1)
    
    #-------------------------#
    # Function to describe pyremethamine
    # length of action from https://onlinelibrary.wiley.com/doi/10.1111/eva.12516
    #------------------------#
    if(drug > 0){
      pyr_length = 3.557-2.586/(1+exp(-8.821+drug))
      if(t<=admin){P <- 0}
      if(t > admin && t <= admin + 1 + pyr_length) {
        P <- mud # only drug action after drug administration and within active time frame
      }
      
      if(t> admin + 1 + pyr_length){P <- 0}
    }
    
    if(drug == 0){P <- 0}
    
    
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
    if(t>alpha+delay && immunity == "ni"){
      S <- exp(-ID + lag1[4])} 
    
    if(t>alpha+delay && immunity != "ni"){
      S <- exp(-ID + lag1[4])
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
    
    if(t<=alpha+delay && immunity == "ni"){
      S <- exp(-ID)} 
    
    if(t<=alpha+delay && immunity != "ni"){
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
      if(t<=alpha+delay){
        Sg <- 0 # not relevent
      } 
      
      if(t>alpha+delay && t<=alpha+alphag+delay){
        Sg <- 0} # does not appear in our equation until first infected RBC burst, which is delay+alpha+alphag
      
      if(t>alpha+alphag+delay){
        Sg <- exp(-ID + lag2[4]) # only due to intrinsic cell death
      }
    }
    
    if(immunity == "tsukushi"){
      if(t<=alpha+delay){
        Sg <- 0 # not relevent
      }
      
      if(t>alpha+delay && t<=alpha+alphag+delay){
        Sg <- 0 # not relevent
      }
      
      if(t>alpha+alphag+delay){
        Sg <- exp(-ID + lag2[4])
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
    
    dS <- S
    
    ## Define the models without lag terms. 
    if(immunity != "tsukushi"){
      dR <- lambda*(1-(R/K))-(mu*R)-(p*R*M)-(p*R*Mg) # change in susceptible RBC
    } 
    
    if(immunity == "tsukushi"){ #Tsukushi exclusive ODEs
      #dR <- lambda*(1-R/K)-mu*R-p*R*M-(mu-log(1-N))*R
      dR <- R1*mu+rho*(R1-R)-(mu-log(1-N))*R-(p*R*M)-(p*R*Mg)
      dI_nolag <- p*R*M-mu*I-(-log(1-N)-log(1-W)-log(1-A))*I-P*I
      #dIg_nolag <- cr(cue_state)*p*R*M-mu*Ig-(-log(1-N))*Ig #assume no targeted clearance
      dIg_nolag <- p*R*Mg-mu*Ig-(-log(1-N)-log(1-W)-log(1-A))*Ig-P*Ig
      #dN <- psin*(I/iota)*(1-N)-(N/phin) # assume Ig does not elicit strong immune response. Not included in cue
      dN <- psin*((I+Ig)/iota)*(1-N)-(N/phin)
      #dW <- psiw*(I/iota)*(1-W)-(W/phiw)
      dW <- psiw*((I+Ig)/iota)*(1-W)-(W/phiw)
      dM_nolag <- (-mum*M)-(p*R*M)
      dMg_nolag <- (-mum*Mg)-(p*R*Mg)
      dG_nolag <- -mug*G
      dID <- mu-log(1-N)-log(1-W)-log(1-A)+P
    }
    
    if(immunity =="kochin"){
      dE <- sigma*I*(1-E)-mue*E # change in innate immune strength
      dI_nolag <- p*R*M-mu*I-gamma*E*I-P*I
      dIg_nolag <- p*R*Mg-mu*Ig-P*Ig
      dM_nolag <- -mum*M-p*R*M
      dMg_nolag <- -mum*Mg-p*R*Mg
      dG_nolag <- -mug*G
      dID <- mu+gamma*E+P
    }
    
    if(immunity == "ni"){
      dI_nolag <- p*R*M-mu*I-P*I # change in infected RBC density
      dIg_nolag <- p*R*Mg-mu*Ig-P*Ig
      dM_nolag <- -mum*M-p*R*M
      dMg_nolag <- -mum*Mg-p*R*Mg
      dG_nolag <- -mug*G
      dID <- mu+P
    }
    
    if(immunity == "i") {
      dI_nolag <- p*R*M-mu*I-(a*I)/(b+I)-P*I # change in infected RBC density with immunity
      dIg_nolag <- p*R*Mg-mu*Ig-P*Ig
      dM_nolag <- -mum*M-p*R*M
      dMg_nolag <- -mum*Mg-p*R*Mg
      dG_nolag <- -mug*G
      dID <- mu+a/(b+I)+P
    } 
    
    if(t<delay){
      dI <- 0
    }
    
    ## Track states in initial cohort of infection
    if(t<=alpha+delay){
      dI <- dI_nolag-pulseBeta*S 
      dM <- dM_nolag+beta*pulseBeta*S # all of them are asexual merozoite
      dMg <- 0 # should have no Mg before day 1
      dIg <- 0 #first wave starts on day alpha
      dG <- 0 # first wave starts on day alpha+alphag
    }
    
    if(t<=alpha+alphag+delay && t>alpha+delay){
      dG <- 0
      dIg <- dIg_nolag
    }
    
    ## Track states after delay 
    if(t>alpha+delay){
      dI <- dI_nolag-p*lag1[1]*lag1[2]*S 
      
      if(log_cue == "log"){
        dM <- dM_nolag+beta*(1-cr(log(cue_lag1)))*p*lag1[1]*lag1[2]*S
        dMg <- dMg_nolag+beta*cr(log(cue_lag1))*p*lag1[1]*lag1[2]*S}
      if(log_cue == "none"){
        dM <- dM_nolag+beta*(1-cr(cue_lag1))*p*lag1[1]*lag1[2]*S 
        dMg <- dMg_nolag+beta*cr(cue_lag1)*p*lag1[1]*lag1[2]*S}
      if(log_cue == "log10"){
        dM <- dM_nolag+beta*(1-cr(log10(cue_lag1)))*p*lag1[1]*lag1[2]*S
        dMg <- dMg_nolag+beta*cr(log10(cue_lag1))*p*lag1[1]*lag1[2]*S}
    }
    
    if(t>alpha+alphag+delay){
      dG <- dG_nolag+p*lag2[1]*lag2[3]*Sg
      dIg <- dIg_nolag-p*lag2[1]*lag2[3]*Sg
    }
    
    ## Return the states. Must be in the same order as states!
    if (immunity == "ni" || immunity == "i") {return(list(c(dR, dM, dMg, dID, dS, dI, dIg, dG, dA)))}
    
    if (immunity == "kochin") {return(list(c(dR, dM, dMg, dID,  dS, dI,dIg, dG, dE, dA)))}
    
    if (immunity == "tsukushi") {return(list(c(dR, dM, dMg, dID, dS, dI, dIg, dG, dN, dW, dA)))}
  }
  #--------------------------#
  # Create event for strain 1 injection (delayed)
  #--------------------------#
  delay_injection <- data.frame(var = "I",
                                time = delay,
                                value = parameters["I0"]*ratio,
                                method = "add")
  
  #-------------------------#
  # Run single-infection model
  #------------------------#
  chabaudi_si.df <- as.data.frame(deSolve::dede(y = state,
                                                times = time_range,
                                                func = chabaudi_si_model_lag,
                                                p = parameters,
                                                method = solver,
                                                events = list(data = delay_injection),
                                                control=list(mxhist = 1e6)))
  
  #-------------------------#
  # Calculate fitness
  #------------------------#
  ## Get Gametocyte density time series data
  gam <- chabaudi_si.df$G
  gam[gam<0] <- 0 # Assign negative gametocyte density to 0
  
  ## Get timeseries interval. Simplify first time after t=0
  int <- 1e-3
  
  ## Define the fitness parameter values
  aval <- -12.69
  bval <- 3.6
  dens <- log10(gam)
  
  ## Calculate the transmission potential at each time t
  tau.ls <- (exp(aval+bval*dens))/(1+exp(aval+bval*dens))
  
  ## Get approximation of cumulative transmission potential
  tau.sum <- sum(tau.ls*int)
  
  # return cumulative transmission potential. Turn negative to maximize
  if(dyn == FALSE){return(tau.sum)} # if running other than ga, set control = list(fnscale = -1)==
  
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
    if(cue != "t"){
      if(lag_deriv == FALSE){
        cue_for_cr.df <- chabaudi_si.df %>% dplyr::mutate(cue_state = eval(parse(text = cue)))
        cue_for_cr <- cue_for_cr.df$cue_state}
      
      if(lag_deriv == TRUE){
        
        if(lag_smooth == 0){
          cue_for_cr.df <- chabaudi_si.df %>% 
            dplyr::mutate(cue_state = eval(parse(text = cue))) %>% 
            dplyr::mutate(cue_dif = (cue_state-dplyr::lag(cue_state))/0.001)
          cue_for_cr <- cue_for_cr.df$cue_dif}
        
        if(lag_smooth >0){
          cue_for_cr.df <- chabaudi_si.df %>% 
            dplyr::mutate(cue_state = eval(parse(text = cue))) %>% 
            dplyr::mutate(cue_dif = (cue_state-dplyr::lag(cue_state, (lag_smooth/0.001)))/lag_smooth)
          cue_for_cr <- cue_for_cr.df$cue_dif}
      }
      
      if(log_cue == "log"){
        cr.ls <- cr(log(cue_for_cr))}
      if(log_cue == "none"){
        cr.ls <- cr(cue_for_cr)}
      if(log_cue == "log10"){
        cr.ls <- cr(log10(cue_for_cr))}
    }
    
    if(cue == "t"){cr.ls <- cr(time_range)}
    
    chabaudi_si.df$cr <- cr.ls
    
    ### processing df for plotting
    #### If no adaptive immunity, filter out adaptive immunity
    if(!adaptive){chabaudi_si.df2 <- chabaudi_si.df %>% 
      dplyr::select(-A) %>% 
      dplyr::mutate(S_new = (S - dplyr::lag(S)*1000)) %>% 
      dplyr::select(-S)}
    
    chabaudi_si.df3 <- chabaudi_si.df2 %>% tidyr::gather(key = "variable", value = "value", -time)
    
    return(chabaudi_si.df3)
  }
}


