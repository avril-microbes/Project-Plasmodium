#-----------------------#
# Single infection model of P. chabaudi without RBC limitation
# Avril Wang
# Last edited 2022-07-11 based on 2022-06-22 iterations of chabaudi_si_clean.R script
# dual cue cable
#-----------------------#

chabaudi_si_clean_R <- function(
  parameters_cr, # input parameters for conversion rate reaction norm
  parameters, # sets of values for parameters in the mosel
  immunity, # immunity selection. Tsukushi's model, saturating immunity, or no immunity possible
  time_range, # time that simulation is ran for
  cue, # cue that parasite use to change conversion rate
  log_cue = "none", # whether to log10 transform cue
  cue_range, # range that cue is active against
  solver = "lsoda", # solver used for numerical integration. "vode" is often faster
  delay = 0, # the time at which merozoite is injected into the host
  drug = "none", # whether drug action is simulated. 
  drug_dose = 0, # drug dosage in mg/kg
  drug_admin = 0, # day at which drug is administered
  cue_b = "none", # if second cue is incorporated
  cue_range_b = "none", # if second cue is used, the range of cue that parasite use to adjust conversion rate to
  log_cue_b = "none",  # whether to log10 transform cue_b
  dyn = FALSE, # whether the function should return simulation dynamics rather than fitness
  neg = FALSE, # set to TRUE if using minimization function
  gam = "te" # set the gam model for dual cue. te = tensor interactions + main effect. ti = tensor interactions only
){
  
  #----------------------#
  # Force argument
  #----------------------#
  force(parameters_cr)
  force(parameters)
  force(immunity)
  force(time_range)
  force(cue)
  force(log_cue)
  force(cue_range)
  force(solver)
  force(delay)
  force(drug)
  force(drug_dose)
  force(drug_admin)
  force(cue_b)
  force(cue_range_b)
  force(log_cue_b)
  force(dyn)
  force(neg)
  force(gam)
  
  #----------------------#
  # Quality checks
  #----------------------#
  ## If dual cue, ensure length is correct
  #if (cue_b != "none" && length(parameters_cr) != 5) {
  #  stop("Must have 5 conversion rate parameters for dual cue!")
  #}
  
  ## Ensure immunity input is correct
  if (immunity != "ni" && immunity != "i" && immunity != "tsukushi"){
    stop("Immunity must be either 'ni', 'i', 'kochin,' or 'tsukushi'")
  }
  
  ## Double checking cues
  ### Ensure that time_range is used as cue_range when t is used
  if(cue == "t" && !isTRUE(all.equal(cue_range, time_range))){
    stop("Time is chosen as cue. Cue_range must equal to time_range")
  }
  
  if(cue_b == "t" && !isTRUE(all.equal(cue_range_b, time_range))){
    stop("Time is chosen as cue)b. Cue_range_b must equal to time_range")
  }
  
  ## Ensure that cue transformation is entered correctly
  if(log_cue != "none" && log_cue != "log10"){
    stop("log_cue must be either 'none' or 'log10'")
  }
  
  if(log_cue_b != "none" && log_cue_b != "log10"){
    stop("log_cue_b must be either 'none' or 'log10'")
  }
  
  ## Checking drug administration conditions
  if(drug > 0 && drug_admin < 0){
    stop("Drug administration date must be above 0!")
  }
  
  ## check ga model
  if(gam != "te" & gam != "ti"){
    stop("Gam model must be either te or ti")
  }
  
  #----------------------#
  # Define initial conditions
  #----------------------#
  ## when no immunity or saturating immunity
  if(immunity == "ni" || immunity == "i"){
    state <- c(R = parameters[["R1"]], # density of RBC
               M = 0, # density of asexually-committed merozoite
               Mg = 0, # density of sexually-committed merozoite
               ID = 0, # overall death rate of iRBC
               I = 0, # asexual iRBC density
               Ig = 0, # sexual iRBC density
               G = 0, # gametocyte density
               cr_t = 0) # conversion rate
  }
  
  ## when using Tsukushi's model of immunity
  if(immunity == "tsukushi"){
    state <- c(R = parameters[["R1"]],
               M = 0,
               Mg = 0,
               ID = 0,
               I = 0,
               Ig = 0,
               G = 0,
               N = 0, # general RBC removal
               W = 0, # targeted RBC removal
               cr_t = 0) 
  }
  
  #----------------------#
  # Describe initial population structure
  # of injected malaria
  #----------------------#
  ## function that takes initial merozoite emergence density (I0), shape parameter
  ## of beta function (sp), and time point (t)
  pulseBeta_fun <- function(I0, sp, t){ 
    res = vector(length = length(t))
    res = I0*(dbeta(t, sp, sp))
    return(res)
  }
  
  #----------------------#
  # define Heaviside transformation
  #----------------------#
  # function that transforms anything that above specificed max value to max. 
  # if value cue_range is below max, does not change the value.
  heaviside_trans <- function(cue_range, max){
    res <- crone::heaviside(cue_range)*(cue_range)+(crone::heaviside(cue_range-max)*(max-cue_range))
    return(res)
  }
  
  #----------------------#
  # Define conversion rate function
  #----------------------#
  ## Single cue conversion rate
  if(cue_b == "none"){
    ### Define dummy data of conversion rate
    dummy_y.vals <- rep(0, length(cue_range)) 
    dummy_cr.data <- as.data.frame(cbind(cue_range, dummy_y.vals))
    
    ### fit basic cubic spline with no internal knots. 
    dummy_cr.mod <- lm(dummy_y.vals ~ splines2::bSpline(x = cue_range, degree = 3))
    dummy_cr.mod$data <- dummy_cr.data
    
    ### assign coefficient to be optimized to the dummy conversion rate function
    dummy_cr.mod$coefficients <- parameters_cr
    
    ### double exponentiation conversion rate to get it between 0 and 1
    cr_fit <- exp(-exp(predict(dummy_cr.mod, newdata = data.frame(cue_range))))
    
    ### fit spline function to predicted conversion rate. This increases processing speed
    cr_fun <- splinefun(cbind(cue_range, cr_fit))
  }
  
  ## Dual cue conversion rate
  ### note that we are using gam with tensor product smoothing k = c(3,3)
  if(cue_b != "none"){
    ## create all combinations of 2 cues
    cr_grid <- expand.grid(cue_range, cue_range_b)
    ## rename
    names(cr_grid) <- c("cue_range", "cue_range_b")
    ## create dummy y
    dummy_y <- runif(length(cue_range_b), 0, 1)
    ## put together df
    dummy_df <- data.frame(cue_range, cue_range_b, dummy_y)
    ## gam model
    if(gam == "te"){ # note that te include tensor full interactions and main effects. Need 9 coefficients
      dummy_cr.mod <- mgcv::gam(dummy_y ~ te(cue_range, cue_range_b, 
                                             k = c(3,3)), 
                                data = dummy_df)
    }
    
    if(gam == "ti"){ # note that ti includes only tensor interactions and NO main effects. Need 5 coefficeints.
      dummy_cr.mod <- mgcv::gam(dummy_y ~ ti(cue_range, cue_range_b, 
                                             k = c(3,3)), 
                                data = dummy_df)
    }
    ## assign parameters
    dummy_cr.mod$coefficients <- parameters_cr
    
    # exponential transformation to limit conversion rate to between 0 and 1
    cr_fun <- function(cue_1, cue_2){
      exp(-exp(mgcv::predict.gam(dummy_cr.mod, 
                                 newdata = data.frame("cue_range" = cue_1,
                                                      "cue_range_b" = cue_2))))
    }
    
  }
  
  #-----------------------------------------#
  #-----------------------------------------#
  # Define dynamics of single infection model
  #-----------------------------------------#
  #-----------------------------------------#
  
  chabaudi_si_dyn <- function(t, state, parameters){
    
    #----------------------#
    # Redefine parameters for cleaner code
    #----------------------#
    ## Rename parameters for cleaner code. With.list not used to speed up computation
    if(immunity == "tsukushi"){
      R1 <- parameters["R1"] # maximum RBC density at homeostasis
      rho <- parameters["rho"] # proportion of RBC deviation from R1 restored/day
      psin <- parameters["psin"] # activation strength for general RBC removal
      psiw <- parameters["psiw"] # activation strength for targeted RBC removal
      phin <- parameters["phin"] # half life for general RBC removal
      phiw <- parameters["phiw"] # half life for targeted RBC removal
      iota <- parameters["iota"] # iaximum iRBC detection limit (higher = lower immunity activations)
    }
    
    mu <- parameters["mu"] # death rate/day of RBC and iRBC
    p <- parameters["p"] # probability of invasion success upon merzoite-RBC contact
    alpha <- parameters["alpha"] # length of asexual development period of iRBC (I)
    alphag <- parameters["alphag"] # length of sexual developement period of iRBC (Ig)
    beta <- parameters["beta"] # number of merozoite produced per iRBC
    mum <- parameters["mum"] # merozoite death rate/day
    mug <- parameters["mug"] # gametocyte death rate/day
    I0 <- parameters["I0"] # initial dosage of iRBC injected
    sp <- parameters["sp"] # shape parameters for initial malaria population structure. 1 for non-synchronous, 100 for high synchrony
    
    if(immunity == "i"){
      a <- parameters["a"] # maximum rate of iRBC removal/day with saturating immunity
      b <- parameters["b"] # iRBC density needed to achieve half maximum iRBC removal rate
    }
    
    if(immunity != "tsukushi"){lambda <- parameters["lambda"]} # maximum RBC replenishment rate
    
    if(drug_dose > 0){mud <- parameters["mud"]} # if drug action is included, add drug induced death rate
    if(drug_dose == 0){mud <- 0} # if no drug is administered, no drug-induced mortality
    
    #----------------------#
    # Redefine states for cleaner code
    #----------------------#
    R <- state["R"]
    I <- state["I"]
    Ig <- state["Ig"]
    ID <- state["ID"]
    S <- state["S"]
    M <- state["M"]
    G <- state["G"]
    Mg <- state["Mg"]
    cr_t <- state["cr_t"]
    if(immunity == "tsukushi"){
      N <- state["N"]
      W <- state["W"]
    }
    
    #-----Defining Pulse beta function based on current time-----#
    pulseBeta <- pulseBeta_fun(I0, sp, t-delay)
    
    #----------------------#
    # Define lag terms
    #----------------------#
    ## lag1 = value 1 day ago
    if(t>alpha+delay){
      lag1 = deSolve::lagvalue(t-alpha)
      
    } 
    
    ## lag2 = value 2 days ago
    if(t>alphag+delay){
      lag2 = deSolve::lagvalue(t-alphag)
    } 
    
    #----------------------#
    # Define cue value
    #----------------------#
    ##------get present cue------##
    ### convert cue to time if time-based conversion rate strategy is used
    if(cue == "t"){
      cue_state <- t
    }
    
    ### get cue_state if it is state-based
    if(cue != "t"){
      cue_state <- state[cue]
      if(cue_b != "none"){cue_state_b <- state[cue_b]}
    }
    
    ##------get lagged cue values------##
    ### if cue is not time-based... get index at which cue is based on
    #### if multiple cues are given, multiple indexes are returned
    if(cue != "t") {
      lag.i <- match(unlist(stringr::str_split(cue, "\\+|\\-|\\*|\\/")), names(state))
      if(cue_b != "none"){lag.i_b <- match(unlist(stringr::str_split(cue_b, "\\+|\\-|\\*|\\/")), names(state))}
    }
    
    ### For simple cues (if it does not contain special characters)
    if(stringr::str_detect(cue, "\\+|\\-|\\*|\\/", negate = TRUE)){
      #### if cue is time-based
      if(t>alpha+delay && cue == "t"){ ##### lag 1 day ago
        cue_lag1 <- t-alpha} 
      
      if(t>alphag+delay && cue == "t") { ##### lag 2 days ago
        cue_lag2 <- t-alphag
      } 
      
      #### if cue is state-based
      if(t>alpha+delay && cue != "t"){ ##### lag 1 day ago
        cue_lag1 <- lag1[lag.i]
        ##### get second cue if second cue is given
        if(cue_b != "none"){cue_lag1_b <- lag1[lag.i_b]}
      }
      
      if(t>alphag+delay && cue != "t") { ##### lag 2 days ago
        cue_lag2 <- lag2[lag.i]
        if(cue_b != "none"){cue_lag2_b <- lag2[lag.i_b]}
      } ### for complex cues that involve additions
    } else{
      if(stringr::str_detect(cue, "\\+")){ 
        #### for cue 1 day ago
        if(t>alpha+delay && cue != "t"){
          cue_lag1 <- lag1[lag.i[1]]+lag1[lag.i[2]] #### add first cue to second cue
          if(cue_b != "none"){cue_lag1_b <- lag1[lag.i_b[1]]+lag1[lag.i_b[2]]}
        }
      }
    }
    
    #------------------#
    # Process cue values
    #------------------#
    if(t>alpha+delay){
      if(log_cue == "log10"){
        cue_lag1_log <- log10(abs(cue_lag1)+5e-324) # log the lagged cue
        ## heaviside transformation
        cue_lag1_p <- heaviside_trans(cue_lag1_log, max(cue_range)) 
      } 
      ## IMPORTANT: none of the cue lags should be naturally negative. adding this prevents
      ## stiff "dipping" of cue from producing NAs. 
      if(log_cue == "none"){
        cue_lag1_p <- heaviside_trans(cue_lag1, max(cue_range))
      } # keep it the same
      
      if(cue_b == "none"){cr <- cr_fun(cue_lag1_p)}
      
      ## if second cue is used
      if(cue_b != "none"){
        if(log_cue_b == "log10"){
          cue_lag1_log_b <- log10(abs(cue_lag1_b)+5e-324) # log the lagged cue
          ## heaviside transformation
          cue_lag1_p_b <- heaviside_trans(cue_lag1_log_b, max(cue_range_b)) 
        } 
        if(log_cue_b == "none"){
          cue_lag1_p_b <- heaviside_trans(cue_lag1_b, max(cue_range_b))
        } 
        # get cr fnction
        cr <- cr_fun(cue_lag1_p, cue_lag1_p_b)
      }
    }
    
    
    #----------------#
    # Drug actions
    #----------------#
    # later
    
    #----------------#
    # Survival functions
    #----------------#
    # before first iRBC maturation, survival function of iRBC
    if(t<=alpha+delay ){
      S <- exp(-ID) ## survival function of asexual iRBC. ID = cumulative hazard rate in 1 day
      Sg <- 0 ## survival function of sexual iRBC. No sexual iRBC before alpha so set to 0
    } 
    
    # between first asexual iRBC burst but before first sexual iRBC burst
    if(t>alpha+delay && t<=alpha+alphag+delay){
      S <- exp(-ID + lag1[4]) # need to account for previous cumulation of cumulative hazard
      Sg <- exp(-ID)} # does not appear in our equation until first infected RBC burst, which is delay+alpha+alphag
    
    # after first sexual iRBC burst
    if(t>alpha+alphag+delay){
      S <- exp(-ID + lag1[4])
      Sg <- exp(-ID + lag2[4]) # sexual iRBC must survive for longer (2 days), hence lag2
    }
    
    #-----------------#
    # Model with lag terms (lose of iRBC that has matured and burst)
    # makes the code cleaner looking
    #-----------------#
    ## Define K, maximum RBC density
    if(immunity != "tsukushi"){K <- lambda*R1/(lambda-mu*R1)}
    
    # asexual merozoite
    dM_nolag <- (-mum*M)-(p*R*M)
    # sexual merozoite
    dMg_nolag <- (-mum*Mg)-(p*R*Mg)
    # gametocyte
    dG_nolag <- -mug*G
    
    # If Tsukushi's model of immunity is used, use the following
    if(immunity == "tsukushi"){
      ## define N and W transformed, 2 variables that might go to negative and give NAN in optimization process. 
      ## do heaviside transformation that maintains N/W at 0.999 if N/W>=1 and N/W at 0 if N/W <0. 
      N_trans <- ((crone::heaviside(N)*N)+crone::heaviside(N-0.999)*(0.999-N))
      W_trans <- ((crone::heaviside(W)*W)+crone::heaviside(W-0.999)*(0.999-W))
      
      ## RBC density. changed to static. R always remain the same
      dR <- 0
      ## asexual iRBC 
      dI_nolag <- (p*R*M)-(mu*I)-((-log(1-N_trans)-log(1-W_trans))*I)
      ## sexual iRBC 
      dIg_nolag <- (p*R*Mg)-(mu*Ig)-((-log(1-N_trans)-log(1-W_trans))*Ig)
      ## indiscriminant RBC removal
      dN <- psin*((I+Ig)/iota)*(1-N_trans)-(N_trans/phin) 
      ## targeted iRBC removal
      dW <- psiw*((I+Ig)/iota)*(1-W_trans)-(W_trans/phiw)
      ## hazard function of iRBC
      dID <- mu-log(1-N_trans)-log(1-W_trans)
    }
    
    # if no immunity is used
    if(immunity == "ni"){
      # RBC density. again changed to static
      dR <- 0
      dI_nolag <- (p*R*M)-(mu*I)
      dIg_nolag <- (p*R*Mg)-(mu*Ig)
      dID <- mu
    }
    
    # if saturating immunity is used
    if(immunity == "i") {
      # RBC density
      dR <- 0
      dI_nolag <- (p*R*M)-(mu*I)-((a*I)/(b+I))
      dIg_nolag <- (p*R*Mg)-(mu*Ig)
      dID <- mu+(a/(b+I))
    }
    
    #--------------#
    # Infection dynamics in the first cohort of 
    # injected parasite
    #--------------#
    # Before delay, no iRBC is produced!
    if(t<delay){
      dI <- 0
    }
    
    # before all first cohort of asexual iRBC bursts (before day 1)
    if(t<=alpha+delay){
      dI <- dI_nolag-(pulseBeta*S) # some initial asexual iRBC burst due to maturation
      dM <- dM_nolag+(beta*pulseBeta*S) # asexual merozoites are produced when asexual iRBC burst
      dMg <- 0 # should have no Mg before day 1
      dIg <- 0 #first wave starts on day alpha
      dG <- 0 # first wave starts on day alpha+alphag
    }
    
    # the period after first production of sexual iRBC but before they burst
    if(t<=alpha+alphag+delay && t>alpha+delay){
      dG <- 0 # no gametocyte production
      dIg <- dIg_nolag # no sexual iRBC death from previous cycle
    }
    
    # after all first cohort of asexual iRBC burst (after day 1)
    if(t>alpha+delay){
      dI <- dI_nolag-(p*lag1[1]*lag1[2]*S) # bursting of iRBC produced alpha days ago
      dM <- dM_nolag+(beta*(1-cr)*p*lag1[1]*lag1[2]*S) # production of asexual merozoite from asexual iRBC burst
      dMg <- dMg_nolag+(beta*cr*p*lag1[1]*lag1[2]*S) # production of sexual merozoite from asexual iRBC burst
    }
    
    # after the first gametocyte production
    if(t>alpha+alphag+delay){
      dG <- dG_nolag+(p*lag2[1]*lag2[3]*Sg) # production of gametocyte from sexual iRBC produced alphag days ago
      dIg <- dIg_nolag-(p*lag2[1]*lag2[3]*Sg) # loss of sexual iRBC
    }
    
    # track cr
    if(t<=alpha+delay){dcr_t <- 0}
    if(t>alpha+delay){dcr_t <- cr}
    
    
    #----------------------#
    # Return the states
    #----------------------#
    if (immunity == "ni" || immunity == "i") {return(list(c(dR, dM, dMg, dID, dI, dIg, dG, dcr_t)))}
    if (immunity == "tsukushi") {return(list(c(dR, dM, dMg, dID, dI, dIg, dG, dN, dW, dcr_t)))}
  } 
  
  #---------------------------------------------------------#
  #---------------------------------------------------------#
  #-------------------End of dynamics function--------------#
  #---------------------------------------------------------#
  #---------------------------------------------------------#
  
  #----------------------#
  # Create injection event
  # needed for delayed infection
  #----------------------#
  delay_injection <- data.frame(var = "I",
                                time = delay,
                                value = parameters["I0"],
                                method = "add")
  
  
  #-------------------------#
  # Run single-infection model
  #------------------------#
  chabaudi_si.df <- as.data.frame(deSolve::dede(y = state,
                                                times = time_range,
                                                func = chabaudi_si_dyn,
                                                p = parameters,
                                                method = solver,
                                                events = list(data = delay_injection),
                                                control=list(mxhist = 1e7)))
  
  #-------------------# 
  # Calculate fitness for optimization
  #-------------------# 
  # Get Gametocyte density time series data
  gam <- chabaudi_si.df$G
  gam[gam<0] <- 0 ## Assign negative gametocyte density to 0. can arise due to stiffness of function
  
  
  
  # Get timeseries interval
  int <- chabaudi_si.df$time[2]
  
  # Define the fitness parameter values
  aval <- -12.69
  bval <- 3.6
  dens <- log10(gam) 
  
  # Calculate the transmission potential at each time point
  tau.ls <- (exp(aval+(bval*dens)))/(1+exp(aval+(bval*dens)))
  
  # Get approximation of cumulative transmission potential
  tau.sum <- sum(tau.ls*int, na.rm = TRUE)
  
  # return cumulative transmission potential. Used for optimization
  if(dyn == FALSE && neg == FALSE){return(tau.sum)} ## 
  if(dyn == FALSE && neg == TRUE){
    tau.sumneg <- tau.sum*-1
    return(tau.sumneg)
  } ## if using minimization algorthm, set neg to TRUE to get maximize
  #----------------------------#
  #----------------------------#
  #-Simulate infection dynamics#
  #-If dyn = TRUE--------------#
  #----------------------------#
  #----------------------------#
  
  if(dyn == TRUE) {
    
    #---------------#
    # Calculate cumulative transmission potential
    #--------------#
    tau_cum.ls <- cumsum(tau.ls*int)
    
    # cbind results
    chabaudi_si.df$tau <- tau.ls
    chabaudi_si.df$tau_cum <- tau_cum.ls
    
    #-------------#
    # calculate cr
    #-------------#
    # time based cue
    if(cue == "t"){
      cr.ls <- cr_fun(time_range)
      chabaudi_si.df$cr <- cr.ls
    }
    
    # state-based cue
    if(cue != "t"){
      chabaudi_si.df <- chabaudi_si.df %>% 
        dplyr::mutate(cr = (cr_t - dplyr::lag(cr_t))*(1/time_range[2]))
    }
    
    # make df long for ease of plotting
    chabaudi_si.df2 <- chabaudi_si.df %>% tidyr::gather(key = "variable", value = "value", -time)
    
    # return
    return(chabaudi_si.df2)
  }
}
