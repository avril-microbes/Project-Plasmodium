#--------------------#
# Coinfection model of malaria chabaudi
#-------------------#
# Used for infection simulation of co-infection model and for optimization of best conversion rate strategy
# Avril Wang
# last edited 2022-03-25
## added heaviside transfromation to constrain cue range
## Note, added "sum" as possibility for cue. This would give us I1+I2+Ig1+Ig2

# code reviewed version of chabaudi_ci_opt_lag
# strain 1 is set as the invading strain whereas strain 2 is set as the residence strain
# all optimization will be performed according to strain 1 performance and subsequent
# infection manipulation (altering injection and timing of invasion) only applies to strain 1
# dual cue not implemented due to potential performance issues

chabaudi_ci_clean <- function(parameters_cr_1, # parameters for strain 1 conversion rate strategy
                              parameters_cr_2, # parameters for strain 2 conversion rate strategy
                              immunity, # mode of immunity. ni, i, or tsukushi available
                              parameters, # parameters for model
                              time_range, # time range at which infection is simulated
                              cue_1, # cue of strain 1
                              cue_2, # cue of strain 2
                              cue_range_1, # cue range of strain 1
                              cue_range_2, # cue range of strain 2
                              log_cue_1 = "none", # whether to log transform cue 1
                              log_cue_2 = "none", # whether to log transform cue 2
                              solver = "lsoda", # solver for numerical integration. Vode often gives faster runs
                              dyn = FALSE, # whether to simulate infection dynamics (false for optimization)
                              delay = 0, # when strain 1 injection occurs days-post infection
                              ratio = 1){ # ratio of # of strain 1 injected compared to strain 2
  
  #-------------------------#
  # Ensure values we inputted 
  # are available in environment
  #------------------------#
  force(parameters_cr_1)
  force(parameters_cr_2)
  force(immunity)
  force(parameters)
  force(time_range)
  force(cue_1)
  force(cue_2)
  force(cue_range_1)
  force(cue_range_2)
  force(solver)
  force(dyn)
  force(log_cue_1)
  force(log_cue_2)
  force(delay)
  force(ratio)
  
  #-------------------------#
  # Define initial conditions
  #------------------------#
  # initial conditions if no immunity or saturating immunity is used
  if(immunity == "ni" || immunity == "i"){
    state <- c(R = parameters[["R1"]], # RBC density
               M1 = 0, # asexual merozoite density of strain 1
               M2 = 0, # asexual merozoite density of strain 2
               Mg1 = 0, # sexual merozoite density of strain 1
               Mg2 = 0, # sexual merozoite density of strain 2
               ID = 0, # death rate of all iRBC (same across strain 1 and strain 2)
               I1 = 0, # initial iRBC density injected for strain 1. set to 0 for delayed infection. If not delayed, I0 parasite inject @ 0.001
               I2 = parameters[["I0"]], # initial iRBC density injected for strain2 (residence)
               Ig1 = 0, # sexual iRBC density of strain 1
               Ig2 = 0, # sexual iRBC density of strain 2
               G1 = 0, # gametocyte density of strain 1
               G2 = 0, # gametocyte density of strain 3
               cr_t1 = 0,
               cr_t2 = 0)  
  } 
  if(immunity == "tsukushi"){ # initial states if tsukushi's mode of immunity is used 
    state <- c(R = parameters[["R1"]],
               M1 = 0,
               M2 = 0,
               Mg1 = 0,
               Mg2 = 0,
               ID = 0,
               I1 = 0, # set to 0 for delayed infection. If not delayed, parasite inject @ 0.001
               I2 = parameters[["I0"]],
               Ig1 = 0,
               Ig2 = 0,
               G1 = 0,
               G2 = 0,
               N = 0, # indiscriminant RBC removal
               W = 0, # targetted iRBC removal
               cr_t1 = 0,
               cr_t2 = 0)
  }
    
    #-------------------------#
    # Ensure inputs are correct
    #------------------------#
    # Ensure immunity input is correct
    if (immunity != "ni" && immunity != "i" && immunity != "tsukushi") {
      stop("Immunity must be either 'ni', 'i', tsukushi'")
    }
    
    # Ensure cue is correct
    ## cue of strain 1
    if (!(unlist(stringr::str_split(cue_1, "\\+|\\-|\\*|\\/")) %in% names(state)) && 
        !((cue_1 %in% names(state)) | cue_1 == "sum") && cue_1 != "t") {
      stop("Cue 1 must be one of the states or time")
    }
    
    ## cue of strain 2
    if (!(unlist(stringr::str_split(cue_2, "\\+|\\-|\\*|\\/")) %in% names(state)) && 
        !((cue_2 %in% names(state)) | cue_2 == "sum") && cue_2 != "t") {
      stop("Cue 2 must be one of the states or time")
    }
    
    # Ensure that time_range is used as cue_range when t is used
    if(cue_1 == "t" && (!isTRUE(all.equal(cue_range_1, time_range)))){
      stop("Time is chosen as cue_1. Cue_range_1 must equal to time_range")
    }
    if(cue_2 == "t" && (!isTRUE(all.equal(cue_range_2, time_range)))){
      stop("Time is chosen as cue_2. Cue_range_2 must equal to time_range")
    }

    # Ensure that cue transformation is entered correctly
    if(log_cue_1 != "none" && log_cue_1 != "log10"){
      stop("log_cue_1 must be either 'none' or 'log10'")
    }
    if(log_cue_2 != "none" && log_cue_2 != "log10"){
      stop("log_cue_2 must be either 'none' or 'log10'")
    }
    
    # Ensure delay is...
    ## not longer then infection period
    if(delay > max(time_range)){
      stop("Delay must be before infection period ends")
    }
    ## not negative
    if(delay < 0){
      stop("Delay must be positive")
    }
    
    # Ensure ratio is within acceptable range (not 0 or negative)
    if(ratio <= 0){
      stop("Ratio must be above 0!")
    }
  
  #-------------------------#
  # Function to describe population 
  # structure of initial inoculum
  #------------------------#
  pulseBeta_fun <- function(I0, sp, t){ 
    res = rep(NA, length(t))
    res = I0*(dbeta(t, sp, sp))
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
  
  #-------------------------#
  # Define conversion rate function. 
  # Simplified to increase performance. 
  #------------------------#
  # Define dummy data of conversion rate
  dummy_y.vals_1 <- rep(0, length(cue_range_1)) 
  dummy_cr.data_1 <- as.data.frame(cbind(cue_range_1, dummy_y.vals_1)) ## make df of dummy data
  
  ## repeat for cue of strain 2
  dummy_y.vals_2 <- rep(0, length(cue_range_2)) 
  dummy_cr.data_2 <- as.data.frame(cbind(cue_range_2, dummy_y.vals_2))
  
  # fit basic cubic spline with no internal knots
  dummy_cr.mod_1 <- lm(dummy_y.vals_1 ~ splines2::bSpline(x = cue_range_1, degree = 3))
  dummy_cr.mod_1$data <- dummy_cr.data_1
  
  ## repeat for strain 2
  dummy_cr.mod_2 <- lm(dummy_y.vals_2 ~ splines2::bSpline(x = cue_range_2, degree = 3))
  dummy_cr.mod_2$data <- dummy_cr.data_2
  
  # Assign coefficient to be optimized to the dummy conversion rate function
  dummy_cr.mod_1$coefficients <- parameters_cr_1
  dummy_cr.mod_2$coefficients <- parameters_cr_2
  
  # use spline function to predict cr. Double exponentitate to get conversion rate to be between 0 and 1
  cr_fit_t_1 <- exp(-exp(predict(dummy_cr.mod_1, newdata = data.frame(cue_range_1))))
  cr_fit_t_2 <- exp(-exp(predict(dummy_cr.mod_2, newdata = data.frame(cue_range_2))))

  # Get spline function where cr ~ cue
  cr_1_fun <- splinefun(cbind(cue_range_1, cr_fit_t_1))
  cr_2_fun <- splinefun(cbind(cue_range_2, cr_fit_t_2))
    
  
  #---------------------------#
  #---------------------------#
  # Define co-infection model
  #---------------------------#
  #---------------------------#
  chabaudi_ci_dyn <- function(t, state, parameters){
    
    #----------------------#
    # Rename parameters for cleaner code
    #----------------------#
    R1 <- parameters["R1"] # maximum RBC density at homeostasis
    mu <- parameters["mu"] # baseline death rate of RBC/iRBC (/day)
    p <- parameters["p"] # probability of merozoite invasion success
    alpha <- parameters["alpha"] # length of asexual iRBC development cycle
    alphag <- parameters["alphag"] # length of sexual iRBC development cycle
    beta <- parameters["beta"] # burst size (number of merozoites produced per asexual iRBC)
    mum <- parameters["mum"] # death rate of merozoite (/day)
    mug <- parameters["mug"] # death rate of gametocyte (/day)
    sp <- parameters["sp"] # shape parameter for beta distribution that describes level of malaria synchrony (1 = asynchronous, 100 = very synchronous)
    I0 <- parameters["I0"] # initial asexual iRBC injected (strain 2). If ratio = 1, also the value for strain 1
    
    # additional parameters if saturating immunity is used
    if(immunity == "i") {
      a <- parameters["a"] # maximum iRBC removal rate /day for saturating immunity (not used if i or tsukushi chosen for immunity)
      b <- parameters["b"] # total iRBC density needed to activate half of maximum iRBC removal 
    }

    # additional parameters if tsukushi's model of immunity is used
    if (immunity == "tsukushi") {
      psin <- parameters["psin"] # activation strength for indiscriminate RBC removal
      psiw <- parameters["psiw"] # activation strength for targeted iRBC removal
      phin <- parameters["phin"] # half life for indiscriminant RBC removal
      phiw <- parameters["phiw"] # half life for targeted iRBC removal
      iota <- parameters["iota"] # cue strength (denominator for total infected iRBC)
      rho <- parameters["rho"] #  proportion of the deviation from the homeostatic equilibrium restored by the host per day 
    }
    
    # only no immunity or saturating immunity uses lambda
    if(immunity != "tsukushi"){lambda <- parameters["lambda"]} # maximum RBC replenishment rate

    #----------------------#
    # Rename states for cleaner code
    #----------------------#
    R <- state["R"] # RBC density
    I1 <- state["I1"] # asexual iRBC density of strain 1
    I2 <- state["I2"] # asexual iRBC density of strain 2
    Ig1 <- state["Ig1"] # sexual iRBC density of strain 1
    Ig2 <- state["Ig2"] # sexual iRBC density of strain 2
    ID <- state["ID"] # death rate of all iRBC for both strain 1 and strain 2
    M1 <- state["M1"] # asexual merozoite density of strain 1
    M2 <- state["M2"] # asexual merozoite density of strain 2
    Mg1 <- state["Mg1"] # sexual merozoite density of strain 1
    Mg2 <- state["Mg2"] # sexual merozoite density of strain 2
    G1 <- state["G1"] # gametocyte density of strain 1
    G2 <- state["G2"]  # gametocyte density of strain 2
    cr_t1 <- state["cr_t1"] # track conversion rate of strain 1
    cr_t2 <- state["cr_t2"] # track conversion rate of strain 2

    if (immunity == "tsukushi"){
      N <- state["N"] # extent of indiscriminant RBC removal
      W <- state["W"] # extent of targetted iRBC removal
    }
    
    #---------------------#
    # Define initial iRBC population structure
    # based on beta distribution
    #---------------------#
    pulseBeta_1 <- pulseBeta_fun(I0*ratio, sp, t-delay) # only strain 1 can experience delay and decrease in injected amount
    pulseBeta_2 <- pulseBeta_fun(I0, sp, t) # strain 2 does not experience delay
    
    #---------------------#
    # Define lag terms
    #---------------------#
    #------strain 1-------# (a for alpha days ago, b = for alphag ago. 1 for strain number)
    if(t>alpha+delay){lag_a_1 = deSolve::lagvalue(t-alpha)} ## lag state for asexual development
    if(t>alphag+delay){lag_b_1 = deSolve::lagvalue(t-alphag)} ## lag state for gametocyte development
    
    #------strain 2-------#
    if(t>alpha){lag_a_2 = deSolve::lagvalue(t-alpha)} # lag state for asexual development
    if(t>alphag){lag_b_2 = deSolve::lagvalue(t-alphag)} # lag state for gametocyte development
    
    #-----Define lag index based on cue-------#
    ## strain 1
    if(cue_1 != "t" && cue_1 != "sum") {
      lag.i_1 <- match(unlist(stringr::str_split(cue_1, "\\+|\\-|\\*|\\/")), names(state)) 
    }
    
    ## strain 2
    if(cue_2 != "t" && cue_2 != "sum") {
      lag.i_2 <- match(unlist(stringr::str_split(cue_2, "\\+|\\-|\\*|\\/")), names(state)) 
    }
    
    #---------------------#
    # Define lagged cue values
    #---------------------#
    # cue_lag_lagextent_strain (a = alpha days ago, b = alphag days ago, 1 = strain 1, 2 = strain 2)
    
    #-----------Simple cue without addition for strain 1------------#
    if(stringr::str_detect(cue_1, "\\+|\\-|\\*|\\/", negate = TRUE)){
      ## time-based cue for strain 1
      if(t>alpha+delay && cue_1 == "t"){ ### cue lagged alpha days ago
        cue_lag_a1 <- t-alpha
      } 
      ## State-based cue for strain 1
      if(t>alpha+delay && cue_1 != "t" && cue_1 != "sum"){ ### cue lagged alpha days ago
        cue_lag_a1 <- lag_a_1[lag.i_1]
      }
      if(t>alpha+delay && cue_1 == "sum"){cue_lag_a1 <- lag_a_1[7]+lag_a_1[8]+lag_a_1[9]+lag_a_1[10]} ## if sum is chosen as cue, add up all iRBC values (I1+I2+Ig1+Ig2)
    } #-----------------------complex cues for addition-based cues-----------------------#
    if(stringr::str_detect(cue_1, "\\+")){ # addition only (currently only supports two cues addition)
        if(t>alpha+delay && cue_1 != "t"){cue_lag_a1 <- lag_a_1[lag.i_1[1]]+lag_a_1[lag.i_1[2]]}
      }

    
    #----------------Simple cue without addition for strain 2-------#
    if(stringr::str_detect(cue_2, "\\+|\\-|\\*|\\/", negate = TRUE)){
      ## time-based cues
      if(t>alpha+delay && cue_2 == "t"){
        cue_lag_a2 <- t-alpha
      } 
      
      ## state-based cues
      if(t>alpha+delay && cue_2 != "t" && cue_2 != "sum"){
        cue_lag_a2 <- lag_a_2[lag.i_2]
      }
      if(t>alpha+delay && cue_2 == "sum"){cue_lag_a2 <- lag_a_2[7]+lag_a_2[8]+lag_a_2[9]+lag_a_2[10]}
      
    } #------------complex cues involving addition of 2 cues-----------------------#
    if(stringr::str_detect(cue_2, "\\+")){ # if it contains plus. strain 1
        if(t>alpha+delay && cue_2 != "t"){cue_lag_a2 <- lag_a_2[lag.i_2[1]]+lag_a_2[lag.i_2[2]]}
      }
    
    
    #------------------#
    # Process cue values
    #------------------#
    # Strain 1
    if(t>alpha+delay){
      if(log_cue_1 == "log10"){cue_lag1_p <- heaviside_trans(log10(abs(cue_lag_a1)+5e-324), max(cue_range_1))} # log the lagged cue
      if(log_cue_1 == "none"){cue_lag1_p <- heaviside_trans(cue_lag_a1, max(cue_range_1))} # keep it the same
      cr_1 <- cr_1_fun(cue_lag1_p)
    }
    
    # Strain 2
    if(t>alpha){
      if(log_cue_2 == "log10"){cue_lag2_p <- heaviside_trans(log10(abs(cue_lag_a2)+5e-324), max(cue_range_2))} # log the lagged cue
      if(log_cue_2 == "none"){cue_lag2_p <- heaviside_trans(cue_lag_a2, max(cue_range_2))} # keep it the same
      cr_2 <- cr_2_fun(cue_lag2_p)
    }
    
    #---------------------#
    # Define survival functions
    #---------------------#
    # S_strain (1 = strain 1, 2 = strain 2)
    
    #-------Survival function of asexual iRBC---------------#
    #-------Strain 1-------#
    ## hazard function is constant (mu) if no immunity is used
    if(immunity == "ni"){
      if(t<=alpha+delay){S_1 <- exp(-mu*t)}
      if(t>alpha+delay){S_1 <- exp(-mu*alpha)}
      } 
    
    ## hazard function varies if tsukushi or saturating immunity is used
    if(immunity != "ni"){
      if(t<=alpha+delay){S_1 <- exp(-ID)}
      if(t>alpha+delay){S_1 <- exp(-ID + lag_a_1[6])} ## subtracts lagged ID from overall hazard rate
    }
    
    #-------Strain 2-------#
    if(immunity == "ni"){
      if(t<=alpha){S_2 <- exp(-mu*t)}
      if(t>alpha){S_2 <- exp(-mu*alpha)} 
    }
    
    if(immunity != "ni"){
      if(t<=alpha){S_2 <- exp(-ID)}
      if(t>alpha){S_2 <- exp(-ID + lag_a_2[6])}
    }
    
    #-------Survival function of sexual iRBC---------------#
    ## before production of sexual iRBC, set to 0
    if(t<=alpha+delay){
      Sg_1 <- 0 ### strain 1
    }
    
    if(t<=alpha){
      Sg_2 <- 0 ### strain 2
    }
    
    ## if saturating immunity or no immunity is used. Sexual iRBC not killed by immunity
    if(immunity != "tsukushi"){
      ## after production of first sexual iRBC but before any of them  burst
      if(t>alpha+delay && t<=alpha+alphag+delay){
        Sg_1 <- exp(-mu*t) ### strain 1
        } 
      
      if(t>alpha && t<=alpha+alphag){
        Sg_2 <- exp(-mu*t) ### strain 2
        } 
      
      ## after first burst of sexual iRBC
      if(t>alpha+alphag+delay){
        Sg_1 <- exp(-mu*alphag) ### strain 1
        }
      
      if(t>alpha+alphag){
        Sg_2 <- exp(-mu*alphag) ### strain 2
        }
    }
    
    ## if tsukushi's mode of immunity is used. Sexual iRBC is killed by immunity
    if(immunity == "tsukushi"){
      if(t>alpha+delay && t<=alpha+alphag+delay){
        Sg_1 <- exp(-ID)
      }
      
      if(t>alpha && t<=alpha+alphag){
        Sg_2 <- exp(-ID)
      }
      
      if(t>alpha+alphag+delay){
        Sg_1 <- exp(-ID + lag_b_1[6])
      }
      
      if(t>alpha+alphag){
        Sg_2 <- exp(-ID + lag_b_2[6])
      }
    }
    
    #---------------------#
    # Define ODEs with no lag terms
    # makes the code cleaner
    #---------------------#
    
    # change in RBC density for no immunity or saturating immunity model
    if(immunity != "tsukushi"){
      dR <- lambda*(1-(R/K))-(mu*R)-(p*R*M1)-(p*R*M2)-(p*R*Mg1)-(p*R*Mg2) 
    } 
    
    # universal no-lagged ODEs
    dM1_nolag <- (-mum*M1)-(p*R*M1) 
    dM2_nolag <- (-mum*M2)-(p*R*M2)
    dMg1_nolag <- (-mum*Mg1)-(p*R*Mg1)
    dMg2_nolag <- (-mum*Mg2)-(p*R*Mg2)
    dG1_nolag <- -mug*G1
    dG2_nolag <- -mug*G2
    
    # tracks cr 
    if(t>alpha+delay){dcr_t1 <- cr_1} else{dcr_t1 <- 0}
    if(t>alpha){dcr_t2 <- cr_2} else{dcr_t2 <- 0}
    
    # ODEs exclusive to tsukushi's model
    if(immunity == "tsukushi"){ #Tsukushi exclusive ODEs
      ## define N and W transformed, 2 variables that might go to negative and give NAN in optimization process. 
      ## do heaviside transformation that maintains N/W at 0.999 if N/W>=1 and N/W at 0 if N/W <0. 
      N_trans <- ((crone::heaviside(N)*N)+crone::heaviside(N-0.999)*(0.999-N))
      W_trans <- ((crone::heaviside(W)*W)+crone::heaviside(W-0.999)*(0.999-W))
      
      dR <- (R1*mu)+(rho*(R1-R))-((mu-log(1-N_trans))*R)-(p*R*M1)-(p*R*M2)-(p*R*Mg1)-(p*R*Mg2)
      dI1_nolag <- (p*R*M1)-(mu*I1)-((-log(1-N_trans)-log(1-W_trans))*I1)
      dI2_nolag <- (p*R*M2)-(mu*I2)-((-log(1-N_trans)-log(1-W_trans))*I2)
      dIg1_nolag <- (p*R*Mg1)-(mu*Ig1)-((-log(1-N_trans)-log(1-W_trans))*Ig1)
      dIg2_nolag <- (p*R*Mg2)-(mu*Ig2)-((-log(1-N_trans)-log(1-W_trans))*Ig2)
      dN <- psin*((I1+I2+Ig1+Ig2)/iota)*(1-N_trans)-(N_trans/phin)
      dW <- psiw*((I1+I2+Ig1+Ig2)/iota)*(1-W_trans)-(W_trans/phiw)
      dID <- mu-log(1-N_trans)-log(1-W_trans)
    }
    
    # ODEs exclusive to no immunity model
    if(immunity == "ni"){
      dI1_nolag <- (p*R*M1)-(mu*I1)
      dI2_nolag <- (p*R*M2)-(mu*I2)
      dIg1_nolag <- (p*R*Mg1)-(mu*Ig1)
      dIg2_nolag <- (p*R*Mg2)-(mu*Ig2)
      dID <- mu
    }
    
    # ODEs exclusive to saturating immunity model
    if(immunity == "i"){
      dI1_nolag <- (p*R*M1)-(mu*I1)-((a*I1)/(b+I1+I2+Ig1+Ig2))
      dI2_nolag <- (p*R*M2)-(mu*I2)-((a*I2)/(b+I1+I2+Ig1+Ig2))
      dIg1_nolag <- (p*R*Mg1)-(mu*Ig1)-((a*Ig1)/(b+I1+I2+Ig1+Ig2))
      dIg2_nolag <- (p*R*Mg2)-(mu*Ig2)-((a*Ig2)/(b+I1+I2+Ig1+Ig2))
      dID <- mu+(a/(b+I1+I2+Ig1+Ig2))
    }
    
    #----------------------------#
    # Tracking variables with lagged states
    #----------------------------#
    #------Before delay----------#
    if(t<delay){dI1 <- 0} # strain 1 not injected yet
    
    #------Before first asexual iRBC burst-------------#
    ## Strain 1
    if(t<=alpha+delay){
      dI1 <- dI1_nolag-(pulseBeta_1*S_1)
      dM1 <- dM1_nolag+(beta*pulseBeta_1*S_1) # all of them are asexual merozoite
      dMg1 <- 0 # should have no Mg before day alpha
      dIg1 <- 0 #first wave starts on day alpha+delay
      dG1 <- 0 # first wave starts on day alpha+alphag+delay
    }
    
    ## Strain 2
    if(t<=alpha){
      dI2 <- dI2_nolag-(pulseBeta_2*S_2) 
      dM2 <- dM2_nolag+(beta*pulseBeta_2*S_2)
      dMg2 <- 0
      dIg2 <- 0
      dG2 <- 0
    }
    
    #------After first asexual iRBC burst but before first sexual iRBC burst-----------#
    ## Strain 1
    if(t>alpha+delay && t<=alpha+alphag+delay){
      dG1 <- 0 # no production of gametocyte yet
      dIg1 <- dIg2_nolag # no death due to mature sexual iRBC bursting
    }
    
    ## Strain 2
    if(t>alpha && t<=alpha+alphag){
      dG2 <- 0 
      dIg2 <- dIg2_nolag 
    }
    
    #-----After first asexual iRBC burst-----#
    ## Strain 1
    if(t>alpha+delay){
      dI1 <- dI1_nolag-(p*lag_a_1[1]*lag_a_1[2]*S_1)
      dM1 <- dM1_nolag+(beta*(1-cr_1)*p*lag_a_1[1]*lag_a_1[2]*S_1)
      dMg1 <- dMg1_nolag+(beta*cr_1*p*lag_a_1[1]*lag_a_1[2]*S_1)

      }
    
    ## Strain 2
    if(t>alpha){
      dI2 <- dI2_nolag-(p*lag_a_2[1]*lag_a_2[3]*S_2)
      dM2 <- dM2_nolag+(beta*(1-cr_2)*p*lag_a_2[1]*lag_a_2[3]*S_2)
      dMg2 <- dMg2_nolag+(beta*cr_2*p*lag_a_2[1]*lag_a_2[3]*S_2)
    }
    
    #-----After first sexual iRBC burst-----#
    ## Strain 1
    if(t>alpha+alphag+delay){
      dG1 <- dG1_nolag+(p*lag_b_1[1]*lag_b_1[4]*Sg_1)
      dIg1 <- dIg1_nolag-(p*lag_b_1[1]*lag_b_1[4]*Sg_1)
    }
    
    ## Strain 2
    if(t>alpha+alphag){
      dG2 <- dG2_nolag+(p*lag_b_2[1]*lag_b_2[5]*Sg_2)
      dIg2 <- dIg2_nolag-(p*lag_b_2[1]*lag_b_2[5]*Sg_2)
    }
    
    #----------------------#
    # Return the states
    #----------------------#
    if (immunity == "ni" || immunity == "i") {return(list(c(dR, dM1, dM2, dMg1, dMg2, dID, dI1, dI2, dIg1, dIg2, dG1, dG2, dcr_t1, dcr_t2)))}
    if (immunity == "tsukushi") {return(list(c(dR, dM1, dM2, dMg1, dMg2, dID, dI1, dI2, dIg1, dIg2, dG1, dG2, dN, dW, dcr_t1, dcr_t2)))}
  }
  
  #----------------------------#
  #----------------------------#
  # End of dynamics function
  #----------------------------#
  #----------------------------#
  
  #--------------------------#
  # Create event for strain 1 injection (delayed)
  #--------------------------#
  delay_injection <- data.frame(var = "I1",
                                time = delay,
                                value = parameters["I0"]*ratio,
                                method = "add")
  
  #-------------------------#
  # Run single-infection model
  #------------------------#
  chabaudi_ci.df <- as.data.frame(deSolve::dede(y = state,
                                                times = time_range,
                                                func = chabaudi_ci_dyn,
                                                p = parameters,
                                                method = solver,
                                                events = list(data = delay_injection),
                                                control=list(mxhist = 1e6)))
  
  #-------------------------#
  # Calculate fitness
  #------------------------#
  # Get Gametocyte density time series data
  gam_1 <- chabaudi_ci.df$G1
  gam_2 <- chabaudi_ci.df$G2
  gam_1[gam_1<0] <- 0 # Assign negative gametocyte density to 0
  gam_2[gam_2<0] <- 0
  
  # Get time series interval
  int <- chabaudi_ci.df$time[2]
  
  # Define the fitness parameter values
  aval <- -12.69
  bval <- 3.6
  dens <- log10(gam_1+gam_2)
  
  # Calculate the transmission potential at each time t
  tau1.ls <- (exp(aval+(bval*dens)))/(1+exp(aval+(bval*dens)))*(gam_1/(gam_1+gam_2))
  tau2.ls <- (exp(aval+(bval*dens)))/(1+exp(aval+(bval*dens)))*(gam_2/(gam_1+gam_2))
  
  # Get approximation of cumulative transmission potential
  tau1.sum <- sum(tau1.ls*int, na.rm = TRUE)
  tau2.sum <- sum(tau2.ls*int, na.rm = TRUE)
  
  # return cumulative transmission potential difference between strain 1 and strain 2
  if(dyn == FALSE){return(tau1.sum-tau2.sum)}
  #-------------------------#
  # Simulating infection dynamics if Dyn == TRUE
  #------------------------# 
  if(dyn == TRUE){
    # calculate cumulative transmission potential gain. NAs or became 0
    tau_cum1.ls <- cumsum(ifelse(is.na(tau1.ls*int), 0, tau1.ls*int)) + tau1.ls*0
    tau_cum2.ls <- cumsum(ifelse(is.na(tau2.ls*int), 0, tau2.ls*int)) + tau2.ls*0
    
    ## cbind results
    chabaudi_ci.df$tau1 <- tau1.ls
    chabaudi_ci.df$tau2 <- tau2.ls
    
    chabaudi_ci.df$tau_cum1 <- tau_cum1.ls
    chabaudi_ci.df$tau_cum2 <- tau_cum2.ls
    
    #------Calculate cr------#
    ## time-based conversion rate
   if(cue_1 == "t"){
    cr.ls1 <- cr_1_fun(time_range)
    chabaudi_ci.df$cr_1 <- cr.ls1}
    
    if(cue_2 == "t"){
      cr.ls2 <- cr_2_fun(time_range)
      chabaudi_ci.df$cr_2 <- cr.ls2}

    ## state-based conversion rate
    if(cue_1 != "t"){
      chabaudi_ci.df <- chabaudi_ci.df %>% 
      dplyr::mutate(cr_1 = (cr_t1 - dplyr::lag(cr_t1))*1000)}
    
   if(cue_2 != "t"){
      chabaudi_ci.df <- chabaudi_ci.df %>% 
        dplyr::mutate(cr_2 = (cr_t2 - dplyr::lag(cr_t2))*1000)
    }
    
    #-----Process for better looking df------#
    ## make longer
    chabaudi_ci.df2 <- chabaudi_ci.df %>% tidyr::gather(key = "variable", value = "value", -time)
    
    ## processing to separate strain from variable
    chabaudi_ci.df3 <- chabaudi_ci.df2 %>% 
      dplyr::mutate(variable_alt = gsub("[[:digit:]]","", variable),
                    strain = gsub("[^[:digit:]]", "", variable))
    
    return(chabaudi_ci.df3)
  }
}
