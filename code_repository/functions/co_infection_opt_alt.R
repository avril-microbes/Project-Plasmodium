# wrapper around co-infection model that will perform "stationary" optimization,
# meaning that the optimal parameter for each iteration is assigned as both the
# the resident and mutant parameter. Mutant parameter is allowed to change
# as to maximize the difference between mutant and resident. L-BFGS-B is used as optimization algorithm.

# note for new iteration, simplifying code but adding second layer: if cannot improve or fitness = 0, repeat optimization with mutant set to 0.5 x4

# By Avril Wang
# last updated 2022-05-05

co_infection_opt_alt <- function(parameters_cr,  # preliminary parameter set
                             limit, # minimum fitness difference between competing strains to break
                             model, # infection model
                             ...){ # additional parameters to be inserted into the infection model
  force(parameters_cr)
  force(limit)
  force(model)
  additional_arg <- list(...) # unpack for do.call
  
  # assign default parameters
  # assign iteration index
  index <- 1
  output_ls <- list()
  lower_repeat <- FALSE # is this a repeat run after lower strain 1 fitness?

  # do-while loop to continue to optimize strain 1 fitness
  repeat{
    ## for first iteration, both strain adopts initial parameter condition assigned
    if(index == 1){
      residence_par <- parameters_cr
      mutant_par <- parameters_cr
      print(residence_par)
      print(mutant_par)
    } else{
      ### if strain 1 is fitter than strain 2, assign both strains to adopt the optimized par
      if(fitness > 0){
        residence_par <- par
        mutant_par <- par
        print("Strain 1 is fitter than strain 2.")
        print(residence_par)
        print(mutant_par)
      }
      ### special case: if fitness difference = 0, need skip out 
      if(fitness == 0){
        if(index > 2){
          residence_par <- par ### assign residence parameter to previous iteratio
          mutant_par <- rep(0.5, 4) #### assign mutant to 0.5 x4 
          print(residence_par)
          print(mutant_par)
        }
        if(index == 2){
          residence_par <- par # special case. If first round of optimization fails, go back to initial parameter set
          mutant_par <- rep(0.5, 4)
          print(residence_par)
          print(mutant_par)
        } 
        lower_repeat <- TRUE
        print("Strain 1 has same fitness than strain 2. Repeat optimization with 0.5x4")
      }
    }
    ## code to execute parallel LGBF-GS optimization
    output <- list() # reset per run
    print(paste("starting iteration", index))
    # add index
    index <- index + 1
    cl <- makeCluster(detectCores()); setDefaultCluster(cl = cl)
    model_output <- do.call(optimParallel::optimParallel, c(list(par = mutant_par, # competing parameter
                                                                 fn = model,
                                                                 control = list(trace = 6, fnscale = -1),
                                                                 parameters_cr_2 = residence_par), # reassign resident parameter to optimal one
                                                            additional_arg))
    
    
    ## save output
    par <- model_output$par
    fitness <- round(model_output$value, 10) # rounding to 10 digits so that almost 0 values becomes 0
    output <- list(index, mutant_par, residence_par, fitness)
    output_ls[[index]] <- output
    print(output)
    
    
    # exit loop IF
    ## previous run converged, strain 1 is fitter than strain 2, and that the fitness difference is minute
    if(fitness < limit && fitness != 0){
      print("Strain 1 difference is minute.")
      stopCluster(cl)
      break
    }
    ## If repeat run after lower strain 1 fitness, break if continue to be 
    if(lower_repeat == TRUE && fitness < limit){
      print("After repetition, either strain 1 is still less fit or that fitness difference is minute.")
      stopCluster(cl)
      break
    }
  }
  
  # final output
  return(output_ls)
  }