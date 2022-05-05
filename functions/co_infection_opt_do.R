# wrapper around co-infection model that will perform "stationary" optimization,
# meaning that the optimal parameter for each iteration is assigned as both the
# the resident and mutant parameter. Mutant parameter is allowed to change
# as to maximize the difference between mutant and resident. DEoptim (GA) is used as optimization algorithm.

# By Avril Wang
# last updated 2022-05-05

co_infection_opt_do <- function(parameters_cr,  # preliminary parameter set
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
  worse_then_previous <- FALSE # is the current run worse then previous?
  
  # do-while loop to continue to optimize strain 1 fitness
  repeat{
    ## for first iteration, both strain adopts initial parameter condition assigned
    if(index == 1){
      residence_par <- parameters_cr
      mutant_par <- parameters_cr
    } else{
      ### if strain 1 is fitter than strain 2, assign both strains to adopt the optimized par
      if(fitness > 0){
        residence_par <- par
        mutant_par <- par
        print("Strain 1 is fitter than strain 2.")
      }
      ### if strain 1 is less fit than strain 2 but optimization did not converge, repeat with 
      ### mutant strain 1 final parameter values. Keep residence parameter as it is
      if(fitness <= 0 && conv != 0){
        if(index > 2){residence_par <- output_ls[[index - 1]][[4]]} ### assign residence parameter to previous iteration
        if(index == 2){residence_par <- parameters_cr} # special case. If first round of optimization fails, go back to initial parameter set
        mutant_par <- par
        print("Did not converge. Strain 1 is less fit than strain 2. Repeat optimization.")
      }
      
      ### if run converged but strain 1 is less fit then strain 2, try again with new initial mutant start point
      if(fitness <= 0 && conv == 0){
        if(index > 2){residence_par <- output_ls[[index - 1]][[4]]}
        if(index == 2){residence_par <- parameters_cr}
        mutant_par <- par
        lower_repeat <- TRUE
        print("Converged. Strain 1 is less fit than strain 2. Repeat optimization.")
      }
    }
    ## code to execute parallel LGBF-GS optimization
    output <- list() # reset per run
    print(paste("starting iteration", index))
    # add index
    index <- index + 1
    model_output <- do.call(DEoptim::DEoptim,
                            c(list(fn = model,
                                   control = list(trace = 1, parallelType = 1),
                                   lower = c(-1, -100, -1000, -5000),
                                   upper = c(1, 100, 1000, 5000),
                                   parameters_cr_2 = residence_par),
                              additional_arg))
    
    
    ## save output
    conv <- model_output$convergence
    par <- model_output$par
    fitness <- model_output$value
    output <- list(index, conv, mutant_par, residence_par, fitness)
    output_ls[[index]] <- output
    print(output)
    
    
    # exit loop IF
    ## previous run converged, strain 1 is fitter than strain 2, and that the fitness difference is minute
    if(conv == 0 && fitness > 0 && fitness < limit){
      print("Strain 1 is fitter than strain 2 but difference is minute.")
      stopCluster(cl)
      break
    }
    ## If repeat run after lower strain 1 fitness, break if continue to be 
    if(conv == 0 && lower_repeat == TRUE && fitness < limit){
      print("After repetition, fitness difference is still minute.")
      stopCluster(cl)
      break
    }
    # if at optimum
    if(conv == 0 && fitness == 0 ){
      print("Reached optimum. Break!")
      stopCluster(cl)
      break
    }
  }
  
  # final output
  return(output_ls)
}