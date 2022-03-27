# wrapper around co-infection model that will perform "stationary" optimization,
# meaning that the optimal parameter for each iteration is assigned as both the
# the resident and mutant parameter. Mutant parameter is allowed to change
# as to maximize the difference between mutant and resident. L-BFGS-B is used as optimization algorithm.

# By Avril Wang
# last updated 2022-03-24

co_infection_opt <- function(parameters_cr,  # preliminary parameter set
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
      if(fitness < 0 && conv != 0){
        if(index > 2){residence_par <- output_ls[[index - 1]][[4]]} ### assign residence parameter to previous iteration
        if(index == 2){residence_par <- parameters_cr} # special case. If first round of optimization fails, go back to initial parameter set
        mutant_par <- par
        print("Did not converge. Strain 1 is less fit than strain 2. Repeat optimization.")
      }
      
      ### if run converged but strain 1 is less fit then strain 2, try again with new initial mutant start point
      if(fitness < 0 && conv == 0){
        if(index > 2){residence_par <- output_ls[[index - 1]][[4]]}
        if(index == 2){residence_par <- parameters_cr}
        mutant_par <- par
        lower_repeat <- TRUE
        print("Converged. Strain 1 is less fit than strain 2. Repeat optimization.")
      }
    }
    ## code to execute parallel LGBF-GS optimization
    output <- list() # reset per run
    index <- index + 1
    print(paste("starting iteration", index))
    cl <- makeCluster(detectCores()); setDefaultCluster(cl = cl)
    model_output <- do.call(optimParallel::optimParallel, c(list(par = mutant_par, # competing parameter
                                                                 fn = model,
                                                                 control = list(trace = 6, fnscale = -1),
                                                                 parameters_cr_2 = residence_par), # reassign resident parameter to optimal one
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
    ## If repeat run after lower strain 1 fitness, if strain 1 is now fitter, break if difference is minute
    if(conv == 0 && fitness < 0 && lower_repeat == TRUE && fitness > 0 && fitness < limit){
      print("Strain 1 was previously less fit than strain 2. Now it is fitter and the difference is minute.")
      stopCluster(cl)
      break
      }
    # If repeat run and there is no improvements or even decrease in fitness, break
    if(conv == 0 && fitness < 0 && lower_repeat == TRUE && fitness < 0 && output_ls[[index]][[5]] <= output_ls[[index - 1]][[5]]){
      worse_then_previous <- TRUE
      print("Strain 1 is less fit than strain 2. Repeating optimization did not improve performance.")
      stopCluster(cl)
      break
      }
  }

  # final output
    return(output_ls)
}