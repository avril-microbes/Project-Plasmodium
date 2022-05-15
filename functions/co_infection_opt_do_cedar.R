# for Computecanada Cedar server
# wrapper around co-infection model that will perform "stationary" optimization,
# meaning that the optimal parameter for each iteration is assigned as both the
# the resident and mutant parameter. Mutant parameter is allowed to change
# as to maximize the difference between mutant and resident. DEoptim (GA) is used as optimization algorithm.

# By Avril Wang
# last updated 2022-05-10

co_infection_opt_do_cedar <- function(parameters_cr,  # preliminary parameter set
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
    } else{
      ### if strain 1 is fitter than strain 2, assign both strains to adopt the optimized par
      if(fitness < 0){
        residence_par <- par
        mutant_par <- par
        print("Strain 1 is fitter than strain 2.")
      }
      ### if strain 1 is less fit than strain 2, repeat optimization with par
      if(fitness > 0){
        if(index > 2){
          residence_par <- output_ls[[index - 1]][[4]] ### assign residence parameter to previous iteratio
          mutant_par <- par #### assign mutant to optimal par
        }
        if(index == 2){
          residence_par <- parameters_cr # special case. If first round of optimization fails, go back to initial parameter set
          mutant_par <- par
        } 
        lower_repeat <- TRUE
        print("Strain 1 is less fit than strain 2. Repeat optimization with 0.5x4")
      }
      ### special case: if fitness difference = 0, need skip out 
      if(fitness == 0){
        if(index > 2){
          residence_par <- output_ls[[index - 1]][[4]] ### assign residence parameter to previous iteratio
          mutant_par <- rep(2, 4) #### assign mutant to 2 x4 
        }
        if(index == 2){
          residence_par <- parameters_cr # special case. If first round of optimization fails, go back to initial parameter set
          mutant_par <- rep(2, 4)
        } 
        lower_repeat <- TRUE
        print("Strain 1 has same fitness than strain 2. Repeat optimization with 2x4")
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
                                   parameters_cr_2 = residence_par,
                                   cluster = cl),
                              additional_arg))
    
    
    ## save output
    par <- model_output$optim$bestmem
    fitness <- model_output$optim$bestval
    output <- list(index, mutant_par, residence_par, fitness)
    output_ls[[index]] <- output
    print(output)
    
    
    # exit loop IF
    ## previous run converged, strain 1 is fitter than strain 2, and that the fitness difference is minute
    if(fitness < 0 && fitness > (limit*-1)){
      print("Strain 1 is fitter than strain 2 but difference is minute.")

      break
    }
    ## If repeat run after lower strain 1 fitness, break if continue to be 
    if(lower_repeat == TRUE && fitness > (limit*-1)){
      print("After repetition, either strain 1 is still less fit or that fitness difference is minute.")

      break
    }
  }
  
  # final output
  return(output_ls)
}