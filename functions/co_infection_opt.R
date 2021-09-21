# wrapper around co-infection model that will perform "stationary" optimization,
# meaning that the optimal parameter for each iteration is assigned as both the
# the resident and invador parameter.

# By Avril Wang

co_infection_opt <- function(parameters_cr,  # preliminary parameter set
                             limit, # minimum fitness difference between competing strains to break
                             model, # infection model
                             ...){ # additional parameters to be inserted into the infection model
  force(parameters_cr)
  force(limit)
  force(model)
  additional_arg <- list(...) # unpack for do.call
  
  # assign iteration index
  index <- 1
  
  # do-while loop to continue to optimize until fitness difference between 
  # strain1 and strain2 is below limit
  repeat{
    ## for first iteration, both strain adopts initial parameter condition assigned
    if(index == 1){
      opt_parm_temp <- parameters_cr
    } else{ # reassign parameters_cr (invader starting strategy) to optimal value after first iteration
            # In the future, we could slightly alter this to prevent trapping into local optimum
      parameters_cr <- opt_parm_temp
    }
    ## code to execute parallel LGBF-GS optimization
    model_output <- do.call(optimParallel::optimParallel, c(list(par = parameters_cr,
                                                                 fn = model,
                                                                 control = list(trace = 6, fnscale = -1),
                                                                 parameters_cr_2 = opt_parm_temp), # reassign resident parameter to optimal one
                                                            additional_arg))
    
    ## save output
    opt_parm_temp <- model_output$par
    opt_value_temp <- model_output$value
    index <- index + 1
    
    # exit loop if limit is reached
    if(opt_value_temp < limit) {
      return(model_output)
      break
    }
  }
}