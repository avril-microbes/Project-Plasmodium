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
  output_ls <- list()
  
  # do-while loop to continue to optimize until fitness difference between 
  # strain1 and strain2 is below limit
  repeat{
    ## for first iteration, both strain adopts initial parameter condition assigned
    if(index == 1){
      opt_parm_temp <- parameters_cr
    } else{ # reassign parameters_cr (invader starting strategy) to optimal value after first iteration
            # In the future, we could slightly alter this to prevent trapping into local optimum
      parameters_cr <- opt_parm_temp # only reassign when strain 1 is better
    }
    ## code to execute parallel LGBF-GS optimization
    output <- list()
    index <- index + 1
    print(paste("starting iteration", index))
    model_output <- do.call(optimParallel::optimParallel, c(list(par = parameters_cr, # competing parameter
                                                                 fn = model,
                                                                 control = list(trace = 6, fnscale = -1),
                                                                 parameters_cr_2 = opt_parm_temp), # reassign resident parameter to optimal one
                                                            additional_arg))
    
    ## save output
    opt_parm_conv <- model_output$convergence
    opt_parm_temp <- model_output$par
    opt_value_temp <- model_output$value
    output <- list(index, opt_parm_temp, opt_value_temp, opt_parm_conv)
    output_ls[[index]] <- output
    print(output)
    
    # exit loop IF
    if(opt_value_temp > 0 && opt_value_temp < limit) {break} # strain 1 is fitter than strain 2, but difference is small
    if(opt_value_temp <= 0 && 
       abs(output_ls[[index]][[3]] - output_ls[[index-1]][[3]]) = limit) 
      {break} # strain 1 is less fit than strain 2 but difference cannot be improved
  }
  opt_parm_conv <- model_output$convergence
  opt_parm_temp <- model_output$par
  opt_value_temp <- model_output$value
  index <- index + 1
  output <- list(index, opt_parm_temp, opt_value_temp, opt_parm_conv)
  output_ls[[index]] <- output
  return(output_ls)
}