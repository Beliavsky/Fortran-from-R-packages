#' optim_SQGDE
#'
#' @description  Runs Stochastic Quasi-Gradient Differential Evolution (SQG-DE; Sala, Baldanzini, and Pierini, 2018) to minimize an objective function f(x). To maximize a function f(x), simply pass g(x)=-f(x) to ObjFun argument.
#' @param ObjFun A scalar-returning function to minimize whose first argument is a real-valued n_params-dimensional vector.
#' @param control_params control parameters for SQG-DE algo. see \code{\link{GetAlgoParams}} function documentation for more details. The only argument you NEED to pass here is n_params.
#' @param ... additional arguments to pass ObjFun.
#' @return list containing solution and it's corresponding weight (i.e. f(solution)).
#' @export
#' @md
#' @examples
#' ##############
#' # Maximum Likelihood Example
#' ##############
#'
#' # simulate from model
#' dataExample=matrix(rnorm(1000,c(-1,1,0,1),c(1,1,1,1)),ncol=4,byrow = TRUE)
#'
#' # list parameter names
#' param_names_example=c("mu_1","mu_2","mu_3","mu_4")
#'
#' # negative log likelihood
#' ExampleObjFun=function(x,data,param_names){
#'   out=0
#'
#'   names(x) <- param_names
#'
#'   # log likelihoods
#'   out=out+sum(dnorm(data[,1],x["mu_1"],sd=1,log=TRUE))
#'   out=out+sum(dnorm(data[,2],x["mu_2"],sd=1,log=TRUE))
#'   out=out+sum(dnorm(data[,3],x["mu_3"],sd=1,log=TRUE))
#'   out=out+sum(dnorm(data[,4],x["mu_4"],sd=1,log=TRUE))
#'
#'   return(out*-1)
#' }
#'
#' ########################
#' # run optimization
#' out <- optim_SQGDE(ObjFun = ExampleObjFun,
#'                    control_params = GetAlgoParams(n_params = length(param_names_example),
#'                                              n_iter = 250,
#'                                               n_particles = 12,
#'                                               n_diff = 2,
#'                                               return_trace = TRUE),
#'                    data = dataExample,
#'                    param_names = param_names_example)
#' old_par <- par() # save graphic state for user
#' # plot particle trajectory
#'
#' par(mfrow=c(2,2))
#' matplot(out$particles_trace[,,1],type='l')
#' matplot(out$particles_trace[,,2],type='l')
#' matplot(out$particles_trace[,,3],type='l')
#' matplot(out$particles_trace[,,4],type='l')
#'
#' #SQG DE solution
#' out$solution
#'
#' #analytic solution
#' apply(dataExample, 2, mean)
#'
#' par(old_par) # restore user graphic state
#'

optim_SQGDE = function(ObjFun, control_params = GetAlgoParams(), ...){

  # create memory structures for storing particle trajectories
  particles = array(NA,
                    dim = c(control_params$n_iters_per_particle,
                            control_params$n_particles,
                            control_params$n_params))
  weights = matrix(Inf,
                   nrow = control_params$n_iters_per_particle,
                   ncol = control_params$n_particles)

  # pop initialization
  message('initalizing population...')
  for(pmem_index in 1:control_params$n_particles){
    count = 0 # establish a count variable to avoid infinite run time
    while(weights[1,pmem_index]==Inf) {
      particles[1, pmem_index, ] = stats::rnorm(control_params$n_params,
                                                control_params$init_center,
                                                control_params$init_sd)

      weights[1, pmem_index] = ObjFun(particles[1, pmem_index, ], ...)

      # catcha NA's and Infinity and assign worst possible value
      if(!is.finite(weights[1, pmem_index])){
        weights[1, pmem_index] = Inf
      }
      count = count + 1
      if(count>control_params$give_up_init){
        stop('population initialization failed.
        inspect objective function or change init_center/init_sd to sample more
             likely parameter values')
      }
    }
    message(paste0(pmem_index, " / ", control_params$n_particles))
  }
  message('population initialization complete  :)')

  # assign adaption scheme
  if(control_params$adapt_scheme=='rand'){
    AdaptSQGDE = SQG_DE_bin_1_rand
  }
  if(control_params$adapt_scheme=='best'){
    AdaptSQGDE = SQG_DE_bin_1_best
  }
  if(control_params$adapt_scheme=='current'){
    AdaptSQGDE = SQG_DE_bin_1_curr
  }

  # cluster initialization
  if(!control_params$parallel_type=='none'){

    message(paste0("initalizing ",
                 control_params$parallel_type, " cluser with ",
                 control_params$n_cores_use, " cores"))

    doParallel::registerDoParallel(control_params$n_cores_use)
    cl_use = parallel::makeCluster(control_params$n_cores_use,
                                   type = control_params$parallel_type)
  }

  message("running SQG-DE...")



  iter_idx=1
  converge_test_passed=FALSE
  for(iter in 1:control_params$n_iter){

    if(control_params$parallel_type=='none'){
      # adapt particles using SQG DE sequentially
      temp=matrix(unlist(lapply(1:control_params$n_particles, AdaptSQGDE,
                                current_params = particles[iter_idx, , ],   # current parameter values (numeric matrix)
                                current_weight = weights[iter_idx, ],  # corresponding weights (numeric vector)
                                objFun = ObjFun,  # objective function (returns scalar)
                                step_size = control_params$step_size,
                                jitter_size = control_params$jitter_size,
                                n_particles = control_params$n_particles,
                                crossover_rate = control_params$crossover_rate,
                                n_diff = control_params$n_diff, ...)),
                  control_params$n_particles,
                  control_params$n_params+1, byrow=TRUE)
    } else {
      # adapt particles using SQG DE in parallel
      temp=matrix(unlist(parallel::parLapplyLB(cl_use, 1:control_params$n_particles, AdaptSQGDE,
                                             current_params = particles[iter_idx, , ],   # current parameter values (numeric matrix)
                                             current_weight = weights[iter_idx, ],  # corresponding weight (numeric vector)
                                             objFun = ObjFun,  # function we want to minimize (returns scalar)
                                             step_size = control_params$step_size,
                                             jitter_size = control_params$jitter_size,
                                             n_particles = control_params$n_particles,
                                             crossover_rate = control_params$crossover_rate,
                                             n_diff = control_params$n_diff, ...)),
                  control_params$n_particles,
                  control_params$n_params+1, byrow=TRUE)

    }
    # update particle after adaption
    weights[iter_idx, ] = temp[, 1]
    particles[iter_idx, , ] = temp[, 2:(control_params$n_params+1)]
    # carry over particles for next iteration
    if(iter<control_params$n_iter){
      weights[iter_idx+1, ] = temp[, 1]
      particles[iter_idx+1, , ] = temp[, 2:(control_params$n_params+1)]
    }

    #####################
    ####### purify
    #####################
    if(iter%%control_params$purify==0){

      if(control_params$parallel_type=='none'){
        temp=matrix(unlist(lapply(1:control_params$n_particles, Purify,
                                  current_params = particles[iter_idx, , ],   # current parameter values (numeric  matrix)
                                  current_weight = weights[iter_idx, ],  # corresponding weights (numeric vector)
                                  objFun = ObjFun,  # objective function (returns scalar)
                                  ...)),
                    nrow = control_params$n_particles,
                    ncol = control_params$n_params+1, byrow=TRUE)
      } else {
        temp=matrix(unlist(parallel::parLapplyLB(cl_use, 1:control_params$n_particles, Purify,
                                               current_params = particles[iter_idx, , ],   # current parameter values (numeric matrix)
                                               current_weight = weights[iter_idx, ],  # corresponding weights (numeric vector)
                                               objFun = ObjFun,  # objective function (returns scalar)
                                               ...)),
                    control_params$n_particles,
                    control_params$n_params+1, byrow=TRUE)
      }

      # update particle after adaption
      weights[iter_idx, ] = temp[, 1]
      particles[iter_idx, , ] = temp[, 2:(control_params$n_params+1)]
      # carry over particles for next iteration
      if(iter<control_params$n_iter){
        weights[iter_idx+1, ] = temp[, 1]
        particles[iter_idx+1, , ] = temp[, 2:(control_params$n_params+1)]
      }
    }

    if(iter %% control_params$stop_check==0){
      # assign convergence method

      if(control_params$converge_crit=='percent'){
        percent_improve=(1-stats::median(weights[iter_idx, ])/stats::median(weights[iter_idx-control_params$stop_check+1, ]))*100
        if(percent_improve<(control_params$stop_tol)){
          message("Convergence criterion met. Stopping optimization early")
          converge_test_passed=TRUE
          break
        }
      }
      if(control_params$converge_crit=='stdev'){
        weight_sd=stats::sd(weights[iter_idx:(iter_idx-control_params$stop_check+1), ])
        if(weight_sd<(control_params$stop_tol)){
          message("Convergence criterion met. Stopping optimization early")
          converge_test_passed=TRUE
          break
        }
      }
    }



    if(iter%%100==0){
      message(paste0('iter ', iter, '/', control_params$n_iter))
    }
    if(iter%%control_params$thin==0 & !(iter==control_params$n_iter)){
      iter_idx = iter_idx+1
    }

  }
  message(paste0('run complete!'))
  # cluster stop
  if(!control_params$parallel_type=='none'){
    parallel::stopCluster(cl = cl_use)
  }
  minIdx = which.min(weights[iter_idx, ])
  minEst = particles[iter_idx, minIdx, ]

  if(control_params$return_trace==TRUE){
    return(list('solution' = minEst,
                'weight' = weights[iter_idx, minIdx],
                'particles_trace' = particles,
                'weights_trace' = weights,
                'converged' = converge_test_passed))
  } else {
    return(list('solution' = minEst,
                'weight' = weights[iter_idx, minIdx],
                'converged' = converge_test_passed))
  }
}
