#' GetAlgoParams
#' @description Get control parameters for optim_SQGDE function.
#' @param n_params The number of parameters estimated/optimized, this integer value NEEDS to be specified.
#' @param n_particles The number of particles (population size), 3*n_params is the default value.
#' @param n_iter The number of iterations to run the algorithm, 1000 is default.
#' @param n_diff The number of mutually exclusive vector pairs to stochastically approximate the gradient.
#' @param crossover_rate A numeric scalar on the interval (0,1]. Determines the probability a parameter on a chain is updated on a given crossover step, sampled from a Bernoulli distribution. The default value is 1.
#' @param init_sd A positive scalar or n_params-dimensional numeric vector, determines the standard deviation of the Gaussian initialization distribution. The default value is 0.01.
#' @param init_center A scalar or n_params-dimensional numeric vector, determines the mean of the Gaussian initialization distribution. The default value is 0.
#' @param n_cores_use An integer specifying the number of cores used when using parallelization. The default value is 1.
#' @param step_size A positive scalar, jump size or "F" in the DE crossover step notation. The default value is 2.38/sqrt(2*n_params).
#' @param jitter_size A positive scalar that determines the jitter (noise) size. Noise is added during adaption step from Uniform(-jitter_size,jitter_size) distribution. 1e-6 is the default value. Set to 0 to turn off jitter.
#' @param parallel_type A string specifying parallelization type. 'none','FORK', or 'PSOCK' are valid values. 'none' is default value. 'FORK' does not work with Windows OS.
#' @param return_trace A boolean, if true, the function returns particle trajectories. This is helpful for assessing convergence or debugging model code. The trace will be an iteration/thin $x$ n_particles $x$ n_params array containing parameter values and an iteration/thin $x$ n_particles array containing particle weights.
#' @param thin A positive integer. Only every 'thin'-th iteration will be stored in memory. The default value is 1. Increasing thin will reduce the memory required when running the algorithim for longer.
#' @param purify A positive integer. On every 'purify'-th iteration the particle weights are recomputed. This is useful if the objective function is stochastic/noisy. If the objective function is deterministic, this computation is redundant. Purify is set to Inf by default, disabling it.
#' @param adapt_scheme A string that must be 'rand','current', or 'best' that determines the DE adaption scheme/strategy. 'rand' uses rand/1/bin DE-like scheme where a random particle and the particle-based quasi-gradient approximation are used to generate proposal updates for a given particle. 'current' uses current/1/bin, and 'best' uses best/1/bin which follow an analogous adaption scheme to rand. 'rand' is the default value.
#' @param give_up_init An integer for how many failed initialization attempts before stopping the optimization routine. 100 is the default value.
#' @param stop_check An integer for how often to check the convergence criterion. The default is 10 iterations.
#' @param stop_tol A convergence metric must be less than value to be labeled as converged. The default is 1e-4.
#' @param converge_crit A string denoting the convergence metric used, valid metrics are 'stdev' (standard deviation of population weight in the last stop_check iterations) and 'percent' (percent improvement in median particle weight in the last stop_check iterations). 'stdev' is the default.
#' @return A list of control parameters for the optim_SQGDE function.
#' @export
GetAlgoParams = function(n_params,
                         n_particles = NULL,
                         n_diff = 2,
                         n_iter = 1000,
                         init_sd = 0.01,
                         init_center = 0,
                         n_cores_use = 1,
                         step_size = NULL,
                         jitter_size = 1e-6,
                         crossover_rate = 1,
                         parallel_type = 'none',
                         return_trace = FALSE,
                         thin = 1,
                         purify = Inf,
                         adapt_scheme = NULL,
                         give_up_init = 100,
                         stop_check = 10,
                         stop_tol = 1e-4,
                         converge_crit = 'stdev'){
  # n_params
  ### catch errors
  n_params = as.integer(n_params)
  if(any(!is.finite(n_params))){
    stop('ERROR: n_params is not finite')
  }  else if( n_params<1 | length(n_params)>1){
    stop('ERROR: n_params must be a postitive integer scalar')
  }

  # n_particles
  ### if null assign default value
  if(is.null(n_particles)){
    n_particles = max(3*n_params,4)
  }
  ### catch errors
  n_particles = as.integer(n_particles)
  if(any(!is.finite(n_particles))){
    stop('ERROR: n_particles is not finite')
  } else if( n_particles<4 | length(n_particles)>1){
    stop('ERROR: n_particles must be a postitive integer scalar, and atleast 4')
  }

  # n_iter
  ### if null assign default value
  if(is.null(n_iter)){
    n_iter = 1000
  }
  ### catch errors
  n_iter = as.integer(n_iter)
  if(any(!is.finite(n_iter))){
    stop('ERROR: n_iter is not finite')
  } else if( n_iter<4 | length(n_iter)>1){
    stop('ERROR: n_iter must be a postitive integer scalar, and atleast 4')
  }

  # init_sd
  init_sd = as.numeric(init_sd)
  if(any(!is.finite(init_sd))){
    stop('ERROR: init_sd is not finite')
  } else if(any(init_sd<= 0 | is.complex(init_sd))){
    stop('ERROR: init_sd must be positive and real-valued')
  } else if(!(length(init_sd) == 1 | length(init_sd) == n_params)){
    stop('ERROR: init_sd vector length must be 1 or n_params')
  }

  # init_center
  init_center = as.numeric(init_center)
  if(any(!is.finite(init_center))){
    stop('ERROR: init_center is not finite')
  } else if(any(is.complex(init_center))){
    stop('ERROR: init_center must be real valued')
  } else if(!(length(init_center) == 1 | length(init_center) == n_params)){
    stop('ERROR: init_center vector length must be 1 or n_params')
  }

  # n_cores_use
  ### assign NULL value default
  if(is.null(n_cores_use)){
    n_cores_use = 1
  }
  ### catch any errors
  n_cores_use = as.integer(n_cores_use)
  if(any(!is.finite(n_cores_use))){
    stop('ERROR: n_cores_use is not finite')
  } else if( n_cores_use<1 | length(n_cores_use)>1){
    stop('ERROR: n_cores_use must be a postitive integer scalar, and atleast 1')
  }


  # step_size
  ### assign NULL value default
  if(is.null(step_size)){
    step_size = 2.38/sqrt(2*n_params) # step size recommend in ter braak's 2006 paper
  }
  ### catch any errors
  if(any(!is.finite(step_size))){
    stop('ERROR: step_size is not finite')
  } else if(any(step_size<= 0 | is.complex(step_size))){
    stop('ERROR: step_size must be positive and real-valued')
  } else if(!(length(step_size) == 1)){
    stop('ERROR: step_size vector length must be 1 ')
  }

  #jitter_size
  ### assign NULL value default
  if(is.null(jitter_size)){
    jitter_size = 1e-6
  }
  ### catch any errors
  if(any(!is.finite(jitter_size))){
    stop('ERROR: jitter_size is not finite')
  } else if(any(jitter_size<= 0 | is.complex(jitter_size))){
    stop('ERROR: jitter_size must be positive and real-valued')
  } else if(!(length(jitter_size) == 1)){
    stop('ERROR: jitter_size vector length must be 1 ')
  }

  # crossover_rate
  ### if null assign default value
  if(any(is.null(crossover_rate))){
    crossover_rate = 1
  }
  ### catch errors
  crossover_rate = as.numeric(crossover_rate)
  if(any(!is.finite(crossover_rate))){
    stop('ERROR: crossover_rate is not finite')
  } else if(any(crossover_rate>1) | any(crossover_rate<= 0) | length(crossover_rate)>1){
    stop('ERROR: crossover_rate must be a numeric scalar on the interval (0,1]')
  } else if(is.complex(crossover_rate)){
    stop('ERROR: crossover_rate cannot be complex')
  }

  #parallel_type
  validParType = c('none','FORK','PSOCK')
  ### assign NULL value default
  if(is.null(parallel_type)){
    parallel_type = 'none'
  }
  ### catch any errors
  if(!parallel_type %in% validParType){
    stop(paste('ERROR: invalid parallel_type.'))
  }

  #converge_crit
  validConType = c('stdev','percent')
  ### assign NULL value default
  if(is.null(converge_crit)){
    converge_crit = 'stdev'
  }
  ### catch any errors
  if(!converge_crit %in% validConType){
    stop(paste('ERROR: invalid converge_crit.'))
  }

  validAdaptType = c('rand','current','best')
  ### assign NULL value default
  if(is.null(adapt_scheme)){
    adapt_scheme = 'rand'
  }
  ### catch any errors
  if(!adapt_scheme %in% validAdaptType){
    stop(paste('ERROR: invalid adaption scheme.'))
  }

  # thin
  ### if null assign default value
  if(is.null(thin)){
    thin = 1
  }
  ### catch errors
  thin = as.integer(thin)
  if(any(!is.finite(thin))){
    stop('ERROR: thin is not finite')
  } else if(any(thin<1) | length(thin)>1){
    stop('ERROR: thin must be a scalar postive integer')
  }

  #number iters stored per particle
  n_iters_per_particle = floor((n_iter)/thin)
  ### catch errors
  if(n_iters_per_particle<1 | (!is.finite(n_iters_per_particle))){
    stop('ERROR: number of stored particle value is negative or non finite.
         n_iters_per_particle = floor((n_iter)/thin)')
  }

  # purify
  if(is.null(purify)){
    purify = Inf
  } else if(is.finite(purify)){
    purify = as.integer(purify)
  }
  ### catch errors
  if(any(purify < 1) | (length(purify)>1)){
    stop('ERROR: purify must be a positive integer or Inf')
  }

  ### catch errors
  n_diff = as.integer(n_diff)
  if(any(!is.finite(n_diff))){
    stop('ERROR: n_diff is not finite')
  } else if(any(n_diff<1) | length(n_diff)>1){
    stop('ERROR: n_diff must be a scalar postive integer')
  } else if(n_diff>(n_particles/2)){
    stop('ERROR: n_diff cannot exceed n_particles/2')
  }

  ##################
  # give_up_init
  if(any(is.null(give_up_init)) | any(is.na(give_up_init))){
    give_up_init = 100
  }
  give_up_int=round(give_up_init)
  ### catch errors
  if(any(give_up_init < 1) | (length(give_up_init)>1)){
    stop('ERROR: give_up_init must be a scalar positive integer')
  }

  ##################
  # stop_check
  if(any(is.null(stop_check)) | any(is.na(stop_check))){
    stop_check = 10
  }
  stop_check=round(stop_check)
  ### catch errors
  if(any(stop_check < 2) | (length(stop_check)>1)){
    stop('ERROR: stop_check must be a scalar positive integer and greater than 2')
  }


  ##################
  # stop_tol
  if(any(is.null(stop_tol)) | any(is.na(stop_tol))){
    stop_tol = 1e-4
  }
  ### catch errors
  if(any(stop_tol < 0) | (length(stop_tol)>1)){
    stop('ERROR: stop_tol must be a scalar positive')
  }

  out = list('n_params' = n_params,
             'n_particles' = n_particles,
             'n_iter' = n_iter,
             'init_sd' = init_sd,
             'init_center' = init_center,
             'n_cores_use' = n_cores_use,
             'step_size' = step_size,
             'crossover_rate' = crossover_rate,
             'jitter_size' = jitter_size,
             'parallel_type' = parallel_type,
             'thin' = thin,
             'purify' = purify,
             'n_iters_per_particle' = n_iters_per_particle,
             'return_trace' = return_trace,
             'n_diff' = n_diff,
             'adapt_scheme' = adapt_scheme,
             'give_up_init'= give_up_init,
             'stop_tol' = stop_tol,
             'stop_check' = stop_check,
             'converge_crit' = converge_crit)

  return(out)
}
