#' SQG_DE_bin_1_curr
#'
#' @param pmem_index  Index of population (particle) member you are you are updating
#' @param current_params Current parameter values for partcle (numeric vector)
#' @param current_weight  weights for current population
#' @param objFun  function we want to minimize
#' @param n_particles number of particles
#' @param crossover_rate i.i.d. bernouli probability for updating a parameter
#' @param n_diff number of particles used to estimate the gradient
#' @param ... additional arguments for objective function
#' @noRd
#' 
SQG_DE_bin_1_curr=function(pmem_index,
                           current_params,
                           current_weight,
                           objFun,
                           step_size = .8,
                           jitter_size = 1e-6,
                           n_particles,
                           crossover_rate = 1,
                           n_diff, ... ){

  # sample parents
  parent_indices = sample(c(1:n_particles)[-pmem_index],size=2*n_diff,replace=F)

  # get statistics about sparticle
  weight_use = current_weight[pmem_index]
  # best_pmem_index = which.max(current_weight)
  params_use = current_params[pmem_index,]
  len_param_use = length(params_use)

  # use binomial to sample which params to update matching crossover rate frequency
  param_idices_bool = stats::rbinom(len_param_use, prob = crossover_rate, size=1)

  # if no pars selected, randomly sample 1 parameter
  if(all(param_idices_bool==0)){
    param_idices_bool[sample(x=1:len_param_use,size=1)] = 1
  }

  # indices of parameters to be updated
  param_indices = seq(1,len_param_use,by=1)[as.logical(param_idices_bool)]


  # get rid of this for loop later
  vec_diff_sum = grad_approx = numeric(length(param_indices))
  for(d in 1:n_diff){
    # difference in vector pairs
    vec_diff_temp = (current_params[parent_indices[d],param_indices] -
                     current_params[parent_indices[d+n_diff],param_indices])

    # sum up all vector differences for normalization step
    vec_diff_sum = vec_diff_sum+vec_diff_temp

    # difference in function values
    weight_diff_temp = (current_weight[parent_indices[d]] -
                        current_weight[parent_indices[d+n_diff]])

    # sum the approximate gradient vectors
    grad_approx = grad_approx + vec_diff_temp*(weight_diff_temp/sqrt(sum(vec_diff_temp^2)))
  }

  # calculate normalization factor for algorithm self scaling
  psi_num = sqrt(sum(vec_diff_sum^2))*(1/n_diff)
  psi_den = sqrt(sum(grad_approx^2))
  psi=psi_num/psi_den


  if(all(is.finite(grad_approx)) & is.finite(psi) & (psi>0)){
    # mate parents for proposal
    params_use[param_indices] = current_params[pmem_index, param_indices] -
      step_size*psi*(grad_approx) + # move in the direction against the gradient
      stats::runif(len_param_use, -jitter_size, jitter_size) # a little noise
  }
  params_use = matrix(params_use,1,len_param_use)

  weight_proposal = NA
  # get weight
  if(all(is.finite(params_use)))weight_proposal = objFun(params_use,...)
  if(is.na(weight_proposal))weight_proposal = Inf

  #negative greedy acceptance rule
  if(weight_proposal < weight_use) {
    current_params[pmem_index,] = params_use
    current_weight[pmem_index] = weight_proposal
  }

  return(c(current_weight[pmem_index],current_params[pmem_index,]))

}
