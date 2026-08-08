#' Purify
#'
#' @param pmem_index  Index of population (particle) member you are you are updating
#' @param current_params Current parameter values for partcle (numeric vector)
#' @param current_weight  weights for current population
#' @param objFun function we want to minimize
#' @param n_particles number of particles
#' @param ... additional arguments for objective function
#' @noRd
#' 
Purify=function(pmem_index,
                current_params,
                current_weight,
                objFun,
                n_particles, ... ){

  # get statistics about particle
  weight_use = current_weight[pmem_index]
  params_use = current_params[pmem_index,]
  len_param_use = length(params_use)

  params_use = matrix(params_use,1,len_param_use)

  # recompute weight
  weight_proposal = objFun(params_use,...)

  if(is.na(weight_proposal))weight_proposal = Inf

  if(is.finite(weight_proposal)) {
    current_params[pmem_index,] = params_use
    current_weight[pmem_index] = weight_proposal
  }

  return(c(current_weight[pmem_index],current_params[pmem_index,]))

}
