# ================================
# Nuisance Parameters
# Tool Function
# ================================

# Returning nuisance parameters data
NuisanceData <- function(nuisPars){
  
  # Evaluation of nuisance parameters
  if(is.null(nuisPars))
    nuisPars <- nuisParsFn() else{
      if(!is.null(nuisPars$mu)){
        nuis.mu <- nuisPars$mu} else{
          nuis.mu <- 0.01}
      if(!is.null(nuisPars$sd)){
        nuis.sd <- nuisPars$sd} else{
          nuis.sd <- 0.05}
      if(!is.null(nuisPars$c)){
        nuis.c <- nuisPars$c} else{
          nuis.c <- 0}
      if(!is.null(nuisPars$alpha)){
        nuis.alpha <- nuisPars$alpha} else{
          nuis.alpha <- 0.1}
      if(!is.null(nuisPars$beta)){
        nuis.beta <- nuisPars$beta} else{
          nuis.beta <- 0.1}
      nuisPars <- nuisParsFn(nuis.mu, nuis.sd, nuis.c, nuis.alpha, nuis.beta)
    }
  
  # Returning nuisance parameters
  return(nuisPars)
}