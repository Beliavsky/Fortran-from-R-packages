## =============================================================================
## steady, -- solves the steady-state condition of
## ordinary differential equation systems
## has similar calling sequence as integration routines from package deSolve
## =============================================================================

steady  <- function (y, time=NULL, func, parms=NULL, method="stode", times = time, ...)  {

  if (!method %in% c("stode", "stodes","runsteady"))
    stop (" 'method' should be one of 'stode', 'stodes', 'runsteady'")   

  if (is.null(times)) {
    if (method %in% c("stode", "stodes")) 
      times <- 0
    else
      times <- c(0,Inf)  
  }  

  if (method=="stode")
    out <- stode(y,times,func,parms=parms,...)  else
  if (method=="stodes")
    out <- stodes(y,times,func,parms=parms,...) else
  if (method=="runsteady")
    out <- runsteady(y,time=times,func,parms=parms,...)
  class(out) <- c("steady","rootSolve","list")    # a steady-state 
  attr(out, "nspec") <- length(y)
  attr(out,"ynames") <- names(y)
  return(out)
}


                                                                                                                                                                                                            