## =============================================================================
## steady.1D -- solves the steady-state condition of
## ordinary differential equation systems resulting from
## multi-component 1-D PDE models
## has a similar calling sequence as ode.1D from package deSolve
## =============================================================================

steady.1D    <- function (y, time=NULL, func, parms=NULL, nspec = NULL,
                       dimens = NULL, names = NULL, method = "stode", 
                       jactype = NULL, cyclicBnd = NULL, bandwidth = 1, times = time, ...) {
  if (any(!is.na(pmatch(names(list(...)), "jacfunc")))) 
    stop ("cannot run steady.1D with jacfunc specified - remove jacfunc from call list")
  if (is.null(dimens) && is.null(nspec)) 
    stop ("cannot run steady.1D: either nspec or dimens should be specified")
  N     <- length(y)
  if (is.null(nspec)  )
    nspec <- N/dimens
  else if (N%%nspec != 0) 
    stop("cannot run steady.1D: nspec is not an integer fraction of number of state variables")
  if (! is.null(names) && length(names) != nspec)
    stop("length of 'names' should equal 'nspec'")
  if (! is.null (cyclicBnd)) method <- "stodes"
  if (method == "stodes" && bandwidth != 1)
    stop ("cannot combine 'method = stodes' with 'bandwidth' not = 1") 
  if (!method %in% c("stode", "stodes","runsteady"))
    stop (" 'method' should be one of 'stode', 'stodes', 'runsteady'")   

  if (is.null(times))
    times <- 0
  if (method == "runsteady" & length(times) == 1)
      times <- c(times, Inf)  
        
  if (nspec == 1 & method == "stode") {
    out <- steady.band(y, times, func, parms, nspec, 
      bandup = nspec*bandwidth, banddown = nspec*bandwidth,...)
  } else if (method=="stodes") {
    if (is.null(jactype)) jactype <- "1D"
    dimens <- N/nspec
    Bnd <- 0
    if (! is.null(cyclicBnd)) {
      if (max(cyclicBnd) > 1 )
        stop ("cannot run steady.1D: cyclicBnd should be NULL or a value not exceeding 1")
      Bnd <-1   
    }

    if (jactype == "1D")
      out <- stodes(y=y,time=times,func=func,parms=parms,
                  nnz=c(nspec,dimens,Bnd),sparsetype=jactype,...)
    else
      out <- stodes(y=y,time=times,func=func,parms=parms,
                    sparsetype=jactype, ...)
         
  } else if (is.character(func)) {
    if (is.null(jactype)) jactype <- "1Dint"
    ii    <- as.vector(t(matrix(ncol=nspec,1:N)))   # from ordering per slice -> per spec
    if (bandwidth != 1)
      stop ("cannot combine DLL with 'bandwidth' not = 1") 
    if (method == "stode")
      out <- stode(y=y[ii],time=times,func=func,parms=parms,
                jactype=jactype,bandup=nspec,banddown=N/nspec,...)                    
    else if (method == "runsteady")
      out <- runsteady (y=y[ii],times = times,func=func,parms=parms,
                jactype=jactype,bandup=nspec,banddown=N/nspec,...)                    
   
    else
      stop ("cannot run steady.1D: method should be 'stode' or 'runsteady' if func is a DLL")
            
    out[[1]][ii] <- out[[1]]
  } else {

  # internal function #
    if (is.null(jactype)) jactype <- "bandint"

    bmodel <- function (time,state,pars,model,...) {
      Modconc <-  model(time,state[ij],pars,...)   # ij: reorder state variables
      c(list(Modconc[[1]][ii]),Modconc[-1])        # ii: reorder rate of change
    }

    ii    <- as.vector(t(matrix(ncol=nspec,1:N)))   # from ordering per slice -> per spec
    ij    <- as.vector(t(matrix(nrow=nspec,1:N)))   # from ordering per spec -> per slice
    
    bmod  <- function(time,state,pars,...) {
      bmodel(time,state,pars,func,...)
    }

    if (method=="stode")
      out <- stode(y[ii],times,func=bmod,parms=parms,
                bandup=nspec*bandwidth,banddown=nspec*bandwidth,
                jactype=jactype,...)    else
      out <- runsteady(y[ii],times=times,func=bmod,parms=parms,
                bandup=nspec*bandwidth,banddown=nspec*bandwidth,
                jactype=jactype,...)

    out[[1]][ii] <- out[[1]]
  }
  class(out) <- c("steady1D","rootSolve","list")    # a steady-state 
    if (! is.null(names)) {
      out[[1]] <- matrix(ncol=nspec,data=out[[1]])
      colnames(out[[1]]) <- names
    }
  if (is.null(dimens))  dimens <- N/nspec
  attr(out,"dimens") <- dimens
  attr(out, "nspec") <- nspec
  attr(out,"ynames") <- names
  return(out)
}
