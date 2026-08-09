model2rjfun <- function(modelformula, pvec, data = NULL, jacobian = TRUE, 
                          testresult = TRUE, ...) {
    # add dots 161127, but they have NOT been tested

    stopifnot(inherits(modelformula, "formula"))

    if (length(modelformula) == 2) { ## Explain this in vignette!?
        residexpr <- modelformula[[2]]
    } else if (length(modelformula) == 3) {
        residexpr <- call("-", modelformula[[3]], modelformula[[2]])
    } else stop("Unrecognized formula")
    
    if (is.null(names(pvec))){
          cat("apparently no names in pvec\n")
	  names(pvec) <- paste0("p_", seq_along(pvec))
    }
    
    if (jacobian)
        residexpr <- deriv(residexpr, names(pvec))
## SHOULD WE TRY:?
##	residexpr <- fnDeriv(residexpr, names(pvec))
##  data <- list(...)	
    if (is.null(data)) {
	  data <- environment(modelformula) # this will handle variables in the parent frame
    }  ## Don't yet handle variable in dot args. But no dot args here. Why?
    else if (is.list(data))
	data <- list2env(data, parent = environment(modelformula))
    else if (!is.environment(data))
        	stop("'data' must be a dataframe, list, or environment")
    rjfun <- function(prm) {
        if (is.null(names(prm))) 
	    names(prm) <- names(pvec)
	  localdata <- list2env(as.list(prm), parent = data)
	  eval(residexpr, envir = localdata)  
      # Saves Jacobian matrix as "gradient" attribute (consistent with deriv())
    }
    
    if (testresult) {
	resids <- rjfun(pvec)
	if (any(bad <- !is.finite(resids))) 
	    stop("residuals contain ", unique(resids[bad]))
	if (jacobian && any(bad <- !is.finite(attr(resids, "gradient"))))
	    stop("Jacobian contains ", unique(attr(resids, "gradient")[bad]))
	rm(resids, bad)  # Don't want to capture these in the environment of rjfun
    }
    rjfun
}

model2ssgrfun <- function(modelformula, pvec, data = NULL, gradient = TRUE, 
                          testresult = TRUE, ...) {
    rjfun <- model2rjfun(modelformula, pvec, data = data, jacobian = gradient, 
                         testresult = testresult, ...)
			
    function(prm) {
    	resids <- rjfun(prm)
	    ss <- as.numeric(crossprod(resids))
	    if (gradient) {
	        jacval <- attr(resids, "gradient")
	        grval <- 2*as.numeric(crossprod(jacval, resids))
	        attr(ss, "gradient") <- grval
	    }
	    attr(ss, "resids") <- resids
	    ss
    }
}

modelexpr <- function(fun) {
    env <- environment(fun)
    if (exists("rjfun", env))
	env <- environment(env$rjfun)
    env$residexpr
}


SSmod2rjfun <- function(modelformula, pvec, data = NULL, jacobian = TRUE, 
                        testresult = TRUE, ...) {
  stopifnot(inherits(modelformula, "formula"))
  
  if (length(modelformula) == 2) { ## Explain this in vignette!?
    residexpr <- modelformula[[2]]
    rhs <- residexpr # for SS case
  } else if (length(modelformula) == 3) {
    residexpr <- call("-", modelformula[[3]], modelformula[[2]])
    rhs <- modelformula[[3]] # for SS case
  } else stop("Unrecognized formula")

  if (is.null(names(pvec))){
    cat("apparently no names in pvec\n")
    names(pvec) <- paste0("p_", seq_along(pvec))
  }
  
  if (jacobian) {
    if (is.null(data)) {
      data <- environment(modelformula) # this will handle variables in the parent frame
    }  ## Don't yet handle variable in dot args. But no dot args here.
    else if (is.list(data))
      data <- list2env(data, parent = environment(modelformula))
    else if (!is.environment(data))
      stop("'data' must be a dataframe, list, or environment")
    rjfun <- function(prm) {
      if (is.null(names(prm))) 
        names(prm) <- names(pvec)
      localdata <- list2env(as.list(prm), parent = data)
      eval(residexpr, envir = localdata)  
      # Saves Jacobian matrix as "gradient" attribute (consistent with deriv())
    }
  }  

  if (testresult) {
    resids <- rjfun(pvec)
    if (any(bad <- !is.finite(resids))) 
      stop("residuals contain ", unique(resids[bad]))
    if (jacobian && any(bad <- !is.finite(attr(resids, "gradient"))))
      stop("Jacobian contains ", unique(attr(resids, "gradient")[bad]))
    rm(resids, bad)  # Don't want to capture these in the environment of rjfun
  }
  rjfun
}

