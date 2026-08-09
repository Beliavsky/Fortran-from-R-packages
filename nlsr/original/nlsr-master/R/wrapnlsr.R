wrapnlsr <- function(formula=NULL, data=NULL, start=NULL, control = NULL, 
     trace = FALSE, subset=NULL, lower = -Inf, upper = Inf, weights = NULL, ...) {
    # A wrapper to call nlxb() and then call nls() with the
    # solution.  The calling sequence matches that of nlsmnq()
    if (is.null(formula)) stop("No formula")
  
  ## NOTE: The next few lines are suggestions by F. Miguez (2022-09-19) mod by JN 2022-9-24
    if (is.null(start)) {
    if (trace) cat("start is NULL -- try getInitial\n")
    issvalues <- try(getInitial(object = formula, data = data), silent = TRUE)
    if(!inherits(issvalues, "try-error")){
      if(is.null(control)){
        control <- list(japprox = "SSJac")
      }else{
        control[["jjaprox"]] <- "SSJac"
      }  
      start <- issvalues
    } 
    }
    if (is.null(start)) stop("No start")
    if (trace) {
       cat("control list for wrapnlsr:\n")
       print(control)
    }
    if (trace && ! is.null(weights)) { 
       cat("weights:"); print(str(weights)) 
    }
    nobs <- dim(data)[1]
    swts <- weights
    if (is.null(weights)) swts <- rep(1,nobs)
    sform<-formula
    if (trace) print(str(sform))
# Ensure we have controls.
    jfn<-NULL # default for nlxb
    ctrl <- nlsr.control()
    if(!missing(control)) {
	control <- as.list(control)
        if (! is.null(control$japprox)) jfn <- control$japprox
	ctrl[names(control)] <- control
    }
    ctrl$japprox<-jfn
    if (is.null(data)) stop("wrapnls() must have 'data' supplied")
    if(! is.null(subset)) {
##        stop("subset NOT used in wrapnlsr(). Please explicitly subset your data.")
       swts[-subset] <- 0 # handle subsetting through weights
       if (trace) { cat("swts:"); print(swts) }
       warning("subset effected by zero weights")
    }
    npar <- length(start)
    if (length(lower) < npar) {
        if (length(lower) == 1) 
            lower <- rep(lower, npar) else stop("lower bounds wrong length")
    }
    if (length(upper) < npar) {
        if (length(upper) == 1) 
            upper <- rep(upper, npar) else stop("upper bounds wrong length")
    }
    if (trace) {
        cat("wrapnlsr call with lower="); print(lower)
        cat("and upper="); print(upper)
    }
    first <- nlxb(formula=sform, start=start, trace = trace, data = data,
        lower = lower, upper = upper, weights=swts, control = ctrl, ...)
    # Should check this has worked, but ...
    if (trace) pshort(first)
    newstart <- first$coefficients
    if (trace) { nvec(newstart) }
    if (all(is.infinite(lower)) && all(is.infinite(upper))) {
        if (trace) cat("nls call with no bounds\n")
        second <- try(do.call("nls", list(formula=sform, data=data, start=newstart, 
		control=ctrl, algorithm=NULL, trace=trace, weights = swts, ...)))
     } else {
        if (trace) cat("Now try nls - bounded\n")
        second <- try(do.call("nls", list(formula=sform, data=data, start=newstart, 
		control=list(), algorithm="port", trace=trace, lower = lower, 
		upper = upper, weights = swts, ...)))

    }
    second
}

   nlsr <- wrapnlsr

## End of file ##

