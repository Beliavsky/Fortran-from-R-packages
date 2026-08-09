rhs <- function(formula) {
  if (length(formula) == 2)
    formula[[2]]
  else
    formula[[3]]
}

nlxb <- function(formula, data=parent.frame(), start, trace = FALSE,  lower = NULL,
                 upper = NULL, weights=NULL, control=list(), ...) {
  
    weights0 <- weights
    
    if("subset" %in% names(list(...)))
        stop("subset NOT used in nlxb(). Please explicitly subset your data.")
    ctrl <- nlsr.control(control)
    if (missing(start)) stop("A start vector MUST be provided")
    pnames <- names(start)
    start <- as.numeric(start) # ensure we convert (e.g., if matrix)
    names(start) <- pnames ## as.numeric strips names, so this is needed
    # bounds
    npar <- length(start)  # number of parameters
    if (is.null(lower)) lower<- -Inf
    if (is.null(upper)) upper<- Inf
    if (length(lower) == 1) 
        lower <- rep(lower, npar) # expand to full dimension
    if (length(upper) == 1) 
        upper <- rep(upper, npar)
# more tests on bounds
    if (length(lower) != npar) 
        stop("Wrong length: lower")
    if (length(upper) != npar) 
        stop("Wrong length: upper")
    if (any(start < lower) || any(start > upper)) {
        if(trace) {  cat("start:"); print(start)
                     cat("\n");cat("lower:"); print(lower)
                     cat("\n"); cat("upper:"); print(upper) 
                  }
        stop("Infeasible start")
    }
    if (trace && ctrl$prtlvl>1) {
        cat("formula: ")
        print(formula)
        cat("lower:"); print(lower); cat("upper:"); print(upper) 
    }
    if (trace && ctrl$prtlvl>2) print(str(ctrl))
# Note spelling of lamda -- a throwback to Ag Can 1974 and way to see if folk are copying code.
# First get all the variable names:
    vn <- all.vars(formula)
    # Then see which ones are parameters (get their positions
    # in the set xx
    pnum <- start
    pnames <- names(pnum)
    bdmsk <- rep(1, npar)  # set all params free for now
    maskidx <- which(lower==upper)
    if (length(maskidx) > 0 && trace) {
        cat("The following parameters are masked:")
        print(pnames[maskidx])
    }
    bdmsk[maskidx] <- 0  # fixed parameters
    if (trace && ctrl$prtlvl>1) { # diagnostic printout
        cat("Finished masks check\n")
        parpos <- match(pnames, vn) # Can we explain the 'why'?
        datvar <- vn[-parpos]  # NOT the parameters
        cat("datvar:")
        print(datvar)
        for (i in 1:length(datvar)) {
            dvn <- datvar[[i]]
            cat("Data variable ", dvn, ":")
            if (is.null(data)) { 
                print(eval(parse(text = dvn)))
            } else {
                print(with(data, eval(parse(text = dvn)))) 
            }
        }
    }
    if (is.null(ctrl$japprox)) {
       trjfn<-model2rjfun(formula, pnum, data=data) 
    } else {
      if(ctrl$japprox == "SSJac") {
         trjfn<-SSmod2rjfun(formula, pnum, data=data)
      } else {
         trjfn<-model2rjfun(formula, pnum, data=data, jacobian=FALSE)
      }
    } 
    if (trace && ctrl$prtlvl>2) { cat("str(trjfn):\n"); print(str(trjfn))}
    if (trace && ctrl$prtlvl>2){ cat("ctrl:\n"); print(ctrl) }
    if (inherits(weights, "formula")) {
      # This is a local function so that it has access to 
      # local variables like data, formula
      weights2fun <- function(weights) {
        if (length(weights) != 2)
          warning("Left hand side of weights formula was ignored.")
        weightsexpr <- rhs(weights)
        # Getting fitted values may be slow, so only get
        # them if necessary
        needfitted <- "fitted" %in% all.vars(weightsexpr)
        if (is.null(data))
          data <- environment(weights)
        else if (is.list(data))
          data <- list2env(data, parent = environment(weights))
        else if (!is.environment(data))
          stop("'data' must be a dataframe, list, or environment")
        function(prm, resids) {
          if (is.null(names(prm))) 
            names(prm) <- names(pnum)
          localdata <- list2env(as.list(prm), parent = data)
          localdata$resid <- resids
          if (needfitted)
            localdata$fitted <- eval(rhs(formula), localdata)
          eval(weightsexpr, localdata)
        }
      } 
      weights <- weights2fun(weights)
    }
    if (is.null(ctrl$japprox)) { # call regularly when no jacobian approx
       resfb <- nlfb(start=pnum, resfn=trjfn, jacfn=trjfn, trace=trace, 
            data=data, lower=lower, upper=upper, 
            weights=weights, ctrlcopy=TRUE, control=ctrl)
    } else { # need approx
      if(!is.character(ctrl$japprox)) stop("Non-character ctrl$japprox")
      if(ctrl$japprox == "SSJac") { # Jacobian from selfStart model
        if (trace) { cat('Using Jacobian code from selfStart model \n') }
       resfb <- nlfb(start=pnum, resfn=trjfn, jacfn=trjfn, trace=trace, 
            data=data, lower=lower, upper=upper, 
            weights=weights, ctrlcopy=TRUE, control=ctrl)
      } else {
      if (trace) { cat('Using approximation ',as.character(ctrl$japprox),' 
		 is.character(ctrl$japprox)=', is.character(ctrl$japprox),'\n') }
       resfb <- nlfb(start=pnum, resfn=trjfn, jacfn=ctrl$japprox, trace=trace, 
            data=data, lower=lower, upper=upper, 
            weights=weights, ctrlcopy=TRUE, control=ctrl)
      }
    }
    resfb$formula <- formula # 190805 to add formula to solution (for predict?)
    resfb$resfn <- trjfn # May be lacking Jacobian
    resfb$data <- substitute(data) # Should there be any ... arguments?
    resfb$weights0 <- weights0
    pnum <- as.vector(resfb$coefficients)
    names(pnum) <- pnames # Make sure names re-attached. Is this needed?
    names(resfb$coefficients) <- pnames
    result <- resfb
    class(result) <- "nlsr" ## needed for print method
    result
}

