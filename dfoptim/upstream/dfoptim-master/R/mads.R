
mads <- (function(par, fn, lower = -Inf, upper = Inf, scale = 1.0, control = list(), ...){
    # par : initial solution, where to start from, x0
    # fn : name of the function to optimize, f(x)
    # lower, upper: optional bounds; can be single value (same for each component in x) or vectors same length of x
    # scaling: multiples of delta for unbalanced variables (usually for unbounded problems); one value for all or one value per component
    # control: see descriptions in mads.control()
   
  mads.unbounded <- function(par, fn, scale = 1.0, control = list(), ...){
    return(mads(par, fn, -Inf, Inf, scale, control))
  }
  
  mads.bounded <- function(par, fn, lower = -Inf, upper = Inf, scale = 1.0, control = list(), ...){
    return(mads(par, fn, lower, upper, scale, control))
  }


     #-- Check initial solution and bounds
    if (!is.numeric(par))
        stop("Argument 'par' must be a numeric vector.", call. = FALSE)
    nvar <- length(par)
    if (nvar == 1)
        stop("For univariate functions use some different method.", call. = FALSE)
    
    #-- Handle bounds
    if(!is.numeric(lower) || !is.numeric(upper))
        stop("Lower and upper limits must be numeric.", call. = FALSE)
    if (length(lower) == 1) lower <- rep(lower, nvar)
    if (length(upper) == 1) upper <- rep(upper, nvar)
    if (!all(lower <= upper))
        stop("All lower limits must be smaller than upper limits.", call. = FALSE)
    if (!all(lower <= par) || !all(par <= upper))
        stop("Infeasible starting values -- check limits.", call. = FALSE)
    if(!is.numeric(scale))
        stop("Scaling factor must be numeric, put 1 if you don't know what to do", call. = FALSE)
    if (length(scale) == 1) scale <- rep(scale, nvar)
    isBounded = TRUE
    offset = rep(0, nvar)
    if(sum(is.finite(lower))==0 & sum(is.finite(upper))==0)
        isBounded=FALSE
    else if(sum(is.finite(lower))<nvar | sum(is.finite(upper))<nvar)
        stop("Bounds must be all finite or all infinite, but not partially finite.", call. = FALSE)
    else
        if(!all(scale<(upper-lower)) ) stop("Scale factor is larger than bound width.", call. = FALSE)
    
    #-- If still alive!
    offset = rep(0, nvar)
    scaling = rep(1,nvar)
    if(isBounded){
        offset = lower
        scaling = (upper-lower)
    }
    
    #-- Control list handling ----------
    param <- list(
      "trace" = TRUE,        # for printing interim information
      "tol" = 0.000001,      # algorithm will stop when mesh size goes smaller than that
      "maxfeval" = 10000,    # set to Inf if no limit wanted
      "maximize" = FALSE,    # set to TRUE for maximization
      "pollStyle" = "lite",  # lite=n+1, full=2n  with n the number of decision variables
      "deltaInit" = 0.01,     # initial mesh size (in the [-1, 1]^n space...)
      "expand" = 4,          # after successful iteration, augment mesh size by this or shrink if unsuccessful
      "lineSearch" = 20,     # max attempts along best direction if poll was successful; set to <1 to disable linesearch
      "seed" = 1138          # randmoness is everywhere, make it reproducible!
    )
    
    namc <- match.arg(names(control), choices=names(param), several.ok=TRUE)
    if (!all(namc %in% names(param))) 
      stop("unknown names in control: ", namc[!(namc %in% names(param))])
    if (!is.null(names(control))) param[namc] <- control
    if(!all(scale>param$tol)) stop("Scaling factor smaller than tolerance, please revise.", call. = FALSE)
    
    #-- Prepare black-box ----------
    set.seed(param$seed)
    fun <- match.fun(fn)
    if(isBounded){
        blackbox <- (function(x, lb=offset, span=scaling) {ifelse(param$maximize,-1,1) * fun(lb+(x+1)*span/2, ...)})
        best_x = 2*(par-offset)/scaling-1  # solver will work in [-1, 1]^n space
    } else {
        blackbox <- (function(x, lb=offset, span=scaling) {ifelse(param$maximize,-1,1) * fun(lb+x*span, ...)})
        best_x = (par-offset)/scaling  # solver will work in [-1, 1]^n space
    }

    #-- Prepare the solver and the inputs
    delta = param$deltaInit  # current mesh size
    pollSize = ifelse(param$pollStyle=="lite", nvar+1, 2*nvar)
    n = 1 # number of objective calls
    fibSeq <- NULL # fibonacci sequence for exploding exploration step sizes (feel lucky? go faster!)
    zoom = param$expand
    if(param$lineSearch>2){
        fibSeq <- c(1,2)
        for(k in 3:param$lineSearch) fibSeq <- c(fibSeq, fibSeq[k-1]+fibSeq[k-2])
    }
    #best_x = (par-offset)/scaling  # solver will work in [-1, 1]^n space
    best_f = blackbox(best_x) # let's start by this!
   
    # Trigger iterations
    goAhead = TRUE
    output <- list(par=best_x, value=best_f, feval=1, convergence=param$deltaInit, iterlog=NULL)
    iterlog <- data.frame(n=1, delta=param$deltaInit, searchSuccess=0, f=best_f)
    while(goAhead){

        foundBetter = FALSE
        searchSuccess = 0
        current_x <- best_x

        # get a pollset and evaluate it
        pollX = .pollSet(best_x, pollSize, delta, scale)
        pollY = apply(pollX, 1, (function(x){ ifelse(all(lower<=x & x<=upper), blackbox(x), NA) }))
        n = n + sum(!is.na(pollY))
        pollY.min = which.min(pollY)
        if(length(pollY.min)>0) {if(pollY[pollY.min] < best_f){
            foundBetter = TRUE
            best_f = pollY[pollY.min]
            best_x = pollX[pollY.min,]
        }}

        # possibly run the line search (stop when decrease stops)
        if(foundBetter & !is.null(fibSeq)){
        	lsearchX = .linesearchSet(best_x, best_x - current_x, delta, fibSeq, scale)
        	for(k in 1:length(fibSeq)) {
        	    lsearchY = ifelse(all(lower<=lsearchX[k,] & lsearchX[k,]<=upper), blackbox(lsearchX[k,]), NA)
        		if(is.numeric(lsearchY)) {
        		    n = n+1
            		if(lsearchY<best_f) {
            			best_f = lsearchY
            			best_x = lsearchX[k,]
            			searchSuccess = fibSeq[k]
            		} else break
        		} else break
        	}
        }

        # update algorithm
        if(foundBetter){
            delta = delta * zoom 
            if(delta>1) { delta = 1 } #prevent useless iteration because delta is too large
        } else{
            delta = delta / sqrt(zoom)
        }
        if(delta<=param$tol | n>=param$maxfeval) { goAhead = FALSE }
        
        iterlog <- rbind(iterlog, c(n, delta, searchSuccess, best_f))
        if(param$trace) { 
        	message(paste("n=", n, "   lsearch=", searchSuccess, "   f=", best_f, "   delta=", delta))
        }

    }
    
    # Terminate solver
    if(isBounded)
        output$par = offset+(best_x+1)*scaling/2
    else
        output$par = offset+best_x*scaling
    output$value = best_f
    output$feval = n
    output$convergence = delta
    output$iterlog= if (param$trace) iterlog else NA
    return(output)
})


#
# --------Internal functions, don't call, don't touch-------------------------------------------------------------------
# 

.pollSet <- (function(center, npoints, meshSize, customScale){
    # Build positive basis and finalize the pollset
    #   - center is where to build around in [-1,1]^n
    #   - npoints is the size of the basis (nvar+1 or 2*nvar)
    #   - delta is the zoom factor (all points within +-delta)
    n = npoints
    m = length(center)
    
    # lower triangular, then shuffle all, then squeeze in [-1, 1]^m (original LT-MADS implementation)
    tempx = matrix(0, nrow=m, ncol=m)
    stepsize = ifelse(meshSize>1, 1, floor(m/(meshSize)))
    for(i in 1:m) {
        # dominant terms on the diagonal
        tempx[i,i] = stepsize * sign(1-2*runif(1))
        # lesser influent terms on lower triangle
        for(j in 1:i-1) { tempx[i,j] = floor(stepsize*(1-2*runif(1))) }
    }
    tempx = apply(tempx,1,function(x) {sample(x)}) #shuffle rows
    tempx = apply(tempx,1,function(x) {sample(x)}) #shuffle columns, remember sample() returns transpose!
    
    # complete to positive basis
    if(n == m+1) {
        # lite poll set, usually works fine, to use in high dimensions
        tempx = meshSize * tempx
        avrg = colMeans(tempx)
        tempx = rbind(tempx, -avrg/sqrt(avrg*avrg))
    } else {
        # exhaustive poll set
        tempx = meshSize * rbind(tempx, -tempx)
    }
    
    # return centered poll set
    return(t(apply(tempx,1,function(x, offset=center, zoom=meshSize*customScale) {zoom*x + offset})))
})

.linesearchSet <- (function(center, slope, delta, increments, customScale){
    # Build series of points along "slope"
    #   - center is where to start from
    #   - slope is the search direction
    #   - delta is the zoom factor 
    #   - increments is an integer vector of increasing step counts, ex. the Fibonacci sequence
    tempx = matrix(0, nrow=length(increments), ncol=length(center))
    for(k in 1:nrow(tempx)) {
        for(j in 1:ncol(tempx)){
            tempx[k,j] = center[j] + delta*customScale[j]*increments[k]*slope[j]
        }
    }
    return((tempx))
})
