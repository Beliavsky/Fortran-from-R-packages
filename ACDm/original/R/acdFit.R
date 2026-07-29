acdFit <- function(durations = NULL, model = "ACD", dist = "exponential", 
                   order = NULL, startPara = NULL,  dailyRestart = FALSE, optimFnc = "optim", 
                   method = "Nelder-Mead", output = TRUE, bootstrapErrors = FALSE, 
                   forceErrExpec = TRUE, fixedParamPos = NULL, bp = NULL,
                   exogenousVariables = NULL, optimFncArgs = list(), control = list()){  
  

# check arguments etc. ----------------------------------------------------
  # checks the 'duration' argument
  vecProvided <- FALSE #indicator of whether 'duration' is a vector
  if(is.data.frame(durations)){   
    if("adjDur" %in% colnames(durations)){
      dur <- as.vector(durations$adjDur)
    } else if("durations" %in% colnames(durations)){
      warning("no 'adjDur' column for diurnally adjusted durations found - used unadjusted durations instead")
      dur <- as.vector(durations$durations)
    } else stop("neither a 'durations' or a 'adjDur' column was found in the data.frame 'durations'")
    if("time" %in% colnames(durations)){
      if(!("POSIXlt" %in% class(durations$time))) durations$time <- as.POSIXlt(durations$time)
      time = durations$time
    } else time <- NULL      
  } else if(is.vector(durations)){
    dur = durations
    vecProvided <- TRUE
    time <- NULL
  } else stop("'durations' must be a data.frame or a vector")
  
  N <- length(dur)
  mean <- mean(dur)
  currentTime <- Sys.time()
  
  #provides the possibility of entering truncated and/or case mismatched arguments:
  model <- match.arg(toupper(model), c("ACD", "LACD1", "LACD2", "EXACD", "AMACD", "ABACD", "AACD", "TACD", "BACD", "BCACD", "SNIACD", "LSNIACD", "TAMACD"))
  dist <- match.arg(tolower(dist), c("exponential", "weibull", "burr", "gengamma", "genf", "qweibull", "mixqwe", "mixqww", "mixinvgauss"))
  
  # --- prepare exogenous variables --------------------------------
  z <- ExoVarNames <- NULL 
  N_exo <- length(exogenousVariables)
  
  if (N_exo != 0) {
    # 1) type check
    if (!is.numeric(exogenousVariables) && !is.character(exogenousVariables)) {
      stop("`exogenousVariables` must be either numeric positions or character names.")
    }
    
    # 2) numeric: integer + in-range
    if (is.numeric(exogenousVariables)) {
      if (any(exogenousVariables %% 1 != 0)) {
        stop("When numeric, `exogenousVariables` must be integer values.")
      }
      if (any(exogenousVariables < 1 | exogenousVariables > ncol(durations))) {
        stop("Numeric `exogenousVariables` must be between 1 and ", ncol(durations), ".")
      }
    } 
    
    # 3) character: must exist in names(durations)
    if (is.character(exogenousVariables)) {
      bad <- setdiff(exogenousVariables, names(durations))
      if (length(bad)) {
        stop("Unknown exogenous variable name(s): ", paste(bad, collapse = ", "))
      }
    }
    
    # 4) extract
    z <- as.matrix(durations[, exogenousVariables, drop = FALSE])
    
    # 5) store the names
    if (is.numeric(exogenousVariables)) {
      ExoVarNames <- names(durations)[exogenousVariables]
    } else {
      ExoVarNames <- exogenousVariables
    }
    
    if (exists("control") && isTRUE(control$use_gradient)) {
      stop("Exogenous variables are not yet supported with `control$use_gradient = TRUE`.")
    }
  }
  # -----------------------------------------------------------------

  # # handles the threshold variable if the model is TACD or TAMACD (for future implementations)
  # if(model %in% c("TACD", "TAMACD")){
  #   if(is.null(thresholdVar) || thresholdVar == "dur") TACDthresholdVar <- c(0, dur[-length(dur)])
  #   else if(thresholdVar == "conditional_dur") TACDthresholdVar <- rep(0, N)
  #   else if(is.vector(thresholdVar)) if(length(thresholdVar) != N) stop("if 'thresholdVar' is  a vector, it must be of the same length as the duration series") 
  #   else stop("incorrect 'thresholdVar' argument")
  # } 
  
  # checks bp (breakpoints):
  if (model %in% c("SNIACD", "LSNIACD", "TACD", "TAMACD")) {
    if (is.null(bp)) {
      bp <- .setBP(1)
    }
    
    if (!is.numeric(bp) || !is.vector(bp)) {
      stop("'bp' must be a numeric vector")
    }
    
    if (any(bp < 0)) {
      stop("'bp' cannot contain negative values")
    }
    
    if (is.unsorted(bp)) {
      bp <- sort(bp)
      warning("'bp' was not sorted; it has been reordered to ascending")
    }
    
    J <- length(bp) + 1
  } else {
    J <- NULL
  }
  
  distCode <- .getDistCode(dist)
  # checks startPara and order input:
  if(length(startPara) != 0){
    if(length(order) == 0) order <- .setOrder(model)
    else .checkOrder(order, model)
    .checkPara(order, startPara, distCode, model, Nexovar = N_exo)
    paraTemp <- .seperateStartPara(startPara, model, distCode, order, N_exo, J)
    distStartPara <- paraTemp$distStartPara
    startPara <- paraTemp$startPara
  }else{
    if(length(order) != 0){
      .checkOrder(order, model)
    } else{
      order <- .setOrder(model)
    }
    startPara <- .setStartPara(model, distCode, mean, order, Nexovar = N_exo, J)
    distStartPara <- startPara$distStartPara
    startPara <- startPara$startPara
  }
  
  # checks the control list arguments:
  con <- list(newDay = 0,
              maxit = 4000,
              xtol_rel = 1.0e-4, 
              n.restarts = 10, 
              n.sim = 100,
              fixed = 1,
              trace = 0,
              B = 999,
              BootRoundTosec = FALSE,
              save_data = TRUE,
              use_gradient = FALSE)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if(length(noNms <- namc[!namc %in% nmsC])) 
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  
  # controls that time was provided if dailyRestart == TRUE, and computes the vector of indices for first observation of a new day
  if(dailyRestart == TRUE && length(time) != 0){
    if(length(con$newDay) == 0 || con$newDay == 0) con$newDay <- .getNewDay(time) #computes the vector of indices for first observation of a news day
  } else if(dailyRestart == TRUE && length(time) == 0 && con$newDay == 0){
    warning("can only use daily restart of the conditional mean if the clocktimes of transactions or the vector newDay (as a control parameter) are provided! Estimation done assuming no daily restart.")
  }
  
  if(con$trace != 0) {
    assign("ACDmOptimTrace", NULL, envir = ACDmGlobalEnv)
    traceMatrix <- NULL
  }
  
  # makes it possible to have fixed parameters when optimizing the log likelihood
  fixedParam <- NULL
  if(length(fixedParamPos) != 0){
    if(length(fixedParamPos) != length(startPara)) stop("the length of 'fixedParamPos' must equal the length of 'startPara'")
    fixedParam <- startPara[fixedParamPos]
    startPara <- startPara[!fixedParamPos]
  }  
  
  # list of all arguments passed to the '.getLLcall' or '.computeLLcpp' function that computes the log likelihood
  getLL.args <- list(
    dur = dur,
    exogenousVar = z,
    model = model,
    order = order,
    mean = mean,
    distCode = distCode,
    newDay = con$newDay,
    returnMu = FALSE,
    breakPoints = bp,
    forceErrExpec = forceErrExpec,
    fixedParam = fixedParam,
    fixedParamPos = fixedParamPos,
    trace = con$trace,
    startType = 2
  )
  if(model %in% c("TACD", "TAMACD", "SNIACD", "LSNIACD")){
    getLL.args <- c(getLL.args, list(J = J))
  } 
    

# optimizes the log likelihood function -----------------------------------
  if(optimFnc == "optim"){   
    # uses the 'control' arguments unless they were given in 'optimFncArgs':
    if(length(optimFncArgs$control$trace) == 0) optimFncArgs$control$trace <- con$trace
    if(length(optimFncArgs$control$maxit) == 0) optimFncArgs$control$maxit <- con$maxit
    
    # boolean vector of what of the possible arguments were entered in 'optimFncArgs':
    possibleArgs <- c("lower", "upper", "control")
    usedArgs <- possibleArgs %in% names(optimFncArgs)
    
    # warns if arguments in 'optimFncArgs' are not used:
    if(length(noNms <- names(optimFncArgs)[!names(optimFncArgs) %in% possibleArgs])) 
      warning("the 'optimFnc' \"", 
              optimFnc,
              "\" does not use the following arguments given in 'optimFncArgs': ",
              paste(noNms, collapse = ", "))
    
    failed <- FALSE
    tryCatch({ #uses the tryCatch to catch the trace path in case of a failure of the optimization function
      if(con$use_gradient){
        
        if(method == "Nelder-Mead") warning("The 'Nelder-Mead' method in optimFnc = 'optim' does not use gradients.")
        
        
        getLL.args$returnIndex <- switch(method,
                                         "Nelder-Mead" = 1, 
                                         "BFGS" = 5, 
                                         "CG" = 5, 
                                         "L-BFGS-B" = 5, 
                                         "SANN" = 1)
        
        fit <- do.call(stats::optim, c(list(par = startPara, fn = .computeLLcpp, gr = .getScoreCpp,
                                            getLL.args = getLL.args, method = method, hessian = TRUE),
                                       list(lower = optimFncArgs$lower,
                                            upper = optimFncArgs$upper,
                                            control = optimFncArgs$control)[usedArgs]))
        
      } else{
        fit <- do.call(stats::optim, c(list(par = startPara, fn = .getLLcall,
                                            getLL.args = getLL.args, method = method, hessian = TRUE),
                                       list(lower = optimFncArgs$lower,
                                            upper = optimFncArgs$upper,
                                            control = optimFncArgs$control)[usedArgs]))
        
      }
    }, error = function(c) {
      
      failed <<- TRUE
      if(con$trace != 0) { #in case the trace option were used and the optimization failed, the trace will be plotted:
        
        numcol <- length(startPara) + length(fixedParam) + 1
        if(length(get("ACDmOptimTrace", envir = ACDmGlobalEnv)) != 0){
          numrow <- ceiling(length(get("ACDmOptimTrace", envir = ACDmGlobalEnv)) / numcol)
          traceMatrix <<- matrix(get("ACDmOptimTrace", envir = ACDmGlobalEnv), 
                                 ncol = numcol, 
                                 byrow = T, 
                                 dimnames = list(1:numrow, c(paste0("para", 1:(numcol - 1)), "LL")))
          .plotTracePath(traceMatrix)
          rm("ACDmOptimTrace", envir = ACDmGlobalEnv)
        }
        
      }
      cat(c$message)
    }
    )
    
    if(failed){ #output in case of optimization failure:
      if(con$trace != 0){
        cat("\n\nacdFit: Oops, seems like the the optimization function failed. The trace path up to the crash was returned\n\n")
        return(traceMatrix)
      } else{
        cat("\n\nError: Oops, seems like the the optimization function failed. Changing the 'optimFnc' or/and its settings, or starting from a diffrent 'startPara' might work. You can also trace the MLE search path by adding the argument 'control = list(trace = 1)'. \n\n")
        return()
      }
    }
    
    
  } else if(optimFnc == "constrOptim"){   
    
    # uses the 'control' arguments unless they were given in 'optimFncArgs':
    if(length(optimFncArgs$control$trace) == 0) optimFncArgs$control$trace <- con$trace
    if(length(optimFncArgs$control$maxit) == 0) optimFncArgs$control$maxit <- con$maxit
    
    if(!all(c("ui", "ci") %in% names(optimFncArgs))) stop("for optimFnc = 'constrOptim', the arguments 'ui' and 'ci' 
                                                   must be present in the list passed to the 'optimFncArgs' argument")
    
    # boolean vector of what of the possible arguments were entered in 'optimFncArgs':
    possibleArgs <- c("mu", "outer.iterations", "outer.eps", "control")
    usedArgs <- possibleArgs %in% names(optimFncArgs)
    
    # warns if arguments in 'optimFncArgs' is not used:
    possibleArgs <- c("ui", "ci", possibleArgs)
    if(length(noNms <- names(optimFncArgs)[!names(optimFncArgs) %in% possibleArgs])) 
      warning("the 'optimFnc' \"", 
              optimFnc,
              "\" does not use the following arguments given in 'optimFncArgs': ",
              paste(noNms, collapse = ", "))
    
    failed <- FALSE
    tryCatch({ #uses the tryCatch to catch the trace path in case of a failure of the optimization function
      if(con$use_gradient){
        
        getLL.args$returnIndex <- switch(method,
                                         "Nelder-Mead" = 1, 
                                         "BFGS" = 5, 
                                         "CG" = 5, 
                                         "L-BFGS-B" = 5, 
                                         "SANN" = 1)
        
        fit <- do.call(stats::constrOptim, c(list(theta = startPara, f = .computeLLcpp, grad = .getScoreCpp,
                                                  getLL.args = getLL.args, method = method, hessian = TRUE,
                                                  ui = optimFncArgs$ui, ci = optimFncArgs$ci),
                                             list(mu = optimFncArgs$mu,
                                                  outer.iterations = optimFncArgs$outer.iterations,
                                                  outer.eps = optimFncArgs$outer.eps,
                                                  control = optimFncArgs$control)[usedArgs]))

        
      } else{
        
        fit <- do.call(stats::constrOptim, c(list(theta = startPara, f = .getLLcall, grad = NULL,
                                                  getLL.args = getLL.args, method = method, hessian = TRUE,
                                                  ui = optimFncArgs$ui, ci = optimFncArgs$ci),
                                             list(mu = optimFncArgs$mu,
                                                  outer.iterations = optimFncArgs$outer.iterations,
                                                  outer.eps = optimFncArgs$outer.eps,
                                                  control = optimFncArgs$control)[usedArgs]))

      }
    }, error = function(c) {
      
      failed <<- TRUE
      if(con$trace != 0) { #in case the trace option were used and the optimization failed, the trace will be plotted:
        
        numcol <- length(startPara) + length(fixedParam) + 1
        if(length(get("ACDmOptimTrace", envir = ACDmGlobalEnv)) != 0){
          numrow <- ceiling(length(get("ACDmOptimTrace", envir = ACDmGlobalEnv)) / numcol)
          traceMatrix <<- matrix(get("ACDmOptimTrace", envir = ACDmGlobalEnv), 
                                 ncol = numcol, 
                                 byrow = T, 
                                 dimnames = list(1:numrow, c(paste0("para", 1:(numcol - 1)), "LL")))
          .plotTracePath(traceMatrix)
          rm("ACDmOptimTrace", envir = ACDmGlobalEnv)
        }
        
      }
      cat(c$message)
    }
    )
    
    if(failed){ #output in case of optimization failure:
      if(con$trace != 0){
        cat("\n\nacdFit: Oops, seems like the the optimization function failed. The trace path up to the crash was returned\n\n")
        return(traceMatrix)
      } else{
        cat("\n\nError: Oops, seems like the the optimization function failed. Changing the 'optimFnc' or/and its settings, or starting from a diffrent 'startPara' might work. You can also trace the MLE search path by adding the argument 'control = list(trace = 1)'. \n\n")
        return()
      }
    }
    
    
  } else if(optimFnc == "nlminb"){
    
    # uses the 'control' arguments unless they were given in 'optimFncArgs':
    if(length(optimFncArgs$control$trace) == 0) optimFncArgs$control$trace <- con$trace
    if(length(optimFncArgs$control$maxit) == 0) optimFncArgs$control$maxit <- con$maxit
    
    # boolean vector of what of the possible arguments were entered in 'optimFncArgs':
    possibleArgs <- c("scale", "lower", "upper", "control")
    usedArgs <- possibleArgs %in% names(optimFncArgs)
    
    # warns if arguments in 'optimFncArgs' is not used:
    if(length(noNms <- names(optimFncArgs)[!names(optimFncArgs) %in% possibleArgs])) 
      warning("the 'optimFnc' \"", 
              optimFnc,
              "\" does not use the following arguments given in 'optimFncArgs': ",
              paste(noNms, collapse = ", "))
    
    failed <- FALSE
    tryCatch({ #uses the tryCatch to catch the trace path in case of a failure of the optimization function
      
      if(con$use_gradient){
        getLL.args$returnIndex <- 5
        
        fit <- do.call(stats::nlminb, c(list(start = startPara, objective = .computeLLcpp, 
                                             gradient = .getScoreCpp, getLL.args = getLL.args),
                                        list(scale = optimFncArgs$scale,
                                             lower = optimFncArgs$lower,
                                             upper = optimFncArgs$upper,
                                             control = optimFncArgs$control)[usedArgs]))
        
      } else{
        fit <- stats::nlminb(start = startPara, objective = .getLLcall, getLL.args = getLL.args,
                             control = list(trace = con$trace, iter.max = con$maxit), lower = -Inf, upper = Inf)
      }
      
      
    }, error = function(c) {
      failed <<- TRUE
      
      if(con$trace != 0) { #in case the trace option were used and the optimization failed, the trace will be plotted:
        
        numcol <- length(startPara) + length(fixedParam) + 1
        if(acdFit != 0){
          numrow <- ceiling(length(get("ACDmOptimTrace", envir = ACDmGlobalEnv)) / numcol)
          traceMatrix <<- matrix(get("ACDmOptimTrace", envir = ACDmGlobalEnv), 
                                 ncol = numcol, 
                                 byrow = T, 
                                 dimnames = list(1:numrow, c(paste0("para", 1:(numcol - 1)), "LL")))
          .plotTracePath(traceMatrix)
          rm("ACDmOptimTrace", envir = ACDmGlobalEnv)
        }
        
        cat(c$message)
        
      }
    }
    )
    
    if(failed){ #output in case of optimization failure:
      if(con$trace != 0){
        cat("\n\nacdFit: Oops, seems like the the optimization function failed. The trace path up to the crash was returned\n\n")
        return(traceMatrix)
      } else{
        cat("\n\nError: Oops, seems like the the optimization function failed. Changing the 'optimFnc' or/and its settings, or starting from a diffrent 'startPara' might work. You can also trace the MLE search path by adding the argument 'control = list(trace = 1)'. \n\n")
        return()
      }
    }
    
    fit <- list(par = fit$par, hessian = NULL, value = fit$objective, convergence = fit$convergence, counts = fit$evaluations[2])
    
  } else if(optimFnc == "solnp"){
    
    # uses the 'control' arguments unless they were given in 'optimFncArgs':
    if(length(optimFncArgs$control$trace) == 0) optimFncArgs$control$trace <- con$trace
    if(length(optimFncArgs$control$maxit) == 0) optimFncArgs$control$maxit <- con$maxit
    
    # boolean vector of what of the possible arguments were entered in 'optimFncArgs':
    possibleArgs <- c("eqfun", "eqB", "ineqfun", "ineqLB", "ineqUB", "LB", "UB", "control")
    usedArgs <- possibleArgs %in% names(optimFncArgs)
    
    # warns if arguments in 'optimFncArgs' is not used:
    if(length(noNms <- names(optimFncArgs)[!names(optimFncArgs) %in% possibleArgs])) 
      warning("the 'optimFnc' \"", 
              optimFnc,
              "\" does not use the following arguments given in 'optimFncArgs': ",
              paste(noNms, collapse = ", "))
    
    failed <- FALSE
    tryCatch({ #uses the tryCatch to catch the trace path in case of a failure of the optimization function
      
      if(con$use_gradient) warning("The optimFnc = 'solnp' optimization does not use analytical gradients.")
      options(warn = -1)
      fit <- do.call(Rsolnp::solnp, c(list(pars = startPara, fun = .getLLcall, getLL.args = getLL.args),
                                      list(eqfun = optimFncArgs$eqfun,
                                           eqB = optimFncArgs$eqB,
                                           ineqfun = optimFncArgs$ineqfun,
                                           ineqLB = optimFncArgs$ineqLB,
                                           ineqUB = optimFncArgs$ineqUB,
                                           LB = optimFncArgs$LB,
                                           UB = optimFncArgs$UB,
                                           control = optimFncArgs$control)[usedArgs]))
      options(warn = 0)
      
      getLL.args$returnIndex <- 5
       
    }, error = function(c) {
      
      failed <<- TRUE
      if(con$trace != 0) { #in case the trace option was used and the optimization failed, the trace will be plotted:
        
        numcol <- length(startPara) + length(fixedParam) + 1
        if(length(get("ACDmOptimTrace", envir = ACDmGlobalEnv)) != 0){
          numrow <- ceiling(length(get("ACDmOptimTrace", envir = ACDmGlobalEnv)) / numcol)
          traceMatrix <<- matrix(get("ACDmOptimTrace", envir = ACDmGlobalEnv), 
                                 ncol = numcol, 
                                 byrow = T, 
                                 dimnames = list(1:numrow, c(paste0("para", 1:(numcol - 1)), "LL")))
          .plotTracePath(traceMatrix)
          rm("ACDmOptimTrace", envir = ACDmGlobalEnv)
        }
        cat(c$message)
        
      }
      
    }
    )  
    
    if(failed){ #output in case of optimization failure:
      if(con$trace != 0){
        cat("\n\nacdFit: Oops, seems like the the optimization function failed. The trace path up to the crash was returned\n\n")
        return(traceMatrix)
      } else{
        cat("\n\nError: Oops, seems like the the optimization function failed. Changing the 'optimFnc' or/and its settings, or starting from a diffrent 'startPara' might work. You can also trace the MLE search path by adding the argument 'control = list(trace = 1)'. \n\n")
        return()
      }
    }
    
    fit <- list(par = fit$pars, hessian = NULL, 
                value = fit$values[length(fit$values)], convergence = fit$convergence, counts = fit$nfuneval)
    
    
  } else if(optimFnc == "optimx"){
    
    #checks if the required package 'optimx' is installed:
    if (!requireNamespace("optimx", quietly = TRUE)) {
      stop("The package 'optimx' needs to be installed. 
           Install it by entering 'install.packages(\"optimx\")' in the console, and try again.",
           call. = FALSE)
    }
    
    failed <- FALSE
    tryCatch({ #uses the tryCatch to catch to still get the trace path in case of a failure of the optimization function
      
      if(con$use_gradient){
        getLL.args$returnIndex <- 5
        fit <- optimx::optimx(startPara, .computeLLcpp, gr = .getScoreCpp,
                              method = method, itnmax = con$maxit, control = list(trace = con$trace, kkt = FALSE), 
                              getLL.args = getLL.args)
        
      } else{fit <- optimx::optimx(startPara, .getLLcall, 
                                   method = method, itnmax = con$maxit, control = list(trace = con$trace, kkt = FALSE), 
                                   getLL.args = getLL.args)
      }
      
      fit <- optimx::optimx(startPara, .getLLcall, 
                            method = method, itnmax = con$maxit, control = list(trace = con$trace, kkt = FALSE), 
                            getLL.args = getLL.args)
      
    }, error = function(c) {
      failed <<- TRUE
      
      if(con$trace != 0) { #in case the trace option were used and the optimization failed, the trace will be plotted:
        
        numcol <- length(startPara) + length(fixedParam) + 1
        if(length(get("ACDmOptimTrace", envir = ACDmGlobalEnv)) != 0){
          numrow <- ceiling(length(get("ACDmOptimTrace", envir = ACDmGlobalEnv)) / numcol)
          traceMatrix <<- matrix(get("ACDmOptimTrace", envir = ACDmGlobalEnv), 
                                 ncol = numcol, 
                                 byrow = T, 
                                 dimnames = list(1:numrow, c(paste0("para", 1:(numcol - 1)), "LL")))
          .plotTracePath(traceMatrix)
          rm("ACDmOptimTrace", envir = ACDmGlobalEnv)
        }
        
        cat(c$message)
        
      }
    }
    )
    
    if(failed){ #output in case of optimization failure:
      if(con$trace != 0){
        cat("\n\nacdFit: Oops, seems like the the optimization function failed. The trace path up to the crash was returned\n\n")
        return(traceMatrix)
      } else{
        cat("\n\nError: Oops, seems like the the optimization function failed. Changing the 'optimFnc' or/and its settings, or starting from a diffrent 'startPara' might work. You can also trace the MLE search path by adding the argument 'control = list(trace = 1)'. \n\n")
        return()
      }
    }
    
    
    fit <- list(par = as.numeric(fit[1, 1:length(startPara)]), hessian = NULL, 
                value = fit$value, convergence = fit$convcode, counts = fit$fevals)
    
  }
  
  #### calculates the hessian numerically:
  tryCatch({
    if(con$use_gradient){
      getLL.args$returnIndex <- 1
      fit$hessian <- numDeriv::hessian(func = .computeLLcpp, x = fit$par, method="Richardson", 
                                       method.args=list(eps=1e-4, d=0.001, zero.tol=sqrt(.Machine$double.eps/7e-7), 
                                                        r=4, v=2, show.details=FALSE), 
                                       getLL.args = getLL.args)
      
    }else{
      fit$hessian <- numDeriv::hessian(func = .getLLcall, x = fit$par, method="Richardson", 
                                       method.args=list(eps=1e-4, d=0.001, zero.tol=sqrt(.Machine$double.eps/7e-7), 
                                                        r=4, v=2, show.details=FALSE), 
                                       getLL.args = getLL.args)
    }
  }, error = function(c) {
    warning("computing the hessian failed: ", c$message)
  })
  
  #returns the muHats (the estimated conditional durations):
  getLL.args$returnMu <- TRUE
  if(con$use_gradient){
    getLL.args$returnIndex <- 2
    getLL.args$trace <- 0
    muTmp <- .computeLLcpp(param = fit$par, getLL.args = getLL.args)     
    mu <- list(mu = muTmp$psi, resi = muTmp$e)
  }else{
    mu <- .getLLcall(param = fit$par, getLL.args = getLL.args)     
  }
  
  #if some parameters were set to be fixed, the fixed and estimated parameters are recombined:
  if(length(fixedParamPos) != 0) parTemp <- .returnfixedPara(fit$par, fixedParam, fixedParamPos)
  else parTemp <- fit$par
  

# computes bootstrap errors -----------------------------------------------
  if(bootstrapErrors){
    if(!con$BootRoundTosec){ #the simulation in the bootstrap won't be rounded to seconds
      
      bootPar <- matrix(nrow = con$B, ncol = length(fit$par))
      i <- 1
      percDone = 5
      failed = 0
      bootConverged <- rep(-99, con$B)
      
      bootStartTime <- Sys.time()
      while(i <= con$B){
        bootDur <- sim_ACD(N, model = model, param = fit$par, order = order, startX = mean, startMu = mean, errors = mu$resi, dist = dist, roundToSec = FALSE)
        getLL.args$dur <- bootDur
        getLL.args$mean <- mean(bootDur)
        getLL.args$returnMu <- FALSE
        bootTemp <- tryCatch(stats::optim(par = fit$par, fn = .getLLcall, getLL.args = getLL.args,  hessian = F, control = list(maxit = con$maxit, trace = con$trace, trace = con$trace)), error = function(e) {NULL})
        bootParTemp <- bootTemp$par 
        
        if(length(bootParTemp) != 0 && all(abs(bootParTemp[-1]) < 1.5) && bootTemp$convergence == 0){
          bootPar[i, ] <- bootParTemp
          bootConverged[i] <- bootTemp$convergence
          i <- i + 1          
        } else failed <- failed + 1
        
        if ((i / con$B) >= .05 && percDone == 5) cat("Estimated time for bootstrap simulation: ", round(difftime(Sys.time(), bootStartTime, units = "secs")*20), " sec \n\nbootstrap % done: \n")
        if ((i / con$B) >= (percDone / 100)) {cat(percDone,"% "); percDone = percDone + 5}
      }      
      cat("\ntime for bootstrap simulation: ", round(difftime(Sys.time(), bootStartTime, units = "secs")), " sec")
      cat("\n\n")
      cat(failed, "of the ", con$B, " bootstrap estimations failed and were resimulated\n")
      bootErr <- sqrt(diag(stats::cov(bootPar)))
      bootCorr <- stats::cor(bootPar)
      bootMean <- apply(bootPar, 2, mean)
    } else{ #the simulation in the bootstraps will be rounded to seconds
      bootPar <- matrix(nrow = con$B, ncol = length(fit$par))
      i <- 1
      percDone = 5
      cat("bootstrap % done: ")
      bootDurTemp <- sim_ACD((con$B*(N+50))*1.3, param = fit$par, model = model, order = order, startX = mean, startMu = mean, errors = mu$resi, roundToSec = FALSE)
      bootDur <- bootDurTemp[bootDurTemp!=0]
      getLL.args$dur <- bootDur
      getLL.args$returnMu <- FALSE
      if(length(bootDur) < N) 
        while(i <= con$B){
          getLL.args$dur <- bootDur[((i-(1+failed))*(N+50)+1):((i+failed)*(N+50))]
          getLL.args$mean <- mean(bootDur)
          bootParTemp <- tryCatch(
            stats::optim(par = fit$par, fn = .getLLcall, getLL.args = getLL.args, 
                         hessian = TRUE, control = list(maxit = con$maxit, trace = con$trace))$par, 
            error = function(e) {NULL}
          )
          if(length(bootParTemp) != 0){
            bootPar[i, ] <- bootParTemp
            i <- i + 1          
          } else failed <- failed + 1
          if ((i / con$B) >= (percDone / 100)) {cat(percDone,"% "); percDone = percDone + 5}
        }      
      cat("\n")
      bootErr <- sqrt(diag(stats::cov(bootPar)))
      bootCorr <- stats::cor(bootPar)
      bootMean <- apply(bootPar, 2, mean)
    }
  }
  

# Robust errors -----------------------------------------------------------
  #if QML (dist = "exponential") and not fixed parameters and not exogenous variables, then robust errors will be computed:
  robustSE <- NULL
  robustCorr <- NULL
  if(dist == "exponential" && length(fixedParamPos) == 0 && length(exogenousVariables) == 0){
    getLL.args$returnIndex <- 4
    getLL.args$trace <- 0
    QLMscore <- .computeLLcpp(fit$par, getLL.args = getLL.args)
    
    QLMscore <- .sumOuter(QLMscore$d_LL_d_theta, QLMscore$d_psi_d_theta, QLMscore$psi)
    
    try({
      sandwich <- solve(QLMscore$A) %*% QLMscore$B %*% solve(QLMscore$A)
      robustSE <- sqrt(diag(sandwich)) 
      robustCorr <- solve(diag(robustSE)) %*% sandwich %*% solve(diag(robustSE))
    })
  }
  # if(model == "ACD" && dist == "exponential" && length(fixedParamPos) == 0 && length(exogenousVariables) == 0){
  #   QLMscore <- .Call("getScoreACDExp",
  #                     as.double(dur),
  #                     as.double(mu$mu),
  #                     as.double(fit$par),                     
  #                     as.integer(order),
  #                     as.integer(0), PACKAGE = "ACDm")
  #   
  #   sandwich <- solve(as.matrix(as.data.frame(QLMscore[3]))) %*% as.matrix(as.data.frame(QLMscore[4])) %*% solve(as.matrix(as.data.frame(QLMscore[3])))
  #   robustSE <- sqrt(diag(sandwich)) 
  #   robustCorr <- solve(diag(robustSE)) %*% sandwich %*% solve(diag(robustSE))
  # }
  
  if(bootstrapErrors) namedParameters <- .getCoef(para = fit$par , model = model, dist = dist, hessian = fit$hessian, 
                                                  order = order, bootError = bootErr, bootCorr = bootCorr, bootMean = bootMean, 
                                                  robustCorr = robustCorr, robustSE = robustSE, fixedParam = fixedParam, 
                                                  fixedParamPos = fixedParamPos, ExoVarNames = ExoVarNames, J = J)
  else namedParameters <- .getCoef(para = fit$par , model = model, dist = dist, hessian = fit$hessian, 
                                   order = order, robustCorr = robustCorr, robustSE = robustSE, fixedParam = fixedParam, 
                                   fixedParamPos = fixedParamPos, ExoVarNames = ExoVarNames, J = J)
  
  
  
  N <- length(dur)
  Npar <- length(fit$par)
  LogLikelihood <- -fit$value
  AIC <- 2 * (Npar - LogLikelihood)
  BIC <- -2 * LogLikelihood + Npar * log(N)
  MSE <- mean((dur-mu$mu)^2)
  GoodnessOfFit <- data.frame("value" = c(LogLikelihood, AIC, BIC, MSE))
  rownames(GoodnessOfFit) <- c("LogLikelihood", "AIC", "BIC", "MSE")
  #computes the value of the unfree distribution parameter as a function of the others (if mean was forced to be 1):
  if(forceErrExpec == 1) forcedDistPara <- .returnFixedMeanPara(distCode, namedParameters$DPar)
  else{ #if forceErrExpec was false (0) this parameter was fixed at 1 in the estimation
    forcedDistPara <- 1
    names(forcedDistPara) <- names(.returnFixedMeanPara(distCode, namedParameters$DPar)) #only to get the name of the parameter
  } 
  
  if(!vecProvided) returnValue <- list(call = match.call(),
                                       durations = durations)
  else returnValue <- list(call = match.call(),
                           durations = data.frame(durations = dur))
  
  returnValue <- append(returnValue, list(
    muHats = mu$mu, 
    residuals = mu$resi,
    model = model,
    order = order,
    distribution = dist, 
    distCode = distCode,
    startPara = startPara,
    mPara = namedParameters$MPar, 
    dPara = namedParameters$DPar,
    exogenousVariables = exogenousVariables,
    breakPoints = bp,
    Npar = Npar[1],
    goodnessOfFit = GoodnessOfFit,
    parameterInference = namedParameters$Inference,
    forcedDistPara = forcedDistPara,
    forceErrExpec = forceErrExpec,
    comments = namedParameters$comment,
    hessian = fit$hessian,
    N = N,
    evals = fit$counts[1],
    convergence = fit$convergence,
    estimationTime = difftime(Sys.time(), currentTime, units = "secs"),
    #description = paste("Estimated at", currentTime, "by user", Sys.info()[["user"]]),
    description = sprintf(
      "Estimated at %s by user %s",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      Sys.info()[["user"]]
    ),
    newDayVector = con$newDay))
  
  
  #if bootstrapp errors: adds the bootstrapp inference
  if(bootstrapErrors) returnValue <- append(returnValue, list(bootstrapEstimates = bootPar, bootConverged = bootConverged, bootErr = bootErr, bootMean = bootMean, bootCorr = namedParameters$bootCorr, bootPar = bootPar))
  #if QML (ACD and exponetial): adds the robust correlation
  if(model == "ACD" && dist == "exponential" && length(fixedParamPos) == 0 && length(exogenousVariables) == 0) returnValue <- append(returnValue, list(robustCorr = namedParameters$robustCorr))

  #plots and append the trace path if the argument 'control = list(trace = 1)' was given
  if(con$trace != 0) {
    numcol <- length(startPara) + length(fixedParam) + 1
    numrow <- ceiling(length(get("ACDmOptimTrace", envir = ACDmGlobalEnv)) / numcol)
    traceMatrix <- matrix(get("ACDmOptimTrace", envir = ACDmGlobalEnv), ncol = numcol, byrow = T, 
                          dimnames = list(1:numrow, c(names(returnValue$mPara), names(returnValue$dPara), "LL")))
    .plotTracePath(traceMatrix)
    rm("ACDmOptimTrace", envir = ACDmGlobalEnv)
    returnValue <- append(returnValue, list(traceMatrix = traceMatrix))
  }
  
  class(returnValue) <-  c("acdFit", class(returnValue))
  if(output) print(returnValue)
  
  if(con$save_data == FALSE){
    
    returnValue$durations <- NULL
    returnValue$muHats <- NULL
    returnValue$residuals <- NULL
    
  }
  
  acdFit <- returnValue
}
