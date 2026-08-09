dex <- function(x, do_substitute = NA, verbose = FALSE) {
  expr <- substitute(x)
  if (is.na(do_substitute) || !do_substitute)
    value <- tryCatch(x, error = function(e) e)
  if (is.na(do_substitute)) {
    if (verbose)
      message("Determining whether to use expression or value...")
    if (inherits(value, "error"))
      do_substitute <- TRUE
    else if (is.character(value))
      do_substitute <- FALSE
    else if (is.name(expr))
      do_substitute <- FALSE
    else
      do_substitute <- TRUE
  } else if (!do_substitute && inherits(value, "error"))
    stop(conditionMessage(value), call. = FALSE)
  if (verbose) {
    if (do_substitute) 
      message("Using expression.")
    else message("Using value.")
  }  
  if (!do_substitute) {
    if (is.character(value)) {
      if (verbose) 
        message("Parsing value.")
      expr <- parse(text=value)
    } else
      expr <- value
  }
  if (is.expression(expr) && length(expr) == 1)
    expr <- expr[[1]]
  else if (is.call(expr) && as.character(expr[[1]]) == "~") {
    if (length(expr) == 3) {
      warning("Left hand side of formula will be ignored.")
      expr <- expr[[3]]
    } else
      expr <- expr[[2]]
  }
  expr
}

sysDerivs <- new.env(parent = emptyenv())

sysSimplifications <- new.env(parent = emptyenv())

newDeriv <- function(expr, deriv, derivEnv = sysDerivs) {
    if (missing(expr))
    	return(ls(derivEnv))
    expr <- substitute(expr) # Get parse tree of expr 
    if (!is.call(expr))
    	stop("expr must be a call to a function")
    fn <- as.character(expr[[1]])
    if (missing(deriv)) 
    	return(derivEnv[[fn]])
    deriv <- substitute(deriv)
    args <- expr[-1]
    argnames <- names(args)
    if (is.null(argnames))
    	argnames <- rep("", length(args))
    required <- which(argnames == "")
    for (i in required)
    	argnames[i] <- as.character(args[[i]])
    value <- list(expr = expr, argnames = argnames, 
                  required = required, deriv = deriv)
    if (!is.null(oldval <- derivEnv[[fn]]) && !identical(value, oldval))
      warning(gettextf("changed derivative for %s", dQuote(fn)))
    assign(fn, value, envir = derivEnv)
    invisible(value) # returns but does not print
} 

newSimplification <- function(expr, test, simplification, 
                              do_eval = FALSE, simpEnv = sysSimplifications) {
    if (missing(expr))
    	return(ls(simpEnv))
    expr <- substitute(expr)
    if (!is.call(expr))
    	stop("expr must be a call to a function")
    fn <- as.character(expr[[1]])
    nargs <- length(expr) - 1L
    simps <- simpEnv[[fn]]
    if (missing(test)) {
        if (nargs <= length(simps))	
    	    return(simps[[nargs]])
    	else
    	    return(NULL)
    }
    test <- substitute(test)
    simplification <- substitute(simplification)
    
    args <- expr[-1]
    if (!is.null(names(args)))
    	stop("expr should not have named arguments")
    if (!all(sapply(args, is.name)))
    	stop("expr should have simple names as arguments")
    argnames <- sapply(args, as.character)
    if (any(duplicated(argnames)))
    	stop("expr names should be unique")
    	
    if (is.null(simps)) simps <- list()
    if (nargs <= length(simps)) 
    	simpn <- simps[[nargs]]
    else
    	simpn <- list()
    simpn <- c(simpn, list(list(expr = expr, argnames = argnames, test = test, 
                                simplification = simplification, do_eval = do_eval)))
    simps[[nargs]] <- simpn
    assign(fn, simps, envir = simpEnv)
}
    	
nlsDeriv <- function(expr, name, derivEnv = sysDerivs, do_substitute = FALSE, verbose = FALSE, ...) {
    Recurse <- function(expr) {
    	if (is.call(expr)) {
    	    if (as.character(expr[[1]]) == "D")
    	    	expr <- nlsDeriv(expr[[2]], name, derivEnv, do_substitute = FALSE, verbose = verbose, ...)
    	    else
    	    	for (i in seq_along(expr)[-1])
    	    	    expr[[i]] <- Recurse(expr[[i]])
    	}
    	expr
    }
    expr <- dex(expr, do_substitute = do_substitute, verbose = verbose)
    if (is.expression(expr))
    	return(as.expression(lapply(expr, nlsDeriv, name = name, derivEnv = derivEnv, do_substitute = FALSE, verbose = verbose, ...)))
    else if (is.numeric(expr) || is.logical(expr))
    	return(0)
    else if (is.call(expr)) {
    	fn <- as.character(expr[[1]])
	if (fn == "expression")
	    return(as.expression(lapply(as.list(expr)[-1], nlsDeriv, name = name, derivEnv = derivEnv, do_substitute = FALSE, verbose = verbose, ...)))
    	model <- derivEnv[[fn]]
    	if (is.null(model))
    	    stop("no derivative known for '", fn, "'")
 	if (missing(name) || verbose) {
 	    message(paste("Expr:", deparse(expr)))
 	    message(if (missing(name)) "Pattern for" else "Using pattern")
 	    message(paste("  ", deparse(model$expr), collapse = "\n"))
 	    message("is")
 	    message(paste("  ", deparse(model$deriv), collapse = "\n"))
 	    if (missing(name))
 	    	return(invisible(NULL))
 	}
        args <- expr[-1]
        argnames <- names(args)
        if (is.null(argnames)) 
            argnames <- rep("", length(args))
        modelnames <- model$argnames
        argnum <- pmatch(argnames, modelnames)
        if (any(bad <- is.na(argnum) & argnames != ""))
            stop("Argument names not matched: ", paste(argnames[bad], collapse = ", "))
        unused <- setdiff(seq_along(modelnames), argnum)
        nonamecount <- sum(is.na(argnum))
        length(unused) <- nonamecount
        argnum[which(is.na(argnum))] <- unused
        default <- setdiff(seq_along(modelnames), argnum)
        if (length(bad <- setdiff(model$required, argnum)))
            stop("Missing required arguments: ", paste(modelnames[bad], collapse = ", "))
        
        # Now do the substitutions
        subst <- list()
        subst[argnum] <- as.list(args)
        subst[default] <- as.list(model$expr[-1])[default]
        names(subst) <- modelnames
        result <- do.call(substitute, list(model$deriv, subst))
        result <- Recurse(result)
        nlsSimplify(result, verbose = verbose, ...)
    } else if (is.name(expr))
        return( as.numeric(as.character(expr) == name) )
}

codeDeriv <- function(expr, namevec, 
       hessian = FALSE, derivEnv = sysDerivs, 
       do_substitute = FALSE, verbose = FALSE, ...) {
  
  # Make sure expr is a single expression object
  expr <- dex(expr, do_substitute = do_substitute, verbose = verbose)
  expr <- as.expression(expr)
  if (length(expr) > 1)
    stop("Only single expressions allowed")
  
  # Set exprs[[1]] to exprs, fill in the rest with the derivatives
  # that were requested
  exprs <- as.list(expr)
  n <- length(namevec)
  length(exprs) <- n + 1L
  for (i in seq_len(n))
    exprs[[i + 1]] <- nlsDeriv(expr[[1]], namevec[i], derivEnv = derivEnv, do_substitute = FALSE, verbose = verbose, ...)
  names(exprs) <- c(".value", namevec)
  
  # If a Hessian is requested, extend exprs to hold its entries
  if (hessian) {
    m <- length(exprs)
	  length(exprs) <- m + n*(n+1)/2
	  for (i in seq_len(n))
	    for (j in seq_len(n-i+1) + i-1) {
		    m <- m + 1
		    exprs[[m]] <- nlsDeriv(exprs[[i + 1]], namevec[j], derivEnv = derivEnv, 
		                    do_substitute = FALSE, verbose = verbose, ...)
	    }
  }
  
  # Find common subexpressions in the result 
  exprs <- as.expression(exprs)
  subexprs <- findSubexprs(exprs)
  m <- length(subexprs)
  final <- subexprs[[m]]
  subexprs[[m]] <- substitute(.value <- expr, list(expr = final[[".value"]]))
  subexprs[[m+1]] <- substitute(.grad <- array(0, c(length(.value), namelen), list(NULL, namevec)),
				  list(namevec = namevec, namelen = length(namevec)))
  if (hessian)
    subexprs[[m+2]] <- substitute(.hessian <- array(0, c(length(.value), namelen, namelen), 
					list(NULL, namevec, namevec)), 
					list(namelen = length(namevec), namevec = namevec))
  m <- length(subexprs)
  for (i in seq_len(n))
    subexprs[[m+i]] <- substitute(.grad[, name] <- expr,
			          list(name = namevec[i], expr = final[[namevec[i]]]))
  m <- length(subexprs)
  h <- 0
  if (hessian) {
    for (i in seq_len(n))
      for (j in seq_len(n-i+1) + i-1) {
		    h <- h + 1
		    if (i == j)
		      subexprs[[m + h]] <- substitute(.hessian[, i, i] <- expr,
				    list(i = namevec[i], expr = final[[1 + n + h]]))
		    else
		      subexprs[[m + h]] <- substitute(.hessian[, i, j] <- .hessian[, j, i] <- expr,
				    list(i = namevec[i], j = namevec[j], expr = final[[1 + n + h]]))
	    }	
	  h <- h + 1
	  subexprs[[m + h]] <- quote(attr(.value, "hessian") <- .hessian)
  }
  m <- length(subexprs)
  subexprs[[m+1]] <- quote(attr(.value, "gradient") <- .grad)
  subexprs[[m+2]] <- quote(.value)
  
  # Return the version that reuses common subexpressions
  subexprs
}

fnDeriv <- function(expr, namevec, args = all.vars(expr), env = environment(expr),
                    do_substitute = FALSE, verbose = FALSE, ...) {
  # Create a dummy function
  fn <- function() NULL
  
  # Set the body to the sequence of computations returned by codeDeriv
  expr <- dex(expr, do_substitute = do_substitute, verbose = verbose)
  body(fn) <- codeDeriv(expr, namevec, do_substitute = FALSE, verbose = verbose, ...)
  
  # Set the args 
  if (is.character(args)) {
    formals <- rep(list(bquote()), length = length(args))
    names(formals) <- args
    args <- formals
  }
  formals(fn) <- args
  
  # Set the environment and return it
  if (is.null(env))
    environment(fn) <- parent.frame()
  else
    environment(fn) <- as.environment(env)
  fn
}

isZERO <- function(x) is.numeric(x) && length(x) == 1 && x == 0

isONE  <- function(x) is.numeric(x) && length(x) == 1 && x == 1

isMINUSONE <- function(x) is.numeric(x) && length(x) == 1 && x == -1

isCALL <- function(x, name) is.call(x) && as.character(x[[1]]) == name


nlsSimplify <- function(expr, simpEnv = sysSimplifications, verbose = FALSE) {
    
    if (is.expression(expr))
    	return(as.expression(lapply(expr, nlsSimplify, simpEnv, verbose = verbose)))    

    if (is.call(expr)) {    	
	for (i in seq_along(expr)[-1])
	    expr[[i]] <- nlsSimplify(expr[[i]], simpEnv, verbose = verbose)
	fn <- as.character(expr[[1]])
	nargs <- length(expr) - 1
	while (!identical(simpEnv, emptyenv())) {
	    simps <- simpEnv[[fn]]
	    if (nargs > length(simps))
		return(expr)
	    simpn <- simps[[nargs]]
	    for (i in seq_along(simpn)) {  	
		argnames <- simpn[[i]]$argnames
		substitutions <- lapply(seq_along(argnames)+1L, function(i) expr[[i]])
		names(substitutions) <- argnames
		test <- simpn[[i]]$test
		if (with(substitutions, eval(test))) {
		    if (verbose) {
	    		message(paste("Simplifying", deparse(expr)))		    	
	    		message("Applying simplification:")
		    	print(simpn[[i]])
		    }
		    simplification <- do.call(substitute, list(simpn[[i]]$simplification, substitutions))
	   	    if (simpn[[i]]$do_eval)
	    		simplification <- eval(simplification)
	    	    return(simplification)
	        }
	     }
	     simpEnv <- parent.env(simpEnv)
	}
    }
    expr
}

findSubexprs <- function(expr, simplify = FALSE, tag = ".expr", verbose = FALSE, ...) {
    digests <- new.env(parent = emptyenv())
    subexprs <- list()
    subcount <- 0
    
    record <- function(index) {
        if (simplify)
	         expr[[index]] <<- subexpr <- nlsSimplify(expr[[index]], 
	                     verbose = verbose, ...)
      	else
     	    subexpr <- expr[[index]]
        	if (is.call(subexpr)) {
	           digest <- digest(subexpr)
	           for (i in seq_along(subexpr))
		         record(c(index,i))
	           prev <- digests[[digest]]
	           if (is.null(prev)) assign(digest, index, envir = digests)
        	   else if (is.numeric(prev))  { # the index where we last saw this
      	        subcount <<- subcount + 1
		            name <- as.name(paste0(tag, subcount))
		            assign(digest, name, envir = digests)
        	   }
	        }
     }
    
     edit <- function(index) {
	     subexpr <- expr[[index]]
	     if (is.call(subexpr)) {
	       digest <- digest(subexpr)	    
	       for (i in seq_along(subexpr)) edit(c(index,i))
    	   prev <- digests[[digest]]
	       if (is.name(prev)) {
		       num <- as.integer(substring(as.character(prev), nchar(tag)+1L))
		       subexprs[[num]] <<- call("<-", prev, expr[[index]])
		       expr[[index]] <<- prev
	       } 
	     }
     }
     for (i in seq_along(expr)) record(i)
     for (i in seq_along(expr)) edit(i)
     result <- quote({})
     result[seq_along(subexprs)+1] <- subexprs
     result[[length(result)+1]] <- expr
     result
}
    
# These are the derivatives supported by deriv()

newDeriv(log(x, base = exp(1)), 
         D(x)/(x*log(base)) - (1/base)*log(x, base)/log(base)*D(base))
newDeriv(exp(x), exp(x)*D(x))
newDeriv(sin(x), cos(x)*D(x))
newDeriv(cos(x), -sin(x)*D(x))
newDeriv(tan(x), 1/cos(x)^2*D(x))
newDeriv(sinh(x), cosh(x)*D(x))
newDeriv(cosh(x), sinh(x)*D(x))
newDeriv(sqrt(x), D(x)/2/sqrt(x))
newDeriv(pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE),
  (if (lower.tail && !log.p) 
  	dnorm((q-mean)/sd)*(D(q)/sd - D(mean)/sd - D(sd)*(q-mean)/sd^2) 
  else if (!lower.tail && !log.p) 
  	-dnorm((q-mean)/sd)*(D(q)/sd - D(mean)/sd - D(sd)*(q-mean)/sd^2) 
  else if (lower.tail && log.p) 
  	(D(q)/sd - D(mean)/sd - D(sd)*(q-mean)/sd^2)*exp(dnorm((q-mean)/sd, log = TRUE) 
  		- pnorm((q-mean)/sd, log = TRUE))
  else if (!lower.tail && log.p)	
  	-(D(q)/sd - D(mean)/sd - D(sd)*(q-mean)/sd^2)*exp(dnorm((q-mean)/sd, log = TRUE) 
  		- pnorm((q-mean)/sd, lower.tail = FALSE, log = TRUE)))
  + stop("cannot take derivative wrt 'lower.tail' or 'log.p'")*(D(lower.tail) + D(log.p)))
newDeriv(dnorm(x, mean = 0, sd = 1, log = FALSE),
  (if (!log) 
  	dnorm((x-mean)/sd)/sd^2*((D(mean)-D(x))*(x-mean)/sd + D(sd)*((x-mean)^2/sd^2 - 1)) 
  else if (log)
        ((D(mean)-D(x))*(x-mean)/sd + D(sd)*((x-mean)^2/sd^2 - 1))/sd) 
  + stop("cannot take derivative wrt 'log'")*D(log))
newDeriv(asin(x), D(x)/sqrt(1+x^2))
newDeriv(acos(x), -D(x)/sqrt(1+x^2))
newDeriv(atan(x), D(x)/(1+x^2))
newDeriv(gamma(x), gamma(x)*digamma(x)*D(x))
newDeriv(lgamma(x), digamma(x)*D(x))
newDeriv(digamma(x), trigamma(x)*D(x))
newDeriv(trigamma(x), psigamma(x, 2L)*D(x))
newDeriv(psigamma(x, deriv = 0L), 
          psigamma(x, deriv + 1L)*D(x) 
        + stop("cannot take derivative wrt 'deriv'")*D(deriv))

newDeriv(x*y, x*D(y) + D(x)*y)
newDeriv(x/y, D(x)/y - x*D(y)/y^2)
newDeriv(x^y, y*x^(y-1)*D(x) + x^y*log(x)*D(y))
newDeriv((x), D(x))
# Need to be careful with unary + or -
newDeriv(`+`(x, y = .MissingVal), if (missing(y)) D(x) else D(x) + D(y))
newDeriv(`-`(x, y = .MissingVal), if (missing(y)) -D(x) else D(x) - D(y))

# These are new

newDeriv(abs(x), sign(x)*D(x))
newDeriv(sign(x), 0)

newDeriv(`~`(x, y = .MissingVal), if (missing(y)) D(x) else D(y))

# Now, the simplifications

newSimplification(+a, TRUE, a)
newSimplification(-a, is.numeric(a), -a, do_eval = TRUE)
newSimplification(-a, isCALL(a, "-") && length(a) == 2, quote(a)[[2]], do_eval = TRUE)

newSimplification(exp(a), isCALL(a, "log") && length(a) == 2, quote(a)[[2]], do_eval = TRUE)
newSimplification(exp(a), is.numeric(a), exp(a), do_eval = TRUE)

newSimplification(log(a), isCALL(a, "exp") && length(a) == 2, quote(a)[[2]], do_eval = TRUE)
newSimplification(log(a), is.numeric(a), log(a), do_eval = TRUE)

newSimplification(!a, isTRUE(a), FALSE)
newSimplification(!a, isFALSE(a), TRUE)

newSimplification((a), TRUE, a)

newSimplification(a + b, isZERO(b), a)
newSimplification(a + b, isZERO(a), b)
newSimplification(a + b, identical(a, b), nlsSimplify(quote(2*a)), do_eval = TRUE)
newSimplification(a + b, is.numeric(a) && is.numeric(b), a+b, do_eval = TRUE)
# Add these to support our error scheme, don't test for stop() everywhere.
newSimplification(a + b, isCALL(a, "stop"), a)
newSimplification(a + b, isCALL(b, "stop"), b)

newSimplification(a - b, isZERO(b), a)
newSimplification(a - b, isZERO(a), nlsSimplify(quote(-b)), do_eval = TRUE)
newSimplification(a - b, identical(a, b), 0)
newSimplification(a - b, is.numeric(a) && is.numeric(b), a - b, do_eval = TRUE)

newSimplification(a * b, isZERO(a), 0)
newSimplification(a * b, isZERO(b), 0)
newSimplification(a * b, isONE(a), b)
newSimplification(a * b, isONE(b), a)
newSimplification(a * b, isMINUSONE(a), nlsSimplify(quote(-b)), do_eval = TRUE)
newSimplification(a * b, isMINUSONE(b), nlsSimplify(quote(-a)), do_eval = TRUE)
newSimplification(a * b, is.numeric(a) && is.numeric(b), a * b, do_eval = TRUE)

newSimplification(a / b, isONE(b), a)
newSimplification(a / b, isMINUSONE(b), nlsSimplify(quote(-a)), do_eval = TRUE)
newSimplification(a / b, isZERO(a), 0)
newSimplification(a / b, is.numeric(a) && is.numeric(b), a / b, do_eval = TRUE)

newSimplification(a ^ b, isONE(b), a)
newSimplification(a ^ b, is.numeric(a) && is.numeric(b), a ^ b, do_eval = TRUE)

newSimplification(log(a, base), isCALL(a, "exp"), nlsSimplify(call("/", quote(a)[[2]], quote(log(base)))), do_eval = TRUE)

newSimplification(a && b, isFALSE(a) || isFALSE(b), FALSE)
newSimplification(a && b, isTRUE(a), b)
newSimplification(a && b, isTRUE(b), a)

newSimplification(a || b, isTRUE(a) || isTRUE(b), TRUE)
newSimplification(a || b, isFALSE(a), b)
newSimplification(a || b, isFALSE(b), a)

newSimplification(if (cond) a, isTRUE(cond), a)
newSimplification(if (cond) a, isFALSE(cond), NULL)

newSimplification(if (cond) a else b, isTRUE(cond), a)
newSimplification(if (cond) a else b, isFALSE(cond), b)
newSimplification(if (cond) a else b, identical(a, b), a)

# This one is used to fix up the unary -
newSimplification(missing(a), identical(a, quote(.MissingVal)), TRUE)
newSimplification(missing(a), !identical(a, quote(.MissingVal)), FALSE)

