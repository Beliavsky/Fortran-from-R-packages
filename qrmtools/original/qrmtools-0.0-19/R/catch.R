### Catching warnings and errors ###############################################

##' @title Catching Results, Warnings and Errors Simultaneously
##' @param expr assignment or function evaluation
##' @return list with components:
##'         'value'  : value of expr or simpleError
##'	    'warning': simpleWarning or NULL
##'         'error'  : simpleError or NULL
##' @author Marius Hofert (see simsalapar)
##' @note https://stat.ethz.ch/pipermail/r-help/2010-December/262626.html
catch <- function(expr)
{
    W <- NULL
    w.handler <- function(w) { # warning handler
	W <<- w
	invokeRestart("muffleWarning")
    }
    res <- list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                            warning = w.handler), warning = W)
    is.err <- is(val <- res$value, "simpleError") # logical indicating an error
    list(value   = if(is.err) NULL else val, # value (or NULL in case of error)
	 warning = res$warning, # warning (or NULL)
	 error   = if(is.err) val else NULL) # error (or NULL if okay)
}
