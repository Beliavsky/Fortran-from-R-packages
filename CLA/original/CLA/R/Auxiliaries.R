## In next version of  package 'sfsmisc' - for now here (not exported)
funEnv <- function(..., envir = NULL, parent = parent.frame(),
                   hash = (...length() > 100), size = max(29L, ...length())) {
    e <- list2env(list(...), envir=envir, parent=parent, hash=hash, size=size)
    for(n in names(e)) ## iff function or formula, set environment to 'e':
	if(is.function(e[[n]]) || (is.call(e[[n]]) &&
				   inherits(e[[n]], "formula")))
	    environment(e[[n]]) <- e
    e
}

if(!is.function(.BaseNamespaceEnv$...length)) # ...length() only in R >= 3.5.0
    ## kludgy substitute, using parent.env() -- but it works (sometimes) in funEnv()
    ...length <- function() eval(quote(length(list(...))), envir = parent.frame())
