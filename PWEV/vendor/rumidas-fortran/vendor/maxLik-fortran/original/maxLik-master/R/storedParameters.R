## Return the stored parameters in a 'maxim' object

storedParameters <- function(x, ...)
    ## stored parameter values at each epoch/iteration
    UseMethod("storedParameters")

storedParameters.maxim <- function(x, ...)
    x$parameterStore
