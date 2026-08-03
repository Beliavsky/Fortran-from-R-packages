## Return the stored values in 'maxim' object

storedValues <- function(x, ...)
    ## stored optimization values at each iteration
    UseMethod("storedValues")

storedValues.maxim <- function(x, ...)
    x$valueStore
