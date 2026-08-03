## Return the objective function, used for optimization

objectiveFn <- function(x, ...)
    UseMethod("objectiveFn")

objectiveFn.maxim <- function(x, ...)
    x$objectiveFn
