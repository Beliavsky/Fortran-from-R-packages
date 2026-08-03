## Return gradient of an object

gradient <- function(x, ...)
    UseMethod("gradient")

gradient.maxim <- function(x, ...)
    x$gradient
