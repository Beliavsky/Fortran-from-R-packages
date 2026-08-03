maxValue <- function(x, ...)
   UseMethod("maxValue")

maxValue.maxim <- function(x, ...)
   x$maximum
