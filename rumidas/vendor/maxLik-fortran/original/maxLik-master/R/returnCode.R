### Returns return code of maxim objects
### This is tells either error, or other cause the iterations ended,
### such as the result converged

returnCode <- function(x, ...)
    UseMethod("returnCode")

returnCode.default <- function(x, ...)
    x$returnCode

returnCode.maxim <- function(x, ...)
    x$code

returnCode.maxLik <- function(x, ...)
    x$code
